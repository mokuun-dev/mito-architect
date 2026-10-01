#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] std::string alignment_as_sa(const ReadRecord &read) {
  const char strand = (read.flags & 0x10U) != 0U ? '-' : '+';
  const auto nm = read.aux_tags.find("NM");
  return read.reference_name + "," + std::to_string(read.reference_start) +
         "," + strand + "," + (read.cigar.empty() ? "*" : read.cigar) + "," +
         std::to_string(static_cast<unsigned int>(read.mapping_quality)) + "," +
         (nm == read.aux_tags.end() ? "0" : nm->second) + ";";
}

[[nodiscard]] std::string alignment_fragment_role(const std::uint16_t flags,
                                                  const bool alignment_input) {
  if (!alignment_input) {
    return "unaligned_read";
  }
  const bool secondary = (flags & 0x100U) != 0U;
  const bool supplementary = (flags & 0x800U) != 0U;
  if (secondary && supplementary) {
    return "secondary_and_supplementary";
  }
  if (secondary) {
    return "secondary";
  }
  if (supplementary) {
    return "supplementary";
  }
  return "primary_candidate";
}

void validate_protocol_tag_value(const std::string &tag,
                                 const std::string &value,
                                 const std::string &read_id) {
  if (value.empty() || value.size() > 256U ||
      std::any_of(value.begin(), value.end(), [](const char value_char) {
        const auto character = static_cast<unsigned char>(value_char);
        return std::iscntrl(character) != 0;
      })) {
    throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                        "invalid protocol tag " + tag + " on read '" + read_id +
                            "'");
  }
}

[[nodiscard]] ProtocolTagValues
protocol_tag_values(const std::vector<ReadRecord> &reads,
                    const std::vector<std::size_t> &indices,
                    const std::string &tag) {
  ProtocolTagValues result;
  if (tag.empty()) {
    return result;
  }
  for (const auto index : indices) {
    const auto found = reads[index].aux_tags.find(tag);
    if (found == reads[index].aux_tags.end()) {
      result.missing = true;
      continue;
    }
    validate_protocol_tag_value(tag, found->second, reads[index].id);
    result.values.insert(found->second);
  }
  return result;
}

[[nodiscard]] std::size_t reference_consuming_length(const ReadRecord &read) {
  std::size_t length = 0U;
  for (const auto &operation : read.cigar_operations) {
    if (operation.code != 'M' && operation.code != '=' &&
        operation.code != 'X' && operation.code != 'D' &&
        operation.code != 'N') {
      continue;
    }
    if (operation.length > std::numeric_limits<std::size_t>::max() - length) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "CIGAR reference span overflows for read '" +
                              read.id + "'");
    }
    length += operation.length;
  }
  return length;
}

void sort_unique_strings(std::vector<std::string> &values) {
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
}

[[nodiscard]] MoleculeAssemblyResult
assemble_molecules(InputRecords &input, const AnalysisConfig &config,
                   const std::size_t reference_length) {
  MoleculeAssemblyResult result;
  if (!input.alignment_input) {
    if (!config.molecule_id_tag.empty() || !config.umi_tag.empty() ||
        !config.duplex_tag.empty()) {
      throw AnalysisError(
          AnalysisErrorCode::invalid_configuration,
          "protocol SAM tags cannot be applied to unaligned FASTQ input");
    }
    result.representatives = std::move(input.reads);
    result.fragments.reserve(result.representatives.size());
    result.molecules.reserve(result.representatives.size());
    std::unordered_map<std::string, std::size_t> name_occurrences;
    name_occurrences.reserve(result.representatives.size());
    for (std::size_t index = 0; index < result.representatives.size();
         ++index) {
      const auto &read = result.representatives[index];
      const auto occurrence = ++name_occurrences[read.id];
      const std::string molecule_id =
          occurrence == 1U ? read.id
                           : read.id + "#record:" + std::to_string(index);
      result.fragments.push_back(
          {AlignmentFragmentId{index}, MoleculeIndex{index}, index, molecule_id,
           alignment_fragment_role(read.flags, false), true, read.flags,
           read.mapping_quality, read.reference_name, read.reference_start,
           read.cigar});
      MoleculeAssemblyEvidence molecule;
      molecule.index = MoleculeIndex{index};
      molecule.id = molecule_id;
      molecule.identity_policy = "fastq_record_proxy";
      molecule.assembly_status = "single_unaligned_fragment";
      molecule.fragment_ids.push_back(AlignmentFragmentId{index});
      molecule.representative_fragment_id = AlignmentFragmentId{index};
      molecule.primary_candidate_count = 1U;
      molecule.source_qnames.push_back(read.id);
      molecule.protocol_flags.emplace_back("READ_PROXY_IDENTITY");
      if (occurrence > 1U) {
        molecule.ambiguous = true;
        molecule.warnings.emplace_back("DUPLICATE_SOURCE_NAME_DISAMBIGUATED");
      }
      result.molecules.push_back(std::move(molecule));
    }
    return result;
  }

  for (auto &record : input.reads) {
    const auto sa = record.aux_tags.find("SA");
    if (sa != record.aux_tags.end()) {
      record.supplementary_alignments =
          parse_supplementary_alignments(sa->second, record.id);
    }
  }

  std::map<std::string, std::vector<std::size_t>> qname_groups;
  for (std::size_t index = 0; index < input.reads.size(); ++index) {
    qname_groups[input.reads[index].id].push_back(index);
  }

  struct AssemblyGroup {
    std::string molecule_id;
    std::vector<std::size_t> indices;
    bool missing_identity_tag = false;
    bool conflicting_identity_tag = false;
    bool partial_identity_tag = false;
  };
  std::map<std::string, AssemblyGroup> groups;
  for (const auto &[qname, indices] : qname_groups) {
    if (config.molecule_id_tag.empty()) {
      auto &group = groups["qname:" + qname];
      group.molecule_id = qname;
      group.indices.insert(group.indices.end(), indices.begin(), indices.end());
      continue;
    }
    const auto tag_values =
        protocol_tag_values(input.reads, indices, config.molecule_id_tag);
    if (tag_values.values.empty()) {
      auto &group = groups["missing:" + qname];
      group.molecule_id = "UNASSIGNED:" + config.molecule_id_tag + ":" + qname;
      group.missing_identity_tag = true;
      group.indices.insert(group.indices.end(), indices.begin(), indices.end());
      continue;
    }
    if (tag_values.values.size() != 1U) {
      auto &group = groups["conflict:" + qname];
      group.molecule_id = "CONFLICT:" + config.molecule_id_tag + ":" + qname;
      group.conflicting_identity_tag = true;
      group.indices.insert(group.indices.end(), indices.begin(), indices.end());
      continue;
    }
    const auto &value = *tag_values.values.begin();
    auto &group = groups["tag:" + value];
    group.molecule_id = config.molecule_id_tag + ":" + value;
    group.partial_identity_tag =
        group.partial_identity_tag || tag_values.missing;
    group.indices.insert(group.indices.end(), indices.begin(), indices.end());
  }

  result.representatives.reserve(groups.size());
  result.fragments.resize(input.reads.size());
  result.molecules.reserve(groups.size());
  for (const auto &[_, group] : groups) {
    const auto &molecule_id = group.molecule_id;
    const auto &indices = group.indices;
    const MoleculeIndex molecule_index{result.molecules.size()};
    std::vector<std::size_t> primary_candidates;
    primary_candidates.reserve(indices.size());
    std::map<std::string, std::size_t> qname_primary_counts;
    std::set<std::string> source_qnames;
    for (const auto index : indices) {
      source_qnames.insert(input.reads[index].id);
      if ((input.reads[index].flags & (0x100U | 0x800U)) == 0U) {
        primary_candidates.push_back(index);
        ++qname_primary_counts[input.reads[index].id];
      } else {
        qname_primary_counts.try_emplace(input.reads[index].id, 0U);
      }
    }
    std::size_t representative_index = primary_candidates.empty()
                                           ? indices.front()
                                           : primary_candidates.front();
    if (!config.molecule_id_tag.empty() && !primary_candidates.empty()) {
      representative_index = *std::max_element(
          primary_candidates.begin(), primary_candidates.end(),
          [&](const auto lhs, const auto rhs) {
            if (input.reads[lhs].mapping_quality !=
                input.reads[rhs].mapping_quality) {
              return input.reads[lhs].mapping_quality <
                     input.reads[rhs].mapping_quality;
            }
            return lhs > rhs;
          });
    }

    MoleculeAssemblyEvidence molecule;
    molecule.index = molecule_index;
    molecule.id = molecule_id;
    molecule.identity_policy = config.molecule_id_tag.empty()
                                   ? "sam_qname"
                                   : "sam_tag:" + config.molecule_id_tag;
    molecule.representative_fragment_id =
        AlignmentFragmentId{representative_index};
    molecule.primary_candidate_count = primary_candidates.size();
    molecule.source_qnames.assign(source_qnames.begin(), source_qnames.end());

    const bool malformed_qname_primary =
        std::any_of(qname_primary_counts.begin(), qname_primary_counts.end(),
                    [](const auto &entry) { return entry.second != 1U; });
    if (config.molecule_id_tag.empty()) {
      molecule.ambiguous = primary_candidates.size() != 1U;
      molecule.analysis_eligible = primary_candidates.size() == 1U;
      molecule.protocol_flags.emplace_back("READ_PROXY_IDENTITY");
      if (primary_candidates.empty()) {
        molecule.assembly_status = "fallback_without_primary";
        molecule.warnings.emplace_back("NO_PRIMARY_ALIGNMENT");
        molecule.exclusion_reasons.emplace_back("NO_PRIMARY_ALIGNMENT");
      } else if (primary_candidates.size() > 1U) {
        molecule.assembly_status = "first_of_multiple_primaries";
        molecule.warnings.emplace_back("MULTIPLE_PRIMARY_ALIGNMENTS");
        molecule.exclusion_reasons.emplace_back("MULTIPLE_PRIMARY_ALIGNMENTS");
      } else {
        molecule.assembly_status = "unique_primary";
      }
    } else {
      molecule.protocol_flags.emplace_back("EXPLICIT_MOLECULE_ID");
      molecule.protocol_metadata["molecule_id_tag"] = config.molecule_id_tag;
      if (!group.missing_identity_tag && !group.conflicting_identity_tag) {
        molecule.protocol_metadata["molecule_id_value"] =
            molecule_id.substr(config.molecule_id_tag.size() + 1U);
      }
      molecule.ambiguous = group.missing_identity_tag ||
                           group.conflicting_identity_tag ||
                           malformed_qname_primary;
      molecule.analysis_eligible = !molecule.ambiguous;
      if (group.missing_identity_tag) {
        molecule.assembly_status = "missing_explicit_molecule_id";
        molecule.warnings.emplace_back("MISSING_MOLECULE_ID_TAG");
        molecule.exclusion_reasons.emplace_back("MISSING_MOLECULE_ID_TAG");
      } else if (group.conflicting_identity_tag) {
        molecule.assembly_status = "conflicting_explicit_molecule_id";
        molecule.warnings.emplace_back("CONFLICTING_MOLECULE_ID_TAG");
        molecule.exclusion_reasons.emplace_back("CONFLICTING_MOLECULE_ID_TAG");
      } else if (malformed_qname_primary) {
        molecule.assembly_status = "invalid_tagged_fragment_group";
        molecule.warnings.emplace_back("INVALID_PRIMARY_COUNT_PER_QNAME");
        molecule.exclusion_reasons.emplace_back(
            "INVALID_PRIMARY_COUNT_PER_QNAME");
      } else if (source_qnames.size() > 1U) {
        molecule.assembly_status = "explicit_tag_multi_qname_group";
      } else {
        molecule.assembly_status = "explicit_tag_single_qname_group";
      }
      if (group.partial_identity_tag) {
        molecule.warnings.emplace_back("PARTIAL_MOLECULE_ID_TAG_INHERITED");
      }
    }

    const auto apply_optional_protocol_tag =
        [&](const std::string &tag, const std::string &name,
            const std::string &recorded_flag) {
          if (tag.empty()) {
            return;
          }
          const auto values = protocol_tag_values(input.reads, indices, tag);
          molecule.protocol_metadata[name + "_tag"] = tag;
          if (values.values.size() == 1U) {
            molecule.protocol_metadata[name + "_value"] =
                *values.values.begin();
            molecule.protocol_flags.push_back(recorded_flag);
            if (values.missing) {
              molecule.warnings.push_back("PARTIAL_" + name + "_TAG_INHERITED");
            }
          } else if (values.values.empty()) {
            molecule.warnings.push_back("MISSING_" + name + "_TAG");
          } else {
            molecule.ambiguous = true;
            molecule.analysis_eligible = false;
            molecule.warnings.push_back("CONFLICTING_" + name + "_TAG_VALUES");
            molecule.exclusion_reasons.push_back("CONFLICTING_" + name +
                                                 "_TAG_VALUES");
          }
        };
    apply_optional_protocol_tag(config.umi_tag, "UMI", "UMI_RECORDED");
    apply_optional_protocol_tag(config.duplex_tag, "DUPLEX",
                                "DUPLEX_METADATA_RECORDED");

    const bool has_eligible_primary = std::any_of(
        primary_candidates.begin(), primary_candidates.end(),
        [&](const auto index) {
          return (input.reads[index].flags & (0x4U | 0x200U | 0x400U)) == 0U;
        });
    if (!primary_candidates.empty() && !has_eligible_primary) {
      molecule.analysis_eligible = false;
      molecule.warnings.emplace_back("NO_ELIGIBLE_PRIMARY_EVIDENCE");
      molecule.exclusion_reasons.emplace_back("NO_ELIGIBLE_PRIMARY_EVIDENCE");
    }

    if (source_qnames.size() > 1U) {
      molecule.protocol_flags.emplace_back("MULTI_QNAME_MOLECULE");
    }
    std::set<std::string> reference_names;
    bool concatemer_candidate = false;
    bool origin_spanning = false;
    for (const auto index : indices) {
      const auto &record = input.reads[index];
      if ((record.flags & 0x100U) != 0U) {
        molecule.protocol_flags.emplace_back("SECONDARY_EVIDENCE_PRESENT");
      }
      if ((record.flags & 0x800U) != 0U) {
        molecule.protocol_flags.emplace_back("SUPPLEMENTARY_EVIDENCE_PRESENT");
      }
      if ((record.flags & 0x400U) != 0U) {
        molecule.protocol_flags.emplace_back("SAM_DUPLICATE_FLAG_PRESENT");
      }
      if (!record.reference_name.empty() && record.reference_name != "*") {
        reference_names.insert(record.reference_name);
      }
      if (reference_length != 0U &&
          record.sequence.size() / 2U >= reference_length) {
        concatemer_candidate = true;
      }
      const auto span = reference_consuming_length(record);
      if (reference_length != 0U && record.reference_start != 0U &&
          span != 0U &&
          span - 1U > reference_length -
                          std::min(reference_length, record.reference_start)) {
        origin_spanning = true;
      }
    }
    if (reference_names.size() > 1U) {
      molecule.protocol_flags.emplace_back("MULTI_REFERENCE_ALIGNMENT");
    }
    if (concatemer_candidate) {
      molecule.protocol_flags.emplace_back("CONCATEMER_LENGTH_CANDIDATE");
    }
    if (origin_spanning) {
      molecule.protocol_flags.emplace_back("ORIGIN_SPANNING_ALIGNMENT");
    }

    molecule.fragment_ids.reserve(indices.size());
    for (const auto index : indices) {
      const auto &record = input.reads[index];
      molecule.fragment_ids.push_back(AlignmentFragmentId{index});
      result.fragments[index] = {
          AlignmentFragmentId{index},
          molecule_index,
          index,
          molecule_id,
          alignment_fragment_role(record.flags, true),
          index == representative_index,
          record.flags,
          record.mapping_quality,
          record.reference_name,
          record.reference_start,
          record.cigar,
      };
    }

    ReadRecord representative = input.reads[representative_index];
    bool has_sa_tag = representative.aux_tags.contains("SA");
    std::string sa =
        has_sa_tag ? representative.aux_tags.at("SA") : std::string{};
    for (const auto index : indices) {
      const bool explicit_fragment_link =
          !config.molecule_id_tag.empty() &&
          (input.reads[index].flags & 0x100U) == 0U;
      if (index == representative_index ||
          ((input.reads[index].flags & 0x800U) == 0U &&
           !explicit_fragment_link) ||
          input.reads[index].reference_name.empty() ||
          input.reads[index].reference_name == "*") {
        continue;
      }
      const auto encoded = alignment_as_sa(input.reads[index]);
      if (encoded == alignment_as_sa(representative)) {
        continue;
      }
      if (sa.find(encoded) == std::string::npos) {
        sa += encoded;
      }
      has_sa_tag = true;
    }
    if (has_sa_tag) {
      representative.aux_tags["SA"] = std::move(sa);
      representative.supplementary_alignments = parse_supplementary_alignments(
          representative.aux_tags.at("SA"), representative.id);
    }
    representative.aux_tags["MA"] = std::to_string(indices.size());
    representative.id = molecule.id;
    sort_unique_strings(molecule.protocol_flags);
    sort_unique_strings(molecule.exclusion_reasons);
    sort_unique_strings(molecule.warnings);
    result.representatives.push_back(std::move(representative));
    result.molecules.push_back(std::move(molecule));
  }
  result.source_alignment_records = std::move(input.reads);
  return result;
}

[[nodiscard]] double mean_quality(const std::string &qualities) {
  if (qualities.empty()) {
    return 0.0;
  }
  const auto total = std::accumulate(
      qualities.begin(), qualities.end(), 0.0, [](double sum, char c) {
        return sum + static_cast<double>(std::max(0, c - 33));
      });
  return total / static_cast<double>(qualities.size());
}

[[nodiscard]] double gc_fraction(std::string_view sequence) {
  if (sequence.empty()) {
    return 0.0;
  }
  std::size_t gc = 0;
  for (const char c : sequence) {
    const char upper =
        static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (upper == 'G' || upper == 'C') {
      ++gc;
    }
  }
  return static_cast<double>(gc) / static_cast<double>(sequence.size());
}

} // namespace mito::detail
