#include "detail/pipeline.hpp"

namespace mito::detail {

void merge_haplogroup_markers(std::vector<ObservedPhyloMutation> &markers,
                              std::vector<ObservedPhyloMutation> additional) {
  markers.insert(markers.end(), std::make_move_iterator(additional.begin()),
                 std::make_move_iterator(additional.end()));
  std::sort(markers.begin(), markers.end(), observed_phylo_mutation_less);
  markers.erase(std::unique(markers.begin(), markers.end(),
                            [](const auto &lhs, const auto &rhs) {
                              return lhs.encoded == rhs.encoded;
                            }),
                markers.end());
}

[[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>>
ranges_from_bitmap(const std::vector<bool> &positions) {
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  std::size_t position = 1U;
  while (position < positions.size()) {
    while (position < positions.size() && !positions[position]) {
      ++position;
    }
    if (position == positions.size()) {
      break;
    }
    const std::size_t start = position;
    while (position + 1U < positions.size() && positions[position + 1U]) {
      ++position;
    }
    ranges.emplace_back(start, position);
    ++position;
  }
  return ranges;
}

[[nodiscard]] std::optional<std::vector<std::pair<std::size_t, std::size_t>>>
phylo_ranges_from_header(const ReadRecord &read) {
  std::istringstream iss(read.id);
  std::string token;
  while (iss >> token) {
    constexpr std::string_view prefix = "phylo_range=";
    if (token.rfind(prefix, 0) != 0) {
      continue;
    }
    const std::string_view encoded =
        std::string_view(token).substr(prefix.size());
    if (encoded.empty()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "empty development phylo_range for read '" + read.id +
                              "'");
    }
    std::vector<bool> positions(
        static_cast<std::size_t>(kDefaultReferenceLength) + 1U, false);
    std::size_t item_start = 0U;
    while (item_start < encoded.size()) {
      const auto item_end = encoded.find(',', item_start);
      const auto item =
          encoded.substr(item_start, item_end == std::string_view::npos
                                         ? encoded.size() - item_start
                                         : item_end - item_start);
      const auto dash = item.find('-');
      const auto start = parse_size(item.substr(0, dash));
      const auto end = dash == std::string_view::npos
                           ? start
                           : parse_size(item.substr(dash + 1U));
      if (!start || !end || *start == 0U || *end < *start ||
          *end > static_cast<std::size_t>(kDefaultReferenceLength)) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "invalid development phylo_range for read '" +
                                read.id + "'");
      }
      for (std::size_t position = *start; position <= *end; ++position) {
        positions[position] = true;
      }
      if (item_end == std::string_view::npos) {
        break;
      }
      item_start = item_end + 1U;
    }
    return ranges_from_bitmap(positions);
  }
  return std::nullopt;
}

void append_callable_position(
    std::vector<std::pair<std::size_t, std::size_t>> &ranges,
    std::optional<std::pair<std::size_t, std::size_t>> &current,
    std::size_t position) {
  if (!current) {
    current = std::pair{position, position};
    return;
  }
  if (position == current->second) {
    return;
  }
  if (current->second < std::numeric_limits<std::size_t>::max() &&
      position == current->second + 1U) {
    current->second = position;
    return;
  }
  ranges.push_back(*current);
  current = std::pair{position, position};
}

[[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>>
normalize_ranges(std::vector<std::pair<std::size_t, std::size_t>> ranges) {
  if (ranges.empty()) {
    return ranges;
  }
  std::sort(ranges.begin(), ranges.end());
  std::vector<std::pair<std::size_t, std::size_t>> normalized;
  normalized.reserve(ranges.size());
  for (const auto &range : ranges) {
    if (normalized.empty() ||
        (normalized.back().second < std::numeric_limits<std::size_t>::max() &&
         range.first > normalized.back().second + 1U)) {
      normalized.push_back(range);
      continue;
    }
    normalized.back().second = std::max(normalized.back().second, range.second);
  }
  return normalized;
}

[[nodiscard]] CallablePhyloRanges
callable_phylo_ranges_from_alignment(const ReadRecord &read,
                                     const std::string &reference,
                                     const AnalysisConfig &config) {
  CallablePhyloRanges result;
  if (reference.empty() || read.reference_start == 0U ||
      read.cigar_operations.empty() ||
      !passes_snp_alignment_filters(read, config)) {
    return result;
  }
  result.known = true;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  std::optional<std::pair<std::size_t, std::size_t>> current;
  std::size_t read_cursor = 0U;
  std::size_t reference_cursor =
      ((read.reference_start - 1U) % reference.size()) + 1U;
  constexpr std::size_t kMaxPhylogeneticIndelLength = 50U;
  for (const auto &operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;
    if (op == 'M' || op == '=' || op == 'X') {
      for (std::size_t offset = 0;
           offset < len && read_cursor + offset < read.sequence.size();
           ++offset) {
        const std::size_t position =
            ((reference_cursor - 1U + offset) % reference.size()) + 1U;
        const char base = static_cast<char>(std::toupper(
            static_cast<unsigned char>(read.sequence[read_cursor + offset])));
        const auto quality = phred_quality_at(read, read_cursor + offset);
        if (is_base(base) && is_base(reference[position - 1U]) && quality &&
            *quality >= config.min_base_quality) {
          append_callable_position(ranges, current, position);
        }
      }
      read_cursor += len;
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'I' || op == 'S') {
      read_cursor += len;
    } else if (op == 'D') {
      if (len <= kMaxPhylogeneticIndelLength) {
        for (std::size_t offset = 0; offset < len; ++offset) {
          append_callable_position(
              ranges, current,
              ((reference_cursor - 1U + offset) % reference.size()) + 1U);
        }
      }
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'N') {
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'H' || op == 'P') {
      // Does not consume query or reference sequence.
    }
  }
  if (current) {
    ranges.push_back(*current);
  }
  result.ranges = normalize_ranges(std::move(ranges));
  return result;
}

[[nodiscard]] const ReadRecord &
source_record_for_fragment(const MoleculeAssemblyResult &assembly,
                           const AlignmentFragmentEvidence &fragment) {
  if (!assembly.source_alignment_records.empty()) {
    if (fragment.source_record_index >=
        assembly.source_alignment_records.size()) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "alignment fragment source record is unresolved");
    }
    return assembly.source_alignment_records[fragment.source_record_index];
  }
  if (fragment.molecule_index.value >= assembly.representatives.size()) {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "FASTQ fragment source record is unresolved");
  }
  return assembly.representatives[fragment.molecule_index.value];
}

[[nodiscard]] std::string
alignment_callability_exclusion(const ReadRecord &read,
                                const AnalysisConfig &config) {
  if ((read.flags & 0x4U) != 0U) {
    return "UNMAPPED_ALIGNMENT";
  }
  if ((read.flags & 0x100U) != 0U) {
    return "SECONDARY_ALIGNMENT";
  }
  if ((read.flags & 0x200U) != 0U) {
    return "QC_FAILED_ALIGNMENT";
  }
  if ((read.flags & 0x400U) != 0U) {
    return "DUPLICATE_ALIGNMENT";
  }
  const auto configured_exclusions =
      static_cast<std::uint16_t>(config.excluded_snp_flags & ~0x800U);
  if ((read.flags & configured_exclusions) != 0U) {
    return "EXCLUDED_ALIGNMENT_FLAGS";
  }
  if (read.mapping_quality < config.min_mapping_quality) {
    return "LOW_MAPPING_QUALITY";
  }
  if (read.reference_start == 0U || read.reference_name.empty() ||
      read.reference_name == "*") {
    return "MISSING_REFERENCE_PLACEMENT";
  }
  if (looks_like_nuclear_contig(read.reference_name)) {
    return "NON_MITOCHONDRIAL_ALIGNMENT";
  }
  if (read.cigar_operations.empty()) {
    return "MISSING_CIGAR";
  }
  if (read.sequence.empty() || read.qualities.empty()) {
    return "MISSING_QUERY_OR_QUALITIES";
  }
  return {};
}

[[nodiscard]] AlignmentCallabilityEvidence
alignment_callability(const AlignmentFragmentEvidence &fragment,
                      const ReadRecord &read, const std::string &reference,
                      const AnalysisConfig &config) {
  AlignmentCallabilityEvidence evidence;
  evidence.alignment_fragment_id = fragment.id;
  const auto excluded = alignment_callability_exclusion(read, config);
  if (!excluded.empty()) {
    evidence.status = excluded;
    return evidence;
  }
  if (reference.empty()) {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "callability requires a non-empty reference");
  }

  evidence.eligible = true;
  evidence.status = "ASSESSED";
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  std::optional<std::pair<std::size_t, std::size_t>> current;
  std::size_t query_cursor = 0U;
  std::size_t reference_cursor =
      ((read.reference_start - 1U) % reference.size()) + 1U;
  for (const auto &operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;
    if (op == 'M' || op == '=' || op == 'X') {
      for (std::size_t offset = 0U; offset < len; ++offset) {
        const auto position =
            ((reference_cursor - 1U + (offset % reference.size())) %
             reference.size()) +
            1U;
        if (query_cursor > std::numeric_limits<std::size_t>::max() - offset) {
          throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                              "CIGAR query coordinate overflow for read '" +
                                  read.id + "'");
        }
        const auto query_index = query_cursor + offset;
        if (query_index >= read.sequence.size() ||
            query_index >= read.qualities.size()) {
          ++evidence
                .reference_exclusion_counts["MISSING_QUERY_BASE_OR_QUALITY"];
          continue;
        }
        const char observed = static_cast<char>(std::toupper(
            static_cast<unsigned char>(read.sequence[query_index])));
        const auto quality = phred_quality_at(read, query_index);
        if (!is_base(reference[position - 1U])) {
          ++evidence.reference_exclusion_counts["NON_CANONICAL_REFERENCE_BASE"];
        } else if (!is_base(observed)) {
          ++evidence.reference_exclusion_counts["NON_CANONICAL_OBSERVED_BASE"];
        } else if (!quality || *quality < config.min_base_quality) {
          ++evidence.reference_exclusion_counts["LOW_OR_INVALID_BASE_QUALITY"];
        } else {
          append_callable_position(ranges, current, position);
        }
      }
      if (len > std::numeric_limits<std::size_t>::max() - query_cursor) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "CIGAR query coordinate overflow for read '" +
                                read.id + "'");
      }
      query_cursor += len;
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'I' || op == 'S') {
      if (len > std::numeric_limits<std::size_t>::max() - query_cursor) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "CIGAR query coordinate overflow for read '" +
                                read.id + "'");
      }
      query_cursor += len;
      if (op == 'I') {
        if (len > std::numeric_limits<std::size_t>::max() -
                      evidence.inserted_query_bases) {
          throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                              "inserted query-base count overflow for read '" +
                                  read.id + "'");
        }
        evidence.inserted_query_bases += len;
        evidence.disrupted_adjacency_anchors.push_back(
            reference_cursor > 1U ? reference_cursor - 1U : reference.size());
      } else {
        if (len > std::numeric_limits<std::size_t>::max() -
                      evidence.soft_clipped_query_bases) {
          throw AnalysisError(
              AnalysisErrorCode::input_parse_failed,
              "soft-clipped query-base count overflow for read '" + read.id +
                  "'");
        }
        evidence.soft_clipped_query_bases += len;
      }
    } else if (op == 'D' || op == 'N') {
      auto &excluded_count =
          evidence.reference_exclusion_counts
              [op == 'D' ? "DELETION_FROM_REFERENCE_PATH" : "REFERENCE_SKIP"];
      if (len > std::numeric_limits<std::size_t>::max() - excluded_count) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "reference exclusion count overflow for read '" +
                                read.id + "'");
      }
      excluded_count += len;
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'H' || op == 'P') {
      // Neither query nor reference is consumed.
    }
  }
  if (current) {
    ranges.push_back(*current);
  }
  evidence.ranges = normalize_ranges(std::move(ranges));
  std::sort(evidence.disrupted_adjacency_anchors.begin(),
            evidence.disrupted_adjacency_anchors.end());
  evidence.disrupted_adjacency_anchors.erase(
      std::unique(evidence.disrupted_adjacency_anchors.begin(),
                  evidence.disrupted_adjacency_anchors.end()),
      evidence.disrupted_adjacency_anchors.end());
  for (const auto &[start, end] : evidence.ranges) {
    evidence.callable_bases += end - start + 1U;
  }
  return evidence;
}

[[nodiscard]] CallabilityResult
build_callability(const MoleculeAssemblyResult &assembly,
                  const std::vector<ReadFeature> &features,
                  const std::string &reference, const AnalysisConfig &config) {
  CallabilityResult result;
  result.molecules.reserve(assembly.molecules.size());
  for (std::size_t molecule_index = 0U;
       molecule_index < assembly.molecules.size(); ++molecule_index) {
    throw_if_cancelled(config);
    const auto &molecule = assembly.molecules[molecule_index];
    MoleculeCallabilityEvidence summary;
    summary.molecule_index = molecule.index;
    if (molecule_index >= features.size()) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "molecule callability feature reference is unresolved");
    }
    if (features[molecule_index].filtered_numt) {
      summary.status = "EXCLUDED_NUMT";
      result.molecules.push_back(std::move(summary));
      continue;
    }

    std::vector<std::pair<std::size_t, std::size_t>> molecule_ranges;
    summary.alignments.reserve(molecule.fragment_ids.size());
    for (const auto fragment_id : molecule.fragment_ids) {
      if (fragment_id.value >= assembly.fragments.size()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "molecule callability fragment reference is unresolved");
      }
      const auto &fragment = assembly.fragments[fragment_id.value];
      const auto &read = source_record_for_fragment(assembly, fragment);
      auto alignment = alignment_callability(fragment, read, reference, config);
      if (alignment.eligible) {
        summary.known = true;
        molecule_ranges.insert(molecule_ranges.end(), alignment.ranges.begin(),
                               alignment.ranges.end());
      }
      summary.alignments.push_back(std::move(alignment));
    }
    summary.ranges = normalize_ranges(std::move(molecule_ranges));
    for (const auto &[start, end] : summary.ranges) {
      summary.callable_bases += end - start + 1U;
    }
    summary.status =
        summary.known
            ? (summary.callable_bases == 0U ? "ASSESSED_EMPTY" : "ASSESSED")
            : "NOT_ASSESSABLE";
    result.molecules.push_back(std::move(summary));
  }
  return result;
}

[[nodiscard]] bool callable_position_in_ranges(
    const std::vector<std::pair<std::size_t, std::size_t>> &ranges,
    const std::size_t position) {
  const auto it =
      std::upper_bound(ranges.begin(), ranges.end(), position,
                       [](const std::size_t value, const auto &range) {
                         return value < range.first;
                       });
  return it != ranges.begin() && position <= std::prev(it)->second;
}

[[nodiscard]] bool
span_in_ranges(const std::vector<std::pair<std::size_t, std::size_t>> &ranges,
               const std::size_t start, const std::size_t end) {
  if (start == 0U || end < start) {
    return false;
  }
  const auto it =
      std::upper_bound(ranges.begin(), ranges.end(), start,
                       [](const std::size_t value, const auto &range) {
                         return value < range.first;
                       });
  return it != ranges.begin() && std::prev(it)->second >= end;
}

[[nodiscard]] bool
reference_adjacency_callable(const MoleculeCallabilityEvidence &callability,
                             const std::size_t left, const std::size_t right) {
  for (const auto &alignment : callability.alignments) {
    if (!alignment.eligible ||
        !callable_position_in_ranges(alignment.ranges, left) ||
        !callable_position_in_ranges(alignment.ranges, right)) {
      continue;
    }
    if ((right == left + 1U || right < left) &&
        !std::binary_search(alignment.disrupted_adjacency_anchors.begin(),
                            alignment.disrupted_adjacency_anchors.end(),
                            left)) {
      return true;
    }
  }
  return false;
}

[[nodiscard]] bool reference_span_with_flanks_callable(
    const MoleculeCallabilityEvidence &callability, const std::size_t start,
    const std::size_t end, const std::size_t reference_length) {
  if (start == 0U || end < start || end > reference_length) {
    return false;
  }
  const auto left = circular_previous_position(start, reference_length);
  const auto right = circular_next_position(end, reference_length);
  for (const auto &alignment : callability.alignments) {
    if (alignment.eligible && span_in_ranges(alignment.ranges, start, end) &&
        callable_position_in_ranges(alignment.ranges, left) &&
        callable_position_in_ranges(alignment.ranges, right)) {
      return true;
    }
  }
  return false;
}

} // namespace mito::detail
