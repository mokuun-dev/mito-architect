#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] std::string escape_json(std::string_view value) {
  std::string out;
  out.reserve(value.size() + 8);
  const auto append_byte_escape = [&out](const unsigned char byte) {
    static constexpr char hex[] = "0123456789abcdef";
    out += "\\u00";
    out.push_back(hex[byte >> 4U]);
    out.push_back(hex[byte & 0x0fU]);
  };
  for (std::size_t index = 0U; index < value.size();) {
    const auto byte = static_cast<unsigned char>(value[index]);
    const char c = static_cast<char>(byte);
    switch (c) {
    case '"':
      out += "\\\"";
      break;
    case '\\':
      out += "\\\\";
      break;
    case '\b':
      out += "\\b";
      break;
    case '\f':
      out += "\\f";
      break;
    case '\n':
      out += "\\n";
      break;
    case '\r':
      out += "\\r";
      break;
    case '\t':
      out += "\\t";
      break;
    default:
      if (byte < 0x20U) {
        std::ostringstream oss;
        oss << "\\u" << std::hex << std::setw(4) << std::setfill('0')
            << static_cast<int>(byte);
        out += oss.str();
      } else if (byte < 0x80U) {
        out += c;
      } else {
        // Preserve valid UTF-8. An invalid byte is represented losslessly as
        // a JSON unicode escape, so a hostile SAM/QNAME cannot make the
        // native JSON invalid for C or Rust callers.
        std::size_t length = 0U;
        const auto continuation = [&](const std::size_t offset) {
          return index + offset < value.size() &&
                 (static_cast<unsigned char>(value[index + offset]) & 0xc0U) ==
                     0x80U;
        };
        if (byte >= 0xc2U && byte <= 0xdfU && continuation(1U)) {
          length = 2U;
        } else if (byte == 0xe0U && continuation(1U) && continuation(2U) &&
                   static_cast<unsigned char>(value[index + 1U]) >= 0xa0U) {
          length = 3U;
        } else if (byte >= 0xe1U && byte <= 0xecU && continuation(1U) &&
                   continuation(2U)) {
          length = 3U;
        } else if (byte == 0xedU && continuation(1U) && continuation(2U) &&
                   static_cast<unsigned char>(value[index + 1U]) <= 0x9fU) {
          length = 3U;
        } else if (byte >= 0xeeU && byte <= 0xefU && continuation(1U) &&
                   continuation(2U)) {
          length = 3U;
        } else if (byte == 0xf0U && continuation(1U) && continuation(2U) &&
                   continuation(3U) &&
                   static_cast<unsigned char>(value[index + 1U]) >= 0x90U) {
          length = 4U;
        } else if (byte >= 0xf1U && byte <= 0xf3U && continuation(1U) &&
                   continuation(2U) && continuation(3U)) {
          length = 4U;
        } else if (byte == 0xf4U && continuation(1U) && continuation(2U) &&
                   continuation(3U) &&
                   static_cast<unsigned char>(value[index + 1U]) <= 0x8fU) {
          length = 4U;
        }
        if (length == 0U) {
          append_byte_escape(byte);
        } else {
          out.append(value.substr(index, length));
          index += length - 1U;
        }
      }
    }
    ++index;
  }
  return out;
}

[[nodiscard]] std::string quoted(std::string_view value) {
  return "\"" + escape_json(value) + "\"";
}

// Prefer the JSON escaper above over std::quoted when argument-dependent lookup
// sees std::string. The latter only escapes quotes/backslashes and is not a
// complete JSON control-character encoder.
[[nodiscard]] std::string quoted(const std::string &value) {
  return quoted(std::string_view(value));
}

[[nodiscard]] std::string quoted(std::string &value) {
  return quoted(std::string_view(value));
}

[[nodiscard]] std::string quoted(const char *value) {
  return quoted(std::string_view(value));
}

class JsonBuffer {
public:
  explicit JsonBuffer(const std::size_t reserve_bytes) {
    buffer_.reserve(reserve_bytes);
  }

  JsonBuffer &operator<<(const char *value) {
    buffer_.append(value);
    return *this;
  }

  JsonBuffer &operator<<(const std::string &value) {
    buffer_.append(value);
    return *this;
  }

  JsonBuffer &operator<<(const std::string_view value) {
    buffer_.append(value);
    return *this;
  }

  JsonBuffer &operator<<(const char value) {
    buffer_.push_back(value);
    return *this;
  }

  JsonBuffer &operator<<(const bool value) {
    buffer_.append(value ? "true" : "false");
    return *this;
  }

  template <typename Integer>
    requires(std::is_integral_v<Integer> && !std::is_same_v<Integer, bool> &&
             !std::is_same_v<Integer, char>)
  JsonBuffer &operator<<(const Integer value) {
    std::array<char, 32> encoded{};
    const auto [end, error] =
        std::to_chars(encoded.data(), encoded.data() + encoded.size(), value);
    if (error != std::errc{}) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "integer JSON serialization failed");
    }
    buffer_.append(encoded.data(), end);
    return *this;
  }

  JsonBuffer &operator<<(const double value) {
    std::array<char, 64> encoded{};
    const auto [end, error] =
        std::to_chars(encoded.data(), encoded.data() + encoded.size(), value,
                      std::chars_format::fixed, 9);
    if (error != std::errc{}) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "floating-point JSON serialization failed");
    }
    buffer_.append(encoded.data(), end);
    return *this;
  }

  [[nodiscard]] std::string take() && { return std::move(buffer_); }

private:
  std::string buffer_;
};

void write_string_array(JsonBuffer &out,
                        const std::vector<std::string> &values) {
  out << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    out << quoted(values[i]);
  }
  out << "]";
}

template <typename Integer>
  requires(std::is_integral_v<Integer>)
void write_integer_array(JsonBuffer &out, const std::vector<Integer> &values) {
  out << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    out << values[i];
  }
  out << "]";
}

[[nodiscard]] bool has_clinical_annotation(const SnpCall &snp) {
  return !snp.clinical_assertions.empty();
}

void write_clinical_annotation(JsonBuffer &out, const SnpCall &snp) {
  out << "{"
      << "\"schema_version\":" << quoted(kClinicalAnnotationSchemaVersion)
      << ","
      << "\"conflict_status\":" << quoted(snp.clinical_conflict_status) << ","
      << "\"consensus_significance\":"
      << quoted(snp.clinical_consensus_significance) << ","
      << "\"pathogenicity\":" << quoted(snp.pathogenicity) << ","
      << "\"phenotype\":" << quoted(snp.phenotype) << ","
      << "\"references\":";
  write_string_array(out, snp.references);
  out << ",\"sources\":";
  write_string_array(out, snp.sources);
  out << ",\"source\":" << quoted("local-cache") << ","
      << "\"assertions\":[";
  for (std::size_t i = 0; i < snp.clinical_assertions.size(); ++i) {
    if (i != 0U) {
      out << ",";
    }
    const auto &assertion = snp.clinical_assertions[i];
    out << "{"
        << "\"source\":" << quoted(assertion.source) << ","
        << "\"assertion_id\":" << quoted(assertion.assertion_id) << ","
        << "\"allele_id\":" << quoted(assertion.allele_id) << ","
        << "\"disease\":" << quoted(assertion.disease) << ","
        << "\"clinical_significance\":"
        << quoted(assertion.clinical_significance) << ","
        << "\"normalized_significance\":"
        << quoted(assertion.normalized_significance) << ","
        << "\"review_status\":" << quoted(assertion.review_status) << ","
        << "\"assertion_date\":" << quoted(assertion.assertion_date) << ","
        << "\"source_url\":" << quoted(assertion.source_url) << ","
        << "\"references\":";
    write_string_array(out, assertion.references);
    out << ",\"resource_version\":" << quoted(assertion.resource_version) << ","
        << "\"retrieved_at\":" << quoted(assertion.retrieved_at) << "}";
  }
  out << "]";
  if (!snp.clinvar_allele_id.empty()) {
    out << ",\"clinvar_allele_id\":" << quoted(snp.clinvar_allele_id);
  }
  if (!snp.mitomap_url.empty()) {
    out << ",\"mitomap_url\":" << quoted(snp.mitomap_url);
  }
  out << "}";
}

void write_quality_summary(JsonBuffer &out, std::size_t count,
                           double quality_sum, std::uint8_t minimum,
                           std::uint8_t maximum) {
  out << "{\"count\":" << count << ",\"mean_phred\":";
  if (count == 0U) {
    out << "null,\"min_phred\":null,\"max_phred\":null}";
    return;
  }
  out << (quality_sum / static_cast<double>(count)) << ","
      << "\"min_phred\":" << static_cast<unsigned int>(minimum) << ","
      << "\"max_phred\":" << static_cast<unsigned int>(maximum) << "}";
}

void write_nullable_mean(JsonBuffer &out, double sum, std::size_t count) {
  if (count == 0U) {
    out << "null";
  } else {
    out << (sum / static_cast<double>(count));
  }
}

void write_string_map(JsonBuffer &out,
                      const std::map<std::string, std::string> &values) {
  out << "{";
  std::size_t index = 0;
  for (const auto &[key, value] : values) {
    if (index++ != 0) {
      out << ",";
    }
    out << quoted(key) << ":" << quoted(value);
  }
  out << "}";
}

[[nodiscard]] std::string
alignment_fragment_id_string(const AlignmentFragmentId id) {
  return "alignment:" + std::to_string(id.value);
}

[[nodiscard]] std::string_view
observation_state_name(const ObservationState state) {
  switch (state) {
  case ObservationState::reference:
    return "REFERENCE";
  case ObservationState::alternate:
    return "ALTERNATE";
  case ObservationState::event_absent:
    return "EVENT_ABSENT";
  case ObservationState::not_callable:
    return "NOT_CALLABLE";
  case ObservationState::low_quality:
    return "LOW_QUALITY";
  case ObservationState::conflict:
    return "CONFLICT";
  }
  throw AnalysisError(AnalysisErrorCode::internal_error,
                      "unknown schema 0.6 observation state");
}

[[nodiscard]] std::string_view architecture_assignment_status_name(
    const ArchitectureAssignmentStatus status) noexcept {
  switch (status) {
  case ArchitectureAssignmentStatus::assigned:
    return "ASSIGNED";
  case ArchitectureAssignmentStatus::ambiguous:
    return "AMBIGUOUS";
  case ArchitectureAssignmentStatus::unassigned:
    return "UNASSIGNED";
  case ArchitectureAssignmentStatus::ineligible:
    return "INELIGIBLE";
  }
  return "INELIGIBLE";
}

[[nodiscard]] bool is_unified_variant_event(const EvidenceEvent &event) {
  return event.type == "SNV" || event.type == "SMALL_INSERTION" ||
         event.type == "SMALL_DELETION";
}

void add_variant_allele_observation(VariantAlleleEvidenceSummary &summary,
                                    const EvidenceObservation &observation) {
  ++summary.count;
  ++summary.strand_depths[observation.strand == '+' ? 0U : 1U];
  summary.mapping_quality_sum +=
      static_cast<double>(observation.mapping_quality);
  summary.mapping_quality_min =
      std::min(summary.mapping_quality_min, observation.mapping_quality);
  summary.mapping_quality_max =
      std::max(summary.mapping_quality_max, observation.mapping_quality);
  if (observation.base_quality) {
    ++summary.base_quality_count;
    summary.base_quality_sum += static_cast<double>(*observation.base_quality);
    summary.base_quality_min =
        std::min(summary.base_quality_min, *observation.base_quality);
    summary.base_quality_max =
        std::max(summary.base_quality_max, *observation.base_quality);
  }
  if (observation.center_proximity) {
    ++summary.read_position_count;
    summary.read_position_sum += *observation.center_proximity;
  }
}

[[nodiscard]] std::vector<UnifiedVariantEvidenceSummary>
summarize_unified_variant_evidence(const SparseEvidenceStore &store) {
  std::vector<UnifiedVariantEvidenceSummary> summaries(store.events.size());
  for (const auto &observation : store.observations) {
    if (observation.event_index >= store.events.size()) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "variant projection references an unknown event");
    }
    const auto &event = store.events[observation.event_index];
    if (!is_unified_variant_event(event)) {
      continue;
    }
    auto &summary = summaries[observation.event_index];
    switch (observation.state) {
    case ObservationState::alternate:
      add_variant_allele_observation(summary.alternate, observation);
      break;
    case ObservationState::reference:
      add_variant_allele_observation(summary.reference, observation);
      break;
    case ObservationState::event_absent:
      ++summary.event_absent;
      if (event.type == "SNV" ||
          observation.evidence_source == "cigar_alternative_small_indel") {
        add_variant_allele_observation(summary.other, observation);
      } else {
        add_variant_allele_observation(summary.reference, observation);
      }
      break;
    case ObservationState::low_quality:
      ++summary.low_quality;
      break;
    case ObservationState::conflict:
      ++summary.conflict;
      break;
    case ObservationState::not_callable:
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "unified variant projection received explicit NOT_CALLABLE");
    }
  }

  std::map<std::pair<std::string, std::size_t>, std::size_t> locus_counts;
  for (const auto &event : store.events) {
    if (is_unified_variant_event(event)) {
      ++locus_counts[{event.type, event.start}];
    }
  }
  for (std::size_t index = 0U; index < store.events.size(); ++index) {
    const auto &event = store.events[index];
    if (is_unified_variant_event(event)) {
      summaries[index].multi_allelic =
          locus_counts[{event.type, event.start}] > 1U;
    }
  }
  return summaries;
}

[[nodiscard]] std::pair<double, double>
wilson_interval(const std::size_t support, const std::size_t depth) {
  if (depth == 0U) {
    return {0.0, 0.0};
  }
  constexpr double z = 1.959963984540054;
  const double n = static_cast<double>(depth);
  const double proportion = static_cast<double>(support) / n;
  const double z2 = z * z;
  const double denominator = 1.0 + z2 / n;
  const double center = (proportion + z2 / (2.0 * n)) / denominator;
  const double margin =
      z *
      std::sqrt((proportion * (1.0 - proportion) / n) + (z2 / (4.0 * n * n))) /
      denominator;
  return {std::max(0.0, center - margin), std::min(1.0, center + margin)};
}

[[nodiscard]] std::size_t
circular_reference_homopolymer_run(const std::string &reference,
                                   const std::size_t position) {
  if (reference.empty() || position == 0U || position > reference.size()) {
    return 0U;
  }
  const char base = reference[position - 1U];
  std::size_t run = 1U;
  for (std::size_t offset = 1U; offset < reference.size(); ++offset) {
    const auto index =
        (position - 1U + reference.size() - offset) % reference.size();
    if (reference[index] != base) {
      break;
    }
    ++run;
  }
  for (std::size_t offset = 1U; run < reference.size(); ++offset) {
    const auto index = (position - 1U + offset) % reference.size();
    if (reference[index] != base) {
      break;
    }
    ++run;
  }
  return run;
}

void write_mapping_quality_summary(
    JsonBuffer &out, const VariantAlleleEvidenceSummary &summary) {
  out << "{\"count\":" << summary.count << ",\"mean\":";
  if (summary.count == 0U) {
    out << "null,\"min\":null,\"max\":null}";
    return;
  }
  out << (summary.mapping_quality_sum / static_cast<double>(summary.count))
      << ",\"min\":" << static_cast<unsigned int>(summary.mapping_quality_min)
      << ",\"max\":" << static_cast<unsigned int>(summary.mapping_quality_max)
      << "}";
}

void write_base_quality_summary(JsonBuffer &out,
                                const VariantAlleleEvidenceSummary &summary) {
  write_quality_summary(out, summary.base_quality_count,
                        summary.base_quality_sum, summary.base_quality_min,
                        summary.base_quality_max);
}

[[nodiscard]] std::vector<std::string> variant_qc_flags(
    const EvidenceEvent &event, const UnifiedVariantEvidenceSummary &summary,
    const std::size_t homopolymer_run, const bool numt_assessable) {
  std::vector<std::string> flags;
  if (!numt_assessable) {
    flags.emplace_back("NUMT_NOT_ASSESSABLE");
  }
  if (!event.absence_assessable) {
    flags.emplace_back("SUPPORT_ONLY_EVENT");
  }
  if (summary.conflict != 0U) {
    flags.emplace_back("CONFLICTING_MOLECULE");
  }
  if (summary.low_quality != 0U) {
    flags.emplace_back("LOW_QUALITY_OBSERVATIONS");
  }
  if (summary.alternate.count != 0U &&
      (summary.alternate.strand_depths[0] == 0U ||
       summary.alternate.strand_depths[1] == 0U)) {
    flags.emplace_back("SINGLE_STRAND_ALT_SUPPORT");
  }
  if (summary.multi_allelic) {
    flags.emplace_back("MULTI_ALLELIC_LOCUS");
  }
  if (homopolymer_run > 1U) {
    flags.emplace_back("HOMOPOLYMER_CONTEXT");
  }
  return flags;
}

[[nodiscard]] std::string
render_json(const std::string &input_path, const std::string &reference_path,
            const AnalysisConfig &config, std::string_view clustering_backend,
            const NumtAssessment &numt_assessment,
            std::size_t input_record_count, const std::string &reference,
            const MoleculeAssemblyResult &assembly,
            const std::vector<ReadRecord> &reads,
            const std::vector<ReadFeature> &features,
            const std::map<std::string, SvCall> &svs,
            const std::map<std::string, ComplexSvCall> &complex_events,
            const CoverageResult &coverage_result,
            const std::map<std::string, SnpAggregate> &variants,
            const SparseEvidenceStore &evidence_store,
            const CallabilityResult &callability,
            const std::vector<PhaseLink> &phase_links,
            const ArchitectureInferenceResult &architecture_inference,
            const std::map<int, ClusterHaplogroupAssignment> &haplogroups,
            const std::vector<ClusterHaplogroupAssignment>
                &architecture_haplogroups,
            const std::vector<ResourceRecord> &resources) {
  throw_if_cancelled(config);
  const auto &coverage = coverage_result.bins;
  std::map<int, std::vector<const ReadFeature *>> clusters;
  for (const auto &feature : features) {
    if (!feature.filtered_numt) {
      clusters[feature.cluster_id].push_back(&feature);
    }
  }

  std::vector<std::array<std::size_t, 5>> molecule_evidence_counts;
  std::vector<std::vector<std::string>> molecule_alternate_event_ids;
  if (config.result_schema == ResultSchema::v0_6) {
    molecule_evidence_counts.resize(assembly.molecules.size());
    molecule_alternate_event_ids.resize(assembly.molecules.size());
    for (const auto &observation : evidence_store.observations) {
      if (observation.molecule_index.value >= assembly.molecules.size() ||
          observation.event_index >= evidence_store.events.size()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "molecule projection received an unresolved observation");
      }
      auto &counts = molecule_evidence_counts[observation.molecule_index.value];
      switch (observation.state) {
      case ObservationState::alternate:
        ++counts[0];
        molecule_alternate_event_ids[observation.molecule_index.value]
            .push_back(evidence_store.events[observation.event_index].id);
        break;
      case ObservationState::reference:
        ++counts[1];
        break;
      case ObservationState::event_absent:
        ++counts[2];
        break;
      case ObservationState::low_quality:
        ++counts[3];
        break;
      case ObservationState::conflict:
        ++counts[4];
        break;
      case ObservationState::not_callable:
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "molecule projection received explicit NOT_CALLABLE");
      }
    }
  }

  const std::size_t read_snp_count =
      std::accumulate(features.begin(), features.end(), std::size_t{0},
                      [](const std::size_t count, const ReadFeature &feature) {
                        return count + feature.snps.size();
                      });
  constexpr std::size_t maximum_initial_reserve = 128U * 1024U * 1024U;
  std::size_t estimated_json_bytes = 64U * 1024U;
  const auto add_reserve_estimate = [&](const std::size_t count,
                                        const std::size_t bytes_each) {
    const std::size_t available =
        maximum_initial_reserve - estimated_json_bytes;
    estimated_json_bytes +=
        std::min(count, available / bytes_each) * bytes_each;
  };
  add_reserve_estimate(reads.size(), 512U);
  add_reserve_estimate(read_snp_count, 48U);
  add_reserve_estimate(variants.size(), 640U);
  JsonBuffer out(estimated_json_bytes);
  out << "{";
  out << "\"metadata\":{";
  out << "\"schema_version\":"
      << quoted(config.result_schema == ResultSchema::v0_6
                    ? kEvidenceGraphSchemaVersion
                    : kResultSchemaVersion)
      << ",";
  out << "\"sv_event_schema_version\":" << quoted(kSvEventSchemaVersion) << ",";
  out << "\"complex_sv_event_schema_version\":"
      << quoted(kComplexSvEventSchemaVersion) << ",";
  out << "\"clinical_annotation_schema_version\":"
      << quoted(kClinicalAnnotationSchemaVersion) << ",";
  out << "\"engine_version\":" << quoted(kEngineVersion) << ",";
  out << "\"sample\":"
      << quoted(config.sample_name.empty()
                    ? std::filesystem::path(input_path).stem().string()
                    : config.sample_name)
      << ",";
  out << "\"input_path\":" << quoted(input_path) << ",";
  out << "\"reference_path\":"
      << quoted(reference_path.empty() ? bundled_rcrs_path() : reference_path)
      << ",";
  out << "\"reference_accession\":"
      << quoted(reference_path.empty() ? "NC_012920.1" : "custom") << ",";
  out << "\"reference_length\":" << reference.size() << ",";
  out << "\"threads\":" << effective_thread_count(config.threads, reads.size())
      << ",";
  out << "\"requested_threads\":" << std::max<std::size_t>(1, config.threads)
      << ",";
  out << "\"calling_parameters\":{";
  out << "\"min_mapping_quality\":"
      << static_cast<unsigned int>(config.min_mapping_quality) << ",";
  out << "\"min_base_quality\":"
      << static_cast<unsigned int>(config.min_base_quality) << ",";
  out << "\"excluded_snp_flags\":" << config.excluded_snp_flags << ",";
  out << "\"numt_threshold\":" << config.numt_threshold << ",";
  out << "\"development_tags_enabled\":"
      << (config.allow_development_tags ? "true" : "false");
  if (config.result_schema == ResultSchema::v0_6) {
    out << ",\"max_evidence_observations\":" << config.max_evidence_observations
        << ","
        << "\"max_phase_links\":" << config.max_phase_links << ","
        << "\"max_phase_work\":" << config.max_phase_work << ","
        << "\"max_phase_molecule_references\":"
        << config.max_phase_molecule_references << ","
        << "\"max_result_bytes\":" << config.max_result_bytes << ","
        << "\"evidence_page_size\":" << config.evidence_page_size << ","
        << "\"min_architecture_molecules\":"
        << config.min_architecture_molecules << ","
        << "\"max_candidate_architectures\":"
        << config.max_candidate_architectures << ","
        << "\"architecture_max_distance\":" << config.architecture_max_distance
        << ","
        << "\"architecture_ambiguity_margin\":"
        << config.architecture_ambiguity_margin << ","
        << "\"architecture_min_overlap_fraction\":"
        << config.architecture_min_overlap_fraction << ","
        << "\"architecture_consensus_fraction\":"
        << config.architecture_consensus_fraction << ","
        << "\"architecture_optional_fraction\":"
        << config.architecture_optional_fraction << ","
        << "\"architecture_stability_replicates\":"
        << config.architecture_stability_replicates << ","
        << "\"architecture_seed\":" << config.architecture_seed << ","
        << "\"molecule_id_tag\":" << quoted(config.molecule_id_tag) << ","
        << "\"umi_tag\":" << quoted(config.umi_tag) << ","
        << "\"duplex_tag\":" << quoted(config.duplex_tag);
  }
  out << "},";
  out << "\"resources\":[";
  for (std::size_t i = 0; i < resources.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    const auto &resource = resources[i];
    out << "{\"name\":" << quoted(resource.name) << ","
        << "\"version\":" << quoted(resource.version) << ","
        << "\"path\":" << quoted(resource.path) << ","
        << "\"sha256\":" << quoted(resource.sha256) << ","
        << "\"source\":" << quoted(resource.source) << ","
        << "\"license\":" << quoted(resource.license) << ","
        << "\"retrieved\":" << quoted(resource.retrieved) << "}";
  }
  out << "],";
  out << "\"algorithm_notes\":["
      << quoted("FASTQ/SAM parser with optional htslib BAM/CRAM reader") << ","
      << quoted("Aligned reads use CIGAR-aware reference comparison for SNP "
                "extraction")
      << ","
      << quoted("NUMT specificity requires competitive "
                "nuclear-plus-mitochondrial alignment")
      << ","
      << quoted("Haplogroups use weighted PhyloTree 17.3 lineage matching")
      << ","
      << quoted("Allele QC metrics are observational and are not unvalidated "
                "hard filters")
      << ","
      << quoted("Clinical schema 1.0 preserves source assertions and "
                "deterministic conflicts")
      << ","
      << quoted("Clinical assertion payloads are emitted once per aggregate "
                "variant")
      << "," << quoted("SV event IDs use canonical mtDNA adjacency schema 1.0")
      << ","
      << quoted("Split and supplementary alignments are reconstructed as "
                "molecule event edges")
      << ","
      << quoted("Multi-junction split paths are coalesced with "
                "strand-invariant exact IDs")
      << ","
      << quoted("Candidate molecular architectures use deterministic "
                "callable-aware positive-event seed assignment; they are "
                "not cell clones")
      << "," << quoted(clustering_backend) << "]";
  out << "},";

  out << "\"filter_stats\":{";
  const std::size_t filtered = static_cast<std::size_t>(
      std::count_if(features.begin(), features.end(),
                    [](const auto &feature) { return feature.filtered_numt; }));
  out << "\"input_reads\":" << reads.size() << ",";
  out << "\"input_alignment_records\":" << input_record_count << ",";
  out << "\"input_molecules\":" << reads.size() << ",";
  out << "\"passed_reads\":" << (features.size() - filtered) << ",";
  out << "\"numt_filtered_reads\":" << filtered << ",";
  if (config.result_schema == ResultSchema::v0_6) {
    std::size_t evidence_eligible = 0U;
    for (std::size_t index = 0U; index < assembly.molecules.size(); ++index) {
      if (assembly.molecules[index].analysis_eligible &&
          index < features.size() && !features[index].filtered_numt) {
        ++evidence_eligible;
      }
    }
    const auto ambiguous = static_cast<std::size_t>(
        std::count_if(assembly.molecules.begin(), assembly.molecules.end(),
                      [](const auto &molecule) { return molecule.ambiguous; }));
    out << "\"evidence_eligible_molecules\":" << evidence_eligible << ","
        << "\"ambiguous_molecules\":" << ambiguous << ",";
  }
  out << "\"numt_threshold\":" << config.numt_threshold << ",";
  out << "\"numt_assessment\":{"
      << "\"mode\":" << quoted(numt_assessment.mode)
      << ",\"nuclear_contigs_present\":"
      << (numt_assessment.nuclear_contigs_present ? "true" : "false")
      << ",\"specificity_assessable\":"
      << (numt_assessment.specificity_assessable ? "true" : "false")
      << ",\"provenance_status\":"
      << quoted(numt_assessment.provenance_status)
      << ",\"provenance_reason\":"
      << quoted(numt_assessment.provenance_reason) << "}";
  out << "},";

  if (config.result_schema == ResultSchema::v0_6) {
    const auto observation_page_count =
        evidence_store.observations.size() / config.evidence_page_size +
        (evidence_store.observations.size() % config.evidence_page_size == 0U
             ? 0U
             : 1U);
    out << "\"evidence_encoding\":{"
        << "\"layout\":" << quoted("paged_columnar_molecule_event") << ","
        << "\"scope\":" << quoted("snv_indel_sv_complex_evidence_rc2") << ","
        << "\"observation_storage\":" << quoted("embedded_columnar_pages")
        << ","
        << "\"missing_pair_state\":" << quoted("NOT_CALLABLE") << ","
        << "\"phase_molecule_policy\":" << quoted("evidence_eligible_only")
        << ",\"phase_molecule_reference\":" << quoted("molecules[].index")
        << ","
        << "\"phase_null_model\":"
        << quoted("independent_marginals_within_jointly_callable") << ","
        << "\"observation_limit\":" << config.max_evidence_observations << ","
        << "\"observation_count\":" << evidence_store.observations.size()
        << ",\"observation_page_size\":" << config.evidence_page_size << ","
        << "\"observation_page_count\":" << observation_page_count << ","
        << "\"phase_link_limit\":" << config.max_phase_links << ","
        << "\"phase_work_limit\":" << config.max_phase_work << ","
        << "\"phase_molecule_reference_limit\":"
        << config.max_phase_molecule_references << ","
        << "\"result_byte_limit\":" << config.max_result_bytes << "},";

    out << "\"alignments\":[";
    for (std::size_t index = 0; index < assembly.fragments.size(); ++index) {
      if ((index & 1023U) == 0U) {
        throw_if_cancelled(config);
      }
      if (index != 0U) {
        out << ",";
      }
      const auto &fragment = assembly.fragments[index];
      const ReadRecord *source_record = nullptr;
      if (!assembly.source_alignment_records.empty() &&
          fragment.source_record_index <
              assembly.source_alignment_records.size()) {
        source_record =
            &assembly.source_alignment_records[fragment.source_record_index];
      } else if (fragment.molecule_index.value <
                 assembly.representatives.size()) {
        source_record =
            &assembly.representatives[fragment.molecule_index.value];
      }
      if (source_record == nullptr) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "schema 0.6 alignment source record is unresolved");
      }
      out << "{\"id\":" << quoted(alignment_fragment_id_string(fragment.id))
          << ","
          << "\"source_record_index\":" << fragment.source_record_index << ","
          << "\"molecule_id\":" << quoted(fragment.molecule_id) << ","
          << "\"molecule_index\":" << fragment.molecule_index.value << ","
          << "\"role\":" << quoted(fragment.role) << ","
          << "\"selected_representative\":"
          << (fragment.selected_representative ? "true" : "false") << ","
          << "\"flags\":" << fragment.flags << ","
          << "\"strand\":" << quoted((fragment.flags & 0x10U) == 0U ? "+" : "-")
          << ","
          << "\"mapping_quality\":"
          << static_cast<unsigned int>(fragment.mapping_quality)
          << ",\"reference_name\":" << quoted(fragment.reference_name) << ","
          << "\"reference_start\":" << fragment.reference_start << ","
          << "\"cigar\":" << quoted(fragment.cigar) << ","
          << "\"query_length\":" << source_record->sequence.size() << ","
          << "\"base_qualities_available\":"
          << (!source_record->qualities.empty() ? "true" : "false") << ","
          << "\"aux_tags\":";
      write_string_map(out, source_record->aux_tags);
      out << "}";
    }
    out << "],";

    out << "\"molecules\":[";
    for (std::size_t index = 0; index < assembly.molecules.size(); ++index) {
      if ((index & 1023U) == 0U) {
        throw_if_cancelled(config);
      }
      if (index != 0U) {
        out << ",";
      }
      const auto &molecule = assembly.molecules[index];
      const auto &molecule_callability = callability.molecules[index];
      const auto &feature = features[index];
      out << "{\"id\":" << quoted(molecule.id) << ","
          << "\"index\":" << molecule.index.value << ","
          << "\"identity_policy\":" << quoted(molecule.identity_policy) << ","
          << "\"assembly_status\":" << quoted(molecule.assembly_status) << ","
          << "\"primary_candidate_count\":" << molecule.primary_candidate_count
          << ","
          << "\"ambiguous\":" << (molecule.ambiguous ? "true" : "false") << ","
          << "\"analysis_eligible\":"
          << (molecule.analysis_eligible ? "true" : "false") << ","
          << "\"evidence_eligible\":"
          << (molecule.analysis_eligible && index < features.size() &&
                      !features[index].filtered_numt
                  ? "true"
                  : "false")
          << ","
          << "\"callability_status\":" << quoted(molecule_callability.status)
          << ","
          << "\"callable_bases\":" << molecule_callability.callable_bases << ","
          << "\"callable_fraction\":"
          << (reference.empty()
                  ? 0.0
                  : static_cast<double>(molecule_callability.callable_bases) /
                        static_cast<double>(reference.size()))
          << ","
          << "\"query_length\":" << feature.length << ","
          << "\"mean_base_quality\":" << feature.mean_quality << ","
          << "\"mapping_quality\":"
          << static_cast<unsigned int>(feature.mapping_quality) << ","
          << "\"numt_score\":" << feature.numt_score << ","
          << "\"numt_evidence\":";
      write_string_array(out, feature.numt_evidence);
      if (index >= architecture_inference.assignments.size()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "molecule architecture assignment projection is incomplete");
      }
      const auto &architecture_assignment =
          architecture_inference.assignments[index];
      out << ",\"cluster_id\":" << feature.cluster_id << ","
          << "\"architecture_assignment\":{\"status\":"
          << quoted(architecture_assignment_status_name(
                 architecture_assignment.status))
          << ",\"architecture_id\":";
      if (architecture_assignment.architecture_index.has_value()) {
        if (*architecture_assignment.architecture_index >=
            architecture_inference.architectures.size()) {
          throw AnalysisError(AnalysisErrorCode::internal_error,
                              "molecule architecture assignment is unresolved");
        }
        out << quoted(
            architecture_inference
                .architectures[*architecture_assignment.architecture_index]
                .id);
      } else {
        out << "null";
      }
      out << ",\"distance\":" << architecture_assignment.distance
          << ",\"score\":" << architecture_assignment.score
          << ",\"overlap_fraction\":"
          << architecture_assignment.overlap_fraction
          << ",\"confidence\":" << architecture_assignment.confidence
          << ",\"candidate_architecture_ids\":[";
      for (std::size_t candidate_index = 0U;
           candidate_index < architecture_assignment.candidate_indices.size();
           ++candidate_index) {
        if (candidate_index != 0U) {
          out << ",";
        }
        const auto architecture_index =
            architecture_assignment.candidate_indices[candidate_index];
        if (architecture_index >= architecture_inference.architectures.size()) {
          throw AnalysisError(AnalysisErrorCode::internal_error,
                              "ambiguous architecture candidate is unresolved");
        }
        out << quoted(
            architecture_inference.architectures[architecture_index].id);
      }
      out << "]},"
          << "\"alternate_event_ids\":";
      write_string_array(out, molecule_alternate_event_ids[index]);
      const auto &evidence_counts = molecule_evidence_counts[index];
      out << ",\"evidence_state_counts\":{\"alternate\":" << evidence_counts[0]
          << ",\"reference\":" << evidence_counts[1]
          << ",\"event_absent\":" << evidence_counts[2]
          << ",\"low_quality\":" << evidence_counts[3]
          << ",\"conflict\":" << evidence_counts[4] << "},"
          << "\"representative_alignment_id\":"
          << quoted(alignment_fragment_id_string(
                 molecule.representative_fragment_id))
          << ","
          << "\"source_qnames\":";
      write_string_array(out, molecule.source_qnames);
      out << ",\"protocol_metadata\":";
      write_string_map(out, molecule.protocol_metadata);
      out << ",\"protocol_flags\":";
      write_string_array(out, molecule.protocol_flags);
      out << ",\"exclusion_reasons\":";
      write_string_array(out, molecule.exclusion_reasons);
      out << ","
          << "\"alignment_ids\":[";
      for (std::size_t fragment_index = 0;
           fragment_index < molecule.fragment_ids.size(); ++fragment_index) {
        if (fragment_index != 0U) {
          out << ",";
        }
        out << quoted(alignment_fragment_id_string(
            molecule.fragment_ids[fragment_index]));
      }
      out << "],\"warnings\":";
      write_string_array(out, molecule.warnings);
      out << "}";
    }
    out << "],";

    out << "\"callability\":[";
    for (std::size_t index = 0U; index < callability.molecules.size();
         ++index) {
      if ((index & 1023U) == 0U) {
        throw_if_cancelled(config);
      }
      if (index != 0U) {
        out << ",";
      }
      const auto &summary = callability.molecules[index];
      out << "{\"molecule_id\":" << quoted(assembly.molecules[index].id) << ","
          << "\"status\":" << quoted(summary.status) << ","
          << "\"known\":" << (summary.known ? "true" : "false") << ","
          << "\"basis\":" << quoted("passing_aligned_reference_bases") << ","
          << "\"callable_bases\":" << summary.callable_bases << ","
          << "\"callable_fraction\":"
          << (reference.empty() ? 0.0
                                : static_cast<double>(summary.callable_bases) /
                                      static_cast<double>(reference.size()))
          << ",\"ranges\":[";
      for (std::size_t range_index = 0U; range_index < summary.ranges.size();
           ++range_index) {
        if (range_index != 0U) {
          out << ",";
        }
        out << "{\"start\":" << summary.ranges[range_index].first
            << ",\"end\":" << summary.ranges[range_index].second << "}";
      }
      out << "],\"alignments\":[";
      for (std::size_t alignment_index = 0U;
           alignment_index < summary.alignments.size(); ++alignment_index) {
        if (alignment_index != 0U) {
          out << ",";
        }
        const auto &alignment = summary.alignments[alignment_index];
        out << "{\"alignment_id\":"
            << quoted(alignment_fragment_id_string(
                   alignment.alignment_fragment_id))
            << ","
            << "\"eligible\":" << (alignment.eligible ? "true" : "false") << ","
            << "\"status\":" << quoted(alignment.status) << ","
            << "\"callable_bases\":" << alignment.callable_bases << ","
            << "\"inserted_query_bases\":" << alignment.inserted_query_bases
            << ","
            << "\"soft_clipped_query_bases\":"
            << alignment.soft_clipped_query_bases << ","
            << "\"reference_exclusion_counts\":{";
        std::size_t exclusion_index = 0U;
        for (const auto &[reason, count] :
             alignment.reference_exclusion_counts) {
          if (exclusion_index++ != 0U) {
            out << ",";
          }
          out << quoted(reason) << ":" << count;
        }
        out << "},\"disrupted_adjacency_anchors\":[";
        for (std::size_t anchor_index = 0U;
             anchor_index < alignment.disrupted_adjacency_anchors.size();
             ++anchor_index) {
          if (anchor_index != 0U) {
            out << ",";
          }
          out << alignment.disrupted_adjacency_anchors[anchor_index];
        }
        out << "],\"ranges\":[";
        for (std::size_t range_index = 0U;
             range_index < alignment.ranges.size(); ++range_index) {
          if (range_index != 0U) {
            out << ",";
          }
          out << "{\"start\":" << alignment.ranges[range_index].first
              << ",\"end\":" << alignment.ranges[range_index].second << "}";
        }
        out << "]}";
      }
      out << "]}";
    }
    out << "],";

    struct SerializedEventCounts {
      std::size_t alternate = 0U;
      std::size_t reference = 0U;
      std::size_t event_absent = 0U;
      std::size_t low_quality = 0U;
      std::size_t conflict = 0U;
    };
    std::vector<SerializedEventCounts> event_counts(
        evidence_store.events.size());
    for (const auto &observation : evidence_store.observations) {
      auto &counts = event_counts[observation.event_index];
      switch (observation.state) {
      case ObservationState::alternate:
        ++counts.alternate;
        break;
      case ObservationState::reference:
        ++counts.reference;
        break;
      case ObservationState::event_absent:
        ++counts.event_absent;
        break;
      case ObservationState::low_quality:
        ++counts.low_quality;
        break;
      case ObservationState::conflict:
        ++counts.conflict;
        break;
      case ObservationState::not_callable:
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 serializer received explicit NOT_CALLABLE");
      }
    }

    out << "\"events\":[";
    for (std::size_t index = 0; index < evidence_store.events.size(); ++index) {
      if ((index & 1023U) == 0U) {
        throw_if_cancelled(config);
      }
      if (index != 0U) {
        out << ",";
      }
      const auto &event = evidence_store.events[index];
      const auto &counts = event_counts[index];
      out << "{\"id\":" << quoted(event.id) << ","
          << "\"index\":" << index << ","
          << "\"type\":" << quoted(event.type) << ",\"start\":";
      if (event.start == 0U) {
        out << "null";
      } else {
        out << event.start;
      }
      out << ",\"end\":";
      if (event.end == 0U) {
        out << "null";
      } else {
        out << event.end;
      }
      out << ",\"length\":" << event.length << ",\"ref\":";
      if (event.reference.empty()) {
        out << "null";
      } else {
        out << quoted(event.reference);
      }
      out << ",\"alt\":";
      if (event.alternate.empty()) {
        out << "null";
      } else {
        out << quoted(event.alternate);
      }
      out << ",\"normalization\":" << quoted(event.normalization) << ","
          << "\"source_projection\":" << quoted(event.source_projection) << ","
          << "\"negative_evidence_rule\":"
          << quoted(event.negative_evidence_rule) << ","
          << "\"assessability\":"
          << quoted(event.absence_assessable ? "REFERENCE_AND_ALTERNATE"
                                             : "ALTERNATE_SUPPORT_ONLY")
          << ",\"component_event_ids\":";
      write_string_array(out, event.component_event_ids);
      out << ",\"supporting_molecule_ids\":";
      write_string_array(out, event.supporting_molecules);
      out << ",\"evidence_counts\":{"
          << "\"alternate\":" << counts.alternate << ","
          << "\"reference\":" << counts.reference << ","
          << "\"event_absent\":" << counts.event_absent << ","
          << "\"callable\":"
          << (counts.alternate + counts.reference + counts.event_absent)
          << ",\"low_quality\":" << counts.low_quality << ","
          << "\"conflict\":" << counts.conflict << "}}";
    }
    out << "],";

    const auto write_observation_column = [&](const std::size_t begin,
                                              const std::size_t end,
                                              const auto &write_value) {
      out << "[";
      for (std::size_t index = begin; index < end; ++index) {
        if (index != begin) {
          out << ",";
        }
        write_value(evidence_store.observations[index]);
      }
      out << "]";
    };
    out << "\"observation_pages\":[";
    for (std::size_t page_index = 0U; page_index < observation_page_count;
         ++page_index) {
      throw_if_cancelled(config);
      if (page_index != 0U) {
        out << ",";
      }
      const auto begin = page_index * config.evidence_page_size;
      const auto end = std::min(evidence_store.observations.size(),
                                begin + config.evidence_page_size);
      out << "{\"index\":" << page_index << ",\"offset\":" << begin
          << ",\"count\":" << (end - begin) << ",\"columns\":{";
      out << "\"molecule_id\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << quoted(assembly.molecules[observation.molecule_index.value].id);
      });
      out << ",\"event_id\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << quoted(evidence_store.events[observation.event_index].id);
      });
      out << ",\"alignment_id\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << quoted(
            alignment_fragment_id_string(observation.alignment_fragment_id));
      });
      out << ",\"state\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << quoted(observation_state_name(observation.state));
      });
      out << ",\"observed_allele\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        if (observation.observed_allele.empty()) {
          out << "null";
        } else {
          out << quoted(observation.observed_allele);
        }
      });
      out << ",\"base_quality\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        if (observation.base_quality) {
          out << static_cast<unsigned int>(*observation.base_quality);
        } else {
          out << "null";
        }
      });
      out << ",\"mapping_quality\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << static_cast<unsigned int>(observation.mapping_quality);
      });
      out << ",\"strand\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << quoted(std::string(1, observation.strand));
      });
      out << ",\"evidence_source\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        out << quoted(observation.evidence_source);
      });
      out << ",\"read_position\":";
      write_observation_column(begin, end, [&](const auto &observation) {
        if (observation.center_proximity) {
          out << *observation.center_proximity;
        } else {
          out << "null";
        }
      });
      out << "}}";
    }
    out << "],\"phase_links\":[";
    for (std::size_t index = 0U; index < phase_links.size(); ++index) {
      if ((index & 1023U) == 0U) {
        throw_if_cancelled(config);
      }
      if (index != 0U) {
        out << ",";
      }
      const auto &link = phase_links[index];
      const auto &event_a = evidence_store.events[link.event_a_index];
      const auto &event_b = evidence_store.events[link.event_b_index];
      if (link.supporting_molecule_indices.size() != link.both_alternate ||
          link.uncertain_molecule_indices.size() != link.jointly_uncertain ||
          std::ranges::any_of(link.supporting_molecule_indices,
                              [&](const auto molecule_index) {
                                return molecule_index >=
                                       assembly.molecules.size();
                              }) ||
          std::ranges::any_of(
              link.uncertain_molecule_indices, [&](const auto molecule_index) {
                return molecule_index >= assembly.molecules.size();
              })) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "phase molecule-index traceability is inconsistent");
      }
      out << "{\"id\":" << quoted("phase:" + event_a.id + "|" + event_b.id)
          << ","
          << "\"event_a_id\":" << quoted(event_a.id) << ","
          << "\"event_b_id\":" << quoted(event_b.id) << ","
          << "\"assessability\":"
          << quoted(link.complete_callability ? "COMPLETE_FOR_BOTH_EVENTS"
                                              : "SUPPORT_CONDITIONED")
          << ",\"jointly_callable\":" << link.jointly_callable << ","
          << "\"jointly_uncertain\":" << link.jointly_uncertain << ","
          << "\"both_alternate\":" << link.both_alternate << ","
          << "\"a_alternate_b_absent\":" << link.a_alternate_b_absent << ","
          << "\"a_absent_b_alternate\":" << link.a_absent_b_alternate << ","
          << "\"neither_alternate\":" << link.neither_alternate << ","
          << "\"co_alternate_fraction\":" << link.co_alternate_fraction << ","
          << "\"co_alternate_ci95_low\":" << link.ci95_low << ","
          << "\"co_alternate_ci95_high\":" << link.ci95_high << ","
          << "\"expected_co_alternate_fraction\":"
          << link.expected_co_alternate_fraction << ","
          << "\"linkage_delta\":" << link.linkage_delta << ","
          << "\"supporting_molecule_indices\":";
      write_integer_array(out, link.supporting_molecule_indices);
      out << ",\"uncertain_molecule_indices\":";
      write_integer_array(out, link.uncertain_molecule_indices);
      out << ",\"qc_flags\":[";
      bool wrote_phase_flag = false;
      if (!link.complete_callability) {
        out << quoted("SUPPORT_CONDITIONED");
        wrote_phase_flag = true;
      }
      if (link.jointly_uncertain != 0U) {
        if (wrote_phase_flag) {
          out << ",";
        }
        out << quoted("UNCERTAIN_COOCCURRENCE_EXCLUDED");
      }
      out << "]}";
    }
    out << "],\"architecture_inference\":{"
        << "\"schema_version\":" << quoted("0.8-rc4-dev1") << ","
        << "\"status\":" << quoted(architecture_inference.status) << ","
        << "\"interpretation\":"
        << quoted("candidate_molecular_architectures_not_cell_clones") << ","
        << "\"method\":"
        << quoted("deterministic_callable_aware_seed_assignment_v1") << ","
        << "\"eligible_molecules\":"
        << architecture_inference.eligible_molecules << ","
        << "\"assigned_molecules\":"
        << architecture_inference.assigned_molecules << ","
        << "\"ambiguous_molecules\":"
        << architecture_inference.ambiguous_molecules << ","
        << "\"unassigned_molecules\":"
        << architecture_inference.unassigned_molecules << ","
        << "\"candidate_count\":" << architecture_inference.architectures.size()
        << ","
        << "\"qc_flags\":";
    write_string_array(out, architecture_inference.qc_flags);
    out << "},\"architectures\":[";
    for (std::size_t architecture_index = 0U;
         architecture_index < architecture_inference.architectures.size();
         ++architecture_index) {
      if (architecture_index != 0U) {
        out << ",";
      }
      const auto &architecture =
          architecture_inference.architectures[architecture_index];
      const ClusterHaplogroupAssignment *architecture_haplogroup =
          architecture_index < architecture_haplogroups.size()
              ? &architecture_haplogroups[architecture_index]
              : nullptr;
      const auto write_event_ids =
          [&](const std::vector<std::size_t> &indices) {
            out << "[";
            for (std::size_t index = 0U; index < indices.size(); ++index) {
              if (index != 0U) {
                out << ",";
              }
              if (indices[index] >= evidence_store.events.size()) {
                throw AnalysisError(
                    AnalysisErrorCode::internal_error,
                    "architecture event projection is unresolved");
              }
              out << quoted(evidence_store.events[indices[index]].id);
            }
            out << "]";
          };
      const auto write_molecule_ids =
          [&](const std::vector<std::size_t> &indices) {
            out << "[";
            for (std::size_t index = 0U; index < indices.size(); ++index) {
              if (index != 0U) {
                out << ",";
              }
              if (indices[index] >= assembly.molecules.size()) {
                throw AnalysisError(
                    AnalysisErrorCode::internal_error,
                    "architecture molecule projection is unresolved");
              }
              out << quoted(assembly.molecules[indices[index]].id);
            }
            out << "]";
          };

      out << "{\"architecture_id\":" << quoted(architecture.id) << ","
          << "\"status\":" << quoted("CANDIDATE") << ","
          << "\"defining_event_signature\":";
      write_event_ids(architecture.defining_event_indices);
      out << ",\"seed_event_signature\":";
      write_event_ids(architecture.seed_signature);
      out << ",\"optional_event_ids\":";
      write_event_ids(architecture.optional_event_indices);
      out << ",\"molecule_count\":"
          << architecture.member_molecule_indices.size() << ","
          << "\"estimated_fraction\":" << architecture.estimated_fraction
          << ",\"confidence_interval\":{\"method\":"
          << quoted("wilson_score_95")
          << ",\"level\":0.95,\"low\":" << architecture.ci95_low
          << ",\"high\":" << architecture.ci95_high << "},"
          << "\"median_callable_fraction\":"
          << architecture.median_callable_fraction << ","
          << "\"assignment_confidence\":{\"mean\":"
          << architecture.mean_assignment_confidence
          << ",\"minimum\":" << architecture.min_assignment_confidence << "},"
          << "\"cluster_stability\":{\"status\":"
          << quoted(config.architecture_stability_replicates == 0U
                        ? "NOT_ESTIMATED"
                        : "RESAMPLING_ESTIMATE")
          << ",\"value\":" << architecture.signature_stability << ",\"method\":"
          << quoted("deterministic_80pct_full_dataset_reinference_v1")
          << ",\"recovery_rate\":" << architecture.recovery_rate
          << ",\"abundance_standard_deviation\":"
          << architecture.abundance_standard_deviation
          << ",\"replicates\":" << config.architecture_stability_replicates
          << ",\"seed\":" << config.architecture_seed << "},"
          << "\"representative_molecule_ids\":";
      write_molecule_ids(architecture.representative_molecule_indices);
      out << ",\"member_molecule_indices\":";
      write_integer_array(out, architecture.member_molecule_indices);
      out << ",\"ambiguous_molecule_ids\":";
      write_molecule_ids(architecture.ambiguous_molecule_indices);
      out << ",\"unassigned_molecule_count\":"
          << architecture_inference.unassigned_molecules << ","
          << "\"haplogroup_evidence\":{\"status\":"
          << quoted(architecture_haplogroup == nullptr
                        ? "NOT_AVAILABLE"
                        : (!architecture_haplogroup->callable_ranges_known
                               ? "CALLABLE_RANGE_UNKNOWN"
                               : (architecture_haplogroup->best == "unassigned"
                                      ? "UNASSIGNED"
                                      : "ASSIGNED")))
          << ",\"resource\":" << quoted("phylotree-rcrs@17.3")
          << ",\"best\":"
          << quoted(architecture_haplogroup == nullptr
                        ? "unassigned"
                        : (!architecture_haplogroup->callable_ranges_known
                               ? "unassigned"
                               : architecture_haplogroup->best))
          << ",\"quality\":"
          << (architecture_haplogroup == nullptr
                  ? 0.0
                  : (architecture_haplogroup->callable_ranges_known
                         ? architecture_haplogroup->quality
                         : 0.0))
          << ",\"callable_ranges_known\":"
          << (architecture_haplogroup != nullptr &&
                      architecture_haplogroup->callable_ranges_known
                  ? "true"
                  : "false")
          << ",\"contamination_warning\":"
          << (architecture_haplogroup != nullptr &&
                      architecture_haplogroup->contamination_warning
                  ? "true"
                  : "false")
          << ",\"observed_markers\":";
      if (architecture_haplogroup == nullptr) {
        out << "[]";
      } else {
        write_string_array(out, architecture_haplogroup->observed_markers);
      }
      out << ",\"callable_ranges\":[";
      if (architecture_haplogroup != nullptr) {
        for (std::size_t range_index = 0U;
             range_index < architecture_haplogroup->callable_ranges.size();
             ++range_index) {
          if (range_index != 0U) {
            out << ",";
          }
          const auto &[start, end] =
              architecture_haplogroup->callable_ranges[range_index];
          out << "{\"start\":" << start << ",\"end\":" << end << "}";
        }
      }
      out << "]},\"structural_path_evidence\":{\"defining_event_ids\":[";
      bool wrote_structural = false;
      for (const auto event_index : architecture.defining_event_indices) {
        const auto &event = evidence_store.events[event_index];
        if (!event.type.starts_with("SV_") && event.type != "COMPLEX_SV_PATH") {
          continue;
        }
        if (wrote_structural) {
          out << ",";
        }
        out << quoted(event.id);
        wrote_structural = true;
      }
      out << "],\"optional_event_ids\":[";
      wrote_structural = false;
      for (const auto event_index : architecture.optional_event_indices) {
        const auto &event = evidence_store.events[event_index];
        if (!event.type.starts_with("SV_") && event.type != "COMPLEX_SV_PATH") {
          continue;
        }
        if (wrote_structural) {
          out << ",";
        }
        out << quoted(event.id);
        wrote_structural = true;
      }
      out << "]},\"event_profile\":[";
      for (std::size_t profile_index = 0U;
           profile_index < architecture.event_profile.size(); ++profile_index) {
        if (profile_index != 0U) {
          out << ",";
        }
        const auto &profile = architecture.event_profile[profile_index];
        out << "{\"event_id\":"
            << quoted(evidence_store.events[profile.event_index].id) << ","
            << "\"alternate\":" << profile.alternate << ","
            << "\"absent\":" << profile.absent << ","
            << "\"uncertain\":" << profile.uncertain << ","
            << "\"not_callable\":" << profile.not_callable << ","
            << "\"alternate_fraction\":" << profile.alternate_fraction
            << ",\"consensus_state\":"
            << quoted(profile.expected_alternate
                          ? "ALTERNATE"
                          : (profile.expected_absent ? "ABSENT" : "OPTIONAL"))
            << "}";
      }
      out << "],\"method\":{\"name\":"
          << quoted("deterministic_callable_aware_seed_assignment_v1") << ","
          << "\"positive_evidence_seeded\":true,"
          << "\"not_callable_is_reference\":false,"
          << "\"minimum_molecules\":" << config.min_architecture_molecules
          << ","
          << "\"maximum_distance\":" << config.architecture_max_distance << ","
          << "\"minimum_overlap_fraction\":"
          << config.architecture_min_overlap_fraction << ","
          << "\"ambiguity_margin\":" << config.architecture_ambiguity_margin
          << "}}";
    }
    out << "],";
  }

  out << "\"genes\":[";
  const auto genes = default_genes();
  for (std::size_t i = 0; i < genes.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    out << "{"
        << "\"name\":" << quoted(genes[i].name) << ","
        << "\"start\":" << genes[i].start << ","
        << "\"end\":" << genes[i].end << ","
        << "\"strand\":" << quoted(genes[i].strand) << ","
        << "\"biotype\":" << quoted(genes[i].biotype) << "}";
  }
  out << "],";

  out << "\"coverage\":[";
  for (std::size_t i = 0; i < coverage.size(); ++i) {
    if (i != 0) {
      out << ",";
    }
    out << "{"
        << "\"start\":" << coverage[i].start << ","
        << "\"end\":" << coverage[i].end << ","
        << "\"depth\":" << coverage[i].depth << "}";
  }
  out << "],";

  const auto bins_gt20 = static_cast<std::size_t>(
      std::count_if(coverage.begin(), coverage.end(),
                    [](const auto &bin) { return bin.depth > 20; }));
  const double pct_bins_gt20 = coverage.empty()
                                   ? 0.0
                                   : (static_cast<double>(bins_gt20) * 100.0) /
                                         static_cast<double>(coverage.size());
  std::map<std::uint8_t, std::size_t> mapq_histogram;
  for (std::size_t i = 0; i < reads.size(); ++i) {
    if (i < features.size() && !features[i].filtered_numt &&
        (reads[i].flags & 0x4U) == 0U && !reads[i].reference_name.empty() &&
        reads[i].reference_name != "*") {
      if (!looks_like_nuclear_contig(reads[i].reference_name)) {
        ++mapq_histogram[reads[i].mapping_quality];
      }
    }
  }
  out << "\"coverage_metrics\":{";
  out << "\"mean_depth\":" << coverage_result.mean_depth << ",";
  out << "\"pct_sites_gt20x\":" << coverage_result.pct_sites_gt20x << ",";
  out << "\"pct_bins_gt20x\":" << pct_bins_gt20 << ",";
  out << "\"max_depth\":" << coverage_result.max_depth << ",";
  out << "\"mapping_quality_histogram\":[";
  std::size_t mapq_index = 0;
  for (const auto &[mapq, count] : mapq_histogram) {
    if (mapq_index++ != 0) {
      out << ",";
    }
    out << "{"
        << "\"mapq\":" << static_cast<int>(mapq) << ","
        << "\"count\":" << count << "}";
  }
  out << "]";
  out << "},";

  out << "\"svs\":[";
  std::size_t sv_index = 0;
  for (const auto &[_, sv] : svs) {
    if (sv_index++ != 0) {
      out << ",";
    }
    out << "{"
        << "\"id\":" << quoted(sv.id) << ",";
    if (config.result_schema == ResultSchema::v0_6) {
      out << "\"event_id\":" << quoted("sv:" + sv.id) << ",";
    }
    out << "\"type\":" << quoted(sv.type) << ","
        << "\"start\":" << sv.start << ","
        << "\"end\":" << sv.end << ","
        << "\"length\":" << sv.length << ","
        << "\"known_event\":" << (sv.known_event ? "true" : "false") << ","
        << "\"evidence_source\":"
        << quoted(sv.evidence_sources.size() == 1U ? sv.evidence_sources.front()
                                                   : "combined")
        << ",\"evidence_sources\":";
    write_string_array(out, sv.evidence_sources);
    out << ",\"segment_count\":" << sv.segment_count << ",";
    if (!sv.orientations.empty()) {
      out << "\"orientation\":"
          << quoted(sv.orientations.size() == 1U ? sv.orientations.front()
                                                 : "mixed")
          << ",\"orientations\":";
      write_string_array(out, sv.orientations);
      out << ",";
    }
    out << "\"supporting_reads\":[";
    for (std::size_t i = 0; i < sv.supporting_reads.size(); ++i) {
      if (i != 0) {
        out << ",";
      }
      out << quoted(sv.supporting_reads[i]);
    }
    out << "]}";
  }
  out << "],";

  out << "\"complex_events\":[";
  std::size_t complex_event_index = 0;
  for (const auto &[_, event] : complex_events) {
    if (complex_event_index++ != 0U) {
      out << ",";
    }
    out << "{"
        << "\"id\":" << quoted(event.id) << ",";
    if (config.result_schema == ResultSchema::v0_6) {
      out << "\"event_id\":" << quoted(event.id) << ",";
    }
    out << "\"junction_count\":" << event.junction_ids.size() << ","
        << "\"segment_count\":" << event.segment_count << ","
        << "\"canonicalization\":" << quoted("strand_invariant_path") << ","
        << "\"junction_ids\":";
    write_string_array(out, event.junction_ids);
    out << ",\"junction_orientations\":";
    write_string_array(out, event.junction_orientations);
    out << ",\"supporting_reads\":";
    write_string_array(out, event.supporting_reads);
    out << "}";
  }
  out << "],";

  out << "\"variants\":[";
  std::size_t variant_index = 0;
  if (config.result_schema == ResultSchema::v0_6) {
    const auto unified_summaries =
        summarize_unified_variant_evidence(evidence_store);
    for (std::size_t event_index = 0U;
         event_index < evidence_store.events.size(); ++event_index) {
      const auto &event = evidence_store.events[event_index];
      if (!is_unified_variant_event(event)) {
        continue;
      }
      if (variant_index++ != 0U) {
        out << ",";
      }
      const auto &summary = unified_summaries[event_index];
      const auto callable_depth = summary.alternate.count +
                                  summary.reference.count + summary.other.count;
      const double heteroplasmy =
          callable_depth == 0U ? 0.0
                               : static_cast<double>(summary.alternate.count) /
                                     static_cast<double>(callable_depth);
      const auto [ci95_low, ci95_high] =
          wilson_interval(summary.alternate.count, callable_depth);
      const auto homopolymer_run =
          circular_reference_homopolymer_run(reference, event.start);
      const auto qc_flags = variant_qc_flags(
          event, summary, homopolymer_run,
          numt_assessment.specificity_assessable);
      const bool circular_origin_event =
          event.type == "SMALL_DELETION" &&
          (event.start == 1U || event.end < event.start);
      const bool vcf_representable = !circular_origin_event;
      const auto vcf_position =
          event.type == "SMALL_DELETION" && event.start > 1U ? event.start - 1U
                                                             : event.start;
      const SnpCall *snp = nullptr;
      if (event.type == "SNV" && event.reference.size() == 1U &&
          event.alternate.size() == 1U) {
        const auto found = variants.find(snp_key(
            event.start, event.reference.front(), event.alternate.front()));
        if (found == variants.end()) {
          throw AnalysisError(
              AnalysisErrorCode::internal_error,
              "schema 0.6 SNV event has no aggregate annotation projection");
        }
        snp = &found->second.call;
      }

      out << "{\"event_id\":" << quoted(event.id) << ","
          << "\"type\":" << quoted(event.type) << ","
          << "\"position\":" << event.start << ","
          << "\"start\":" << event.start << ","
          << "\"end\":" << event.end << ","
          << "\"length\":" << event.length << ","
          << "\"ref\":" << quoted(event.reference) << ","
          << "\"alt\":" << quoted(event.alternate) << ","
          << "\"normalization\":" << quoted(event.normalization) << ","
          << "\"negative_evidence_rule\":"
          << quoted(event.negative_evidence_rule) << ","
          << "\"assessability\":"
          << quoted(event.absence_assessable ? "REFERENCE_AND_ALTERNATE"
                                             : "ALTERNATE_SUPPORT_ONLY")
          << ",\"vcf_position\":" << vcf_position << ","
          << "\"vcf_representable\":" << (vcf_representable ? "true" : "false")
          << ","
          << "\"alt_depth\":" << summary.alternate.count << ","
          << "\"ref_depth\":" << summary.reference.count << ","
          << "\"other_depth\":" << summary.other.count << ","
          << "\"event_absent_depth\":" << summary.event_absent << ","
          << "\"low_quality_depth\":" << summary.low_quality << ","
          << "\"conflict_depth\":" << summary.conflict << ","
          << "\"callable_depth\":" << callable_depth << ","
          << "\"heteroplasmy\":" << heteroplasmy << ","
          << "\"ci95_low\":" << ci95_low << ","
          << "\"ci95_high\":" << ci95_high << ","
          << "\"filter_status\":" << quoted("NOT_CALIBRATED") << ","
          << "\"qc_flags\":";
      write_string_array(out, qc_flags);
      out << ",\"numt_assessability\":"
          << quoted(numt_assessment.specificity_assessable ? "ASSESSABLE"
                                                           : "NOT_ASSESSABLE")
          << ",\"multi_allelic\":" << (summary.multi_allelic ? "true" : "false")
          << ","
          << "\"homopolymer_context\":{\"reference_base\":";
      if (event.start == 0U || event.start > reference.size()) {
        out << "null";
      } else {
        out << quoted(std::string(1U, reference[event.start - 1U]));
      }
      out << ",\"run_length\":" << homopolymer_run << "},"
          << "\"molecule_support\":{\"alternate\":" << summary.alternate.count
          << ",\"reference\":" << summary.reference.count
          << ",\"other\":" << summary.other.count
          << ",\"callable\":" << callable_depth
          << "},\"strand_support\":{\"alt_forward\":"
          << summary.alternate.strand_depths[0]
          << ",\"alt_reverse\":" << summary.alternate.strand_depths[1]
          << ",\"ref_forward\":" << summary.reference.strand_depths[0]
          << ",\"ref_reverse\":" << summary.reference.strand_depths[1]
          << ",\"other_forward\":" << summary.other.strand_depths[0]
          << ",\"other_reverse\":" << summary.other.strand_depths[1]
          << "},\"strand_bias_delta\":";
      if (summary.alternate.count == 0U || summary.reference.count == 0U) {
        out << "null";
      } else {
        const double alternate_forward_fraction =
            static_cast<double>(summary.alternate.strand_depths[0]) /
            static_cast<double>(summary.alternate.count);
        const double reference_forward_fraction =
            static_cast<double>(summary.reference.strand_depths[0]) /
            static_cast<double>(summary.reference.count);
        out << std::abs(alternate_forward_fraction -
                        reference_forward_fraction);
      }
      out << ",\"allele_quality\":{\"alternate\":";
      write_base_quality_summary(out, summary.alternate);
      out << ",\"reference\":";
      write_base_quality_summary(out, summary.reference);
      out << ",\"other\":";
      write_base_quality_summary(out, summary.other);
      out << "},\"mapping_quality\":{\"alternate\":";
      write_mapping_quality_summary(out, summary.alternate);
      out << ",\"reference\":";
      write_mapping_quality_summary(out, summary.reference);
      out << ",\"other\":";
      write_mapping_quality_summary(out, summary.other);
      out << "},\"read_position\":{\"definition\":"
          << quoted("normalized_center_proximity") << ",\"alternate_mean\":";
      write_nullable_mean(out, summary.alternate.read_position_sum,
                          summary.alternate.read_position_count);
      out << ",\"reference_mean\":";
      write_nullable_mean(out, summary.reference.read_position_sum,
                          summary.reference.read_position_count);
      out << ",\"other_mean\":";
      write_nullable_mean(out, summary.other.read_position_sum,
                          summary.other.read_position_count);
      out << ",\"bias_delta\":";
      if (summary.alternate.read_position_count == 0U ||
          summary.reference.read_position_count == 0U) {
        out << "null";
      } else {
        const double alternate_mean =
            summary.alternate.read_position_sum /
            static_cast<double>(summary.alternate.read_position_count);
        const double reference_mean =
            summary.reference.read_position_sum /
            static_cast<double>(summary.reference.read_position_count);
        out << std::abs(alternate_mean - reference_mean);
      }
      out << "},\"supporting_molecule_ids\":";
      write_string_array(out, event.supporting_molecules);
      out << ",\"supporting_reads\":";
      write_string_array(out, event.supporting_molecules);
      if (snp != nullptr && !snp->gene.empty()) {
        out << ",\"gene\":" << quoted(snp->gene);
      }
      if (snp != nullptr && !snp->consequence.empty()) {
        out << ",\"consequence\":" << quoted(snp->consequence);
      }
      if (snp != nullptr && !snp->protein.empty()) {
        out << ",\"protein\":" << quoted(snp->protein);
      }
      if (snp != nullptr && !snp->residue.empty()) {
        out << ",\"residue\":" << quoted(snp->residue);
      }
      if (snp != nullptr && has_clinical_annotation(*snp)) {
        out << ",\"annotation\":";
        write_clinical_annotation(out, *snp);
      }
      if (snp != nullptr && !snp->structure_id.empty()) {
        out << ",\"structure\":{\"structure_id\":" << quoted(snp->structure_id);
        if (!snp->structure_chain.empty()) {
          out << ",\"chain\":" << quoted(snp->structure_chain);
        }
        if (snp->structure_residue != 0) {
          out << ",\"residue_index\":" << snp->structure_residue;
        }
        if (!snp->structure_complex.empty()) {
          out << ",\"complex\":" << quoted(snp->structure_complex);
        }
        out << "}";
      }
      out << "}";
    }
  } else {
    for (const auto &[_, aggregate] : variants) {
      if (variant_index++ != 0) {
        out << ",";
      }
      const auto &snp = aggregate.call;
      out << "{\"position\":" << snp.position << ","
          << "\"ref\":" << quoted(std::string(1, snp.reference)) << ","
          << "\"alt\":" << quoted(std::string(1, snp.alternate)) << ","
          << "\"alt_depth\":" << aggregate.alternate_depth << ","
          << "\"ref_depth\":" << aggregate.reference_depth << ","
          << "\"other_depth\":" << aggregate.other_depth << ","
          << "\"callable_depth\":" << aggregate.callable_depth << ","
          << "\"heteroplasmy\":" << aggregate.heteroplasmy << ","
          << "\"ci95_low\":" << aggregate.ci95_low << ","
          << "\"ci95_high\":" << aggregate.ci95_high << ","
          << "\"molecule_support\":{"
          << "\"alternate\":" << aggregate.alternate_depth << ","
          << "\"reference\":" << aggregate.reference_depth << ","
          << "\"other\":" << aggregate.other_depth << ","
          << "\"callable\":" << aggregate.callable_depth << "},"
          << "\"strand_support\":{"
          << "\"alt_forward\":" << aggregate.alternate_strand_depths[0] << ","
          << "\"alt_reverse\":" << aggregate.alternate_strand_depths[1] << ","
          << "\"ref_forward\":" << aggregate.reference_strand_depths[0] << ","
          << "\"ref_reverse\":" << aggregate.reference_strand_depths[1] << ","
          << "\"other_forward\":" << aggregate.other_strand_depths[0] << ","
          << "\"other_reverse\":" << aggregate.other_strand_depths[1] << "},"
          << "\"strand_bias_delta\":";
      if (aggregate.alternate_depth == 0U || aggregate.reference_depth == 0U) {
        out << "null";
      } else {
        const double alternate_forward_fraction =
            static_cast<double>(aggregate.alternate_strand_depths[0]) /
            static_cast<double>(aggregate.alternate_depth);
        const double reference_forward_fraction =
            static_cast<double>(aggregate.reference_strand_depths[0]) /
            static_cast<double>(aggregate.reference_depth);
        out << std::abs(alternate_forward_fraction -
                        reference_forward_fraction);
      }
      out << ",\"allele_quality\":{\"alternate\":";
      write_quality_summary(
          out, aggregate.alternate_depth, aggregate.alternate_quality_sum,
          aggregate.alternate_quality_min, aggregate.alternate_quality_max);
      out << ",\"reference\":";
      write_quality_summary(
          out, aggregate.reference_depth, aggregate.reference_quality_sum,
          aggregate.reference_quality_min, aggregate.reference_quality_max);
      out << ",\"other\":";
      write_quality_summary(
          out, aggregate.other_depth, aggregate.other_quality_sum,
          aggregate.other_quality_min, aggregate.other_quality_max);
      out << "},\"read_position\":{"
          << "\"definition\":" << quoted("normalized_center_proximity") << ","
          << "\"alternate_mean\":";
      write_nullable_mean(out, aggregate.alternate_read_position_sum,
                          aggregate.alternate_depth);
      out << ",\"reference_mean\":";
      write_nullable_mean(out, aggregate.reference_read_position_sum,
                          aggregate.reference_depth);
      out << ",\"other_mean\":";
      write_nullable_mean(out, aggregate.other_read_position_sum,
                          aggregate.other_depth);
      out << ",\"bias_delta\":";
      if (aggregate.alternate_depth == 0U || aggregate.reference_depth == 0U) {
        out << "null";
      } else {
        const double alternate_mean =
            aggregate.alternate_read_position_sum /
            static_cast<double>(aggregate.alternate_depth);
        const double reference_mean =
            aggregate.reference_read_position_sum /
            static_cast<double>(aggregate.reference_depth);
        out << std::abs(alternate_mean - reference_mean);
      }
      out << "},\"supporting_reads\":";
      write_string_array(out, aggregate.supporting_reads);
      if (!snp.gene.empty()) {
        out << ",\"gene\":" << quoted(snp.gene);
      }
      if (!snp.consequence.empty()) {
        out << ",\"consequence\":" << quoted(snp.consequence);
      }
      if (!snp.protein.empty()) {
        out << ",\"protein\":" << quoted(snp.protein);
      }
      if (!snp.residue.empty()) {
        out << ",\"residue\":" << quoted(snp.residue);
      }
      if (has_clinical_annotation(snp)) {
        out << ",\"annotation\":";
        write_clinical_annotation(out, snp);
      }
      if (!snp.structure_id.empty()) {
        out << ",\"structure\":{"
            << "\"structure_id\":" << quoted(snp.structure_id);
        if (!snp.structure_chain.empty()) {
          out << ",\"chain\":" << quoted(snp.structure_chain);
        }
        if (snp.structure_residue != 0) {
          out << ",\"residue_index\":" << snp.structure_residue;
        }
        if (!snp.structure_complex.empty()) {
          out << ",\"complex\":" << quoted(snp.structure_complex);
        }
        out << "}";
      }
      out << "}";
    }
  }
  out << "],";

  out << "\"clusters\":[";
  std::size_t cluster_index = 0;
  for (const auto &[cluster_id, members] : clusters) {
    if (cluster_index++ != 0) {
      out << ",";
    }
    std::map<std::string, std::size_t> sv_counts;
    std::map<std::string, std::size_t> complex_event_counts;
    std::vector<std::string> read_ids;
    for (const auto *member : members) {
      read_ids.push_back(member->id);
      for (const auto &sv_id : member->sv_ids) {
        ++sv_counts[sv_id];
      }
      for (const auto &complex_event_id : member->complex_event_ids) {
        ++complex_event_counts[complex_event_id];
      }
    }
    const auto assignment_it = haplogroups.find(cluster_id);
    const ClusterHaplogroupAssignment *assignment =
        assignment_it == haplogroups.end() ? nullptr : &assignment_it->second;
    out << "{"
        << "\"id\":" << cluster_id << ","
        << "\"label\":"
        << quoted(cluster_id < 0 ? "Outliers"
                                 : "Cluster " + std::to_string(cluster_id + 1))
        << ","
        << "\"haplogroup\":"
        << quoted(assignment == nullptr ? "unassigned" : assignment->best)
        << ","
        << "\"haplogroup_assignment\":{";
    out << "\"resource\":" << quoted("phylotree-rcrs@17.3") << ","
        << "\"quality\":" << (assignment == nullptr ? 0.0 : assignment->quality)
        << ","
        << "\"contamination_warning\":"
        << (assignment != nullptr && assignment->contamination_warning
                ? "true"
                : "false")
        << ","
        << "\"observed_markers\":";
    if (assignment == nullptr) {
      out << "[]";
    } else {
      write_string_array(out, assignment->observed_markers);
    }
    out << ",\"callable_ranges\":[";
    if (assignment != nullptr) {
      for (std::size_t range_index = 0;
           range_index < assignment->callable_ranges.size(); ++range_index) {
        if (range_index != 0U) {
          out << ",";
        }
        const auto &[start, end] = assignment->callable_ranges[range_index];
        out << "{\"start\":" << start << ",\"end\":" << end << "}";
      }
    }
    out << "],\"candidates\":[";
    if (assignment != nullptr) {
      for (std::size_t candidate_index = 0;
           candidate_index < assignment->candidates.size(); ++candidate_index) {
        if (candidate_index != 0) {
          out << ",";
        }
        const auto &candidate = assignment->candidates[candidate_index];
        out << "{\"name\":" << quoted(candidate.name)
            << ",\"score\":" << candidate.score << ",\"matched\":";
        write_string_array(out, candidate.matched);
        out << ",\"missing\":";
        write_string_array(out, candidate.missing);
        out << ",\"extra\":";
        write_string_array(out, candidate.extra);
        out << "}";
      }
    }
    out << "]},"
        << "\"size\":" << members.size() << ","
        << "\"outlier\":" << (cluster_id < 0 ? "true" : "false") << ","
        << "\"consensus_haplotype\":"
        << quoted(cluster_id < 0
                      ? "noise"
                      : "feature-consensus-" + std::to_string(cluster_id))
        << ","
        << "\"reads\":[";
    for (std::size_t i = 0; i < read_ids.size(); ++i) {
      if (i != 0) {
        out << ",";
      }
      out << quoted(read_ids[i]);
    }
    out << "],\"sv_signature\":[";
    std::size_t signature_index = 0;
    for (const auto &[sv_id, count] : sv_counts) {
      if (signature_index++ != 0) {
        out << ",";
      }
      out << "{"
          << "\"sv_id\":" << quoted(sv_id) << ","
          << "\"support\":" << count << "}";
    }
    out << "],\"complex_event_signature\":[";
    std::size_t complex_signature_index = 0;
    for (const auto &[event_id, count] : complex_event_counts) {
      if (complex_signature_index++ != 0U) {
        out << ",";
      }
      out << "{"
          << "\"event_id\":" << quoted(event_id) << ","
          << "\"support\":" << count << "}";
    }
    out << "]}";
  }
  out << "],";

  out << "\"reads\":[";
  for (std::size_t i = 0; i < features.size(); ++i) {
    if ((i & 1023U) == 0U) {
      throw_if_cancelled(config);
    }
    if (i != 0) {
      out << ",";
    }
    const auto &feature = features[i];
    out << "{"
        << "\"id\":" << quoted(feature.id) << ","
        << "\"length\":" << feature.length << ","
        << "\"mean_quality\":" << feature.mean_quality << ","
        << "\"numt_score\":" << feature.numt_score << ","
        << "\"filtered_numt\":" << (feature.filtered_numt ? "true" : "false")
        << ","
        << "\"numt_evidence\":";
    write_string_array(out, feature.numt_evidence);
    out << ","
        << "\"cluster_id\":" << feature.cluster_id << ","
        << "\"outlier\":" << (feature.outlier ? "true" : "false") << ","
        << "\"mapping_quality\":" << static_cast<int>(feature.mapping_quality)
        << ","
        << "\"flags\":" << feature.flags << ",";
    if (!feature.reference_name.empty()) {
      out << "\"reference_name\":" << quoted(feature.reference_name) << ",";
    }
    if (!feature.aux_tags.empty()) {
      out << "\"aux_tags\":";
      write_string_map(out, feature.aux_tags);
      out << ",";
    }
    out << "\"snps\":[";
    for (std::size_t snp_index = 0; snp_index < feature.snps.size();
         ++snp_index) {
      if (snp_index != 0) {
        out << ",";
      }
      const auto &snp = feature.snps[snp_index];
      out << "{"
          << "\"position\":" << snp.position << ","
          << "\"ref\":" << quoted(std::string(1, snp.reference)) << ","
          << "\"alt\":" << quoted(std::string(1, snp.alternate));
      out << "}";
    }
    out << "],\"haplogroup_markers\":[";
    for (std::size_t marker_index = 0;
         marker_index < feature.haplogroup_markers.size(); ++marker_index) {
      if (marker_index != 0U) {
        out << ",";
      }
      out << quoted(feature.haplogroup_markers[marker_index].encoded);
    }
    out << "],\"haplogroup_range_known\":"
        << (feature.haplogroup_range_known ? "true" : "false")
        << ",\"haplogroup_callable_ranges\":[";
    for (std::size_t range_index = 0;
         range_index < feature.haplogroup_ranges.size(); ++range_index) {
      if (range_index != 0U) {
        out << ",";
      }
      const auto &[start, end] = feature.haplogroup_ranges[range_index];
      out << "{\"start\":" << start << ",\"end\":" << end << "}";
    }
    out << "],\"sv_ids\":[";
    for (std::size_t sv_id_index = 0; sv_id_index < feature.sv_ids.size();
         ++sv_id_index) {
      if (sv_id_index != 0) {
        out << ",";
      }
      out << quoted(feature.sv_ids[sv_id_index]);
    }
    out << "],\"complex_event_ids\":";
    write_string_array(out, feature.complex_event_ids);
    out << "}";
  }
  out << "]";
  out << "}";
  return std::move(out).take();
}

} // namespace mito::detail
