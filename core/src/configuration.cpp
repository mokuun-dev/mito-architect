#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] bool valid_sam_tag_name(const std::string &tag) {
  return tag.empty() ||
         (tag.size() == 2U &&
          std::isalpha(static_cast<unsigned char>(tag[0])) != 0 &&
          std::isalnum(static_cast<unsigned char>(tag[1])) != 0);
}

void validate_config(const AnalysisConfig &config) {
  if (!std::isfinite(config.cluster_epsilon) || config.cluster_epsilon < 0.0 ||
      config.cluster_epsilon > 1.0) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "cluster_epsilon must be finite and between 0 and 1");
  }
  if (!std::isfinite(config.numt_threshold) || config.numt_threshold < 0.0 ||
      config.numt_threshold > 1.0) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "numt_threshold must be finite and between 0 and 1");
  }
  if (config.min_cluster_size == 0) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "min_cluster_size must be at least 1");
  }
  if (config.sv_min_length == 0) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "sv_min_length must be at least 1");
  }
  if (config.min_base_quality > 93U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "min_base_quality must not exceed Phred 93");
  }
  if (config.max_evidence_observations == 0U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "max_evidence_observations must be at least 1");
  }
  if (config.max_phase_links == 0U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "max_phase_links must be at least 1");
  }
  if (config.max_phase_work == 0U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "max_phase_work must be at least 1");
  }
  if (config.max_phase_molecule_references == 0U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "max_phase_molecule_references must be at least 1");
  }
  if (config.max_result_bytes == 0U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "max_result_bytes must be at least 1");
  }
  if (config.min_architecture_molecules < 2U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "min_architecture_molecules must be at least 2");
  }
  if (config.max_candidate_architectures == 0U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "max_candidate_architectures must be at least 1");
  }
  if (config.max_candidate_architectures > 4096U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "max_candidate_architectures must not exceed 4096");
  }
  if (config.architecture_stability_replicates > 10'000U) {
    throw AnalysisError(
        AnalysisErrorCode::invalid_configuration,
        "architecture_stability_replicates must not exceed 10000");
  }
  const auto require_unit_interval = [](const double value, const char *name) {
    if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
      throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                          std::string(name) +
                              " must be finite and between 0 and 1");
    }
  };
  require_unit_interval(config.architecture_max_distance,
                        "architecture_max_distance");
  require_unit_interval(config.architecture_ambiguity_margin,
                        "architecture_ambiguity_margin");
  require_unit_interval(config.architecture_min_overlap_fraction,
                        "architecture_min_overlap_fraction");
  require_unit_interval(config.architecture_consensus_fraction,
                        "architecture_consensus_fraction");
  require_unit_interval(config.architecture_optional_fraction,
                        "architecture_optional_fraction");
  if (config.architecture_consensus_fraction <= 0.5) {
    throw AnalysisError(
        AnalysisErrorCode::invalid_configuration,
        "architecture_consensus_fraction must be greater than 0.5");
  }
  if (config.architecture_optional_fraction >
      config.architecture_consensus_fraction) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "architecture_optional_fraction must not exceed "
                        "architecture_consensus_fraction");
  }
  if (config.evidence_page_size == 0U ||
      config.evidence_page_size > 1'000'000U) {
    throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                        "evidence_page_size must be between 1 and 1000000");
  }
  for (const auto *tag :
       {&config.molecule_id_tag, &config.umi_tag, &config.duplex_tag}) {
    if (!valid_sam_tag_name(*tag)) {
      throw AnalysisError(AnalysisErrorCode::invalid_configuration,
                          "protocol SAM tags must match [A-Za-z][A-Za-z0-9]");
    }
  }
  if (config.result_schema != ResultSchema::v0_6 &&
      (!config.molecule_id_tag.empty() || !config.umi_tag.empty() ||
       !config.duplex_tag.empty())) {
    throw AnalysisError(
        AnalysisErrorCode::invalid_configuration,
        "protocol molecule tags require the schema 0.6 evidence graph");
  }
}

} // namespace mito::detail
