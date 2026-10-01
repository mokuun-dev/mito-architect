#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] std::vector<ClusterHaplogroupAssignment>
assign_architecture_haplogroups(
    const std::vector<ReadFeature> &features,
    const ArchitectureInferenceResult &architecture_inference,
    const std::vector<HaplogroupDefinition> &definitions,
    const std::unordered_map<std::string, double> &weights,
    const std::vector<PhyloAlignmentRule> &alignment_rules,
    const AnalysisConfig &config) {
  std::vector<ReadFeature> projected;
  std::size_t projected_count = 0U;
  for (const auto &architecture : architecture_inference.architectures) {
    projected_count += architecture.member_molecule_indices.size();
  }
  projected.reserve(projected_count);

  for (std::size_t architecture_index = 0U;
       architecture_index < architecture_inference.architectures.size();
       ++architecture_index) {
    if (architecture_index >
        static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                          "architecture haplogroup index exceeds int range");
    }
    for (const auto molecule_index :
         architecture_inference.architectures[architecture_index]
             .member_molecule_indices) {
      if (molecule_index >= features.size()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "architecture haplogroup molecule is unresolved");
      }
      auto feature = features[molecule_index];
      feature.cluster_id = static_cast<int>(architecture_index);
      feature.outlier = false;
      projected.push_back(std::move(feature));
    }
  }

  const auto ranked = assign_haplogroups(projected, definitions, weights,
                                         alignment_rules, config);
  std::vector<ClusterHaplogroupAssignment> result(
      architecture_inference.architectures.size());
  for (std::size_t index = 0U; index < result.size(); ++index) {
    if (const auto found = ranked.find(static_cast<int>(index));
        found != ranked.end()) {
      result[index] = found->second;
    }
  }
  return result;
}

} // namespace mito::detail
