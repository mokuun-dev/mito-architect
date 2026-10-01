#include "detail/pipeline.hpp"

namespace mito::detail {
namespace {

[[nodiscard]] std::uint64_t stable_molecule_hash(std::string_view value) {
  std::uint64_t hash = 1469598103934665603ULL;
  for (const unsigned char byte : value) {
    hash ^= byte;
    hash *= 1099511628211ULL;
  }
  return hash;
}

struct RunningVariance {
  std::size_t count = 0U;
  double mean = 0.0;
  double sum_squared_delta = 0.0;

  void add(const double value) {
    ++count;
    const double delta = value - mean;
    mean += delta / static_cast<double>(count);
    sum_squared_delta += delta * (value - mean);
  }

  [[nodiscard]] double standard_deviation() const {
    return count < 2U ? 0.0
                      : std::sqrt(sum_squared_delta /
                                  static_cast<double>(count - 1U));
  }
};

struct Match {
  double similarity = 0.0;
  std::size_t original = 0U;
  std::size_t resampled = 0U;
};

} // namespace

void evaluate_architecture_stability(
    ArchitectureInferenceResult &result, const SparseEvidenceStore &store,
    const MoleculeAssemblyResult &assembly,
    const CallabilityResult &callability,
    const std::vector<ReadFeature> &features,
    const std::size_t reference_length, const AnalysisConfig &config) {
  if (config.architecture_stability_replicates == 0U ||
      result.architectures.empty()) {
    return;
  }
  std::vector<RunningVariance> abundances(result.architectures.size());
  std::vector<double> signature_sums(result.architectures.size(), 0.0);
  std::vector<std::size_t> recovered(result.architectures.size(), 0U);
  AnalysisConfig resample_config = config;
  resample_config.architecture_stability_replicates = 0U;

  for (std::size_t replicate = 0U;
       replicate < config.architecture_stability_replicates; ++replicate) {
    throw_if_cancelled(config);
    auto resampled_assembly = assembly;
    for (std::size_t index = 0U; index < resampled_assembly.molecules.size();
         ++index) {
      const auto selector = mix_architecture_seed(
          config.architecture_seed ^
          (static_cast<std::uint64_t>(replicate) << 32U) ^
          stable_molecule_hash(resampled_assembly.molecules[index].id));
      if (selector % 100U >= 80U) {
        resampled_assembly.molecules[index].analysis_eligible = false;
      }
    }
    const auto resampled = infer_architectures(
        store, resampled_assembly, callability, features, reference_length,
        resample_config);
    std::vector<Match> candidates;
    for (std::size_t original = 0U; original < result.architectures.size();
         ++original) {
      for (std::size_t sample = 0U; sample < resampled.architectures.size();
           ++sample) {
        candidates.push_back(
            {weighted_signature_jaccard(
                 result.architectures[original].defining_event_indices,
                 resampled.architectures[sample].defining_event_indices,
                 store),
             original, sample});
      }
    }
    std::sort(candidates.begin(), candidates.end(), [](const auto &lhs,
                                                        const auto &rhs) {
      if (lhs.similarity != rhs.similarity) {
        return lhs.similarity > rhs.similarity;
      }
      if (lhs.original != rhs.original) {
        return lhs.original < rhs.original;
      }
      return lhs.resampled < rhs.resampled;
    });
    std::vector<std::optional<std::size_t>> matches(result.architectures.size());
    std::vector<bool> sample_used(resampled.architectures.size(), false);
    for (const auto &candidate : candidates) {
      if (candidate.similarity < 0.5 || matches[candidate.original].has_value() ||
          sample_used[candidate.resampled]) {
        continue;
      }
      matches[candidate.original] = candidate.resampled;
      sample_used[candidate.resampled] = true;
      signature_sums[candidate.original] += candidate.similarity;
      ++recovered[candidate.original];
    }
    for (std::size_t original = 0U; original < result.architectures.size();
         ++original) {
      const double abundance = matches[original].has_value()
                                   ? resampled.architectures[*matches[original]]
                                         .estimated_fraction
                                   : 0.0;
      abundances[original].add(abundance);
    }
  }
  const auto replicates = static_cast<double>(config.architecture_stability_replicates);
  for (std::size_t index = 0U; index < result.architectures.size(); ++index) {
    result.architectures[index].signature_stability = signature_sums[index] / replicates;
    result.architectures[index].recovery_rate =
        static_cast<double>(recovered[index]) / replicates;
    result.architectures[index].abundance_standard_deviation =
        abundances[index].standard_deviation();
  }
}

} // namespace mito::detail
