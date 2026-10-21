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
  const std::size_t replicates = config.architecture_stability_replicates;
  AnalysisConfig resample_config = config;
  resample_config.architecture_stability_replicates = 0U;

  // infer_architectures reads only molecule identity and eligibility from the
  // assembly. Resample a slim copy instead of cloning every read record.
  MoleculeAssemblyResult slim_assembly;
  slim_assembly.molecules.resize(assembly.molecules.size());
  std::vector<std::uint64_t> molecule_hashes(assembly.molecules.size());
  for (std::size_t index = 0U; index < assembly.molecules.size(); ++index) {
    const auto &molecule = assembly.molecules[index];
    auto &slim = slim_assembly.molecules[index];
    slim.index = molecule.index;
    slim.id = molecule.id;
    slim.analysis_eligible = molecule.analysis_eligible;
    molecule_hashes[index] = stable_molecule_hash(molecule.id);
  }

  // Per replicate and original architecture: matched similarity and the
  // matched candidate's abundance. Replicates run in parallel; their results
  // are folded below in replicate order so floating-point sums are unchanged.
  using ReplicateMatches =
      std::vector<std::optional<std::pair<double, double>>>;
  std::vector<ReplicateMatches> replicate_matches(replicates);
  const auto run_replicate = [&](const std::size_t replicate, std::size_t) {
    auto resampled_assembly = slim_assembly;
    for (std::size_t index = 0U; index < resampled_assembly.molecules.size();
         ++index) {
      const auto selector =
          mix_architecture_seed(config.architecture_seed ^
                                (static_cast<std::uint64_t>(replicate) << 32U) ^
                                molecule_hashes[index]);
      if (selector % 100U >= 80U) {
        resampled_assembly.molecules[index].analysis_eligible = false;
      }
    }
    const auto resampled =
        infer_architectures(store, resampled_assembly, callability, features,
                            reference_length, resample_config);
    std::vector<Match> candidates;
    for (std::size_t original = 0U; original < result.architectures.size();
         ++original) {
      for (std::size_t sample = 0U; sample < resampled.architectures.size();
           ++sample) {
        candidates.push_back(
            {weighted_signature_jaccard(
                 result.architectures[original].defining_event_indices,
                 resampled.architectures[sample].defining_event_indices, store),
             original, sample});
      }
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const auto &lhs, const auto &rhs) {
                if (lhs.similarity != rhs.similarity) {
                  return lhs.similarity > rhs.similarity;
                }
                if (lhs.original != rhs.original) {
                  return lhs.original < rhs.original;
                }
                return lhs.resampled < rhs.resampled;
              });
    ReplicateMatches matches(result.architectures.size());
    std::vector<bool> sample_used(resampled.architectures.size(), false);
    for (const auto &candidate : candidates) {
      if (candidate.similarity < 0.5 ||
          matches[candidate.original].has_value() ||
          sample_used[candidate.resampled]) {
        continue;
      }
      matches[candidate.original] = std::pair{
          candidate.similarity,
          resampled.architectures[candidate.resampled].estimated_fraction};
      sample_used[candidate.resampled] = true;
    }
    replicate_matches[replicate] = std::move(matches);
  };

  parallel_for_index(replicates, config, 1U, run_replicate);

  std::vector<RunningVariance> abundances(result.architectures.size());
  std::vector<double> signature_sums(result.architectures.size(), 0.0);
  std::vector<std::size_t> recovered(result.architectures.size(), 0U);
  for (const auto &matches : replicate_matches) {
    for (std::size_t original = 0U; original < result.architectures.size();
         ++original) {
      if (matches[original].has_value()) {
        signature_sums[original] += matches[original]->first;
        ++recovered[original];
        abundances[original].add(matches[original]->second);
      } else {
        abundances[original].add(0.0);
      }
    }
  }
  const auto replicate_count = static_cast<double>(replicates);
  for (std::size_t index = 0U; index < result.architectures.size(); ++index) {
    result.architectures[index].signature_stability =
        signature_sums[index] / replicate_count;
    result.architectures[index].recovery_rate =
        static_cast<double>(recovered[index]) / replicate_count;
    result.architectures[index].abundance_standard_deviation =
        abundances[index].standard_deviation();
  }
}

} // namespace mito::detail
