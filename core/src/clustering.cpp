#include "detail/pipeline.hpp"
#ifdef MITO_HAS_HDBSCAN_CPP
#include <Hdbscan/hdbscan.hpp>
#endif

namespace mito::detail {

[[nodiscard]] FeatureTokens feature_tokens(const ReadFeature &feature) {
  FeatureTokens tokens;
  tokens.reserve(feature.snps.size() + feature.haplogroup_markers.size() +
                 feature.sv_ids.size() + feature.complex_event_ids.size() + 1U);
  for (const auto &snp : feature.snps) {
    tokens.push_back("snp:" +
                     snp_key(snp.position, snp.reference, snp.alternate));
  }
  for (const auto &marker : feature.haplogroup_markers) {
    if (marker.kind != PhyloMutationKind::substitution) {
      tokens.push_back("phylo_indel:" + marker.encoded);
    }
  }
  for (const auto &sv_id : feature.sv_ids) {
    tokens.push_back("sv:" + sv_id);
  }
  for (const auto &complex_event_id : feature.complex_event_ids) {
    tokens.push_back("complex_sv:" + complex_event_id);
  }
  if (tokens.empty()) {
    tokens.push_back("length:" + std::to_string(feature.length / 100U));
  }
  std::sort(tokens.begin(), tokens.end());
  tokens.erase(std::unique(tokens.begin(), tokens.end()), tokens.end());
  return tokens;
}

[[nodiscard]] double token_distance(const FeatureTokens &lhs,
                                    const FeatureTokens &rhs) {
  std::size_t intersection = 0;
  auto lhs_it = lhs.begin();
  auto rhs_it = rhs.begin();
  while (lhs_it != lhs.end() && rhs_it != rhs.end()) {
    if (*lhs_it == *rhs_it) {
      ++intersection;
      ++lhs_it;
      ++rhs_it;
    } else if (*lhs_it < *rhs_it) {
      ++lhs_it;
    } else {
      ++rhs_it;
    }
  }

  const std::size_t union_size = lhs.size() + rhs.size() - intersection;
  if (union_size == 0) {
    return 0.0;
  }
  return 1.0 -
         (static_cast<double>(intersection) / static_cast<double>(union_size));
}

[[nodiscard]] InvertedTokenIndex
build_inverted_token_index(const std::vector<TokenProfile> &profiles) {
  InvertedTokenIndex index;
  for (std::size_t profile_index = 0; profile_index < profiles.size();
       ++profile_index) {
    for (const auto &token : profiles[profile_index].tokens) {
      index[token].push_back(profile_index);
    }
  }
  return index;
}

[[nodiscard]] Neighborhood
region_query(const std::vector<TokenProfile> &profiles,
             const InvertedTokenIndex &inverted_index,
             std::vector<std::size_t> &candidate_marks, std::size_t &mark_epoch,
             std::size_t index, double epsilon) {
  Neighborhood neighborhood;
  std::vector<std::size_t> candidates;
  if (epsilon >= 1.0) {
    candidates.resize(profiles.size());
    std::iota(candidates.begin(), candidates.end(), 0U);
  } else {
    if (mark_epoch == std::numeric_limits<std::size_t>::max()) {
      std::fill(candidate_marks.begin(), candidate_marks.end(), 0U);
      mark_epoch = 0;
    }
    ++mark_epoch;
    for (const auto &token : profiles[index].tokens) {
      const auto posting = inverted_index.find(token);
      if (posting == inverted_index.end()) {
        continue;
      }
      for (const auto candidate : posting->second) {
        if (candidate_marks[candidate] != mark_epoch) {
          candidate_marks[candidate] = mark_epoch;
          candidates.push_back(candidate);
        }
      }
    }
    std::sort(candidates.begin(), candidates.end());
  }

  for (const auto candidate : candidates) {
    const auto lhs_size = profiles[index].tokens.size();
    const auto rhs_size = profiles[candidate].tokens.size();
    const double maximum_similarity =
        static_cast<double>(std::min(lhs_size, rhs_size)) /
        static_cast<double>(std::max(lhs_size, rhs_size));
    if (maximum_similarity < 1.0 - epsilon) {
      continue;
    }
    if (token_distance(profiles[index].tokens, profiles[candidate].tokens) <=
        epsilon) {
      neighborhood.profiles.push_back(candidate);
      neighborhood.read_count += profiles[candidate].feature_indices.size();
    }
  }
  return neighborhood;
}

void append_new_neighbors(std::vector<std::size_t> &queue,
                          std::vector<unsigned char> &queued,
                          const std::vector<std::size_t> &candidates) {
  for (const std::size_t candidate : candidates) {
    if (candidate >= queued.size() || queued[candidate] != 0U) {
      continue;
    }
    queued[candidate] = 1U;
    queue.push_back(candidate);
  }
}

void assign_dbscan_clusters(std::vector<ReadFeature> &features, double epsilon,
                            std::size_t min_cluster_size,
                            const AnalysisConfig &config) {
  std::vector<std::size_t> active_indices;
  active_indices.reserve(features.size());
  for (std::size_t i = 0; i < features.size(); ++i) {
    if (!features[i].filtered_numt) {
      active_indices.push_back(i);
      features[i].cluster_id = -1;
      features[i].outlier = false;
    }
  }

  std::vector<TokenProfile> profiles;
  profiles.reserve(active_indices.size());
  {
    std::map<FeatureTokens, std::size_t> profile_lookup;
    for (const auto index : active_indices) {
      auto tokens = feature_tokens(features[index]);
      const auto existing = profile_lookup.find(tokens);
      if (existing != profile_lookup.end()) {
        profiles[existing->second].feature_indices.push_back(index);
        continue;
      }
      const auto profile_index = profiles.size();
      profile_lookup.emplace(tokens, profile_index);
      profiles.push_back({std::move(tokens), {index}});
    }
  }

  constexpr int kUnvisited = -2;
  std::vector<int> labels(profiles.size(), kUnvisited);
  const auto inverted_index = build_inverted_token_index(profiles);
  std::vector<std::size_t> candidate_marks(profiles.size(), 0U);
  std::size_t mark_epoch = 0;
  const std::size_t min_points = std::max<std::size_t>(1, min_cluster_size);
  int next_cluster_id = 0;

  for (std::size_t point = 0; point < profiles.size(); ++point) {
    throw_if_cancelled(config);
    if (labels[point] != kUnvisited) {
      continue;
    }

    auto neighborhood = region_query(profiles, inverted_index, candidate_marks,
                                     mark_epoch, point, epsilon);
    if (neighborhood.read_count < min_points) {
      labels[point] = -1;
      continue;
    }

    const int cluster_id = next_cluster_id++;
    labels[point] = cluster_id;
    auto &neighbors = neighborhood.profiles;
    std::vector<unsigned char> queued(profiles.size(), 0U);
    for (const auto neighbor : neighbors) {
      if (neighbor < queued.size()) {
        queued[neighbor] = 1U;
      }
    }
    for (std::size_t cursor = 0; cursor < neighbors.size(); ++cursor) {
      throw_if_cancelled(config);
      const std::size_t neighbor = neighbors[cursor];
      if (labels[neighbor] == -1) {
        labels[neighbor] = cluster_id;
      }
      if (labels[neighbor] != kUnvisited) {
        continue;
      }
      labels[neighbor] = cluster_id;
      auto expanded = region_query(profiles, inverted_index, candidate_marks,
                                   mark_epoch, neighbor, epsilon);
      if (expanded.read_count >= min_points) {
        append_new_neighbors(neighbors, queued, expanded.profiles);
      }
    }
  }

  for (std::size_t i = 0; i < profiles.size(); ++i) {
    for (const auto feature_index : profiles[i].feature_indices) {
      auto &feature = features[feature_index];
      feature.cluster_id = labels[i];
      feature.outlier = labels[i] < 0;
    }
  }
}

#ifdef MITO_HAS_HDBSCAN_CPP
[[nodiscard]] bool assign_hdbscan_clusters(std::vector<ReadFeature> &features,
                                           std::size_t min_cluster_size,
                                           const AnalysisConfig &config) {
  std::vector<std::size_t> active_indices;
  active_indices.reserve(features.size());
  for (std::size_t i = 0; i < features.size(); ++i) {
    if (!features[i].filtered_numt) {
      active_indices.push_back(i);
    }
  }
  if (active_indices.empty()) {
    return true;
  }

  std::vector<FeatureTokens> token_sets;
  token_sets.reserve(active_indices.size());
  std::map<std::string, std::size_t> token_columns;
  for (const auto index : active_indices) {
    token_sets.push_back(feature_tokens(features[index]));
    for (const auto &token : token_sets.back()) {
      if (!token_columns.contains(token)) {
        token_columns.emplace(token, token_columns.size());
      }
    }
  }

  std::vector<std::vector<double>> dataset(
      token_sets.size(),
      std::vector<double>(std::max<std::size_t>(1, token_columns.size()), 0.0));
  for (std::size_t row = 0; row < token_sets.size(); ++row) {
    for (const auto &token : token_sets[row]) {
      dataset[row][token_columns.at(token)] = 1.0;
    }
  }

  try {
    throw_if_cancelled(config);
    Hdbscan hdbscan("");
    hdbscan.dataset = std::move(dataset);
    const auto min_points =
        static_cast<int>(std::max<std::size_t>(1, min_cluster_size));
    hdbscan.execute(min_points, min_points, "Euclidean");
    if (hdbscan.normalizedLabels_.size() != active_indices.size()) {
      return false;
    }
    for (std::size_t i = 0; i < active_indices.size(); ++i) {
      const int label = hdbscan.normalizedLabels_[i];
      auto &feature = features[active_indices[i]];
      feature.cluster_id = label <= 0 ? -1 : label - 1;
      feature.outlier = label <= 0;
    }
    return true;
  } catch (const std::bad_alloc &) {
    throw;
  } catch (...) {
    throw_if_cancelled(config);
    return false;
  }
}
#endif

} // namespace mito::detail
