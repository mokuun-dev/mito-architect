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

[[nodiscard]] double token_distance(const TokenRanks &lhs,
                                    const TokenRanks &rhs) {
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
build_inverted_token_index(const std::vector<TokenProfile> &profiles,
                           const std::size_t token_count) {
  InvertedTokenIndex index(token_count);
  for (std::size_t profile_index = 0; profile_index < profiles.size();
       ++profile_index) {
    const auto &profile = profiles[profile_index];
    for (std::size_t i = 0; i < profile.prefix_length; ++i) {
      index[profile.tokens[i]].push_back(
          {profile_index, static_cast<std::uint32_t>(i)});
    }
  }
  return index;
}

namespace {

// Smallest overlap two profiles of the given sizes need for Jaccard similarity
// >= 1 - epsilon (J >= t  <=>  overlap >= t / (1 + t) * (|x| + |y|)). The
// bound is lowered slightly, so rounding can only keep extra candidates.
[[nodiscard]] std::size_t required_jaccard_overlap(const std::size_t lhs_size,
                                                   const std::size_t rhs_size,
                                                   const double epsilon) {
  const double similarity = (1.0 - epsilon) - 1.0e-9;
  if (similarity <= 0.0) {
    return 1U;
  }
  const double bound = similarity / (1.0 + similarity) *
                       static_cast<double>(lhs_size + rhs_size);
  return std::max<std::size_t>(1U, static_cast<std::size_t>(std::ceil(bound)));
}

} // namespace

[[nodiscard]] Neighborhood
region_query(const std::vector<TokenProfile> &profiles,
             const InvertedTokenIndex &inverted_index,
             RegionQueryScratch &scratch, std::size_t index, double epsilon) {
  Neighborhood neighborhood;
  std::vector<std::size_t> candidates;
  if (epsilon >= 1.0) {
    candidates.resize(profiles.size());
    std::iota(candidates.begin(), candidates.end(), 0U);
  } else {
    if (scratch.epoch == std::numeric_limits<std::size_t>::max()) {
      std::fill(scratch.marks.begin(), scratch.marks.end(), 0U);
      scratch.epoch = 0;
    }
    const auto epoch = ++scratch.epoch;
    const auto &profile = profiles[index];
    const auto size = profile.tokens.size();
    // Prefix filtering with positional pruning (PPJoin): a candidate is
    // dropped once the overlap it could still reach falls below the required
    // overlap. Common tokens before the current one are all in both indexed
    // prefixes, so overlaps[] counts them exactly and the bound is safe.
    for (std::size_t i = 0; i < profile.prefix_length; ++i) {
      for (const auto &posting : inverted_index[profile.tokens[i]]) {
        const auto candidate = posting.profile;
        if (scratch.marks[candidate] != epoch) {
          scratch.marks[candidate] = epoch;
          scratch.overlaps[candidate] = 0U;
          scratch.pruned[candidate] = 0U;
          candidates.push_back(candidate);
        } else if (scratch.pruned[candidate] != 0U) {
          continue;
        }
        const auto candidate_size = profiles[candidate].tokens.size();
        const auto remaining =
            std::min(size - i - 1U, candidate_size - posting.position - 1U);
        if (scratch.overlaps[candidate] + 1U + remaining <
            required_jaccard_overlap(size, candidate_size, epsilon)) {
          scratch.pruned[candidate] = 1U;
        } else {
          ++scratch.overlaps[candidate];
        }
      }
    }
    std::erase_if(candidates, [&](const std::size_t candidate) {
      return scratch.pruned[candidate] != 0U;
    });
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

namespace {

// Number of leading (rarest) tokens that must be indexed so that every pair
// with Jaccard similarity >= 1 - epsilon shares at least one indexed token
// (prefix-filtering lemma). The similarity bound is lowered slightly so that
// floating-point rounding in token_distance can only enlarge the candidate
// set; every candidate is still verified with the exact distance.
[[nodiscard]] std::size_t jaccard_prefix_length(const std::size_t token_count,
                                                const double epsilon) {
  if (token_count == 0U || epsilon >= 1.0) {
    return token_count;
  }
  const double minimum_similarity = (1.0 - epsilon) - 1.0e-9;
  std::size_t required_overlap = 1U;
  if (minimum_similarity > 0.0) {
    required_overlap = std::max<std::size_t>(
        1U, static_cast<std::size_t>(std::ceil(
                minimum_similarity * static_cast<double>(token_count))));
  }
  required_overlap = std::min(required_overlap, token_count);
  return token_count - required_overlap + 1U;
}

} // namespace

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

  // Intern string tokens; identical token sets collapse into one profile in
  // first-seen order, exactly as with string keys.
  std::vector<TokenProfile> profiles;
  profiles.reserve(active_indices.size());
  std::size_t token_count = 0;
  {
    std::unordered_map<std::string, std::uint32_t> token_ids;
    std::map<TokenRanks, std::size_t> profile_lookup;
    for (const auto index : active_indices) {
      TokenRanks ids;
      for (auto &token : feature_tokens(features[index])) {
        const auto next_id = static_cast<std::uint32_t>(token_ids.size());
        if (token_ids.size() >= std::numeric_limits<std::uint32_t>::max()) {
          throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                              "clustering token dictionary overflow");
        }
        ids.push_back(
            token_ids.try_emplace(std::move(token), next_id).first->second);
      }
      std::sort(ids.begin(), ids.end());
      const auto existing = profile_lookup.find(ids);
      if (existing != profile_lookup.end()) {
        profiles[existing->second].feature_indices.push_back(index);
        continue;
      }
      const auto profile_index = profiles.size();
      profile_lookup.emplace(ids, profile_index);
      profiles.push_back({std::move(ids), 0U, {index}});
    }
    token_count = token_ids.size();
  }

  // Global token order for prefix filtering: rare tokens first, ties by id.
  {
    std::vector<std::size_t> frequency(token_count, 0U);
    for (const auto &profile : profiles) {
      for (const auto id : profile.tokens) {
        ++frequency[id];
      }
    }
    std::vector<std::uint32_t> order(token_count);
    std::iota(order.begin(), order.end(), 0U);
    std::sort(order.begin(), order.end(),
              [&](const std::uint32_t lhs, const std::uint32_t rhs) {
                return frequency[lhs] != frequency[rhs]
                           ? frequency[lhs] < frequency[rhs]
                           : lhs < rhs;
              });
    std::vector<std::uint32_t> rank(token_count);
    for (std::size_t position = 0; position < order.size(); ++position) {
      rank[order[position]] = static_cast<std::uint32_t>(position);
    }
    for (auto &profile : profiles) {
      for (auto &id : profile.tokens) {
        id = rank[id];
      }
      std::sort(profile.tokens.begin(), profile.tokens.end());
      profile.prefix_length =
          jaccard_prefix_length(profile.tokens.size(), epsilon);
    }
  }

  constexpr int kUnvisited = -2;
  std::vector<int> labels(profiles.size(), kUnvisited);
  const auto inverted_index = build_inverted_token_index(profiles, token_count);
  RegionQueryScratch scratch{std::vector<std::size_t>(profiles.size(), 0U),
                             std::vector<std::size_t>(profiles.size(), 0U),
                             std::vector<unsigned char>(profiles.size(), 0U),
                             0U};
  const std::size_t min_points = std::max<std::size_t>(1, min_cluster_size);
  int next_cluster_id = 0;
  // Reset after each cluster through the neighbor list, which holds every
  // queued profile, instead of allocating a fresh array per cluster.
  std::vector<unsigned char> queued(profiles.size(), 0U);

  for (std::size_t point = 0; point < profiles.size(); ++point) {
    throw_if_cancelled(config);
    if (labels[point] != kUnvisited) {
      continue;
    }

    auto neighborhood =
        region_query(profiles, inverted_index, scratch, point, epsilon);
    if (neighborhood.read_count < min_points) {
      labels[point] = -1;
      continue;
    }

    const int cluster_id = next_cluster_id++;
    labels[point] = cluster_id;
    auto &neighbors = neighborhood.profiles;
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
      auto expanded =
          region_query(profiles, inverted_index, scratch, neighbor, epsilon);
      if (expanded.read_count >= min_points) {
        append_new_neighbors(neighbors, queued, expanded.profiles);
      }
    }
    for (const auto neighbor : neighbors) {
      queued[neighbor] = 0U;
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
