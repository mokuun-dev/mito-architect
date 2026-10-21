#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] std::vector<HaplogroupDefinition> load_haplogroups() {
  std::ifstream input(phylotree_path());
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "could not open PhyloTree resource: " +
                            phylotree_path());
  }
  std::vector<HaplogroupDefinition> definitions;
  std::vector<HaplogroupDefinition> stack;
  std::string line;
  std::size_t mutation_count = 0;
  std::size_t source_order = 0;
  while (std::getline(input, line)) {
    const auto opening = line.find("<haplogroup name=\"");
    if (opening != std::string::npos) {
      const auto name_start =
          opening + std::string_view("<haplogroup name=\"").size();
      const auto name_end = line.find('"', name_start);
      if (name_end == std::string::npos) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "malformed haplogroup name in PhyloTree resource");
      }
      HaplogroupDefinition node;
      node.name = line.substr(name_start, name_end - name_start);
      node.source_order = source_order++;
      if (!stack.empty()) {
        node.mutations = stack.back().mutations;
      }
      stack.push_back(std::move(node));
    }

    const auto poly_start = line.find("<poly>");
    if (poly_start != std::string::npos) {
      if (stack.empty()) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "PhyloTree mutation is outside a haplogroup");
      }
      const auto value_start = poly_start + std::string_view("<poly>").size();
      const auto value_end = line.find("</poly>", value_start);
      if (value_end == std::string::npos) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "unterminated PhyloTree mutation element");
      }
      const auto encoded = trim_copy(
          std::string_view(line).substr(value_start, value_end - value_start));
      const auto mutation = parse_phylo_mutation(encoded);
      if (!mutation) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "unsupported PhyloTree mutation token: " + encoded);
      }
      ++mutation_count;
      auto mutation_it = std::lower_bound(
          stack.back().mutations.begin(), stack.back().mutations.end(),
          mutation->locus, [](const auto &item, const auto &locus) {
            return item.locus < locus;
          });
      const bool same_locus = mutation_it != stack.back().mutations.end() &&
                              mutation_it->locus == mutation->locus;
      if (mutation->backmutation) {
        if (same_locus) {
          stack.back().mutations.erase(mutation_it);
        }
      } else if (same_locus) {
        *mutation_it = *mutation;
      } else {
        stack.back().mutations.insert(mutation_it, *mutation);
      }
    }

    if (line.find("</haplogroup>") != std::string::npos) {
      if (stack.empty()) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "unbalanced PhyloTree haplogroup closing tag");
      }
      definitions.push_back(stack.back());
      stack.pop_back();
    }
  }
  if (!stack.empty() || definitions.size() < 5000U || mutation_count == 0U) {
    throw AnalysisError(AnalysisErrorCode::resource_invalid,
                        "PhyloTree resource is incomplete or unbalanced");
  }
  return definitions;
}

[[nodiscard]] std::unordered_map<std::string, double> load_phylo_weights() {
  std::ifstream input(phylotree_weights_path());
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "could not open PhyloTree weights: " +
                            phylotree_weights_path());
  }
  std::unordered_map<std::string, double> weights;
  std::string line;
  while (std::getline(input, line)) {
    const auto fields = split_tab(line);
    if (fields.size() < 2) {
      continue;
    }
    try {
      const double weight = std::stod(fields[1]);
      if (std::isfinite(weight) && weight > 0.0) {
        weights[fields[0]] = weight;
      }
    } catch (const std::exception &) {
      // Ignore a malformed optional weight and use the deterministic default.
    }
  }
  return weights;
}

[[nodiscard]] double
mutation_weight(const std::unordered_map<std::string, double> &weights,
                const std::string &mutation) {
  const auto it = weights.find(mutation);
  return it == weights.end() ? 1.0 : it->second;
}

[[nodiscard]] std::string macro_haplogroup(std::string_view name) {
  if (name.empty()) {
    return {};
  }
  std::string macro(1, name.front());
  if (name.front() == 'L' && name.size() > 1 &&
      std::isdigit(static_cast<unsigned char>(name[1])) != 0) {
    macro.push_back(name[1]);
  }
  return macro;
}

[[nodiscard]] bool
phylo_mutation_matches(const PhyloMutation &expected,
                       const ObservedPhyloMutation &observed) {
  if (expected.locus.position != observed.position ||
      expected.kind != observed.kind) {
    return false;
  }
  if (expected.kind == PhyloMutationKind::substitution) {
    return expected.alternate == observed.alternate;
  }
  if (expected.kind == PhyloMutationKind::deletion) {
    return true;
  }
  if (expected.inserted_bases.empty() || observed.inserted_bases.empty()) {
    return false;
  }
  if (expected.locus.insertion_index == "X") {
    for (std::size_t i = 0; i < observed.inserted_bases.size(); ++i) {
      if (observed.inserted_bases[i] !=
          expected.inserted_bases[i % expected.inserted_bases.size()]) {
        return false;
      }
    }
    return true;
  }

  const auto expected_index = parse_size(expected.locus.insertion_index);
  const auto observed_index = parse_size(observed.insertion_index);
  if (!expected_index || !observed_index || *expected_index < *observed_index) {
    return false;
  }
  const std::size_t offset = *expected_index - *observed_index;
  return offset <= observed.inserted_bases.size() &&
         expected.inserted_bases.size() <=
             observed.inserted_bases.size() - offset &&
         std::equal(expected.inserted_bases.begin(),
                    expected.inserted_bases.end(),
                    observed.inserted_bases.begin() +
                        static_cast<std::ptrdiff_t>(offset));
}

[[nodiscard]] bool
observed_contains_mutation(const std::vector<ObservedPhyloMutation> &observed,
                           const PhyloMutation &expected) {
  auto candidate = std::lower_bound(observed.begin(), observed.end(),
                                    expected.locus.position,
                                    [](const auto &item, std::size_t position) {
                                      return item.position < position;
                                    });
  while (candidate != observed.end() &&
         candidate->position == expected.locus.position) {
    if (phylo_mutation_matches(expected, *candidate)) {
      return true;
    }
    ++candidate;
  }
  return false;
}

[[nodiscard]] bool
definition_contains_observed(const std::vector<PhyloMutation> &expected,
                             const ObservedPhyloMutation &observed) {
  auto candidate =
      std::lower_bound(expected.begin(), expected.end(), observed.position,
                       [](const auto &item, std::size_t position) {
                         return item.locus.position < position;
                       });
  while (candidate != expected.end() &&
         candidate->locus.position == observed.position) {
    if (phylo_mutation_matches(*candidate, observed)) {
      return true;
    }
    ++candidate;
  }
  return false;
}

[[nodiscard]] bool position_in_ranges(
    const std::vector<std::pair<std::size_t, std::size_t>> &ranges,
    std::size_t position) {
  const auto candidate =
      std::lower_bound(ranges.begin(), ranges.end(), position,
                       [](const auto &range, std::size_t value) {
                         return range.second < value;
                       });
  return candidate != ranges.end() && candidate->first <= position;
}

// Visits mutations whose position is marked callable. Callable ranges are
// normalized (sorted, disjoint), so this visits exactly the mutations that a
// per-range scan would, in the same order, without one binary search per
// range for every candidate definition.
template <typename Visitor>
void for_each_callable_mutation(const std::vector<PhyloMutation> &mutations,
                                const std::vector<char> &callable_positions,
                                Visitor &&visitor) {
  for (const auto &mutation : mutations) {
    const auto position = mutation.locus.position;
    if (position < callable_positions.size() &&
        callable_positions[position] != 0) {
      visitor(mutation);
    }
  }
}

[[nodiscard]] CallablePhyloRanges
majority_callable_ranges(const std::vector<const ReadFeature *> &members) {
  CallablePhyloRanges result;
  if (members.size() == 1U) {
    result.known = members.front()->haplogroup_range_known;
    if (result.known) {
      result.ranges = members.front()->haplogroup_ranges;
    } else {
      result.ranges.emplace_back(
          1U, static_cast<std::size_t>(kDefaultReferenceLength));
    }
    return result;
  }
  const std::size_t known_range_count = static_cast<std::size_t>(
      std::count_if(members.begin(), members.end(), [](const auto *item) {
        return item->haplogroup_range_known;
      }));
  result.known = known_range_count * 2U > members.size();
  if (!result.known) {
    result.ranges.emplace_back(
        1U, static_cast<std::size_t>(kDefaultReferenceLength));
    return result;
  }

  struct RangeEvent {
    std::size_t position = 0U;
    std::int64_t delta = 0;
  };
  std::vector<RangeEvent> events;
  std::size_t range_count = 0U;
  for (const auto *member : members) {
    const auto member_ranges =
        member->haplogroup_range_known ? member->haplogroup_ranges.size() : 0U;
    if (member_ranges > std::numeric_limits<std::size_t>::max() - range_count) {
      throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                          "haplogroup callable range count overflow");
    }
    range_count += member_ranges;
  }
  if (range_count > std::numeric_limits<std::size_t>::max() / 2U) {
    throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                        "haplogroup callable event count overflow");
  }
  events.reserve(range_count * 2U);
  for (const auto *member : members) {
    if (!member->haplogroup_range_known) {
      continue;
    }
    for (const auto &[start, end] : member->haplogroup_ranges) {
      if (start == 0U || end < start ||
          end > static_cast<std::size_t>(kDefaultReferenceLength)) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "haplogroup callable range is outside the reference");
      }
      events.push_back({start, 1});
      events.push_back({end + 1U, -1});
    }
  }
  std::sort(events.begin(), events.end(), [](const auto &lhs, const auto &rhs) {
    return std::tie(lhs.position, lhs.delta) <
           std::tie(rhs.position, rhs.delta);
  });

  const auto threshold = static_cast<std::int64_t>(members.size() / 2U);
  std::int64_t active = 0;
  std::size_t cursor = 1U;
  std::size_t event_index = 0U;
  while (cursor <= static_cast<std::size_t>(kDefaultReferenceLength)) {
    while (event_index < events.size() &&
           events[event_index].position == cursor) {
      active += events[event_index].delta;
      ++event_index;
    }
    if (active < 0) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "haplogroup callable range coverage underflow");
    }
    const std::size_t next =
        event_index < events.size()
            ? std::min(events[event_index].position,
                       static_cast<std::size_t>(kDefaultReferenceLength) + 1U)
            : static_cast<std::size_t>(kDefaultReferenceLength) + 1U;
    if (active > threshold && next > cursor) {
      result.ranges.emplace_back(cursor, next - 1U);
    }
    if (next <= cursor) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "haplogroup callable range events are not ordered");
    }
    cursor = next;
  }
  result.ranges = normalize_ranges(std::move(result.ranges));
  return result;
}

namespace {

// Positions and weights of every definition's mutations in flat arrays, so
// scoring thousands of candidate definitions per cluster streams compact data
// instead of the full PhyloMutation records.
struct CompactDefinitions {
  std::vector<std::size_t> offsets;
  std::vector<std::size_t> positions;
  std::vector<double> weights;

  explicit CompactDefinitions(
      const std::vector<HaplogroupDefinition> &definitions) {
    offsets.reserve(definitions.size() + 1U);
    offsets.push_back(0U);
    for (const auto &definition : definitions) {
      for (const auto &mutation : definition.mutations) {
        positions.push_back(mutation.locus.position);
        weights.push_back(mutation.weight);
      }
      offsets.push_back(positions.size());
    }
  }
};

struct ClusterHaplogroupInput {
  int cluster_id = -1;
  /** Schema 0.6 cluster without known callable ranges: result unavailable. */
  bool unavailable = false;
  std::vector<ObservedPhyloMutation> observed;
  CallablePhyloRanges callable;
  std::string signature;
};

// Majority markers of one cluster, restricted to its callable ranges and
// corrected by the PhyloTree alignment rules.
[[nodiscard]] ClusterHaplogroupInput prepare_cluster_haplogroup_input(
    const int cluster_id, const std::vector<const ReadFeature *> &members,
    const std::vector<PhyloAlignmentRule> &alignment_rules,
    const AnalysisConfig &config) {
  ClusterHaplogroupInput input;
  input.cluster_id = cluster_id;
  struct ObservationCount {
    ObservedPhyloMutation mutation;
    std::size_t count = 0;
  };
  std::unordered_map<std::string_view, ObservationCount> observed_counts;
  for (const auto *member : members) {
    std::unordered_set<std::string_view> molecule_mutations;
    molecule_mutations.reserve(member->haplogroup_markers.size());
    for (const auto &mutation : member->haplogroup_markers) {
      if (!molecule_mutations.insert(mutation.encoded).second) {
        continue;
      }
      auto [count, inserted] = observed_counts.try_emplace(
          mutation.encoded, ObservationCount{mutation, 0U});
      static_cast<void>(inserted);
      ++count->second.count;
    }
  }
  auto &observed = input.observed;
  for (const auto &[_, counted] : observed_counts) {
    if (counted.count * 2U <= members.size()) {
      continue;
    }
    observed.push_back(counted.mutation);
  }
  std::sort(observed.begin(), observed.end(), observed_phylo_mutation_less);

  input.callable = majority_callable_ranges(members);
  const auto &callable = input.callable;
  if (!callable.known && config.result_schema == ResultSchema::v0_6) {
    // A marker without the interval in which it was assessed cannot support
    // a haplogroup conclusion: missing diagnostic markers may simply be
    // outside the molecule.  Keep the result explicitly unavailable.
    input.unavailable = true;
    return input;
  }
  if (callable.known) {
    observed.erase(std::remove_if(observed.begin(), observed.end(),
                                  [&](const auto &mutation) {
                                    return !position_in_ranges(
                                        callable.ranges, mutation.position);
                                  }),
                   observed.end());
  }
  apply_phylo_alignment_rules(observed, alignment_rules);
  if (callable.known) {
    observed.erase(std::remove_if(observed.begin(), observed.end(),
                                  [&](const auto &mutation) {
                                    return !position_in_ranges(
                                        callable.ranges, mutation.position);
                                  }),
                   observed.end());
  }
  for (const auto &mutation : observed) {
    input.signature.append(mutation.encoded);
    input.signature.push_back('\0');
  }
  for (const auto &[start, end] : callable.ranges) {
    input.signature.append("range:");
    input.signature.append(std::to_string(start));
    input.signature.push_back('-');
    input.signature.append(std::to_string(end));
    input.signature.push_back('\0');
  }
  return input;
}

// Scores PhyloTree definitions against one cluster's observed markers.
[[nodiscard]] ClusterHaplogroupAssignment score_cluster_haplogroups(
    const ClusterHaplogroupInput &input,
    const std::vector<HaplogroupDefinition> &definitions,
    const CompactDefinitions &compact,
    const std::vector<std::vector<std::size_t>> &definitions_by_position,
    const std::unordered_map<std::string, double> &weights) {
  const auto &observed = input.observed;
  const auto &callable_ranges = input.callable.ranges;
  constexpr auto kPositions = static_cast<std::size_t>(kDefaultReferenceLength);
  std::vector<char> callable_positions(kPositions + 1U, 0);
  for (const auto &[start, end] : callable_ranges) {
    const auto last = std::min(end, kPositions);
    for (std::size_t position = start; position <= last; ++position) {
      callable_positions[position] = 1;
    }
  }
  // First observed marker at each PhyloTree position; observed is sorted by
  // position, so this replaces a binary search per expected mutation.
  constexpr std::size_t kNoMarker = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> first_observed(kPositions + 1U, kNoMarker);
  for (std::size_t index = observed.size(); index-- > 0U;) {
    if (observed[index].position <= kPositions) {
      first_observed[observed[index].position] = index;
    }
  }
  const auto contains = [&](const PhyloMutation &expected) {
    const auto position = expected.locus.position;
    if (position > kPositions) {
      return false;
    }
    for (auto index = first_observed[position];
         index < observed.size() && observed[index].position == position;
         ++index) {
      if (phylo_mutation_matches(expected, observed[index])) {
        return true;
      }
    }
    return false;
  };

  struct RankedDefinition {
    const HaplogroupDefinition *definition = nullptr;
    double score = 0.0;
    std::size_t matched_count = 0;
  };
  double observed_weight = 0.0;
  for (const auto &mutation : observed) {
    observed_weight += mutation_weight(weights, mutation.encoded);
  }
  std::vector<char> candidate_seen(definitions.size(), 0);
  std::vector<std::size_t> candidate_indices;
  for (const auto &mutation : observed) {
    if (mutation.position >= definitions_by_position.size()) {
      continue;
    }
    for (const auto definition_index :
         definitions_by_position[mutation.position]) {
      if (candidate_seen[definition_index] == 0) {
        candidate_seen[definition_index] = 1;
        candidate_indices.push_back(definition_index);
      }
    }
  }
  std::vector<RankedDefinition> ranked;
  ranked.reserve(candidate_indices.size());
  for (const auto definition_index : candidate_indices) {
    const auto &definition = definitions[definition_index];
    double expected_weight = 0.0;
    double matched_weight = 0.0;
    std::size_t matched_count = 0;
    // Same visit order and summation order as for_each_callable_mutation.
    const auto begin = compact.offsets[definition_index];
    const auto end = compact.offsets[definition_index + 1U];
    for (auto entry = begin; entry < end; ++entry) {
      const auto position = compact.positions[entry];
      if (position >= callable_positions.size() ||
          callable_positions[position] == 0) {
        continue;
      }
      const auto weight = compact.weights[entry];
      expected_weight += weight;
      if (first_observed[position] != kNoMarker &&
          contains(definition.mutations[entry - begin])) {
        matched_weight += weight;
        ++matched_count;
      }
    }
    double score = 0.0;
    if (expected_weight == 0.0 && observed_weight == 0.0) {
      score = 100.0;
    } else if (expected_weight > 0.0 && observed_weight > 0.0) {
      score = 50.0 * ((matched_weight / expected_weight) +
                      (matched_weight / observed_weight));
    }
    if (matched_count > 0 && score > 0.0) {
      ranked.push_back({&definition, score, matched_count});
    }
  }
  const auto ranked_less = [](const auto &lhs, const auto &rhs) {
    if (lhs.score != rhs.score) {
      return lhs.score > rhs.score;
    }
    if (lhs.definition->source_order != rhs.definition->source_order) {
      return lhs.definition->source_order < rhs.definition->source_order;
    }
    return lhs.definition->name < rhs.definition->name;
  };
  const auto keep = std::min<std::size_t>(3, ranked.size());
  std::partial_sort(ranked.begin(),
                    ranked.begin() + static_cast<std::ptrdiff_t>(keep),
                    ranked.end(), ranked_less);
  ClusterHaplogroupAssignment assignment;
  assignment.callable_ranges_known = input.callable.known;
  assignment.callable_ranges = callable_ranges;
  assignment.observed_markers.reserve(observed.size());
  for (const auto &mutation : observed) {
    assignment.observed_markers.push_back(mutation.encoded);
  }
  for (std::size_t i = 0; i < keep; ++i) {
    HaplogroupCandidate candidate;
    candidate.name = ranked[i].definition->name;
    candidate.score = ranked[i].score;
    for_each_callable_mutation(
        ranked[i].definition->mutations, callable_positions,
        [&](const auto &mutation) {
          if (contains(mutation)) {
            candidate.matched.push_back(mutation.encoded);
          } else {
            candidate.missing.push_back(mutation.encoded);
          }
        });
    for (const auto &mutation : observed) {
      if (!definition_contains_observed(ranked[i].definition->mutations,
                                        mutation)) {
        candidate.extra.push_back(mutation.encoded);
      }
    }
    assignment.candidates.push_back(std::move(candidate));
  }
  if (!assignment.candidates.empty()) {
    assignment.best = assignment.candidates.front().name;
    assignment.quality = assignment.candidates.front().score;
  }
  if (assignment.candidates.size() > 1 &&
      assignment.candidates[1].score > 0.0 &&
      assignment.candidates.front().score - assignment.candidates[1].score <
          3.0 &&
      macro_haplogroup(assignment.candidates.front().name) !=
          macro_haplogroup(assignment.candidates[1].name)) {
    assignment.contamination_warning = true;
  }
  return assignment;
}

} // namespace

[[nodiscard]] std::map<int, ClusterHaplogroupAssignment>
assign_haplogroups(const std::vector<ReadFeature> &features,
                   const std::vector<HaplogroupDefinition> &definitions,
                   const std::unordered_map<std::string, double> &weights,
                   const std::vector<PhyloAlignmentRule> &alignment_rules,
                   const AnalysisConfig &config) {
  std::map<int, std::vector<const ReadFeature *>> clusters;
  for (const auto &feature : features) {
    if (!feature.filtered_numt) {
      clusters[feature.cluster_id].push_back(&feature);
    }
  }
  std::vector<std::pair<int, const std::vector<const ReadFeature *> *>>
      cluster_list;
  cluster_list.reserve(clusters.size());
  for (const auto &[cluster_id, members] : clusters) {
    cluster_list.emplace_back(cluster_id, &members);
  }

  std::vector<std::vector<std::size_t>> definitions_by_position(
      static_cast<std::size_t>(kDefaultReferenceLength) + 1U);
  for (std::size_t definition_index = 0; definition_index < definitions.size();
       ++definition_index) {
    for (const auto &mutation : definitions[definition_index].mutations) {
      if (mutation.locus.position < definitions_by_position.size()) {
        definitions_by_position[mutation.locus.position].push_back(
            definition_index);
      }
    }
  }

  // Clusters are independent: prepare them in parallel, score each distinct
  // marker/range signature once (in parallel), then assemble in cluster
  // order. Identical signatures receive identical assignments, as with the
  // former sequential cache.
  std::vector<ClusterHaplogroupInput> inputs(cluster_list.size());
  parallel_for_index(cluster_list.size(), config, 8U,
                     [&](const std::size_t index, std::size_t) {
                       inputs[index] = prepare_cluster_haplogroup_input(
                           cluster_list[index].first,
                           *cluster_list[index].second, alignment_rules,
                           config);
                     });

  std::unordered_map<std::string_view, std::size_t> signature_jobs;
  std::vector<std::size_t> job_inputs;
  std::vector<std::size_t> input_jobs(inputs.size(), 0U);
  for (std::size_t index = 0; index < inputs.size(); ++index) {
    if (inputs[index].unavailable) {
      continue;
    }
    const auto [job, inserted] =
        signature_jobs.try_emplace(inputs[index].signature, job_inputs.size());
    if (inserted) {
      job_inputs.push_back(index);
    }
    input_jobs[index] = job->second;
  }
  const CompactDefinitions compact(definitions);
  std::vector<ClusterHaplogroupAssignment> job_results(job_inputs.size());
  parallel_for_index(job_inputs.size(), config, 1U,
                     [&](const std::size_t job, std::size_t) {
                       job_results[job] = score_cluster_haplogroups(
                           inputs[job_inputs[job]], definitions, compact,
                           definitions_by_position, weights);
                     });

  std::map<int, ClusterHaplogroupAssignment> assignments;
  for (std::size_t index = 0; index < inputs.size(); ++index) {
    if (inputs[index].unavailable) {
      ClusterHaplogroupAssignment assignment;
      assignment.callable_ranges_known = false;
      assignments.emplace(inputs[index].cluster_id, std::move(assignment));
      continue;
    }
    assignments.emplace(inputs[index].cluster_id,
                        job_results[input_jobs[index]]);
  }
  return assignments;
}

} // namespace mito::detail
