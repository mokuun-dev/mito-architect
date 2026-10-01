#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] bool
is_architecture_callable_state(const ObservationState state) noexcept {
  return state == ObservationState::reference ||
         state == ObservationState::alternate ||
         state == ObservationState::event_absent;
}

[[nodiscard]] double
architecture_event_weight(const EvidenceEvent &event) noexcept {
  if (event.type == "COMPLEX_SV_PATH") {
    return 3.0;
  }
  if (event.type.starts_with("SV_")) {
    return 2.0;
  }
  if (event.type == "SMALL_INSERTION" || event.type == "SMALL_DELETION") {
    return 1.5;
  }
  return 1.0;
}

[[nodiscard]] std::vector<MoleculeArchitectureEvidence>
build_molecule_architecture_evidence(const SparseEvidenceStore &store,
                                     const MoleculeAssemblyResult &assembly,
                                     const CallabilityResult &callability,
                                     const std::vector<ReadFeature> &features,
                                     const std::size_t reference_length) {
  std::vector<MoleculeArchitectureEvidence> evidence(assembly.molecules.size());
  for (std::size_t index = 0U; index < evidence.size(); ++index) {
    evidence[index].eligible = assembly.molecules[index].analysis_eligible &&
                               index < features.size() &&
                               !features[index].filtered_numt &&
                               index < callability.molecules.size() &&
                               callability.molecules[index].known;
    if (evidence[index].eligible && reference_length != 0U) {
      evidence[index].callable_fraction =
          static_cast<double>(callability.molecules[index].callable_bases) /
          static_cast<double>(reference_length);
    }
  }
  for (const auto &observation : store.observations) {
    const auto molecule_index = observation.molecule_index.value;
    if (molecule_index >= evidence.size() ||
        observation.event_index >= store.events.size()) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "architecture inference received an unresolved observation");
    }
    if (!evidence[molecule_index].eligible) {
      continue;
    }
    evidence[molecule_index].observations.emplace_back(observation.event_index,
                                                       observation.state);
    if (observation.state == ObservationState::alternate) {
      evidence[molecule_index].alternate_signature.push_back(
          observation.event_index);
    }
  }
  for (std::size_t index = 0U; index < evidence.size(); ++index) {
    auto &molecule = evidence[index];
    std::sort(
        molecule.observations.begin(), molecule.observations.end(),
        [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });
    if (std::adjacent_find(molecule.observations.begin(),
                           molecule.observations.end(),
                           [](const auto &lhs, const auto &rhs) {
                             return lhs.first == rhs.first;
                           }) != molecule.observations.end()) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "architecture inference found duplicate molecule/event evidence");
    }
    std::sort(molecule.alternate_signature.begin(),
              molecule.alternate_signature.end());
    // A directly observed deletion can be informative even when every aligned
    // reference base fails the base-quality threshold.  Zero callable bases
    // without positive event evidence, however, cannot support an assignment.
    if (molecule.eligible &&
        callability.molecules[index].callable_bases == 0U &&
        molecule.alternate_signature.empty()) {
      molecule.eligible = false;
    }
  }
  return evidence;
}

[[nodiscard]] std::vector<ArchitectureEventProfile>
build_architecture_event_profile(
    const std::vector<std::size_t> &members,
    const std::vector<MoleculeArchitectureEvidence> &molecule_evidence,
    const AnalysisConfig &config, const std::size_t minimum_callable_support) {
  struct Counts {
    std::size_t alternate = 0U;
    std::size_t absent = 0U;
    std::size_t uncertain = 0U;
  };
  std::map<std::size_t, Counts> counts;
  for (const auto molecule_index : members) {
    if (molecule_index >= molecule_evidence.size()) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "architecture member index is out of range");
    }
    for (const auto &[event_index, state] :
         molecule_evidence[molecule_index].observations) {
      auto &entry = counts[event_index];
      if (state == ObservationState::alternate) {
        ++entry.alternate;
      } else if (state == ObservationState::reference ||
                 state == ObservationState::event_absent) {
        ++entry.absent;
      } else {
        ++entry.uncertain;
      }
    }
  }

  std::vector<ArchitectureEventProfile> profile;
  profile.reserve(counts.size());
  for (const auto &[event_index, count] : counts) {
    const auto callable = count.alternate + count.absent;
    if (callable < minimum_callable_support) {
      continue;
    }
    ArchitectureEventProfile entry;
    entry.event_index = event_index;
    entry.alternate = count.alternate;
    entry.absent = count.absent;
    entry.uncertain = count.uncertain;
    entry.not_callable = members.size() - callable - count.uncertain;
    entry.alternate_fraction =
        static_cast<double>(count.alternate) / static_cast<double>(callable);
    entry.expected_alternate =
        entry.alternate_fraction >= config.architecture_consensus_fraction;
    entry.expected_absent = entry.alternate_fraction <=
                            (1.0 - config.architecture_consensus_fraction);
    profile.push_back(entry);
  }
  return profile;
}

void update_architecture_signatures(CandidateArchitecture &architecture,
                                    const AnalysisConfig &config) {
  architecture.defining_event_indices.clear();
  architecture.optional_event_indices.clear();
  for (const auto &event : architecture.event_profile) {
    if (event.expected_alternate) {
      architecture.defining_event_indices.push_back(event.event_index);
    } else if (event.alternate != 0U &&
               event.alternate_fraction >=
                   config.architecture_optional_fraction) {
      architecture.optional_event_indices.push_back(event.event_index);
    }
  }
}

[[nodiscard]] std::string
architecture_id_for_signature(const std::vector<std::size_t> &signature,
                              const SparseEvidenceStore &store) {
  if (signature.empty()) {
    return "architecture:v1:reference-like";
  }
  std::string id = "architecture:v1:";
  for (const auto event_index : signature) {
    if (event_index >= store.events.size()) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "architecture signature event is unresolved");
    }
    const auto &event_id = store.events[event_index].id;
    id += std::to_string(event_id.size());
    id.push_back('#');
    id += event_id;
  }
  return id;
}

[[nodiscard]] ArchitectureScore
score_architecture_candidate(const std::size_t architecture_index,
                             const CandidateArchitecture &architecture,
                             const MoleculeArchitectureEvidence &molecule,
                             const SparseEvidenceStore &store) {
  ArchitectureScore result;
  result.architecture_index = architecture_index;
  if (architecture.event_profile.empty()) {
    if (architecture.seed_signature.empty() &&
        molecule.alternate_signature.empty()) {
      result.distance = 0.0;
      result.score = 0.0;
      result.overlap_fraction = 1.0;
    }
    return result;
  }

  double candidate_weight = 0.0;
  double jointly_callable_weight = 0.0;
  double compared_weight = 0.0;
  double disagreement_weight = 0.0;
  for (const auto &profile : architecture.event_profile) {
    // Mixed, low-frequency events describe within-architecture noise.  They
    // remain in the report, but must not turn a shared consensus signature
    // into a rejection solely because each long read has different noise.
    if (!profile.expected_alternate && !profile.expected_absent) {
      continue;
    }
    if (profile.event_index >= store.events.size()) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "architecture profile event is unresolved");
    }
    const auto weight =
        architecture_event_weight(store.events[profile.event_index]);
    candidate_weight += weight;
    const auto found = std::lower_bound(
        molecule.observations.begin(), molecule.observations.end(),
        profile.event_index,
        [](const auto &observation, const std::size_t event_index) {
          return observation.first < event_index;
        });
    if (found == molecule.observations.end() ||
        found->first != profile.event_index ||
        !is_architecture_callable_state(found->second)) {
      continue;
    }
    jointly_callable_weight += weight;
    compared_weight += weight;
    const bool alternate = found->second == ObservationState::alternate;
    disagreement_weight +=
        weight * (alternate ? 1.0 - profile.alternate_fraction
                            : profile.alternate_fraction);
  }

  for (const auto event_index : molecule.alternate_signature) {
    const auto profile = std::lower_bound(
        architecture.event_profile.begin(), architecture.event_profile.end(),
        event_index, [](const auto &entry, const std::size_t target) {
          return entry.event_index < target;
        });
    if (profile != architecture.event_profile.end() &&
        profile->event_index == event_index) {
      continue;
    }
    if (event_index >= store.events.size()) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "molecule architecture event is unresolved");
    }
    const auto weight = architecture_event_weight(store.events[event_index]);
    compared_weight += weight;
    disagreement_weight += weight;
  }

  if (compared_weight == 0.0 || candidate_weight == 0.0) {
    return result;
  }
  result.distance = disagreement_weight / compared_weight;
  result.overlap_fraction =
      std::min(1.0, jointly_callable_weight / candidate_weight);
  constexpr double missing_profile_penalty = 0.25;
  result.score = result.distance +
                 missing_profile_penalty * (1.0 - result.overlap_fraction);
  return result;
}

[[nodiscard]] std::uint64_t mix_architecture_seed(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

[[nodiscard]] double
weighted_signature_jaccard(const std::vector<std::size_t> &lhs,
                           const std::vector<std::size_t> &rhs,
                           const SparseEvidenceStore &store) {
  std::size_t left = 0U;
  std::size_t right = 0U;
  double intersection = 0.0;
  double union_weight = 0.0;
  while (left < lhs.size() || right < rhs.size()) {
    if (right == rhs.size() || (left < lhs.size() && lhs[left] < rhs[right])) {
      union_weight += architecture_event_weight(store.events[lhs[left]]);
      ++left;
    } else if (left == lhs.size() || rhs[right] < lhs[left]) {
      union_weight += architecture_event_weight(store.events[rhs[right]]);
      ++right;
    } else {
      const auto weight = architecture_event_weight(store.events[lhs[left]]);
      intersection += weight;
      union_weight += weight;
      ++left;
      ++right;
    }
  }
  return union_weight == 0.0 ? 1.0 : intersection / union_weight;
}

void finalize_architecture_statistics(
    ArchitectureInferenceResult &result,
    const std::vector<MoleculeArchitectureEvidence> &molecule_evidence,
    const MoleculeAssemblyResult &assembly, const SparseEvidenceStore &store,
    const CallabilityResult &callability, const std::size_t reference_length,
    const AnalysisConfig &config) {
  (void)callability;
  (void)store;
  constexpr std::size_t representative_limit = 5U;
  for (std::size_t architecture_index = 0U;
       architecture_index < result.architectures.size(); ++architecture_index) {
    auto &architecture = result.architectures[architecture_index];
    architecture.event_profile = build_architecture_event_profile(
        architecture.member_molecule_indices, molecule_evidence, config,
        config.min_architecture_molecules);
    update_architecture_signatures(architecture, config);
    const auto member_count = architecture.member_molecule_indices.size();
    architecture.estimated_fraction =
        result.eligible_molecules == 0U
            ? 0.0
            : static_cast<double>(member_count) /
                  static_cast<double>(result.eligible_molecules);
    std::tie(architecture.ci95_low, architecture.ci95_high) =
        wilson_interval(member_count, result.eligible_molecules);

    std::vector<double> callable_fractions;
    callable_fractions.reserve(member_count);
    double confidence_sum = 0.0;
    architecture.min_assignment_confidence = 1.0;
    for (const auto molecule_index : architecture.member_molecule_indices) {
      callable_fractions.push_back(
          molecule_evidence[molecule_index].callable_fraction);
      confidence_sum += result.assignments[molecule_index].confidence;
      architecture.min_assignment_confidence =
          std::min(architecture.min_assignment_confidence,
                   result.assignments[molecule_index].confidence);
    }
    std::sort(callable_fractions.begin(), callable_fractions.end());
    if (!callable_fractions.empty()) {
      const auto middle = callable_fractions.size() / 2U;
      architecture.median_callable_fraction =
          callable_fractions.size() % 2U == 0U
              ? (callable_fractions[middle - 1U] + callable_fractions[middle]) /
                    2.0
              : callable_fractions[middle];
      architecture.mean_assignment_confidence =
          confidence_sum / static_cast<double>(member_count);
    } else {
      architecture.min_assignment_confidence = 0.0;
    }

    architecture.representative_molecule_indices =
        architecture.member_molecule_indices;
    std::sort(architecture.representative_molecule_indices.begin(),
              architecture.representative_molecule_indices.end(),
              [&](const std::size_t lhs, const std::size_t rhs) {
                const auto lhs_confidence = result.assignments[lhs].confidence;
                const auto rhs_confidence = result.assignments[rhs].confidence;
                if (lhs_confidence != rhs_confidence) {
                  return lhs_confidence > rhs_confidence;
                }
                return assembly.molecules[lhs].id < assembly.molecules[rhs].id;
              });
    if (architecture.representative_molecule_indices.size() >
        representative_limit) {
      architecture.representative_molecule_indices.resize(representative_limit);
    }

    architecture.signature_stability = 0.0;
    architecture.recovery_rate = 0.0;
    architecture.abundance_standard_deviation = 0.0;
    (void)reference_length;
  }
}

[[nodiscard]] ArchitectureInferenceResult infer_architectures(
    const SparseEvidenceStore &store, const MoleculeAssemblyResult &assembly,
    const CallabilityResult &callability,
    const std::vector<ReadFeature> &features,
    const std::size_t reference_length, const AnalysisConfig &config) {
  ArchitectureInferenceResult result;
  result.assignments.resize(assembly.molecules.size());
  const auto molecule_evidence = build_molecule_architecture_evidence(
      store, assembly, callability, features, reference_length);
  std::map<std::vector<std::size_t>, std::vector<std::size_t>> seed_groups;
  for (std::size_t molecule_index = 0U;
       molecule_index < molecule_evidence.size(); ++molecule_index) {
    if (!molecule_evidence[molecule_index].eligible) {
      continue;
    }
    ++result.eligible_molecules;
    seed_groups[molecule_evidence[molecule_index].alternate_signature]
        .push_back(molecule_index);
  }
  for (auto &[signature, members] : seed_groups) {
    if (members.size() < config.min_architecture_molecules) {
      continue;
    }
    if (result.architectures.size() >= config.max_candidate_architectures) {
      throw AnalysisError(
          AnalysisErrorCode::resource_exhausted,
          "candidate architecture limit exceeded "
          "(max_candidate_architectures=" +
              std::to_string(config.max_candidate_architectures) + ")");
    }
    CandidateArchitecture architecture;
    architecture.id = architecture_id_for_signature(signature, store);
    architecture.seed_signature = signature;
    architecture.member_molecule_indices = std::move(members);
    result.architectures.push_back(std::move(architecture));
  }

  // Exact signatures are a useful high-specificity seed when repeated.  With
  // noisy long reads, however, a real shared event can be obscured by a
  // different one-off event in every molecule.  Only when no exact signature
  // reached independent support, fall back to supported singleton events.
  // This preserves the stricter seeds whenever the data provide them.
  if (result.architectures.empty()) {
    std::map<std::size_t, std::vector<std::size_t>> event_support;
    for (std::size_t molecule_index = 0U;
         molecule_index < molecule_evidence.size(); ++molecule_index) {
      if (!molecule_evidence[molecule_index].eligible) {
        continue;
      }
      for (const auto event_index :
           molecule_evidence[molecule_index].alternate_signature) {
        event_support[event_index].push_back(molecule_index);
      }
    }
    for (auto &[event_index, members] : event_support) {
      if (members.size() < config.min_architecture_molecules) {
        continue;
      }
      if (result.architectures.size() >= config.max_candidate_architectures) {
        throw AnalysisError(
            AnalysisErrorCode::resource_exhausted,
            "candidate architecture limit exceeded "
            "(max_candidate_architectures=" +
                std::to_string(config.max_candidate_architectures) + ")");
      }
      CandidateArchitecture architecture;
      architecture.seed_signature = {event_index};
      architecture.id =
          architecture_id_for_signature(architecture.seed_signature, store);
      architecture.member_molecule_indices = std::move(members);
      result.architectures.push_back(std::move(architecture));
    }
  }

  if (result.eligible_molecules < config.min_architecture_molecules) {
    result.status = "INSUFFICIENT_ELIGIBLE_MOLECULES";
    result.unassigned_molecules = result.eligible_molecules;
    result.qc_flags.emplace_back("INSUFFICIENT_INDEPENDENT_MOLECULE_SUPPORT");
    for (std::size_t index = 0U; index < molecule_evidence.size(); ++index) {
      if (molecule_evidence[index].eligible) {
        result.assignments[index].status =
            ArchitectureAssignmentStatus::unassigned;
      }
    }
    return result;
  }
  if (result.architectures.empty()) {
    result.status = "NO_SUPPORTED_SEEDS";
    result.unassigned_molecules = result.eligible_molecules;
    result.qc_flags.emplace_back("NO_SIGNATURE_REACHED_MINIMUM_SUPPORT");
    for (std::size_t index = 0U; index < molecule_evidence.size(); ++index) {
      if (molecule_evidence[index].eligible) {
        result.assignments[index].status =
            ArchitectureAssignmentStatus::unassigned;
      }
    }
    return result;
  }

  constexpr std::size_t maximum_iterations = 8U;
  bool converged = false;
  for (std::size_t iteration = 0U; iteration < maximum_iterations;
       ++iteration) {
    throw_if_cancelled(config);
    for (auto &architecture : result.architectures) {
      architecture.event_profile = build_architecture_event_profile(
          architecture.member_molecule_indices, molecule_evidence, config,
          config.min_architecture_molecules);
      update_architecture_signatures(architecture, config);
      architecture.ambiguous_molecule_indices.clear();
    }

    std::vector<ArchitectureAssignment> next_assignments(
        assembly.molecules.size());
    std::vector<std::vector<std::size_t>> next_members(
        result.architectures.size());
    std::vector<std::vector<std::size_t>> next_ambiguous(
        result.architectures.size());
    for (std::size_t molecule_index = 0U;
         molecule_index < molecule_evidence.size(); ++molecule_index) {
      if (!molecule_evidence[molecule_index].eligible) {
        continue;
      }
      auto &assignment = next_assignments[molecule_index];
      assignment.status = ArchitectureAssignmentStatus::unassigned;
      std::vector<ArchitectureScore> scores;
      scores.reserve(result.architectures.size());
      for (std::size_t architecture_index = 0U;
           architecture_index < result.architectures.size();
           ++architecture_index) {
        scores.push_back(score_architecture_candidate(
            architecture_index, result.architectures[architecture_index],
            molecule_evidence[molecule_index], store));
      }
      std::sort(scores.begin(), scores.end(),
                [](const auto &lhs, const auto &rhs) {
                  if (lhs.score != rhs.score) {
                    return lhs.score < rhs.score;
                  }
                  return lhs.architecture_index < rhs.architecture_index;
                });
      const auto passes = [&](const ArchitectureScore &score) {
        return score.distance <= config.architecture_max_distance &&
               score.overlap_fraction >=
                   config.architecture_min_overlap_fraction;
      };
      if (scores.empty() || !passes(scores.front())) {
        continue;
      }
      const auto &best = scores.front();
      assignment.distance = best.distance;
      assignment.score = best.score;
      assignment.overlap_fraction = best.overlap_fraction;
      const bool ambiguous =
          scores.size() > 1U && passes(scores[1]) &&
          scores[1].score - best.score <= config.architecture_ambiguity_margin;
      if (ambiguous) {
        assignment.status = ArchitectureAssignmentStatus::ambiguous;
        for (const auto &score : scores) {
          if (!passes(score) ||
              score.score - best.score > config.architecture_ambiguity_margin) {
            break;
          }
          assignment.candidate_indices.push_back(score.architecture_index);
          next_ambiguous[score.architecture_index].push_back(molecule_index);
        }
        continue;
      }
      assignment.status = ArchitectureAssignmentStatus::assigned;
      assignment.architecture_index = best.architecture_index;
      const auto separation =
          scores.size() == 1U
              ? 1.0
              : std::clamp(
                    (scores[1].score - best.score) /
                        std::max(config.architecture_ambiguity_margin, 1.0e-12),
                    0.0, 1.0);
      assignment.confidence = std::clamp(
          (1.0 - best.distance) * best.overlap_fraction * separation, 0.0, 1.0);
      next_members[best.architecture_index].push_back(molecule_index);
    }

    std::vector<std::vector<std::size_t>> previous_memberships;
    previous_memberships.reserve(result.architectures.size());
    for (const auto &architecture : result.architectures) {
      previous_memberships.push_back(architecture.member_molecule_indices);
    }
    std::vector<CandidateArchitecture> retained;
    retained.reserve(result.architectures.size());
    std::vector<std::size_t> old_to_new(
        result.architectures.size(), std::numeric_limits<std::size_t>::max());
    for (std::size_t index = 0U; index < result.architectures.size(); ++index) {
      if (next_members[index].size() < config.min_architecture_molecules) {
        continue;
      }
      old_to_new[index] = retained.size();
      auto architecture = std::move(result.architectures[index]);
      architecture.member_molecule_indices = std::move(next_members[index]);
      architecture.ambiguous_molecule_indices =
          std::move(next_ambiguous[index]);
      retained.push_back(std::move(architecture));
    }
    if (retained.empty()) {
      result.architectures.clear();
      result.assignments.assign(assembly.molecules.size(), {});
      for (std::size_t index = 0U; index < molecule_evidence.size(); ++index) {
        if (molecule_evidence[index].eligible) {
          result.assignments[index].status =
              ArchitectureAssignmentStatus::unassigned;
        }
      }
      result.status = "NO_STABLE_ARCHITECTURES";
      result.unassigned_molecules = result.eligible_molecules;
      result.qc_flags.emplace_back(
          "ASSIGNMENTS_DID_NOT_RETAIN_MINIMUM_SUPPORT");
      return result;
    }

    for (auto &assignment : next_assignments) {
      if (assignment.architecture_index.has_value()) {
        const auto remapped = old_to_new[*assignment.architecture_index];
        if (remapped == std::numeric_limits<std::size_t>::max()) {
          assignment.status = ArchitectureAssignmentStatus::unassigned;
          assignment.architecture_index.reset();
          assignment.confidence = 0.0;
        } else {
          assignment.architecture_index = remapped;
        }
      }
      std::vector<std::size_t> remapped_candidates;
      for (const auto candidate : assignment.candidate_indices) {
        if (old_to_new[candidate] != std::numeric_limits<std::size_t>::max()) {
          remapped_candidates.push_back(old_to_new[candidate]);
        }
      }
      assignment.candidate_indices = std::move(remapped_candidates);
      if (assignment.status == ArchitectureAssignmentStatus::ambiguous &&
          assignment.candidate_indices.size() < 2U) {
        assignment.status = ArchitectureAssignmentStatus::unassigned;
      }
    }

    bool same_memberships = retained.size() == previous_memberships.size();
    if (same_memberships) {
      for (std::size_t index = 0U; index < retained.size(); ++index) {
        if (retained[index].member_molecule_indices !=
            previous_memberships[index]) {
          same_memberships = false;
          break;
        }
      }
    }
    result.architectures = std::move(retained);
    result.assignments = std::move(next_assignments);
    if (same_memberships) {
      converged = true;
      break;
    }
  }

  if (!converged) {
    result.qc_flags.emplace_back("MAXIMUM_ASSIGNMENT_ITERATIONS_REACHED");
  }
  result.assigned_molecules = 0U;
  result.ambiguous_molecules = 0U;
  result.unassigned_molecules = 0U;
  for (std::size_t index = 0U; index < result.assignments.size(); ++index) {
    if (!molecule_evidence[index].eligible) {
      result.assignments[index].status =
          ArchitectureAssignmentStatus::ineligible;
      continue;
    }
    switch (result.assignments[index].status) {
    case ArchitectureAssignmentStatus::assigned:
      ++result.assigned_molecules;
      break;
    case ArchitectureAssignmentStatus::ambiguous:
      ++result.ambiguous_molecules;
      break;
    case ArchitectureAssignmentStatus::unassigned:
      ++result.unassigned_molecules;
      break;
    case ArchitectureAssignmentStatus::ineligible:
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "eligible molecule retained an ineligible architecture status");
    }
  }
  result.status =
      converged ? "CANDIDATES_INFERRED" : "CANDIDATES_INFERRED_WITH_WARNING";
  finalize_architecture_statistics(result, molecule_evidence, assembly, store,
                                   callability, reference_length, config);
  evaluate_architecture_stability(result, store, assembly, callability,
                                  features, reference_length, config);
  return result;
}

void validate_architecture_inference(const ArchitectureInferenceResult &result,
                                     const MoleculeAssemblyResult &assembly) {
  if (result.assignments.size() != assembly.molecules.size()) {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "architecture assignments do not cover all molecules");
  }
  std::set<std::string> architecture_ids;
  std::vector<std::optional<std::size_t>> declared_owner(
      assembly.molecules.size());
  for (std::size_t architecture_index = 0U;
       architecture_index < result.architectures.size(); ++architecture_index) {
    const auto &architecture = result.architectures[architecture_index];
    if (architecture.id.empty() ||
        !architecture_ids.insert(architecture.id).second) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "candidate architecture IDs are not unique");
    }
    if (architecture.member_molecule_indices.size() < 2U ||
        !std::is_sorted(architecture.member_molecule_indices.begin(),
                        architecture.member_molecule_indices.end()) ||
        !std::is_sorted(architecture.ambiguous_molecule_indices.begin(),
                        architecture.ambiguous_molecule_indices.end())) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "candidate architecture membership is unsupported or unordered");
    }
    for (const auto molecule_index : architecture.member_molecule_indices) {
      if (molecule_index >= declared_owner.size() ||
          declared_owner[molecule_index].has_value()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "candidate architectures contain an unresolved or "
                            "duplicate member");
      }
      declared_owner[molecule_index] = architecture_index;
    }
    for (const auto molecule_index :
         architecture.representative_molecule_indices) {
      if (molecule_index >= declared_owner.size() ||
          std::find(architecture.member_molecule_indices.begin(),
                    architecture.member_molecule_indices.end(),
                    molecule_index) ==
              architecture.member_molecule_indices.end()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "architecture representative is not an assigned member");
      }
    }
  }

  std::size_t eligible = 0U;
  std::size_t assigned = 0U;
  std::size_t ambiguous = 0U;
  std::size_t unassigned = 0U;
  for (std::size_t molecule_index = 0U;
       molecule_index < result.assignments.size(); ++molecule_index) {
    const auto &assignment = result.assignments[molecule_index];
    switch (assignment.status) {
    case ArchitectureAssignmentStatus::assigned:
      ++eligible;
      ++assigned;
      if (!assignment.architecture_index.has_value() ||
          *assignment.architecture_index >= result.architectures.size() ||
          declared_owner[molecule_index] != assignment.architecture_index) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "assigned molecule is not owned by its candidate architecture");
      }
      break;
    case ArchitectureAssignmentStatus::ambiguous:
      ++eligible;
      ++ambiguous;
      if (assignment.architecture_index.has_value() ||
          assignment.candidate_indices.size() < 2U) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "ambiguous architecture assignment lacks multiple candidates");
      }
      for (const auto candidate : assignment.candidate_indices) {
        if (candidate >= result.architectures.size() ||
            !std::binary_search(result.architectures[candidate]
                                    .ambiguous_molecule_indices.begin(),
                                result.architectures[candidate]
                                    .ambiguous_molecule_indices.end(),
                                molecule_index)) {
          throw AnalysisError(
              AnalysisErrorCode::internal_error,
              "ambiguous molecule is absent from a candidate projection");
        }
      }
      break;
    case ArchitectureAssignmentStatus::unassigned:
      ++eligible;
      ++unassigned;
      if (assignment.architecture_index.has_value() ||
          declared_owner[molecule_index].has_value()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "unassigned molecule is retained by an architecture");
      }
      break;
    case ArchitectureAssignmentStatus::ineligible:
      if (assignment.architecture_index.has_value() ||
          declared_owner[molecule_index].has_value()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "ineligible molecule is retained by an architecture");
      }
      break;
    }
  }
  if (eligible != result.eligible_molecules ||
      assigned != result.assigned_molecules ||
      ambiguous != result.ambiguous_molecules ||
      unassigned != result.unassigned_molecules) {
    throw AnalysisError(
        AnalysisErrorCode::internal_error,
        "architecture assignment summary does not partition molecule states");
  }
}

} // namespace mito::detail
