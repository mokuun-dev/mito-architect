#include "detail/pipeline.hpp"

namespace mito::detail {

void append_sparse_observation(SparseEvidenceStore &store,
                               EvidenceObservation observation,
                               const AnalysisConfig &config) {
  if (store.observations.size() >= config.max_evidence_observations) {
    throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                        "schema 0.6 sparse observation limit exceeded "
                        "(max_evidence_observations=" +
                            std::to_string(config.max_evidence_observations) +
                            ")");
  }
  observation.id = store.observations.size();
  store.observations.push_back(std::move(observation));
}

[[nodiscard]] std::unordered_map<std::string, std::size_t>
molecule_index_by_id(const MoleculeAssemblyResult &assembly) {
  std::unordered_map<std::string, std::size_t> indices;
  indices.reserve(assembly.molecules.size());
  for (std::size_t index = 0U; index < assembly.molecules.size(); ++index) {
    if (!indices.emplace(assembly.molecules[index].id, index).second) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "schema 0.6 molecule ID is not unique");
    }
  }
  return indices;
}

[[nodiscard]] std::size_t
circular_previous_position(const std::size_t position,
                           const std::size_t reference_length) {
  return position > 1U ? position - 1U : reference_length;
}

[[nodiscard]] std::size_t
circular_next_position(const std::size_t position,
                       const std::size_t reference_length) {
  return position < reference_length ? position + 1U : 1U;
}

void append_small_indel_evidence(SparseEvidenceStore &store,
                                 const MoleculeAssemblyResult &assembly,
                                 const std::vector<ReadFeature> &features,
                                 const CallabilityResult &callability,
                                 const std::string &reference,
                                 const AnalysisConfig &config) {
  struct Support {
    std::size_t molecule_index = 0U;
    AlignmentFragmentId alignment_fragment_id;
    std::uint8_t mapping_quality = 0U;
    char strand = '+';
    NormalizedSmallIndel call;
  };
  struct Aggregate {
    NormalizedSmallIndel representative;
    std::map<std::size_t, Support> supports;
  };
  std::map<std::string, Aggregate> aggregates;
  AnalysisConfig fragment_config = config;
  fragment_config.excluded_snp_flags =
      static_cast<std::uint16_t>(fragment_config.excluded_snp_flags & ~0x800U);
  for (std::size_t molecule_index = 0U;
       molecule_index < assembly.molecules.size(); ++molecule_index) {
    throw_if_cancelled(config);
    if (molecule_index >= features.size() ||
        features[molecule_index].filtered_numt ||
        !assembly.molecules[molecule_index].analysis_eligible) {
      continue;
    }
    for (const auto fragment_id :
         assembly.molecules[molecule_index].fragment_ids) {
      const auto &fragment = assembly.fragments[fragment_id.value];
      const auto &read = source_record_for_fragment(assembly, fragment);
      auto calls =
          call_small_indels_from_alignment(read, fragment_config, reference);
      for (auto &call : calls) {
        auto [it, inserted] = aggregates.try_emplace(call.id);
        if (inserted) {
          it->second.representative = call;
        } else if (it->second.representative.type != call.type ||
                   it->second.representative.start != call.start ||
                   it->second.representative.end != call.end ||
                   it->second.representative.reference != call.reference ||
                   it->second.representative.alternate != call.alternate) {
          throw AnalysisError(AnalysisErrorCode::internal_error,
                              "normalized small-indel ID collision for '" +
                                  call.id + "'");
        }
        Support candidate{molecule_index, fragment_id, read.mapping_quality,
                          (read.flags & 0x10U) == 0U ? '+' : '-',
                          std::move(call)};
        const auto support = it->second.supports.find(molecule_index);
        if (support == it->second.supports.end() ||
            candidate.call.base_quality.value_or(0U) >
                support->second.call.base_quality.value_or(0U) ||
            (candidate.call.base_quality == support->second.call.base_quality &&
             candidate.alignment_fragment_id <
                 support->second.alignment_fragment_id)) {
          it->second.supports.insert_or_assign(molecule_index,
                                               std::move(candidate));
        }
      }
    }
  }

  using IndelLocus = std::pair<std::string, std::size_t>;
  std::map<IndelLocus,
           std::vector<std::pair<const std::string *, const Aggregate *>>>
      aggregates_by_locus;
  for (const auto &[event_id, aggregate] : aggregates) {
    aggregates_by_locus[{aggregate.representative.type,
                         aggregate.representative.start}]
        .emplace_back(&event_id, &aggregate);
  }

  for (auto &[event_id, aggregate] : aggregates) {
    const auto event_index = store.events.size();
    EvidenceEvent event;
    event.id = aggregate.representative.id;
    event.type = aggregate.representative.type;
    event.start = aggregate.representative.start;
    event.end = aggregate.representative.end;
    event.length = aggregate.representative.length;
    event.reference = aggregate.representative.reference;
    event.alternate = aggregate.representative.alternate;
    event.normalization = aggregate.representative.normalization;
    event.source_projection = "small_indels";
    event.absence_assessable =
        event.type == "SMALL_INSERTION" || event.end >= event.start;
    event.negative_evidence_rule =
        event.type == "SMALL_INSERTION"
            ? "same_fragment_callable_reference_adjacency"
            : "same_fragment_callable_deleted_span_with_flanks";
    for (const auto &[molecule_index, support] : aggregate.supports) {
      (void)support;
      event.supporting_molecules.push_back(
          assembly.molecules[molecule_index].id);
    }
    std::sort(event.supporting_molecules.begin(),
              event.supporting_molecules.end());
    event.supporting_molecules.erase(
        std::unique(event.supporting_molecules.begin(),
                    event.supporting_molecules.end()),
        event.supporting_molecules.end());
    store.events.push_back(event);

    std::vector<const Support *> support_by_molecule(assembly.molecules.size(),
                                                     nullptr);
    for (const auto &[molecule_index, support] : aggregate.supports) {
      support_by_molecule[molecule_index] = &support;
    }
    for (std::size_t molecule_index = 0U;
         molecule_index < assembly.molecules.size(); ++molecule_index) {
      if (molecule_index >= features.size() ||
          features[molecule_index].filtered_numt ||
          !assembly.molecules[molecule_index].analysis_eligible) {
        continue;
      }
      const auto &molecule = assembly.molecules[molecule_index];
      EvidenceObservation observation;
      observation.molecule_index = molecule.index;
      observation.event_index = event_index;
      observation.alignment_fragment_id = molecule.representative_fragment_id;
      const auto &read = assembly.representatives[molecule_index];
      observation.mapping_quality = read.mapping_quality;
      observation.strand = (read.flags & 0x10U) == 0U ? '+' : '-';
      if (const auto *support = support_by_molecule[molecule_index]) {
        observation.alignment_fragment_id = support->alignment_fragment_id;
        observation.mapping_quality = support->mapping_quality;
        observation.strand = support->strand;
        observation.state = ObservationState::alternate;
        observation.observed_allele = support->call.alternate;
        observation.base_quality = support->call.base_quality;
        observation.evidence_source = "cigar_small_indel";
        append_sparse_observation(store, std::move(observation), config);
        continue;
      }
      const Support *alternative_support = nullptr;
      const auto locus = aggregates_by_locus.find({event.type, event.start});
      if (locus == aggregates_by_locus.end()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "small-indel locus index is unresolved");
      }
      for (const auto &[alternative_id, alternative] : locus->second) {
        if (*alternative_id == event_id) {
          continue;
        }
        const auto found = alternative->supports.find(molecule_index);
        if (found != alternative->supports.end()) {
          alternative_support = &found->second;
          break;
        }
      }
      if (alternative_support != nullptr) {
        observation.alignment_fragment_id =
            alternative_support->alignment_fragment_id;
        observation.mapping_quality = alternative_support->mapping_quality;
        observation.strand = alternative_support->strand;
        observation.state = ObservationState::event_absent;
        observation.observed_allele = alternative_support->call.alternate;
        observation.base_quality = alternative_support->call.base_quality;
        observation.evidence_source = "cigar_alternative_small_indel";
        append_sparse_observation(store, std::move(observation), config);
        continue;
      }
      bool absent_is_callable = false;
      if (event.type == "SMALL_INSERTION") {
        absent_is_callable = reference_adjacency_callable(
            callability.molecules[molecule_index], event.start,
            circular_next_position(event.start, reference.size()));
      } else if (event.end >= event.start) {
        absent_is_callable = reference_span_with_flanks_callable(
            callability.molecules[molecule_index], event.start, event.end,
            reference.size());
      }
      if (absent_is_callable) {
        observation.state = ObservationState::event_absent;
        observation.observed_allele = event.reference;
        observation.evidence_source = "callable_reference_path";
        append_sparse_observation(store, std::move(observation), config);
      }
    }
  }
}

[[nodiscard]] std::string uppercase_ascii(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](const unsigned char character) {
                   return static_cast<char>(std::toupper(character));
                 });
  return value;
}

void append_structural_evidence(
    SparseEvidenceStore &store, const MoleculeAssemblyResult &assembly,
    const std::vector<ReadFeature> &features,
    const CallabilityResult &callability,
    const std::map<std::string, SvCall> &svs,
    const std::map<std::string, ComplexSvCall> &complex_events,
    const std::size_t reference_length, const AnalysisConfig &config) {
  const auto molecule_indices = molecule_index_by_id(assembly);
  for (const auto &[_, sv] : svs) {
    const auto event_index = store.events.size();
    EvidenceEvent event;
    event.id = "sv:" + sv.id;
    event.type = "SV_" + uppercase_ascii(sv.type);
    event.start = sv.start;
    event.end = sv.end;
    event.length = sv.length;
    event.normalization = "canonical_mtdna_adjacency_v1";
    event.source_projection = "svs";
    event.absence_assessable =
        sv.type == "insertion" || (sv.type == "deletion" && sv.end >= sv.start);
    event.negative_evidence_rule =
        sv.type == "insertion"
            ? "same_fragment_callable_reference_adjacency"
            : (sv.type == "deletion" && sv.end >= sv.start
                   ? "same_fragment_callable_deleted_span_with_flanks"
                   : "support_only_no_negative_inference");
    for (const auto &molecule_id : sv.supporting_reads) {
      const auto found = molecule_indices.find(molecule_id);
      if (found == molecule_indices.end()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "SV support references an unknown molecule '" +
                                molecule_id + "'");
      }
      const auto molecule_index = found->second;
      if (molecule_index < features.size() &&
          !features[molecule_index].filtered_numt &&
          assembly.molecules[molecule_index].analysis_eligible) {
        event.supporting_molecules.push_back(molecule_id);
      }
    }
    std::sort(event.supporting_molecules.begin(),
              event.supporting_molecules.end());
    event.supporting_molecules.erase(
        std::unique(event.supporting_molecules.begin(),
                    event.supporting_molecules.end()),
        event.supporting_molecules.end());
    if (event.supporting_molecules.empty()) {
      continue;
    }
    store.events.push_back(event);

    std::vector<bool> supports(assembly.molecules.size(), false);
    for (const auto &molecule_id : event.supporting_molecules) {
      const auto found = molecule_indices.find(molecule_id);
      if (found == molecule_indices.end()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "SV support references an unknown molecule '" +
                                molecule_id + "'");
      }
      supports[found->second] = true;
    }
    for (std::size_t molecule_index = 0U;
         molecule_index < assembly.molecules.size(); ++molecule_index) {
      if (molecule_index >= features.size() ||
          features[molecule_index].filtered_numt ||
          !assembly.molecules[molecule_index].analysis_eligible) {
        continue;
      }
      const auto &molecule = assembly.molecules[molecule_index];
      EvidenceObservation observation;
      observation.molecule_index = molecule.index;
      observation.event_index = event_index;
      observation.alignment_fragment_id = molecule.representative_fragment_id;
      const auto &read = assembly.representatives[molecule_index];
      observation.mapping_quality = read.mapping_quality;
      observation.strand = (read.flags & 0x10U) == 0U ? '+' : '-';
      if (supports[molecule_index]) {
        observation.state = ObservationState::alternate;
        observation.observed_allele = "<" + uppercase_ascii(sv.type) + ">";
        observation.evidence_source = "normalized_sv_support";
        append_sparse_observation(store, std::move(observation), config);
        continue;
      }

      bool absent_is_callable = false;
      if (sv.type == "insertion") {
        absent_is_callable = reference_adjacency_callable(
            callability.molecules[molecule_index], sv.start,
            circular_next_position(sv.start, reference_length));
      } else if (sv.type == "deletion" && sv.end >= sv.start) {
        absent_is_callable = reference_span_with_flanks_callable(
            callability.molecules[molecule_index], sv.start, sv.end,
            reference_length);
      }
      if (absent_is_callable) {
        observation.state = ObservationState::event_absent;
        observation.observed_allele = "REFERENCE_ADJACENCY";
        observation.evidence_source = "callable_reference_path";
        append_sparse_observation(store, std::move(observation), config);
      }
    }
  }

  for (const auto &[_, complex] : complex_events) {
    const auto event_index = store.events.size();
    EvidenceEvent event;
    event.id = complex.id;
    event.type = "COMPLEX_SV_PATH";
    event.normalization = "strand_invariant_ordered_path_v1";
    event.source_projection = "complex_events";
    event.negative_evidence_rule = "support_only_no_negative_inference";
    event.component_event_ids.reserve(complex.junction_ids.size());
    for (const auto &junction_id : complex.junction_ids) {
      event.component_event_ids.push_back("sv:" + junction_id);
    }
    for (const auto &molecule_id : complex.supporting_reads) {
      const auto found = molecule_indices.find(molecule_id);
      if (found == molecule_indices.end()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "complex-event support references an unknown molecule '" +
                molecule_id + "'");
      }
      const auto molecule_index = found->second;
      if (molecule_index < features.size() &&
          !features[molecule_index].filtered_numt &&
          assembly.molecules[molecule_index].analysis_eligible) {
        event.supporting_molecules.push_back(molecule_id);
      }
    }
    std::sort(event.supporting_molecules.begin(),
              event.supporting_molecules.end());
    event.supporting_molecules.erase(
        std::unique(event.supporting_molecules.begin(),
                    event.supporting_molecules.end()),
        event.supporting_molecules.end());
    if (event.supporting_molecules.empty()) {
      continue;
    }
    store.events.push_back(event);
    for (const auto &molecule_id : event.supporting_molecules) {
      const auto found = molecule_indices.find(molecule_id);
      if (found == molecule_indices.end()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "complex-event support references an unknown molecule '" +
                molecule_id + "'");
      }
      const auto molecule_index = found->second;
      const auto &molecule = assembly.molecules[molecule_index];
      EvidenceObservation observation;
      observation.molecule_index = molecule.index;
      observation.event_index = event_index;
      observation.alignment_fragment_id = molecule.representative_fragment_id;
      observation.state = ObservationState::alternate;
      observation.observed_allele = "COMPLEX_PATH";
      observation.evidence_source = "split_alignment_path";
      const auto &read = assembly.representatives[molecule_index];
      observation.mapping_quality = read.mapping_quality;
      observation.strand = (read.flags & 0x10U) == 0U ? '+' : '-';
      append_sparse_observation(store, std::move(observation), config);
    }
  }
}

void validate_unified_evidence_graph(const MoleculeAssemblyResult &assembly,
                                     const CallabilityResult &callability,
                                     const SparseEvidenceStore &store,
                                     const std::size_t reference_length) {
  if (callability.molecules.size() != assembly.molecules.size()) {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "schema 0.6 callability/molecule cardinality mismatch");
  }
  for (const auto &molecule : assembly.molecules) {
    if (molecule.source_qnames.empty() || molecule.identity_policy.empty() ||
        (molecule.analysis_eligible && !molecule.exclusion_reasons.empty()) ||
        (!molecule.analysis_eligible && molecule.exclusion_reasons.empty())) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "schema 0.6 molecule assembly contract is invalid");
    }
  }
  for (std::size_t index = 0U; index < callability.molecules.size(); ++index) {
    const auto &summary = callability.molecules[index];
    if (summary.molecule_index.value != index) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 callability molecule reference is unresolved");
    }
    std::size_t callable_bases = 0U;
    std::size_t previous_end = 0U;
    for (const auto &[start, end] : summary.ranges) {
      if (start == 0U || end < start || end > reference_length ||
          (previous_end != 0U && start <= previous_end)) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 callable ranges are invalid or overlapping");
      }
      callable_bases += end - start + 1U;
      previous_end = end;
    }
    if (callable_bases != summary.callable_bases ||
        (!summary.known && summary.callable_bases != 0U)) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "schema 0.6 callable base count invariant failed");
    }
  }

  std::set<std::string> event_ids;
  for (const auto &event : store.events) {
    if (event.id.empty() || !event_ids.insert(event.id).second) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 normalized event IDs are empty or duplicated");
    }
  }
  std::vector<std::vector<std::string>> projected_support(store.events.size());
  std::set<std::pair<std::size_t, std::size_t>> molecule_event_pairs;
  for (std::size_t index = 0U; index < store.observations.size(); ++index) {
    const auto &observation = store.observations[index];
    if (observation.id != index ||
        observation.molecule_index.value >= assembly.molecules.size() ||
        observation.event_index >= store.events.size() ||
        observation.alignment_fragment_id.value >= assembly.fragments.size()) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 unified observation reference is unresolved");
    }
    if (assembly.fragments[observation.alignment_fragment_id.value]
            .molecule_index != observation.molecule_index) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "schema 0.6 unified observation alignment belongs to "
                          "another molecule");
    }
    if (observation.state == ObservationState::not_callable) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "schema 0.6 sparse store materialized NOT_CALLABLE");
    }
    if (observation.evidence_source.empty()) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "schema 0.6 observation evidence source is empty");
    }
    if (!molecule_event_pairs
             .emplace(observation.molecule_index.value, observation.event_index)
             .second) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 unified molecule/event observation is not unique");
    }
    if (observation.state == ObservationState::alternate) {
      projected_support[observation.event_index].push_back(
          assembly.molecules[observation.molecule_index.value].id);
    }
  }
  for (std::size_t event_index = 0U; event_index < store.events.size();
       ++event_index) {
    const auto &event = store.events[event_index];
    if (event.negative_evidence_rule.empty() ||
        (event.absence_assessable == (event.negative_evidence_rule ==
                                      "support_only_no_negative_inference"))) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 event negative-evidence contract is invalid for " +
              event.id);
    }
    auto &support = projected_support[event_index];
    std::sort(support.begin(), support.end());
    support.erase(std::unique(support.begin(), support.end()), support.end());
    if (support != event.supporting_molecules) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 event/support projection invariant failed for " +
              event.id);
    }
  }
}

} // namespace mito::detail
