#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] EvidenceAggregationResult
aggregate_snps(const MoleculeAssemblyResult &assembly,
               const std::vector<ReadFeature> &features,
               const AnalysisConfig &config, std::size_t reference_length,
               const std::string &reference,
               const std::map<std::string, SnpCall> &clinical_annotations) {
  EvidenceAggregationResult result;
  auto &aggregates = result.variants;
  const auto &reads = assembly.representatives;
  const bool capture_sparse_store = config.result_schema == ResultSchema::v0_6;
  std::unordered_set<std::size_t> target_positions;
  for (std::size_t feature_index = 0U; feature_index < features.size();
       ++feature_index) {
    const auto &feature = features[feature_index];
    if (feature.filtered_numt ||
        (capture_sparse_store &&
         (feature_index >= assembly.molecules.size() ||
          !assembly.molecules[feature_index].analysis_eligible))) {
      continue;
    }
    for (const auto &snp : feature.snps) {
      const auto key = snp_key(snp.position, snp.reference, snp.alternate);
      auto [it, inserted] = aggregates.try_emplace(key);
      if (inserted) {
        it->second.call = snp;
      }
      target_positions.insert(snp.position);
    }
  }
  AnalysisConfig fragment_config = config;
  fragment_config.excluded_snp_flags =
      static_cast<std::uint16_t>(fragment_config.excluded_snp_flags & ~0x800U);
  if (capture_sparse_store) {
    for (std::size_t molecule_index = 0U;
         molecule_index < assembly.molecules.size(); ++molecule_index) {
      if (molecule_index >= features.size() ||
          features[molecule_index].filtered_numt ||
          !assembly.molecules[molecule_index].analysis_eligible) {
        continue;
      }
      for (const auto fragment_id :
           assembly.molecules[molecule_index].fragment_ids) {
        const auto &fragment = assembly.fragments[fragment_id.value];
        const auto &read = source_record_for_fragment(assembly, fragment);
        auto fragment_snps =
            call_snps_from_alignment(read, reference, fragment_config);
        apply_clinical_annotations(fragment_snps, clinical_annotations);
        for (const auto &snp : fragment_snps) {
          const auto key = snp_key(snp.position, snp.reference, snp.alternate);
          auto [it, inserted] = aggregates.try_emplace(key);
          if (inserted) {
            it->second.call = snp;
          }
          target_positions.insert(snp.position);
        }
      }
    }
  }

  struct EventAggregateBinding {
    SnpAggregate *aggregate = nullptr;
    std::size_t event_index = 0;
  };
  std::unordered_map<std::size_t, std::vector<EventAggregateBinding>>
      aggregates_by_position;
  aggregates_by_position.reserve(target_positions.size());
  if (capture_sparse_store) {
    result.store.events.reserve(aggregates.size());
  }
  for (auto &[key, aggregate] : aggregates) {
    const std::size_t event_index =
        capture_sparse_store ? result.store.events.size() : 0U;
    if (capture_sparse_store) {
      EvidenceEvent event;
      event.id = "snv:" + key;
      event.type = "SNV";
      event.start = aggregate.call.position;
      event.end = aggregate.call.position;
      event.length = 1U;
      event.reference.assign(1U, aggregate.call.reference);
      event.alternate.assign(1U, aggregate.call.alternate);
      event.normalization = "rcrs_circular_snv_v1";
      event.source_projection = "variants";
      event.negative_evidence_rule = "callable_base_allele";
      event.absence_assessable = true;
      result.store.events.push_back(std::move(event));
    }
    aggregates_by_position[aggregate.call.position].push_back(
        {&aggregate, event_index});
  }

  std::unordered_map<std::size_t, LocusAlleleAccumulator> allele_evidence;
  allele_evidence.reserve(target_positions.size());
  for (std::size_t read_index = 0; read_index < reads.size(); ++read_index) {
    throw_if_cancelled(config);
    if (read_index >= features.size() || features[read_index].filtered_numt ||
        (capture_sparse_store &&
         !assembly.molecules[read_index].analysis_eligible)) {
      continue;
    }
    std::map<std::size_t, MoleculeAlleleObservation> molecule_alleles;
    std::vector<std::pair<const ReadRecord *, AlignmentFragmentId>>
        evidence_records;
    if (capture_sparse_store) {
      const auto &molecule = assembly.molecules[read_index];
      evidence_records.reserve(molecule.fragment_ids.size());
      for (const auto fragment_id : molecule.fragment_ids) {
        const auto &fragment = assembly.fragments[fragment_id.value];
        evidence_records.emplace_back(
            &source_record_for_fragment(assembly, fragment), fragment.id);
      }
    } else {
      evidence_records.emplace_back(
          &reads[read_index],
          assembly.molecules[read_index].representative_fragment_id);
    }

    for (const auto &[read_pointer, fragment_id] : evidence_records) {
      const auto &read = *read_pointer;
      const auto &effective_config =
          capture_sparse_store ? fragment_config : config;
      if (!passes_snp_alignment_filters(read, effective_config) ||
          read.cigar_operations.empty() || reference_length == 0U) {
        continue;
      }
      std::size_t query_cursor = 0U;
      std::size_t reference_cursor =
          ((read.reference_start - 1U) % reference_length) + 1U;
      for (const auto &operation : read.cigar_operations) {
        const auto len = operation.length;
        const char op = operation.code;
        if (op == 'M' || op == '=' || op == 'X') {
          for (std::size_t offset = 0U;
               offset < len && query_cursor + offset < read.sequence.size();
               ++offset) {
            const std::size_t position =
                ((reference_cursor - 1U + offset) % reference_length) + 1U;
            if (!target_positions.contains(position)) {
              continue;
            }
            auto &molecule_observation = molecule_alleles[position];
            if (!molecule_observation.covered ||
                fragment_id < molecule_observation.alignment_fragment_id) {
              molecule_observation.alignment_fragment_id = fragment_id;
              molecule_observation.mapping_quality = read.mapping_quality;
              molecule_observation.strand =
                  (read.flags & 0x10U) == 0U ? '+' : '-';
            }
            molecule_observation.covered = true;
            const auto quality = phred_quality_at(read, query_cursor + offset);
            const char base =
                static_cast<char>(std::toupper(static_cast<unsigned char>(
                    read.sequence[query_cursor + offset])));
            if (!quality || *quality < config.min_base_quality ||
                !base_index(base)) {
              continue;
            }
            const auto query_index = query_cursor + offset;
            const double center_proximity =
                read.sequence.size() <= 1U
                    ? 0.0
                    : (2.0 * static_cast<double>(
                                 std::min(query_index, read.sequence.size() -
                                                           1U - query_index))) /
                          static_cast<double>(read.sequence.size() - 1U);
            if (!molecule_observation.has_passing_observation) {
              molecule_observation.base = base;
              molecule_observation.quality = *quality;
              molecule_observation.center_proximity = center_proximity;
              molecule_observation.query_index = query_index;
              molecule_observation.alignment_fragment_id = fragment_id;
              molecule_observation.mapping_quality = read.mapping_quality;
              molecule_observation.strand =
                  (read.flags & 0x10U) == 0U ? '+' : '-';
              molecule_observation.has_passing_observation = true;
            } else if (!molecule_observation.conflicted) {
              if (molecule_observation.base != base) {
                molecule_observation.conflicted = true;
                molecule_observation.base = 'N';
              } else if (*quality > molecule_observation.quality ||
                         (*quality == molecule_observation.quality &&
                          (fragment_id <
                               molecule_observation.alignment_fragment_id ||
                           (fragment_id ==
                                molecule_observation.alignment_fragment_id &&
                            query_index < molecule_observation.query_index)))) {
                molecule_observation.quality = *quality;
                molecule_observation.center_proximity = center_proximity;
                molecule_observation.query_index = query_index;
                molecule_observation.alignment_fragment_id = fragment_id;
                molecule_observation.mapping_quality = read.mapping_quality;
                molecule_observation.strand =
                    (read.flags & 0x10U) == 0U ? '+' : '-';
              }
            }
          }
          query_cursor += len;
          reference_cursor =
              ((reference_cursor - 1U + (len % reference_length)) %
               reference_length) +
              1U;
        } else if (op == 'I' || op == 'S') {
          query_cursor += len;
        } else if (op == 'D' || op == 'N') {
          reference_cursor =
              ((reference_cursor - 1U + (len % reference_length)) %
               reference_length) +
              1U;
        }
      }
    }

    for (const auto &[position, observation] : molecule_alleles) {
      const auto position_aggregates = aggregates_by_position.find(position);
      if (position_aggregates == aggregates_by_position.end()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "target SNP position has no normalized event");
      }
      if (capture_sparse_store) {
        if (read_index >= assembly.molecules.size()) {
          throw AnalysisError(
              AnalysisErrorCode::internal_error,
              "molecule assembly and representative order diverged");
        }
        const auto &molecule = assembly.molecules[read_index];
        for (const auto &binding : position_aggregates->second) {
          if (result.store.observations.size() >=
              config.max_evidence_observations) {
            throw AnalysisError(
                AnalysisErrorCode::resource_exhausted,
                "schema 0.6 sparse observation limit exceeded "
                "(max_evidence_observations=" +
                    std::to_string(config.max_evidence_observations) + ")");
          }
          EvidenceObservation stored;
          stored.id = result.store.observations.size();
          stored.molecule_index = molecule.index;
          stored.event_index = binding.event_index;
          stored.alignment_fragment_id = observation.alignment_fragment_id;
          stored.mapping_quality = observation.mapping_quality;
          stored.strand = observation.strand;
          stored.evidence_source = "aligned_base";
          if (observation.conflicted) {
            stored.state = ObservationState::conflict;
          } else if (!observation.has_passing_observation) {
            stored.state = ObservationState::low_quality;
          } else {
            stored.observed_allele.assign(1U, observation.base);
            stored.base_quality = observation.quality;
            stored.center_proximity = observation.center_proximity;
            if (observation.base == binding.aggregate->call.alternate) {
              stored.state = ObservationState::alternate;
            } else if (observation.base == binding.aggregate->call.reference) {
              stored.state = ObservationState::reference;
            } else {
              stored.state = ObservationState::event_absent;
            }
          }
          result.store.observations.push_back(std::move(stored));
        }
      }
      if (observation.conflicted || !observation.has_passing_observation) {
        continue;
      }
      if (const auto index = base_index(observation.base)) {
        const std::size_t strand_index = observation.strand == '+' ? 0U : 1U;
        auto &evidence = allele_evidence[position];
        ++evidence.depths[*index];
        ++evidence.strand_depths[*index][strand_index];
        evidence.quality_sums[*index] +=
            static_cast<double>(observation.quality);
        evidence.quality_mins[*index] =
            std::min(evidence.quality_mins[*index], observation.quality);
        evidence.quality_maxes[*index] =
            std::max(evidence.quality_maxes[*index], observation.quality);
        evidence.read_position_sums[*index] += observation.center_proximity;
        for (const auto &binding : position_aggregates->second) {
          if (observation.base == binding.aggregate->call.alternate) {
            binding.aggregate->supporting_reads.push_back(
                features[read_index].id);
          }
        }
      }
    }
  }

  constexpr double z = 1.959963984540054;
  for (auto &[_, aggregate] : aggregates) {
    auto &support = aggregate.supporting_reads;
    std::sort(support.begin(), support.end());
    support.erase(std::unique(support.begin(), support.end()), support.end());
    const auto evidence_it = allele_evidence.find(aggregate.call.position);
    if (evidence_it == allele_evidence.end()) {
      continue;
    }
    const auto &evidence = evidence_it->second;
    const auto &depths = evidence.depths;
    aggregate.callable_depth =
        std::accumulate(depths.begin(), depths.end(), std::size_t{0});
    const auto alternate = base_index(aggregate.call.alternate);
    const auto reference_allele = base_index(aggregate.call.reference);
    if (!alternate || !reference_allele) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "aggregate SNP contains a non-canonical allele");
    }
    aggregate.alternate_depth = depths[*alternate];
    aggregate.reference_depth = depths[*reference_allele];
    aggregate.other_depth = aggregate.callable_depth -
                            aggregate.alternate_depth -
                            aggregate.reference_depth;
    aggregate.alternate_strand_depths = evidence.strand_depths[*alternate];
    aggregate.reference_strand_depths =
        evidence.strand_depths[*reference_allele];
    aggregate.alternate_quality_sum = evidence.quality_sums[*alternate];
    aggregate.reference_quality_sum = evidence.quality_sums[*reference_allele];
    aggregate.alternate_quality_min = evidence.quality_mins[*alternate];
    aggregate.reference_quality_min = evidence.quality_mins[*reference_allele];
    aggregate.alternate_quality_max = evidence.quality_maxes[*alternate];
    aggregate.reference_quality_max = evidence.quality_maxes[*reference_allele];
    aggregate.alternate_read_position_sum =
        evidence.read_position_sums[*alternate];
    aggregate.reference_read_position_sum =
        evidence.read_position_sums[*reference_allele];
    for (std::size_t index = 0U; index < depths.size(); ++index) {
      if (index == *alternate || index == *reference_allele) {
        continue;
      }
      for (std::size_t strand = 0U; strand < 2U; ++strand) {
        aggregate.other_strand_depths[strand] +=
            evidence.strand_depths[index][strand];
      }
      aggregate.other_quality_sum += evidence.quality_sums[index];
      aggregate.other_read_position_sum += evidence.read_position_sums[index];
      if (depths[index] != 0U) {
        aggregate.other_quality_min =
            std::min(aggregate.other_quality_min, evidence.quality_mins[index]);
        aggregate.other_quality_max = std::max(aggregate.other_quality_max,
                                               evidence.quality_maxes[index]);
      }
    }
    if (aggregate.callable_depth == 0) {
      continue;
    }
    const double n = static_cast<double>(aggregate.callable_depth);
    const double proportion =
        static_cast<double>(aggregate.alternate_depth) / n;
    aggregate.heteroplasmy = proportion;
    const double z2 = z * z;
    const double denominator = 1.0 + z2 / n;
    const double center = (proportion + z2 / (2.0 * n)) / denominator;
    const double margin = z *
                          std::sqrt((proportion * (1.0 - proportion) / n) +
                                    (z2 / (4.0 * n * n))) /
                          denominator;
    aggregate.ci95_low = std::max(0.0, center - margin);
    aggregate.ci95_high = std::min(1.0, center + margin);
  }

  if (capture_sparse_store) {
    if (assembly.molecules.size() != reads.size()) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 molecule/representative cardinality mismatch");
    }
    std::vector<std::size_t> fragment_reference_counts(
        assembly.fragments.size(), 0U);
    std::set<std::string> molecule_ids;
    for (std::size_t molecule_index = 0;
         molecule_index < assembly.molecules.size(); ++molecule_index) {
      const auto &molecule = assembly.molecules[molecule_index];
      if (molecule.index.value != molecule_index ||
          !molecule_ids.insert(molecule.id).second) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 molecule IDs are not unique and contiguous");
      }
      if (molecule.representative_fragment_id.value >=
          assembly.fragments.size()) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 representative fragment reference is unresolved");
      }
      for (const auto fragment_id : molecule.fragment_ids) {
        if (fragment_id.value >= assembly.fragments.size() ||
            assembly.fragments[fragment_id.value].molecule_index !=
                molecule.index) {
          throw AnalysisError(
              AnalysisErrorCode::internal_error,
              "schema 0.6 molecule fragment reference is unresolved");
        }
        ++fragment_reference_counts[fragment_id.value];
      }
    }
    if (std::any_of(fragment_reference_counts.begin(),
                    fragment_reference_counts.end(),
                    [](const std::size_t count) { return count != 1U; })) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "schema 0.6 alignment fragment must resolve to exactly one molecule");
    }

    std::set<std::string> event_ids;
    for (std::size_t event_index = 0; event_index < result.store.events.size();
         ++event_index) {
      if (!event_ids.insert(result.store.events[event_index].id).second) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "schema 0.6 normalized event IDs are not unique");
      }
    }

    struct ProjectionCounts {
      std::size_t alternate = 0;
      std::size_t reference = 0;
      std::size_t other = 0;
      std::vector<std::string> supporting_molecules;
    };
    std::vector<ProjectionCounts> projected(result.store.events.size());
    std::set<std::pair<std::size_t, std::size_t>> molecule_event_pairs;
    for (std::size_t observation_index = 0;
         observation_index < result.store.observations.size();
         ++observation_index) {
      const auto &observation = result.store.observations[observation_index];
      if (observation.id != observation_index ||
          observation.molecule_index.value >= assembly.molecules.size() ||
          observation.event_index >= result.store.events.size() ||
          observation.alignment_fragment_id.value >=
              assembly.fragments.size()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "schema 0.6 observation reference is unresolved");
      }
      if (assembly.fragments[observation.alignment_fragment_id.value]
              .molecule_index != observation.molecule_index) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 observation fragment belongs to another molecule");
      }
      if (!molecule_event_pairs
               .emplace(observation.molecule_index.value,
                        observation.event_index)
               .second) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 molecule/event observation is not unique");
      }
      auto &counts = projected[observation.event_index];
      switch (observation.state) {
      case ObservationState::alternate:
        ++counts.alternate;
        counts.supporting_molecules.push_back(
            assembly.molecules[observation.molecule_index.value].id);
        break;
      case ObservationState::reference:
        ++counts.reference;
        break;
      case ObservationState::event_absent:
        ++counts.other;
        break;
      case ObservationState::low_quality:
      case ObservationState::conflict:
        break;
      case ObservationState::not_callable:
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "sparse schema 0.6 store must encode NOT_CALLABLE implicitly");
      }
    }

    for (std::size_t event_index = 0; event_index < result.store.events.size();
         ++event_index) {
      const auto &event = result.store.events[event_index];
      const auto aggregate_it = aggregates.find(snp_key(
          event.start, event.reference.front(), event.alternate.front()));
      if (aggregate_it == aggregates.end()) {
        throw AnalysisError(AnalysisErrorCode::internal_error,
                            "schema 0.6 event has no variant projection");
      }
      const auto &aggregate = aggregate_it->second;
      auto &counts = projected[event_index];
      std::sort(counts.supporting_molecules.begin(),
                counts.supporting_molecules.end());
      result.store.events[event_index].supporting_molecules =
          counts.supporting_molecules;
      if (counts.alternate != aggregate.alternate_depth ||
          counts.reference != aggregate.reference_depth ||
          counts.other != aggregate.other_depth ||
          counts.alternate + counts.reference + counts.other !=
              aggregate.callable_depth ||
          counts.supporting_molecules != aggregate.supporting_reads) {
        throw AnalysisError(
            AnalysisErrorCode::internal_error,
            "schema 0.6 evidence/variant projection invariant failed for " +
                event.id);
      }
    }
  }
  return result;
}

} // namespace mito::detail
