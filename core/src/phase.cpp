#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] bool is_callable_phase_state(const ObservationState state) {
  return state == ObservationState::reference ||
         state == ObservationState::alternate ||
         state == ObservationState::event_absent;
}

[[nodiscard]] std::vector<PhaseLink>
build_phase_links(const SparseEvidenceStore &store,
                  const MoleculeAssemblyResult &assembly,
                  const AnalysisConfig &config) {
  throw_if_cancelled(config);
  // Each unit corresponds to a candidate or a molecule/pair evaluation. This
  // turns the worst-case quadratic phase projection into an explicit contract.
  std::size_t work = 0;
  std::size_t molecule_references = 0;
  const auto consume_work = [&]() {
    if (work >= config.max_phase_work) {
      throw AnalysisError(
          AnalysisErrorCode::resource_exhausted,
          "schema 0.6 phase-work limit exceeded (max_phase_work=" +
              std::to_string(config.max_phase_work) + ")");
    }
    ++work;
    if ((work & 1023U) == 0U) {
      throw_if_cancelled(config);
    }
  };
  const auto retain_molecule_reference = [&]() {
    if (molecule_references >= config.max_phase_molecule_references) {
      throw AnalysisError(
          AnalysisErrorCode::resource_exhausted,
          "schema 0.6 phase molecule-reference limit exceeded "
          "(max_phase_molecule_references=" +
              std::to_string(config.max_phase_molecule_references) + ")");
    }
    ++molecule_references;
  };
  std::vector<std::vector<const EvidenceObservation *>> by_molecule(
      assembly.molecules.size());
  for (const auto &observation : store.observations) {
    consume_work();
    if (assembly.molecules[observation.molecule_index.value]
            .analysis_eligible) {
      by_molecule[observation.molecule_index.value].push_back(&observation);
    }
  }
  for (auto &observations : by_molecule) {
    consume_work();
    std::sort(observations.begin(), observations.end(),
              [](const auto *lhs, const auto *rhs) {
                return lhs->event_index < rhs->event_index;
              });
  }

  using Pair = std::pair<std::size_t, std::size_t>;
  // Keep only one ordered map. The old implementation kept a set, another map
  // containing the same keys, and a neighbour index before aggregation.
  std::map<Pair, PhaseLink> projected;
  for (const auto &observations : by_molecule) {
    for (std::size_t first = 0U; first < observations.size(); ++first) {
      if (observations[first]->state != ObservationState::alternate) {
        continue;
      }
      for (std::size_t second = 0U; second < observations.size(); ++second) {
        consume_work();
        if (first == second) {
          continue;
        }
        // A pair whose two endpoints are alternate would otherwise be seen
        // twice during candidate discovery.
        if (observations[second]->state == ObservationState::alternate &&
            first > second) {
          continue;
        }
        const Pair pair{std::min(observations[first]->event_index,
                                 observations[second]->event_index),
                        std::max(observations[first]->event_index,
                                 observations[second]->event_index)};
        if (pair.first == pair.second) {
          continue;
        }
        auto found = projected.find(pair);
        if (found == projected.end() &&
            projected.size() >= config.max_phase_links) {
          throw AnalysisError(
              AnalysisErrorCode::resource_exhausted,
              "schema 0.6 phase-link limit exceeded (max_phase_links=" +
                  std::to_string(config.max_phase_links) + ")");
        }
        if (found == projected.end()) {
          PhaseLink link;
          link.event_a_index = pair.first;
          link.event_b_index = pair.second;
          link.complete_callability =
              store.events[pair.first].absence_assessable &&
              store.events[pair.second].absence_assessable;
          projected.emplace(pair, std::move(link));
        }
      }
    }
  }

  // A second, direct molecule-pair pass is necessary so links discovered in a
  // later molecule still include earlier REF/REF and uncertain observations.
  // It avoids the previous per-observation neighbour scan and binary search.
  for (std::size_t molecule_index = 0U; molecule_index < by_molecule.size();
       ++molecule_index) {
    const auto &observations = by_molecule[molecule_index];
    for (std::size_t first = 0U; first < observations.size(); ++first) {
      for (std::size_t second = first + 1U; second < observations.size();
           ++second) {
        consume_work();
        const auto *observation_a = observations[first];
        const auto *observation_b = observations[second];
        const Pair pair{observation_a->event_index, observation_b->event_index};
        const auto found = projected.find(pair);
        if (found == projected.end()) {
          continue;
        }
        auto &link = found->second;
        const auto state_a = observation_a->state;
        const auto state_b = observation_b->state;
        if (!is_callable_phase_state(state_a) ||
            !is_callable_phase_state(state_b)) {
          ++link.jointly_uncertain;
          retain_molecule_reference();
          link.uncertain_molecule_indices.push_back(molecule_index);
          continue;
        }
        ++link.jointly_callable;
        const bool alternate_a = state_a == ObservationState::alternate;
        const bool alternate_b = state_b == ObservationState::alternate;
        if (alternate_a && alternate_b) {
          ++link.both_alternate;
          retain_molecule_reference();
          link.supporting_molecule_indices.push_back(molecule_index);
        } else if (alternate_a) {
          ++link.a_alternate_b_absent;
        } else if (alternate_b) {
          ++link.a_absent_b_alternate;
        } else {
          ++link.neither_alternate;
        }
      }
    }
  }

  constexpr double z = 1.959963984540054;
  std::vector<PhaseLink> links;
  links.reserve(projected.size());
  for (auto &[_, link] : projected) {
    consume_work();
    if (link.jointly_callable != 0U) {
      const double n = static_cast<double>(link.jointly_callable);
      const double both = static_cast<double>(link.both_alternate);
      link.co_alternate_fraction = both / n;
      const double z2 = z * z;
      const double denominator = 1.0 + z2 / n;
      const double center =
          (link.co_alternate_fraction + z2 / (2.0 * n)) / denominator;
      const double margin = z *
                            std::sqrt((link.co_alternate_fraction *
                                       (1.0 - link.co_alternate_fraction) / n) +
                                      (z2 / (4.0 * n * n))) /
                            denominator;
      link.ci95_low = std::max(0.0, center - margin);
      link.ci95_high = std::min(1.0, center + margin);
      const double p_a =
          static_cast<double>(link.both_alternate + link.a_alternate_b_absent) /
          n;
      const double p_b =
          static_cast<double>(link.both_alternate + link.a_absent_b_alternate) /
          n;
      link.expected_co_alternate_fraction = p_a * p_b;
      link.linkage_delta =
          link.co_alternate_fraction - link.expected_co_alternate_fraction;
    }
    links.push_back(link);
  }
  return links;
}

} // namespace mito::detail
