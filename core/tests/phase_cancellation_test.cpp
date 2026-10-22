#include "detail/pipeline.hpp"
#include <iostream>

int main() {
  using namespace mito;
  using namespace mito::detail;
  SparseEvidenceStore store;
  MoleculeAssemblyResult assembly;
  assembly.molecules.resize(1);
  for (std::size_t i = 0; i < 256; ++i) {
    store.events.emplace_back();
    EvidenceObservation observation;
    observation.event_index = i;
    observation.state = ObservationState::alternate;
    store.observations.push_back(observation);
  }
  AnalysisConfig config;
  unsigned polls = 0;
  config.should_cancel = [&]() { return ++polls >= 8; };
  try {
    (void)build_phase_links(store, assembly, config);
    std::cerr << "phase construction ignored cancellation\n";
    return 1;
  } catch (const AnalysisError &error) {
    if (error.code() != AnalysisErrorCode::analysis_cancelled)
      return 2;
  }
  config.should_cancel = {};
  config.max_phase_work = 1024;
  config.max_phase_links = 100'000;
  try {
    (void)build_phase_links(store, assembly, config);
    std::cerr << "phase construction ignored work budget\n";
    return 6;
  } catch (const AnalysisError &error) {
    if (error.code() != AnalysisErrorCode::resource_exhausted)
      return 7;
  }
  config.max_phase_work = 20'000'000;
  config.max_phase_links = 10;
  try {
    (void)build_phase_links(store, assembly, config);
    std::cerr << "phase construction ignored pair budget\n";
    return 3;
  } catch (const AnalysisError &error) {
    if (error.code() != AnalysisErrorCode::resource_exhausted)
      return 4;
  }
  // Recovery preserves the complete contingency table for a small valid graph.
  store.events.resize(2);
  store.observations.resize(2);
  const auto links = build_phase_links(store, assembly, config);
  if (links.size() != 1 || links[0].both_alternate != 1 ||
      links[0].jointly_callable != 1 ||
      links[0].supporting_molecule_indices.size() != 1)
    return 5;
  // Per-link molecule traceability has an independent limit.
  assembly.molecules.resize(2);
  store.observations.resize(4);
  store.observations[0].molecule_index.value = 0;
  store.observations[1].molecule_index.value = 0;
  store.observations[2].molecule_index.value = 1;
  store.observations[3].molecule_index.value = 1;
  store.observations[2].event_index = 0;
  store.observations[3].event_index = 1;
  store.observations[2].state = ObservationState::alternate;
  store.observations[3].state = ObservationState::alternate;
  config.max_phase_links = 10;
  config.max_phase_molecule_references = 1;
  try {
    (void)build_phase_links(store, assembly, config);
    std::cerr << "phase construction ignored molecule-reference budget\n";
    return 8;
  } catch (const AnalysisError &error) {
    if (error.code() != AnalysisErrorCode::resource_exhausted)
      return 9;
  }
  return 0;
}
