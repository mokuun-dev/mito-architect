#include "detail/pipeline.hpp"
#include "detail/precision.hpp"

namespace mito::detail {

using AnalysisClock = std::chrono::steady_clock;

[[nodiscard]] std::uint64_t
elapsed_microseconds(const AnalysisClock::time_point start) {
  return static_cast<std::uint64_t>(
      std::chrono::duration_cast<std::chrono::microseconds>(
          AnalysisClock::now() - start)
          .count());
}

[[nodiscard]] std::string analyze_impl(const std::string &input_path,
                                       const std::string &reference_path,
                                       const AnalysisConfig &config,
                                       AnalysisPhaseTimings *timings) {
  if (timings != nullptr) {
    *timings = {};
  }
  const auto total_start = AnalysisClock::now();
  validate_config(config);
  throw_if_cancelled(config);

  const auto reference_start = AnalysisClock::now();
  const auto reference = load_reference(reference_path);
  if (timings != nullptr) {
    timings->reference_load_us = elapsed_microseconds(reference_start);
  }
  throw_if_cancelled(config);

  const auto input_start = AnalysisClock::now();
  auto input =
      read_input(input_path, config, reference,
                 reference_path.empty() ? bundled_rcrs_path() : reference_path);
  if (input.reads.empty()) {
    throw AnalysisError(AnalysisErrorCode::input_empty,
                        "input contains no read records: " + input_path);
  }
  const auto numt_assessment = assess_numt_provenance(input_path, input);
  const std::size_t input_record_count = input.reads.size();
  auto assembly = assemble_molecules(input, config, reference.size());
  const auto &reads = assembly.representatives;
  if (timings != nullptr) {
    timings->input_ingest_us = elapsed_microseconds(input_start);
  }
  throw_if_cancelled(config);

  std::vector<ReadFeature> features;
  features.reserve(reads.size());
  std::map<std::string, SvCall> svs;
  std::map<std::string, ComplexSvCall> complex_events;
  const auto resource_start = AnalysisClock::now();
  const auto clinical_annotations = load_clinical_annotations();
  static const auto base_resources = load_resource_manifest();
  auto resources = base_resources;
  apply_runtime_resource_overrides(resources);
  static const auto haplogroup_weights = load_phylo_weights();
  static const auto haplogroup_alignment_rules = load_phylo_alignment_rules();
  static const auto haplogroup_definitions = [] {
    auto definitions = load_haplogroups();
    for (auto &definition : definitions) {
      for (auto &mutation : definition.mutations) {
        mutation.weight = mutation_weight(haplogroup_weights, mutation.encoded);
      }
    }
    return definitions;
  }();
  if (timings != nullptr) {
    timings->resource_load_us = elapsed_microseconds(resource_start);
  }

  const auto extraction_start = AnalysisClock::now();
  auto extraction_results =
      extract_features_parallel(reads, reference, config, clinical_annotations);
  if (timings != nullptr) {
    timings->feature_extraction_us = elapsed_microseconds(extraction_start);
  }

  const auto event_merge_start = AnalysisClock::now();
  for (auto &result : extraction_results) {
    throw_if_cancelled(config);
    for (auto &sv : result.svs) {
      merge_sv(svs, std::move(sv));
    }
    for (auto &complex_event : result.complex_events) {
      merge_complex_sv(complex_events, std::move(complex_event));
    }
    features.push_back(std::move(result.feature));
  }
  finalize_merged_svs(svs, complex_events);
  if (timings != nullptr) {
    timings->event_merge_us = elapsed_microseconds(event_merge_start);
  }

  const auto precision_start = AnalysisClock::now();
  PrecisionCallResult precision;
  if (input.alignment_input) {
    precision = call_precision_variants(assembly, features, reference, config);
  }
  // Noisy long reads carry thousands of per-read errors. Above the configured
  // size the read, cluster and event projections keep only supported events.
  const bool compact_result = precision.assessed &&
                              config.compact_result_min_molecules != 0U &&
                              assembly.molecules.size() >= config.compact_result_min_molecules;
  if (compact_result) {
    apply_compact_projection(precision, reference, features, svs, complex_events);
  }
  if (timings != nullptr) {
    timings->precision_calling_us = elapsed_microseconds(precision_start);
  }
  throw_if_cancelled(config);

  std::string_view clustering_backend =
      "Read clusters use DBSCAN over SNP/SV feature tokens";
  const auto clustering_start = AnalysisClock::now();
  if (compact_result) {
    // Token clustering needs each read to span most of its lineage's markers. The
    // mixture model treats an uncovered site as missing, so partial reads phase.
    assign_haplotype_clusters(precision, features);
    clustering_backend =
        "Read clusters are haplotypes from a Bernoulli mixture over supported "
        "heteroplasmic SNVs";
  } else {
#ifdef MITO_HAS_HDBSCAN_CPP
    if (assign_hdbscan_clusters(features, config.min_cluster_size, config)) {
      clustering_backend =
          "Read clusters use HDBSCAN-C++ over dense SNP/SV feature vectors";
    } else {
      assign_dbscan_clusters(features, config.cluster_epsilon,
                             config.min_cluster_size, config);
      clustering_backend =
          "Read clusters use DBSCAN fallback after HDBSCAN-C++ adapter failure";
    }
#else
    assign_dbscan_clusters(features, config.cluster_epsilon,
                           config.min_cluster_size, config);
#endif
  }
  if (timings != nullptr) {
    timings->clustering_us = elapsed_microseconds(clustering_start);
  }

  throw_if_cancelled(config);
  const auto aggregation_start = AnalysisClock::now();
  auto coverage = compute_coverage(reads, features, config, reference.size());
  if (compact_result) {
    coverage = coverage_from_molecule_depth(
        precision.depth_by_position,
        coverage.bins.empty() ? std::size_t{180} : coverage.bins.size());
  }
  auto evidence = aggregate_snps(assembly, features, config, reference.size(),
                                 reference, clinical_annotations);
  CallabilityResult callability;
  std::vector<PhaseLink> phase_links;
  ArchitectureInferenceResult architecture_inference;
  if (config.result_schema == ResultSchema::v0_6) {
    callability = build_callability(assembly, features, reference, config);
    append_small_indel_evidence(evidence.store, assembly, features, callability,
                                reference, config);
    append_structural_evidence(evidence.store, assembly, features, callability,
                               svs, complex_events, reference.size(), config);
    validate_unified_evidence_graph(assembly, callability, evidence.store,
                                    reference.size());
    const auto phase_start = AnalysisClock::now();
    phase_links = build_phase_links(evidence.store, assembly, config);
    if (timings != nullptr) {
      timings->phase_link_build_us = elapsed_microseconds(phase_start);
    }
  }
  if (timings != nullptr) {
    timings->evidence_aggregation_us = elapsed_microseconds(aggregation_start);
  }

  const auto architecture_start = AnalysisClock::now();
  if (config.result_schema == ResultSchema::v0_6) {
    architecture_inference =
        infer_architectures(evidence.store, assembly, callability, features,
                            reference.size(), config);
    validate_architecture_inference(architecture_inference, assembly);
  }
  if (timings != nullptr) {
    timings->architecture_inference_us =
        elapsed_microseconds(architecture_start);
  }

  const auto haplogroup_start = AnalysisClock::now();
  const auto haplogroups =
      assign_haplogroups(features, haplogroup_definitions, haplogroup_weights,
                         haplogroup_alignment_rules, config);
  const auto architecture_haplogroups =
      config.result_schema == ResultSchema::v0_6
          ? assign_architecture_haplogroups(
                features, architecture_inference, haplogroup_definitions,
                haplogroup_weights, haplogroup_alignment_rules, config)
          : std::vector<ClusterHaplogroupAssignment>{};
  if (timings != nullptr) {
    timings->haplogroup_assignment_us = elapsed_microseconds(haplogroup_start);
  }

  if (precision.assessed) {
    const auto consensus_feature =
        consensus_haplogroup_feature(precision, reference, assembly);
    const auto consensus_assignment = assign_haplogroups(
        {consensus_feature}, haplogroup_definitions, haplogroup_weights,
        haplogroup_alignment_rules, config);
    if (const auto found = consensus_assignment.find(0);
        found != consensus_assignment.end()) {
      precision.consensus_haplogroup = found->second;
    }
    std::vector<ReadFeature> haplotype_features;
    haplotype_features.reserve(precision.haplotypes.size());
    for (const auto &haplotype : precision.haplotypes) {
      haplotype_features.push_back(haplotype.feature);
    }
    const auto haplotype_assignments = assign_haplogroups(
        haplotype_features, haplogroup_definitions, haplogroup_weights,
        haplogroup_alignment_rules, config);
    for (auto &haplotype : precision.haplotypes) {
      if (const auto found = haplotype_assignments.find(haplotype.id);
          found != haplotype_assignments.end()) {
        haplotype.haplogroup = found->second;
      }
    }
  }

  const auto serialization_start = AnalysisClock::now();
  auto json =
      render_json(input_path, reference_path, config, clustering_backend,
                  numt_assessment, input_record_count,
                  reference, assembly, reads, features, svs, complex_events,
                  coverage, evidence.variants, evidence.store, callability,
                  phase_links, architecture_inference, haplogroups,
                  architecture_haplogroups, resources);
  // The precision section is additive, so it is spliced before the final brace
  // and neither schema's existing fields change.
  if (!json.empty() && json.back() == '}') {
    json.pop_back();
    json += ",\"precision_calls\":";
    json += render_precision_calls_json(precision, config);
    json.push_back('}');
  }
  if (json.size() > config.max_result_bytes) {
    throw AnalysisError(
        AnalysisErrorCode::resource_exhausted,
        "result JSON limit exceeded (max_result_bytes=" +
            std::to_string(config.max_result_bytes) + ")");
  }
  if (timings != nullptr) {
    timings->serialization_us = elapsed_microseconds(serialization_start);
    timings->total_us = elapsed_microseconds(total_start);
  }
  return json;
}

} // namespace mito::detail
