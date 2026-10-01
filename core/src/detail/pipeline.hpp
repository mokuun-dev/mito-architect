#pragma once
#include "evidence.hpp"
#include "helpers.hpp"
#include "reads.hpp"
#include "resources.hpp"
#include "variants.hpp"

namespace mito::detail {

// serialization
[[nodiscard]] std::string escape_json(std::string_view value);
[[nodiscard]] std::string quoted(std::string_view value);
// Prefer the JSON escaper above over std::quoted when argument-dependent lookup
// sees std::string. The latter only escapes quotes/backslashes and is not a
// complete JSON control-character encoder.
[[nodiscard]] std::string quoted(const std::string &value);
[[nodiscard]] std::string quoted(std::string &value);
[[nodiscard]] std::string quoted(const char *value);
[[nodiscard]] bool has_clinical_annotation(const SnpCall &snp);
[[nodiscard]] std::string
alignment_fragment_id_string(const AlignmentFragmentId id);
[[nodiscard]] std::string_view
observation_state_name(const ObservationState state);
[[nodiscard]] std::string_view architecture_assignment_status_name(
    const ArchitectureAssignmentStatus status) noexcept;
[[nodiscard]] bool is_unified_variant_event(const EvidenceEvent &event);
void add_variant_allele_observation(VariantAlleleEvidenceSummary &summary,
                                    const EvidenceObservation &observation);
[[nodiscard]] std::vector<UnifiedVariantEvidenceSummary>
summarize_unified_variant_evidence(const SparseEvidenceStore &store);
[[nodiscard]] std::pair<double, double>
wilson_interval(const std::size_t support, const std::size_t depth);
[[nodiscard]] std::size_t
circular_reference_homopolymer_run(const std::string &reference,
                                   const std::size_t position);
[[nodiscard]] std::vector<std::string>
variant_qc_flags(const EvidenceEvent &event,
                 const UnifiedVariantEvidenceSummary &summary,
                 const std::size_t homopolymer_run, const bool numt_assessable);
[[nodiscard]] std::string
render_json(const std::string &input_path, const std::string &reference_path,
            const AnalysisConfig &config, std::string_view clustering_backend,
            const NumtAssessment &numt_assessment,
            std::size_t input_record_count, const std::string &reference,
            const MoleculeAssemblyResult &assembly,
            const std::vector<ReadRecord> &reads,
            const std::vector<ReadFeature> &features,
            const std::map<std::string, SvCall> &svs,
            const std::map<std::string, ComplexSvCall> &complex_events,
            const CoverageResult &coverage_result,
            const std::map<std::string, SnpAggregate> &variants,
            const SparseEvidenceStore &evidence_store,
            const CallabilityResult &callability,
            const std::vector<PhaseLink> &phase_links,
            const ArchitectureInferenceResult &architecture_inference,
            const std::map<int, ClusterHaplogroupAssignment> &haplogroups,
            const std::vector<ClusterHaplogroupAssignment>
                &architecture_haplogroups,
            const std::vector<ResourceRecord> &resources);

// resources
[[nodiscard]] std::vector<std::string> split_tab(std::string_view line);
[[nodiscard]] std::string trim_copy(std::string_view value);
void strip_trailing_carriage_return(std::string &value);
[[nodiscard]] std::vector<std::string> split_semicolon(std::string_view line);
[[nodiscard]] std::string parse_reference_fasta(std::istream &in,
                                                const std::string &source);
[[nodiscard]] std::string bundled_rcrs_path();
[[nodiscard]] std::string bundled_clinical_annotations_path();
[[nodiscard]] std::string clinical_annotations_path();
[[nodiscard]] std::string phylotree_path();
[[nodiscard]] std::string phylotree_weights_path();
[[nodiscard]] std::string phylotree_alignment_rules_path();
[[nodiscard]] std::string resource_manifest_path();
[[nodiscard]] std::vector<ResourceRecord> load_resource_manifest();
void apply_runtime_resource_overrides(std::vector<ResourceRecord> &records);
void throw_if_cancelled(const AnalysisConfig &config);
[[nodiscard]] std::string load_reference(const std::string &reference_path);
[[nodiscard]] std::vector<GeneAnnotation> default_genes();

// configuration
[[nodiscard]] bool valid_sam_tag_name(const std::string &tag);
void validate_config(const AnalysisConfig &config);

// haplogroup_markers
[[nodiscard]] std::optional<std::size_t> parse_size(std::string_view value);
[[nodiscard]] std::optional<PhyloMutation>
parse_phylo_mutation(std::string_view raw_token);
[[nodiscard]] std::vector<CigarOperation> parse_cigar(std::string_view cigar,
                                                      std::string_view read_id);
[[nodiscard]] std::vector<std::string>
split_alignment_rule_tokens(std::string_view value);
[[nodiscard]] bool is_reference_restore_token(std::string_view token);
[[nodiscard]] std::vector<PhyloAlignmentRule> load_phylo_alignment_rules();
void apply_phylo_alignment_rules(std::vector<ObservedPhyloMutation> &observed,
                                 const std::vector<PhyloAlignmentRule> &rules);
[[nodiscard]] std::vector<ObservedPhyloMutation>
phylo_tags_from_header(const ReadRecord &read);
[[nodiscard]] std::vector<ObservedPhyloMutation>
call_phylo_indels_from_alignment(const ReadRecord &read,
                                 const AnalysisConfig &config,
                                 const std::string &reference);

// input
[[nodiscard]] std::vector<SupplementaryAlignment>
parse_supplementary_alignments(std::string_view encoded_tag,
                               std::string_view read_id);
[[nodiscard]] std::map<std::string, std::string>
parse_sam_aux_tags(const std::vector<std::string> &fields,
                   std::size_t line_number);
[[nodiscard]] std::string snp_key(std::size_t position, char reference,
                                  char alternate);
[[nodiscard]] InputRecords parse_fastq(std::istream &input,
                                       const AnalysisConfig &config);
[[nodiscard]] InputRecords parse_sam(std::istream &input,
                                     const AnalysisConfig &config);
[[nodiscard]] bool
has_extension(const std::filesystem::path &path,
              std::initializer_list<std::string_view> extensions);
[[nodiscard]] InputRecords read_input(const std::string &input_path,
                                      const AnalysisConfig &config);
[[nodiscard]] NumtAssessment assess_numt_provenance(
    const std::string &input_path, const InputRecords &input);

// clinical
[[nodiscard]] std::string lowercase_token(std::string_view value);
[[nodiscard]] std::string
normalize_clinical_significance(std::string_view significance);
[[nodiscard]] std::string
clinical_significance_group(std::string_view normalized);
[[nodiscard]] std::string join_values(const std::vector<std::string> &values,
                                      std::string_view delimiter);
void finalize_clinical_annotation(SnpCall &call);
[[nodiscard]] std::map<std::string, SnpCall> load_clinical_annotations();
void apply_clinical_annotations(
    std::vector<SnpCall> &snps,
    const std::map<std::string, SnpCall> &annotations);

// molecules
[[nodiscard]] std::string alignment_as_sa(const ReadRecord &read);
[[nodiscard]] std::string alignment_fragment_role(const std::uint16_t flags,
                                                  const bool alignment_input);
void validate_protocol_tag_value(const std::string &tag,
                                 const std::string &value,
                                 const std::string &read_id);
[[nodiscard]] ProtocolTagValues
protocol_tag_values(const std::vector<ReadRecord> &reads,
                    const std::vector<std::size_t> &indices,
                    const std::string &tag);
[[nodiscard]] std::size_t reference_consuming_length(const ReadRecord &read);
void sort_unique_strings(std::vector<std::string> &values);
[[nodiscard]] MoleculeAssemblyResult
assemble_molecules(InputRecords &input, const AnalysisConfig &config,
                   const std::size_t reference_length);
[[nodiscard]] double mean_quality(const std::string &qualities);
[[nodiscard]] double gc_fraction(std::string_view sequence);

// snps
[[nodiscard]] bool looks_like_nuclear_contig(std::string_view reference_name);
[[nodiscard]] bool has_nuclear_supplementary_alignment(const ReadRecord &read);
[[nodiscard]] std::vector<std::string>
numt_evidence(const ReadRecord &read, std::size_t reference_length);
[[nodiscard]] double numt_score(const std::vector<std::string> &evidence);
[[nodiscard]] bool is_base(char base);
[[nodiscard]] bool passes_snp_alignment_filters(const ReadRecord &read,
                                                const AnalysisConfig &config);
[[nodiscard]] std::optional<std::uint8_t>
phred_quality_at(const ReadRecord &read, std::size_t query_index);
[[nodiscard]] std::optional<std::size_t> base_index(char base);
[[nodiscard]] std::vector<SnpCall>
call_snps_from_alignment(const ReadRecord &read, const std::string &reference,
                         const AnalysisConfig &config);
[[nodiscard]] std::vector<SnpCall>
call_snps_from_reference_span(const ReadRecord &read,
                              const std::string &reference,
                              const AnalysisConfig &config);
[[nodiscard]] std::optional<SnpCall> parse_snp_tag(std::string_view token);
[[nodiscard]] std::vector<SnpCall> snp_tags_from_header(const ReadRecord &read);
void merge_snps(std::vector<SnpCall> &snps, std::vector<SnpCall> tagged_snps);
[[nodiscard]] bool
observed_phylo_mutation_less(const ObservedPhyloMutation &lhs,
                             const ObservedPhyloMutation &rhs);
[[nodiscard]] ObservedPhyloMutation
observed_phylo_mutation(const PhyloMutation &mutation);

// indels
[[nodiscard]] std::vector<NormalizedSmallIndel>
call_small_indels_from_alignment(const ReadRecord &read,
                                 const AnalysisConfig &config,
                                 const std::string &reference);

// callability
void merge_haplogroup_markers(std::vector<ObservedPhyloMutation> &markers,
                              std::vector<ObservedPhyloMutation> additional);
[[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>>
ranges_from_bitmap(const std::vector<bool> &positions);
[[nodiscard]] std::optional<std::vector<std::pair<std::size_t, std::size_t>>>
phylo_ranges_from_header(const ReadRecord &read);
void append_callable_position(
    std::vector<std::pair<std::size_t, std::size_t>> &ranges,
    std::optional<std::pair<std::size_t, std::size_t>> &current,
    std::size_t position);
[[nodiscard]] std::vector<std::pair<std::size_t, std::size_t>>
normalize_ranges(std::vector<std::pair<std::size_t, std::size_t>> ranges);
[[nodiscard]] CallablePhyloRanges
callable_phylo_ranges_from_alignment(const ReadRecord &read,
                                     const std::string &reference,
                                     const AnalysisConfig &config);
[[nodiscard]] const ReadRecord &
source_record_for_fragment(const MoleculeAssemblyResult &assembly,
                           const AlignmentFragmentEvidence &fragment);
[[nodiscard]] std::string
alignment_callability_exclusion(const ReadRecord &read,
                                const AnalysisConfig &config);
[[nodiscard]] AlignmentCallabilityEvidence
alignment_callability(const AlignmentFragmentEvidence &fragment,
                      const ReadRecord &read, const std::string &reference,
                      const AnalysisConfig &config);
[[nodiscard]] CallabilityResult
build_callability(const MoleculeAssemblyResult &assembly,
                  const std::vector<ReadFeature> &features,
                  const std::string &reference, const AnalysisConfig &config);
[[nodiscard]] bool callable_position_in_ranges(
    const std::vector<std::pair<std::size_t, std::size_t>> &ranges,
    const std::size_t position);
[[nodiscard]] bool
span_in_ranges(const std::vector<std::pair<std::size_t, std::size_t>> &ranges,
               const std::size_t start, const std::size_t end);
[[nodiscard]] bool
reference_adjacency_callable(const MoleculeCallabilityEvidence &callability,
                             const std::size_t left, const std::size_t right);
[[nodiscard]] bool reference_span_with_flanks_callable(
    const MoleculeCallabilityEvidence &callability, const std::size_t start,
    const std::size_t end, const std::size_t reference_length);

// structural
[[nodiscard]] std::optional<SvCall> parse_sv_tag(std::string_view read_id,
                                                 std::string_view token);
[[nodiscard]] std::vector<SvCall> sv_tags_from_header(const ReadRecord &read);
[[nodiscard]] std::vector<SvCall> parse_cigar_svs(const ReadRecord &read,
                                                  std::size_t sv_min_length,
                                                  std::size_t reference_length);
[[nodiscard]] bool is_mitochondrial_contig(std::string_view name);
[[nodiscard]] std::optional<AlignmentSegment>
make_alignment_segment(std::string reference_name, std::size_t reference_start,
                       char strand, std::uint8_t mapping_quality,
                       const std::vector<CigarOperation> &operations);
[[nodiscard]] std::vector<AlignmentSegment>
alignment_segments(const ReadRecord &read);
[[nodiscard]] std::vector<SvCall>
split_alignment_svs(const ReadRecord &read, std::size_t reference_length,
                    std::size_t sv_min_length);
[[nodiscard]] char flip_strand(char strand);
[[nodiscard]] std::string
reverse_complement_orientation(std::string_view orientation);
[[nodiscard]] CanonicalJunctionPath
canonical_junction_path(const std::vector<SvCall> &junctions);
[[nodiscard]] std::optional<ComplexSvCall>
coalesce_complex_sv(std::string_view read_id,
                    const std::vector<SvCall> &junctions);
void merge_sv(std::map<std::string, SvCall> &svs, SvCall sv);
void merge_complex_sv(std::map<std::string, ComplexSvCall> &events,
                      ComplexSvCall event);

// features
[[nodiscard]] std::size_t effective_thread_count(std::size_t requested_threads,
                                                 std::size_t work_items);
[[nodiscard]] FeatureExtractionResult extract_read_feature(
    const ReadRecord &read, const std::string &reference,
    const AnalysisConfig &config,
    const std::map<std::string, SnpCall> &clinical_annotations);
[[nodiscard]] std::vector<FeatureExtractionResult> extract_features_parallel(
    const std::vector<ReadRecord> &reads, const std::string &reference,
    const AnalysisConfig &config,
    const std::map<std::string, SnpCall> &clinical_annotations);

// clustering
[[nodiscard]] FeatureTokens feature_tokens(const ReadFeature &feature);
[[nodiscard]] double token_distance(const FeatureTokens &lhs,
                                    const FeatureTokens &rhs);
[[nodiscard]] InvertedTokenIndex
build_inverted_token_index(const std::vector<TokenProfile> &profiles);
[[nodiscard]] Neighborhood
region_query(const std::vector<TokenProfile> &profiles,
             const InvertedTokenIndex &inverted_index,
             std::vector<std::size_t> &candidate_marks, std::size_t &mark_epoch,
             std::size_t index, double epsilon);
void append_new_neighbors(std::vector<std::size_t> &queue,
                          std::vector<unsigned char> &queued,
                          const std::vector<std::size_t> &candidates);
void assign_dbscan_clusters(std::vector<ReadFeature> &features, double epsilon,
                            std::size_t min_cluster_size,
                            const AnalysisConfig &config);
[[nodiscard]] bool assign_hdbscan_clusters(std::vector<ReadFeature> &features,
                                           std::size_t min_cluster_size,
                                           const AnalysisConfig &config);

// coverage
void add_circular_coverage_span(std::vector<std::int64_t> &difference,
                                std::uint64_t &uniform_depth, std::size_t start,
                                std::size_t length);
[[nodiscard]] CoverageResult
compute_coverage(const std::vector<ReadRecord> &reads,
                 const std::vector<ReadFeature> &features,
                 const AnalysisConfig &config, std::size_t reference_length,
                 std::size_t bin_count = 180);

// snp_evidence
[[nodiscard]] EvidenceAggregationResult
aggregate_snps(const MoleculeAssemblyResult &assembly,
               const std::vector<ReadFeature> &features,
               const AnalysisConfig &config, std::size_t reference_length,
               const std::string &reference,
               const std::map<std::string, SnpCall> &clinical_annotations);

// event_evidence
void append_sparse_observation(SparseEvidenceStore &store,
                               EvidenceObservation observation,
                               const AnalysisConfig &config);
[[nodiscard]] std::unordered_map<std::string, std::size_t>
molecule_index_by_id(const MoleculeAssemblyResult &assembly);
[[nodiscard]] std::size_t
circular_previous_position(const std::size_t position,
                           const std::size_t reference_length);
[[nodiscard]] std::size_t
circular_next_position(const std::size_t position,
                       const std::size_t reference_length);
void append_small_indel_evidence(SparseEvidenceStore &store,
                                 const MoleculeAssemblyResult &assembly,
                                 const std::vector<ReadFeature> &features,
                                 const CallabilityResult &callability,
                                 const std::string &reference,
                                 const AnalysisConfig &config);
[[nodiscard]] std::string uppercase_ascii(std::string value);
void append_structural_evidence(
    SparseEvidenceStore &store, const MoleculeAssemblyResult &assembly,
    const std::vector<ReadFeature> &features,
    const CallabilityResult &callability,
    const std::map<std::string, SvCall> &svs,
    const std::map<std::string, ComplexSvCall> &complex_events,
    const std::size_t reference_length, const AnalysisConfig &config);
void validate_unified_evidence_graph(const MoleculeAssemblyResult &assembly,
                                     const CallabilityResult &callability,
                                     const SparseEvidenceStore &store,
                                     const std::size_t reference_length);

// phase
[[nodiscard]] bool is_callable_phase_state(const ObservationState state);
[[nodiscard]] std::vector<PhaseLink>
build_phase_links(const SparseEvidenceStore &store,
                  const MoleculeAssemblyResult &assembly,
                  const AnalysisConfig &config);

// architectures
[[nodiscard]] bool
is_architecture_callable_state(const ObservationState state) noexcept;
[[nodiscard]] double
architecture_event_weight(const EvidenceEvent &event) noexcept;
[[nodiscard]] std::vector<MoleculeArchitectureEvidence>
build_molecule_architecture_evidence(const SparseEvidenceStore &store,
                                     const MoleculeAssemblyResult &assembly,
                                     const CallabilityResult &callability,
                                     const std::vector<ReadFeature> &features,
                                     const std::size_t reference_length);
[[nodiscard]] std::vector<ArchitectureEventProfile>
build_architecture_event_profile(
    const std::vector<std::size_t> &members,
    const std::vector<MoleculeArchitectureEvidence> &molecule_evidence,
    const AnalysisConfig &config, const std::size_t minimum_callable_support);
void update_architecture_signatures(CandidateArchitecture &architecture,
                                    const AnalysisConfig &config);
[[nodiscard]] std::string
architecture_id_for_signature(const std::vector<std::size_t> &signature,
                              const SparseEvidenceStore &store);
[[nodiscard]] ArchitectureScore
score_architecture_candidate(const std::size_t architecture_index,
                             const CandidateArchitecture &architecture,
                             const MoleculeArchitectureEvidence &molecule,
                             const SparseEvidenceStore &store);
[[nodiscard]] std::uint64_t mix_architecture_seed(std::uint64_t value);
[[nodiscard]] double
weighted_signature_jaccard(const std::vector<std::size_t> &lhs,
                           const std::vector<std::size_t> &rhs,
                           const SparseEvidenceStore &store);
void finalize_architecture_statistics(
    ArchitectureInferenceResult &result,
    const std::vector<MoleculeArchitectureEvidence> &molecule_evidence,
    const MoleculeAssemblyResult &assembly, const SparseEvidenceStore &store,
    const CallabilityResult &callability, const std::size_t reference_length,
    const AnalysisConfig &config);
void evaluate_architecture_stability(
    ArchitectureInferenceResult &result, const SparseEvidenceStore &store,
    const MoleculeAssemblyResult &assembly,
    const CallabilityResult &callability,
    const std::vector<ReadFeature> &features,
    std::size_t reference_length, const AnalysisConfig &config);
[[nodiscard]] ArchitectureInferenceResult infer_architectures(
    const SparseEvidenceStore &store, const MoleculeAssemblyResult &assembly,
    const CallabilityResult &callability,
    const std::vector<ReadFeature> &features,
    const std::size_t reference_length, const AnalysisConfig &config);
void validate_architecture_inference(const ArchitectureInferenceResult &result,
                                     const MoleculeAssemblyResult &assembly);

// haplogroups
[[nodiscard]] std::vector<HaplogroupDefinition> load_haplogroups();
[[nodiscard]] std::unordered_map<std::string, double> load_phylo_weights();
[[nodiscard]] double
mutation_weight(const std::unordered_map<std::string, double> &weights,
                const std::string &mutation);
[[nodiscard]] std::string macro_haplogroup(std::string_view name);
[[nodiscard]] bool
phylo_mutation_matches(const PhyloMutation &expected,
                       const ObservedPhyloMutation &observed);
[[nodiscard]] bool
observed_contains_mutation(const std::vector<ObservedPhyloMutation> &observed,
                           const PhyloMutation &expected);
[[nodiscard]] bool
definition_contains_observed(const std::vector<PhyloMutation> &expected,
                             const ObservedPhyloMutation &observed);
[[nodiscard]] bool position_in_ranges(
    const std::vector<std::pair<std::size_t, std::size_t>> &ranges,
    std::size_t position);
[[nodiscard]] CallablePhyloRanges
majority_callable_ranges(const std::vector<const ReadFeature *> &members);
[[nodiscard]] std::map<int, ClusterHaplogroupAssignment>
assign_haplogroups(const std::vector<ReadFeature> &features,
                   const std::vector<HaplogroupDefinition> &definitions,
                   const std::unordered_map<std::string, double> &weights,
                   const std::vector<PhyloAlignmentRule> &alignment_rules,
                   const AnalysisConfig &config);
[[nodiscard]] std::vector<ClusterHaplogroupAssignment>
assign_architecture_haplogroups(
    const std::vector<ReadFeature> &features,
    const ArchitectureInferenceResult &architecture_inference,
    const std::vector<HaplogroupDefinition> &definitions,
    const std::unordered_map<std::string, double> &weights,
    const std::vector<PhyloAlignmentRule> &alignment_rules,
    const AnalysisConfig &config);

// pipeline
[[nodiscard]] std::string analyze_impl(const std::string &input_path,
                                       const std::string &reference_path,
                                       const AnalysisConfig &config,
                                       AnalysisPhaseTimings *timings);

} // namespace mito::detail
