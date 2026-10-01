import type { ArchitectureAssignment } from './architectures';

export type ObservationState =
  | 'REFERENCE'
  | 'ALTERNATE'
  | 'EVENT_ABSENT'
  | 'NOT_CALLABLE'
  | 'LOW_QUALITY'
  | 'CONFLICT';

/** Original SAM/BAM/CRAM record retained before molecule assembly. */

export interface AlignmentFragment {
  id: string;
  source_record_index: number;
  molecule_id: string;
  molecule_index: number;
  role:
    | 'primary_candidate'
    | 'secondary'
    | 'supplementary'
    | 'secondary_and_supplementary'
    | 'unaligned_read'
    | string;
  selected_representative: boolean;
  flags: number;
  strand: '+' | '-' | string;
  mapping_quality: number;
  reference_name: string;
  reference_start: number;
  cigar: string;
  query_length: number;
  base_qualities_available: boolean;
  aux_tags: Record<string, string>;
}

/** Deterministic result of the configured fragment-to-molecule policy. */

export interface AssembledMolecule {
  id: string;
  index: number;
  identity_policy: 'sam_qname' | 'fastq_record_proxy' | string;
  assembly_status:
    | 'unique_primary'
    | 'fallback_without_primary'
    | 'first_of_multiple_primaries'
    | 'single_unaligned_fragment'
    | string;
  primary_candidate_count: number;
  ambiguous: boolean;
  analysis_eligible: boolean;
  evidence_eligible: boolean;
  callability_status: string;
  callable_bases: number;
  callable_fraction: number;
  query_length: number;
  mean_base_quality: number;
  mapping_quality: number;
  numt_score: number;
  numt_evidence: string[];
  cluster_id: number;
  architecture_assignment?: ArchitectureAssignment;
  alternate_event_ids: string[];
  evidence_state_counts: {
    alternate: number;
    reference: number;
    event_absent: number;
    low_quality: number;
    conflict: number;
  };
  representative_alignment_id: string;
  source_qnames: string[];
  protocol_metadata: Record<string, string>;
  protocol_flags: string[];
  exclusion_reasons: string[];
  alignment_ids: string[];
  warnings: string[];
}

export interface ReferenceRange {
  start: number;
  end: number;
}

export interface AlignmentCallability {
  alignment_id: string;
  eligible: boolean;
  status: string;
  callable_bases: number;
  inserted_query_bases: number;
  soft_clipped_query_bases: number;
  reference_exclusion_counts: Record<string, number>;
  disrupted_adjacency_anchors: number[];
  ranges: ReferenceRange[];
}

/** Per-molecule base-callability projection; absence is never inferred from a coverage gap. */

export interface MoleculeCallability {
  molecule_id: string;
  status: string;
  known: boolean;
  basis: 'passing_aligned_reference_bases' | string;
  callable_bases: number;
  callable_fraction: number;
  ranges: ReferenceRange[];
  alignments: AlignmentCallability[];
}

export interface EvidenceCounts {
  alternate: number;
  reference: number;
  event_absent: number;
  callable: number;
  low_quality: number;
  conflict: number;
}

/** Normalized schema 0.6 SNV, small-indel, SV, or complex-path event. */

export interface EvidenceEvent {
  id: string;
  index: number;
  type:
    | 'SNV'
    | 'SMALL_INSERTION'
    | 'SMALL_DELETION'
    | 'COMPLEX_SV_PATH'
    | `SV_${string}`
    | string;
  start: number | null;
  end: number | null;
  length: number;
  ref: string | null;
  alt: string | null;
  normalization: string;
  source_projection: 'variants' | 'small_indels' | 'svs' | 'complex_events' | string;
  negative_evidence_rule:
    | 'callable_base_allele'
    | 'same_fragment_callable_reference_adjacency'
    | 'same_fragment_callable_deleted_span_with_flanks'
    | 'support_only_no_negative_inference'
    | string;
  assessability: 'REFERENCE_AND_ALTERNATE' | 'ALTERNATE_SUPPORT_ONLY' | string;
  component_event_ids: string[];
  supporting_molecule_ids: string[];
  evidence_counts: EvidenceCounts;
}

/** Sparse molecule/event evidence; an absent pair is explicitly NOT_CALLABLE. */

export interface EvidenceObservation {
  id: string;
  molecule_id: string;
  event_id: string;
  alignment_id: string;
  state: ObservationState;
  observed_allele: string | null;
  base_quality: number | null;
  mapping_quality: number;
  strand: '+' | '-' | string;
  evidence_source: string;
  read_position: number | null;
}

/** Columnar page used to keep observation keys and browser allocations bounded. */

export interface EvidenceObservationPage {
  index: number;
  offset: number;
  count: number;
  columns: {
    molecule_id: string[];
    event_id: string[];
    alignment_id: string[];
    state: ObservationState[];
    observed_allele: Array<string | null>;
    base_quality: Array<number | null>;
    mapping_quality: number[];
    strand: string[];
    evidence_source: string[];
    read_position: Array<number | null>;
  };
}

export interface EvidenceEncoding {
  layout: 'paged_columnar_molecule_event' | string;
  scope: 'snv_indel_sv_complex_evidence_rc2' | string;
  observation_storage: 'embedded_columnar_pages' | 'remote_http_pages' | string;
  /** Server transport projection; canonical CLI JSON keeps embedded pages. */
  observation_page_endpoint?: string;
  /** Exact bounded server-side filters over every immutable observation page. */
  observation_search_endpoint?: string;
  missing_pair_state: 'NOT_CALLABLE';
  phase_molecule_policy: 'evidence_eligible_only' | string;
  phase_molecule_reference: 'molecules[].index' | string;
  phase_null_model: 'independent_marginals_within_jointly_callable' | string;
  observation_limit: number;
  observation_count: number;
  observation_page_size: number;
  observation_page_count: number;
  phase_link_limit: number;
}

/** Callable-aware pairwise physical linkage projection. */

export interface PhaseLink {
  id: string;
  event_a_id: string;
  event_b_id: string;
  assessability: 'COMPLETE_FOR_BOTH_EVENTS' | 'SUPPORT_CONDITIONED' | string;
  jointly_callable: number;
  jointly_uncertain: number;
  both_alternate: number;
  a_alternate_b_absent: number;
  a_absent_b_alternate: number;
  neither_alternate: number;
  co_alternate_fraction: number;
  co_alternate_ci95_low: number;
  co_alternate_ci95_high: number;
  expected_co_alternate_fraction: number;
  linkage_delta: number;
  /** Stable references into molecules[].index; avoids repeating long IDs per pair. */
  supporting_molecule_indices: number[];
  uncertain_molecule_indices: number[];
  qc_flags: string[];
}
