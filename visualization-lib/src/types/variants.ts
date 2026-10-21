import type { EvidenceEvent } from './evidence';

export interface StructuralVariant {
  id: string;
  /** Schema 0.6 unified evidence-graph reference. */
  event_id?: string;
  type: string;
  start: number;
  end: number;
  length: number;
  known_event: boolean;
  supporting_reads: string[];
  /** @deprecated Use evidence_sources; retained for result-schema 0.4 readers. */
  evidence_source?: 'cigar' | 'split_alignment' | 'development_tag' | 'combined' | string;
  /** Sorted unique provenance categories merged under the canonical event ID. */
  evidence_sources?: Array<'cigar' | 'split_alignment' | 'development_tag' | string>;
  /** @deprecated Use orientations; retained for result-schema 0.4 readers. */
  orientation?: string;
  /** Sorted unique strand transitions observed for split-alignment evidence. */
  orientations?: string[];
  segment_count?: number;
  annotation?: ClinicalAnnotation;
}

/** Strand-invariant ordered path of two or more junctions observed on one molecule. */

export interface ComplexStructuralEvent {
  id: string;
  /** Schema 0.6 unified evidence-graph reference. */
  event_id?: string;
  junction_count: number;
  segment_count: number;
  canonicalization: 'strand_invariant_path' | string;
  /** Canonical SV edge IDs in molecule traversal order, normalized against path reversal. */
  junction_ids: string[];
  /** Strand transition paired by index with each canonical junction ID. */
  junction_orientations: string[];
  supporting_reads: string[];
}

/** Optional protein-structure mapping for non-synonymous mtDNA variants. */

export interface ProteinStructureMapping {
  structure_id?: string;
  chain?: string;
  residue_index?: number;
  residue_label?: string;
  complex?: string;
}

/** Single-nucleotide variant call and optional annotation payload. */

export interface SnpCall {
  position: number;
  ref: string;
  alt: string;
  gene?: string;
  consequence?: string;
  protein?: string;
  residue?: string;
  annotation?: ClinicalAnnotation;
  structure?: ProteinStructureMapping;
}

/** Unified schema 0.6 SNV/small-indel projection derived from the evidence graph. */

export interface AggregateVariant extends SnpCall {
  /** Schema 0.6 normalized event reference. */
  event_id?: string;
  type?: 'SNV' | 'SMALL_INSERTION' | 'SMALL_DELETION' | string;
  start?: number;
  end?: number;
  length?: number;
  normalization?: string;
  negative_evidence_rule?: EvidenceEvent['negative_evidence_rule'];
  assessability?: EvidenceEvent['assessability'];
  vcf_position?: number;
  vcf_representable?: boolean;
  alt_depth: number;
  ref_depth: number;
  other_depth: number;
  event_absent_depth?: number;
  low_quality_depth?: number;
  conflict_depth?: number;
  callable_depth: number;
  heteroplasmy: number;
  ci95_low: number;
  ci95_high: number;
  supporting_reads: string[];
  supporting_molecule_ids?: string[];
  filter_status?: 'NOT_CALIBRATED' | string;
  qc_flags?: string[];
  numt_assessability?: 'ASSESSABLE' | 'NOT_ASSESSABLE' | string;
  multi_allelic?: boolean;
  homopolymer_context?: {
    reference_base: string | null;
    run_length: number;
  };
  molecule_support?: {
    alternate: number;
    reference: number;
    other: number;
    callable: number;
  };
  strand_support?: {
    alt_forward: number;
    alt_reverse: number;
    ref_forward: number;
    ref_reverse: number;
    other_forward: number;
    other_reverse: number;
  };
  /** Absolute alternate-vs-reference forward-strand fraction difference; observational only. */
  strand_bias_delta?: number | null;
  allele_quality?: {
    alternate: AlleleQualitySummary;
    reference: AlleleQualitySummary;
    other: AlleleQualitySummary;
  };
  mapping_quality?: {
    alternate: NumericQualitySummary;
    reference: NumericQualitySummary;
    other: NumericQualitySummary;
  };
  read_position?: {
    definition: 'normalized_center_proximity' | string;
    alternate_mean: number | null;
    reference_mean: number | null;
    other_mean: number | null;
    bias_delta: number | null;
  };
}

export interface NumericQualitySummary {
  count: number;
  mean: number | null;
  min: number | null;
  max: number | null;
}

export interface AlleleQualitySummary {
  count: number;
  mean_phred: number | null;
  min_phred: number | null;
  max_phred: number | null;
}

/** One source-preserving clinical record; empty optional source fields are not inferred. */

export interface ClinicalAssertion {
  source: string;
  assertion_id: string;
  allele_id: string;
  disease: string;
  clinical_significance: string;
  normalized_significance: string;
  review_status: string;
  assertion_date: string;
  source_url: string;
  references: string[];
  resource_version: string;
  retrieved_at: string;
}

/** Deterministic summary plus lossless source assertions from the local clinical cache. */

export interface ClinicalAnnotation {
  schema_version?: string;
  phenotype?: string;
  /** @deprecated Prefer consensus_significance and inspect assertions/conflict_status. */
  pathogenicity?: string;
  consensus_significance?: string;
  conflict_status?: 'single_assertion' | 'consistent' | 'conflicting' | string;
  assertions?: ClinicalAssertion[];
  references?: string[];
  source?: 'MITOMAP' | 'ClinVar' | 'local-cache' | string;
  sources?: string[];
  clinvar_allele_id?: string;
  mitomap_url?: string;
}

/** Summary statistics for the coverage layer. */
