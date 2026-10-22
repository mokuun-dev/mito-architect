
export type LayerName = 'genes' | 'coverage' | 'clusters' | 'svs' | 'snps';

/** Run-level metadata emitted by the native analysis core. */

export interface MitoMetadata {
  schema_version?: string;
  sv_event_schema_version?: string;
  complex_sv_event_schema_version?: string;
  clinical_annotation_schema_version?: string;
  engine_version: string;
  sample: string;
  input_path?: string;
  reference_path?: string;
  reference_accession?: string;
  reference_length: number;
  threads?: number;
  algorithm_notes?: string[];
  calling_parameters?: {
    min_mapping_quality: number;
    min_base_quality: number;
    excluded_snp_flags: number;
    numt_threshold?: number;
    development_tags_enabled?: boolean;
    max_evidence_observations?: number;
    max_phase_links?: number;
    evidence_page_size?: number;
    min_architecture_molecules?: number;
    max_candidate_architectures?: number;
    architecture_max_distance?: number;
    architecture_ambiguity_margin?: number;
    architecture_min_overlap_fraction?: number;
    architecture_consensus_fraction?: number;
    architecture_optional_fraction?: number;
    architecture_stability_replicates?: number;
    architecture_seed?: number;
    molecule_id_tag?: string;
    umi_tag?: string;
    duplex_tag?: string;
  };
  resources?: AnalysisResource[];
}

export interface AnalysisResource {
  name: string;
  version: string;
  path: string;
  sha256: string;
  source: string;
  license: string;
  retrieved: string;
}

/** NUMT and input-read accounting. */

export interface FilterStats {
  input_reads: number;
  passed_reads: number;
  numt_filtered_reads: number;
  numt_threshold: number;
  input_alignment_records?: number;
  input_molecules?: number;
  evidence_eligible_molecules?: number;
  ambiguous_molecules?: number;
  numt_assessment?: {
    mode: 'competitive_alignment' | 'mt_only_or_unknown' | 'unaligned_fastq' | string;
    nuclear_contigs_present: boolean;
    specificity_assessable: boolean;
    provenance_status?: 'VERIFIED' | 'NOT_VERIFIED' | 'NOT_APPLICABLE' | string;
    provenance_reason?: string;
  };
}

/** rCRS gene interval used by circular tracks and gene search. */

export interface GeneAnnotation {
  name: string;
  start: number;
  end: number;
  strand: '+' | '-' | string;
  biotype: string;
}

/** Binned coverage depth over inclusive rCRS coordinates. */

export interface CoverageBin {
  start: number;
  end: number;
  depth: number;
}

/** Large deletion/insertion call with read-level support. */
