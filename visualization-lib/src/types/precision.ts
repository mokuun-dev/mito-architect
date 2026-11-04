import type { ProteinStructureMapping } from './variants';

/** Functional effect of a variant in one rCRS feature. */
export interface PrecisionEffect {
  gene?: string;
  biotype: string;
  consequence: string;
  hgvs_p?: string;
  codons?: string;
  residue?: number;
  grantham?: number;
  structure?: ProteinStructureMapping;
}

export interface PrecisionRepeat {
  homopolymer_run: number;
  unit: number;
  span: number;
}

export interface PrecisionStrands {
  alt_forward: number;
  alt_reverse: number;
  depth_forward: number;
  depth_reverse: number;
}

interface PrecisionCallBase {
  id: string;
  position: number;
  status: 'PASS' | 'FILTERED';
  filters: string[];
  flags: string[];
  fraction: number;
  ci95_low: number;
  ci95_high: number;
  depth: number;
  strands: PrecisionStrands;
  p_value: number;
  q_value: number;
  background_rate: number;
  repeat: PrecisionRepeat;
  effects: PrecisionEffect[];
}

export interface PrecisionSnv extends PrecisionCallBase {
  ref: string;
  alt: string;
  alt_depth: number;
  alt_mean_quality: number;
}

export interface PrecisionIndel extends PrecisionCallBase {
  type: 'deletion' | 'insertion';
  /** Leftmost equivalent placement; position is the same as VCF after its anchor. */
  position_3prime: number;
  length: number;
  allele: string;
  vcf_position: number;
  vcf_ref: string;
  vcf_alt: string;
  hgvs: string;
  raw_alt_depth: number;
  effective_alt_depth: number;
  effective_depth: number;
  lod: number;
  rescored: boolean;
  supporting_molecule_count: number;
  supporting_molecules: string[];
}

export interface PrecisionDeletionImpact {
  genes_lost: string[];
  genes_partial: string[];
  trnas_affected: number;
  fusion?: string;
  heavy_strand_origin_lost: boolean;
  light_strand_origin_lost: boolean;
  major_arc: boolean;
  contact_zone: boolean;
}

export interface PrecisionDeletion {
  id: string;
  /** First and last deleted base under the HGVS 3' rule. */
  start: number;
  end: number;
  left_aligned_start: number;
  left_aligned_end: number;
  length: number;
  microhomology: string;
  status: 'PASS' | 'FILTERED';
  filters: string[];
  flags: string[];
  fraction: number | null;
  ci95_low: number;
  ci95_high: number;
  junction_molecules: number;
  junction_forward: number;
  junction_reverse: number;
  wild_type_left: number;
  wild_type_right: number;
  breakpoint_variants: number;
  depth_ratio: number | null;
  p_value: number;
  q_value: number;
  known_common_deletion: boolean;
  impact: PrecisionDeletionImpact;
  supporting_molecule_count: number;
  supporting_molecules: string[];
}

export interface PrecisionHaplotype {
  id: number;
  fraction: number;
  molecules: number;
  defining_variants: string[];
  haplogroup?: string;
  haplogroup_quality?: number;
}

export interface PrecisionDiagnostics {
  molecules_used: number;
  molecules_excluded_divergent: number;
  reference_positions_covered: number;
  median_depth: number;
  mean_substitution_rate: number;
  median_molecule_divergence: number;
  divergence_limit: number;
  snv_hypotheses: number;
  indel_candidates: number;
  indels_rescored: number;
  /** 99th percentile of one-base indel rates outside repeats in this sample. */
  indel_noise_p99: number;
  junction_molecules: number;
  singleton_junctions: number;
  homopolymer_deletion_rate: Array<{ run: number; rate: number }>;
  homopolymer_insertion_rate: Array<{ run: number; rate: number }>;
}

/** Statistically filtered calls, calibrated on the sample's own error profile. */
export interface PrecisionCalls {
  schema_version: string;
  assessed: boolean;
  reason?: string;
  method?: {
    false_discovery_rate: number;
    snv_model: string;
    indel_model: string;
    deletion_model: string;
  };
  summary?: { snvs_pass: number; indels_pass: number; deletions_pass: number };
  diagnostics?: PrecisionDiagnostics;
  consensus_haplogroup?: {
    best: string;
    quality: number;
    observed_markers: number;
    candidates: Array<{ name: string; score: number }>;
  };
  haplotypes?: PrecisionHaplotype[];
  /** Everything outside these 1-based ranges is NOT_CALLABLE. */
  callable_ranges?: Array<[number, number]>;
  snvs?: PrecisionSnv[];
  indels?: PrecisionIndel[];
  deletions?: PrecisionDeletion[];
}
