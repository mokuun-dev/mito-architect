import type { SnpCall } from './variants';

export interface CoverageMetrics {
  mean_depth: number;
  pct_sites_gt20x?: number;
  pct_bins_gt20x?: number;
  max_depth: number;
  mapping_quality_histogram: Array<{ mapq: number; count: number }>;
}

export interface ReadFeature {
  id: string;
  length: number;
  mean_quality: number;
  numt_score: number;
  filtered_numt: boolean;
  numt_evidence?: string[];
  cluster_id: number;
  mapping_quality?: number;
  flags?: number;
  reference_name?: string;
  aux_tags?: Record<string, string>;
  outlier?: boolean;
  snps: SnpCall[];
  haplogroup_markers?: string[];
  haplogroup_range_known?: boolean;
  haplogroup_callable_ranges?: Array<{ start: number; end: number }>;
  sv_ids: string[];
  complex_event_ids?: string[];
}

export interface ClusterSummary {
  id: number;
  label: string;
  haplogroup?: string;
  size: number;
  consensus_haplotype: string;
  reads: string[];
  sv_signature: Array<{ sv_id: string; support: number }>;
  complex_event_signature?: Array<{ event_id: string; support: number }>;
  haplogroup_assignment?: {
    resource: string;
    quality: number;
    contamination_warning: boolean;
    observed_markers?: string[];
    callable_ranges?: Array<{ start: number; end: number }>;
    candidates: Array<{
      name: string;
      score: number;
      matched: string[];
      missing: string[];
      extra: string[];
    }>;
  };
}
