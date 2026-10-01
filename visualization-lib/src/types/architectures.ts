
export interface ArchitectureAssignment {
  status: 'ASSIGNED' | 'AMBIGUOUS' | 'UNASSIGNED' | 'INELIGIBLE' | string;
  architecture_id: string | null;
  distance: number;
  score: number;
  overlap_fraction: number;
  confidence: number;
  candidate_architecture_ids: string[];
}

export interface ArchitectureEventProfile {
  event_id: string;
  alternate: number;
  absent: number;
  uncertain: number;
  not_callable: number;
  alternate_fraction: number;
  consensus_state: 'ALTERNATE' | 'ABSENT' | 'OPTIONAL' | string;
}

export interface MolecularArchitecture {
  architecture_id: string;
  status: 'CANDIDATE' | string;
  defining_event_signature: string[];
  seed_event_signature: string[];
  optional_event_ids: string[];
  molecule_count: number;
  estimated_fraction: number;
  confidence_interval: {
    method: 'wilson_score_95' | string;
    level: number;
    low: number;
    high: number;
  };
  median_callable_fraction: number;
  assignment_confidence: {
    mean: number;
    minimum: number;
  };
  cluster_stability: {
    status: 'RESAMPLING_ESTIMATE' | 'NOT_ESTIMATED' | string;
    value: number;
    method: string;
    recovery_rate: number;
    abundance_standard_deviation: number;
    replicates: number;
    seed: number;
  };
  representative_molecule_ids: string[];
  member_molecule_indices: number[];
  ambiguous_molecule_ids: string[];
  unassigned_molecule_count: number;
  haplogroup_evidence: {
    status: string;
    resource: string;
    best: string;
    quality: number;
    callable_ranges_known: boolean;
    contamination_warning: boolean;
    observed_markers: string[];
    callable_ranges: Array<{ start: number; end: number }>;
  };
  structural_path_evidence: {
    defining_event_ids: string[];
    optional_event_ids: string[];
  };
  event_profile: ArchitectureEventProfile[];
  method: {
    name: string;
    positive_evidence_seeded: boolean;
    not_callable_is_reference: boolean;
    minimum_molecules: number;
    maximum_distance: number;
    minimum_overlap_fraction: number;
    ambiguity_margin: number;
  };
}

export interface ArchitectureInferenceSummary {
  schema_version: string;
  status: string;
  interpretation: 'candidate_molecular_architectures_not_cell_clones' | string;
  method: string;
  eligible_molecules: number;
  assigned_molecules: number;
  ambiguous_molecules: number;
  unassigned_molecules: number;
  candidate_count: number;
  qc_flags: string[];
}

/** Full analysis payload shared by core, CLI, server, report, and UI. */
