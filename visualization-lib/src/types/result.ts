import type { LayerName, MitoMetadata, FilterStats, GeneAnnotation, CoverageBin } from './metadata';
import type { StructuralVariant, ComplexStructuralEvent, AggregateVariant } from './variants';
import type { CoverageMetrics, ReadFeature, ClusterSummary } from './reads';
import type { AlignmentFragment, AssembledMolecule, MoleculeCallability, EvidenceEvent, EvidenceObservation, EvidenceObservationPage, EvidenceEncoding, PhaseLink } from './evidence';
import type { MolecularArchitecture, ArchitectureInferenceSummary } from './architectures';

export interface MitoAnalysisData {
  metadata: MitoMetadata;
  filter_stats: FilterStats;
  coverage_metrics?: CoverageMetrics;
  genes: GeneAnnotation[];
  coverage: CoverageBin[];
  svs: StructuralVariant[];
  complex_events?: ComplexStructuralEvent[];
  variants?: AggregateVariant[];
  clusters: ClusterSummary[];
  reads: ReadFeature[];
  evidence_encoding?: EvidenceEncoding;
  alignments?: AlignmentFragment[];
  molecules?: AssembledMolecule[];
  callability?: MoleculeCallability[];
  events?: EvidenceEvent[];
  observations?: EvidenceObservation[];
  observation_pages?: EvidenceObservationPage[];
  phase_links?: PhaseLink[];
  architecture_inference?: ArchitectureInferenceSummary;
  architectures?: MolecularArchitecture[];
}

export interface MitoCircosOptions {
  width?: number;
  height?: number;
  theme?: 'dark' | 'light';
  showControls?: boolean;
  activeLayers?: Partial<Record<LayerName, boolean>>;
  palette?: string[];
}

export interface LayerState {
  genes: boolean;
  coverage: boolean;
  clusters: boolean;
  svs: boolean;
  snps: boolean;
}
