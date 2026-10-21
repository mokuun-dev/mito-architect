import { MitoCircos } from './MitoCircos';

export { MitoCircos };
export { circularVariants, legacyEventId } from './circular/dataAdapter';
export { inclusiveIntervalToSpans, positionToAngle } from './circular/coordinates';
export type {
  AggregateVariant,
  AlignmentCallability,
  AlignmentFragment,
  AlleleQualitySummary,
  ArchitectureAssignment,
  ArchitectureEventProfile,
  ArchitectureInferenceSummary,
  AssembledMolecule,
  ClinicalAnnotation,
  ClinicalAssertion,
  ClusterSummary,
  ComplexStructuralEvent,
  CoverageBin,
  CoverageMetrics,
  FilterStats,
  EvidenceEncoding,
  EvidenceCounts,
  EvidenceEvent,
  EvidenceObservation,
  GeneAnnotation,
  LayerName,
  LayerState,
  MitoAnalysisData,
  MitoCircosOptions,
  MitoMetadata,
  MoleculeCallability,
  MolecularArchitecture,
  ObservationState,
  PhaseLink,
  ProteinStructureMapping,
  ReadFeature,
  ReferenceRange,
  SnpCall,
  StructuralVariant
} from './types';

declare global {
  interface Window {
    MitoCircos?: typeof MitoCircos;
  }
}

if (typeof window !== 'undefined') {
  window.MitoCircos = MitoCircos;
}
