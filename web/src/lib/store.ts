import type { LayerName, MitoAnalysisData } from '@mito-architect/visualization-lib';
import { create } from 'zustand';

type LayerMap = Record<LayerName, boolean>;
export type AnalysisModule = 'overview' | 'coverage' | 'variants' | 'evidence' | 'architectures' | 'rearrangements' | 'haplogroups' | 'clinical' | 'protein' | 'method';
export type GenomeRegion = { start: number; end: number; label?: string };

interface MitoState {
  selectedFile?: File;
  selectedEvidenceGraph: boolean;
  jobId?: string;
  data?: MitoAnalysisData;
  activeLayers: LayerMap;
  selectedCluster?: number;
  selectedSvId?: string;
  selectedEventId?: string;
  selectedMoleculeId?: string;
  selectedPhaseId?: string;
  selectedGene?: string;
  selectedRegion?: GenomeRegion;
  activeModule: AnalysisModule;
  geneQuery: string;
  minQuality: number;
  setSelectedFile: (file?: File) => void;
  setSelectedEvidenceGraph: (enabled: boolean) => void;
  setJobId: (jobId?: string) => void;
  setData: (data?: MitoAnalysisData) => void;
  setLayer: (layer: LayerName, active: boolean) => void;
  setSelectedCluster: (clusterId?: number) => void;
  setSelectedSvId: (svId?: string) => void;
  setSelectedEventId: (eventId?: string) => void;
  setSelectedMoleculeId: (moleculeId?: string) => void;
  setSelectedPhaseId: (phaseId?: string) => void;
  setSelectedGene: (gene?: string) => void;
  setSelectedRegion: (region?: GenomeRegion) => void;
  setActiveModule: (module: AnalysisModule) => void;
  clearInspection: () => void;
  setGeneQuery: (query: string) => void;
  setMinQuality: (quality: number) => void;
}

export const useMitoStore = create<MitoState>((set) => ({
  selectedEvidenceGraph: false,
  activeLayers: {
    genes: true,
    coverage: true,
    clusters: false,
    svs: true,
    snps: true
  },
  geneQuery: '',
  minQuality: 0,
  activeModule: 'overview',
  setSelectedFile: (selectedFile) => set({ selectedFile }),
  setSelectedEvidenceGraph: (selectedEvidenceGraph) => set({ selectedEvidenceGraph }),
  setJobId: (jobId) => set({ jobId }),
  setData: (data) => set({
    data,
    selectedCluster: undefined,
    selectedSvId: undefined,
    selectedEventId: undefined,
    selectedMoleculeId: undefined,
    selectedPhaseId: undefined,
    selectedGene: undefined,
    selectedRegion: undefined,
    geneQuery: '',
    minQuality: 0
  }),
  setLayer: (layer, active) =>
    set((state) => ({
      activeLayers: {
        ...state.activeLayers,
        [layer]: active
      }
    })),
  setSelectedCluster: (selectedCluster) => set({ selectedCluster, selectedSvId: undefined, selectedEventId: undefined, selectedMoleculeId: undefined, selectedGene: undefined }),
  setSelectedSvId: (selectedSvId) => set({ selectedSvId, selectedEventId: undefined, selectedCluster: undefined, selectedMoleculeId: undefined, selectedGene: undefined }),
  setSelectedEventId: (selectedEventId) => set({ selectedEventId, selectedSvId: undefined, selectedCluster: undefined, selectedMoleculeId: undefined, selectedGene: undefined }),
  setSelectedMoleculeId: (selectedMoleculeId) => set({ selectedMoleculeId, selectedEventId: undefined, selectedSvId: undefined, selectedGene: undefined }),
  setSelectedPhaseId: (selectedPhaseId) => set({ selectedPhaseId }),
  setSelectedGene: (selectedGene) => set({ selectedGene, selectedEventId: undefined, selectedSvId: undefined, selectedCluster: undefined, selectedMoleculeId: undefined }),
  setSelectedRegion: (selectedRegion) => set({ selectedRegion }),
  setActiveModule: (activeModule) => set({ activeModule }),
  clearInspection: () => set({
    selectedCluster: undefined,
    selectedSvId: undefined,
    selectedEventId: undefined,
    selectedMoleculeId: undefined,
    selectedPhaseId: undefined,
    selectedGene: undefined
  }),
  setGeneQuery: (geneQuery) => set({ geneQuery }),
  setMinQuality: (minQuality) => set({ minQuality })
}));
