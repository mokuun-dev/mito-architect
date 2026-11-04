import type { MitoAnalysisData } from '@mito-architect/visualization-lib';

export interface TruthScore {
  tp: number;
  fp: number;
  fn: number;
  precision: number;
  recall: number;
  f1: number;
  fraction_mae: number | null;
}

export interface DemoSampleEntry {
  id: string;
  title: string;
  summary: string;
  shows: string[];
  accessions: string[];
  simulated: boolean;
  platform: string;
  file: string;
  sha256: string;
  haplogroup?: string | null;
  haplogroup_quality?: number | null;
  haplotypes: number;
  median_depth: number;
  calls: { snvs_pass: number; indels_pass: number; deletions_pass: number };
  truth_check?: { snv: TruthScore; snv_easy: TruthScore; indel: TruthScore; indel_easy: TruthScore; large_deletion: TruthScore };
  source?: string;
}

export interface DemoIndex {
  version: number;
  sources: { genomes: string; reads: string; spiked_variants: string };
  samples: DemoSampleEntry[];
}

const PLATFORM_LABELS: Record<string, string> = {
  hifi: 'Simulated PacBio HiFi',
  ont_r10: 'Simulated nanopore R10',
  ont_r9: 'Simulated nanopore R9',
  ont_public: 'Public nanopore run'
};

export function platformLabel(platform: string): string {
  return PLATFORM_LABELS[platform] ?? platform;
}

async function fetchJson<T>(path: string, signal?: AbortSignal): Promise<T> {
  const response = await fetch(path, { signal });
  if (!response.ok) {
    throw new Error(`The demo panel file ${path} could not be loaded (HTTP ${response.status}). Generate it with scripts/mtsim/build_demo_panel.py.`);
  }
  return (await response.json()) as T;
}

export function fetchDemoIndex(signal?: AbortSignal): Promise<DemoIndex> {
  return fetchJson<DemoIndex>('/demo/index.json', signal);
}

/** Requested sample if the panel has it, otherwise the first one. */
export function resolveDemoSample(index: DemoIndex, requested: string | null): DemoSampleEntry | undefined {
  return index.samples.find((sample) => sample.id === requested) ?? index.samples[0];
}

export function fetchDemoSample(entry: DemoSampleEntry, signal?: AbortSignal): Promise<MitoAnalysisData> {
  // The name comes from the panel index, never directly from the URL.
  return fetchJson<MitoAnalysisData>(`/demo/${encodeURIComponent(entry.file)}`, signal);
}
