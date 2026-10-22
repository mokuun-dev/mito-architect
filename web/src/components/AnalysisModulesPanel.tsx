import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { useMemo } from 'react';
import ClusterMutationPanel from './ClusterMutationPanel';
import { ArchitecturePanel } from './analysis/ArchitecturePanel';
import { ClinicalPanel } from './analysis/ClinicalPanel';
import { CoveragePanel } from './analysis/CoveragePanel';
import { EvidencePanel } from './analysis/EvidencePanel';
import { HaplogroupPanel } from './analysis/HaplogroupPanel';
import { ProteinPanel } from './analysis/ProteinPanel';
import { RearrangementsPanel } from './analysis/RearrangementsPanel';
import { summarizeVariants } from './analysis/summaries';
import { summarizeCoverage } from '../lib/coverageModel';
import { useMitoStore, type AnalysisModule } from '../lib/store';
import { MODULES } from './workspace/ModuleNavigation';
import { MethodSources } from './workspace/MethodSources';

const DESCRIPTIONS: Partial<Record<AnalysisModule, string>> = {
  coverage: 'Read depth across the reference and how it is distributed.',
  evidence: 'Molecules, per-position observations, and phase links behind each event.',
  architectures: 'Candidate molecular architectures inferred from co-occurring events.',
  rearrangements: 'Structural events, complex junction paths, and their supporting molecules.',
  haplogroups: 'Haplogroup assignments per feature cluster, with ranked alternatives.',
  clinical: 'Source-preserving clinical assertions for observed variants.',
  protein: 'Protein-level context and 3D structure mapping for coding variants.',
  method: 'Engine, parameters, and versioned resources used for this result.'
};

/** Detail modules other than the overview; the overview has its own layout. */
export default function AnalysisModulesPanel({ data, jobId }: { data: MitoAnalysisData; jobId?: string }) {
  const active = useMitoStore((state) => state.activeModule);
  const passedReads = useMemo(() => data.reads.filter((read) => !read.filtered_numt), [data.reads]);
  const variants = useMemo(() => summarizeVariants(passedReads, data), [passedReads, data]);
  const coverage = useMemo(() => summarizeCoverage(data), [data]);
  const label = MODULES.find((module) => module.id === active)?.label ?? 'Analysis';
  return (
    <section className="card" aria-labelledby="analysis-module-title">
      <div className="card-header">
        <div>
          <h2 id="analysis-module-title" className="card-title">{label}</h2>
          {DESCRIPTIONS[active] && <p className="card-subtitle">{DESCRIPTIONS[active]}</p>}
        </div>
      </div>
      <div className="min-w-0 p-5">
        {active === 'coverage' && <CoveragePanel data={data} metrics={coverage} />}
        {active === 'evidence' && <EvidencePanel data={data} jobId={jobId} />}
        {active === 'architectures' && (
          <div className="grid gap-5">
            <ArchitecturePanel data={data} />
            <details className="rounded-lg border border-line">
              <summary className="cursor-pointer px-4 py-3 text-sm font-semibold">Baseline feature clusters</summary>
              <div className="border-t border-line p-4">
                <p className="mb-4 text-[13px] text-muted">These clusters differ from the candidate molecular architectures and do not denote cellular clones.</p>
                <ClusterMutationPanel data={data} />
              </div>
            </details>
          </div>
        )}
        {active === 'rearrangements' && <RearrangementsPanel data={data} />}
        {active === 'haplogroups' && <HaplogroupPanel data={data} />}
        {active === 'clinical' && <ClinicalPanel variants={variants} />}
        {active === 'protein' && <ProteinPanel variants={variants} />}
        {active === 'method' && <MethodSources data={data} />}
      </div>
    </section>
  );
}
