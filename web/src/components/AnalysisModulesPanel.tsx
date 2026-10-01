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
import VariantEvidenceTable from './VariantEvidenceTable';
import { summarizeCoverage } from '../lib/coverageModel';
import { useMitoStore } from '../lib/store';
import { MODULES } from './workspace/ModuleNavigation';
import { MethodSources } from './workspace/MethodSources';

export default function AnalysisModulesPanel({ data, jobId }: { data: MitoAnalysisData; jobId?: string }) {
  const active = useMitoStore((state) => state.activeModule);
  const passedReads = useMemo(() => data.reads.filter((read) => !read.filtered_numt), [data.reads]);
  const variants = useMemo(() => summarizeVariants(passedReads, data), [passedReads, data]);
  const coverage = useMemo(() => summarizeCoverage(data), [data]);
  const label = MODULES.find((module) => module.id === active)?.label ?? 'Анализ';
  return <section className="workspace-section min-w-0" aria-labelledby="analysis-module-title">
    <div className="border-b border-line pb-4">
      <div className="text-xs text-muted">Раздел анализа</div>
      <h2 id="analysis-module-title" className="mt-1 text-lg font-semibold">{label}</h2>
    </div>
    <div className="min-w-0 pt-5">
      {active === 'coverage' && <CoveragePanel data={data} metrics={coverage} />}
      {active === 'variants' && <VariantEvidenceTable data={data} />}
      {active === 'evidence' && <EvidencePanel data={data} jobId={jobId} />}
      {active === 'architectures' && <div className="grid gap-5"><ArchitecturePanel data={data} /><details className="border-t border-line pt-4"><summary className="cursor-pointer text-sm font-semibold text-aqua">Базовые группы по признакам</summary><p className="mt-2 text-sm text-muted">Эти группы отличаются от кандидатов молекулярных архитектур и не обозначают клеточные клоны.</p><div className="mt-4"><ClusterMutationPanel data={data} /></div></details></div>}
      {active === 'rearrangements' && <RearrangementsPanel data={data} />}
      {active === 'haplogroups' && <HaplogroupPanel data={data} />}
      {active === 'clinical' && <ClinicalPanel variants={variants} />}
      {active === 'protein' && <ProteinPanel variants={variants} />}
      {active === 'method' && <MethodSources data={data} />}
    </div>
  </section>;
}
