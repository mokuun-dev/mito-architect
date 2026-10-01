import type { GeneAnnotation, MitoAnalysisData } from '@mito-architect/visualization-lib';
import { circularVariants } from '@mito-architect/visualization-lib';
import { useQuery } from '@tanstack/react-query';
import { useCallback, useEffect, useRef } from 'react';
import { useParams, useSearchParams } from 'react-router-dom';
import AnalysisModulesPanel from '../components/AnalysisModulesPanel';
import CircosWrapper from '../components/CircosWrapper';
import { GenomeLegend } from '../components/genome/GenomeLegend';
import { GenomeToolbar } from '../components/genome/GenomeToolbar';
import { LinearGenomeView } from '../components/genome/LinearGenomeView';
import { SelectionInspector } from '../components/genome/SelectionInspector';
import { isAnalysisModule, ModuleNavigation } from '../components/workspace/ModuleNavigation';
import { SampleHeader } from '../components/workspace/SampleHeader';
import { getResult } from '../lib/api';
import { summarizeCoverage } from '../lib/coverageModel';
import { demoData } from '../lib/demoData';
import { useMitoStore } from '../lib/store';

export default function ResultPage() {
  const { jobId } = useParams<{ jobId: string }>();
  const isDemo = jobId === 'demo';
  const activeLayers = useMitoStore((state) => state.activeLayers);
  const activeModule = useMitoStore((state) => state.activeModule);
  const selectedRegion = useMitoStore((state) => state.selectedRegion);
  const setActiveModule = useMitoStore((state) => state.setActiveModule);
  const setData = useMitoStore((state) => state.setData);
  const setSelectedCluster = useMitoStore((state) => state.setSelectedCluster);
  const setSelectedSvId = useMitoStore((state) => state.setSelectedSvId);
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);
  const setSelectedGene = useMitoStore((state) => state.setSelectedGene);
  const setSelectedRegion = useMitoStore((state) => state.setSelectedRegion);
  const [searchParams, setSearchParams] = useSearchParams();
  const lastUrlModule = useRef<string | null | undefined>(undefined);
  const suppressUrlWrite = useRef(false);

  const result = useQuery({
    queryKey: ['result', jobId],
    queryFn: ({ signal }) => getResult(jobId!, signal),
    enabled: Boolean(jobId) && !isDemo
  });

  useEffect(() => { setData(isDemo ? demoData : result.data); }, [isDemo, result.data, setData]);

  // Read an incoming URL before writing store changes back to browser history.
  useEffect(() => {
    const urlModule = searchParams.get('module');
    if (urlModule !== lastUrlModule.current) {
      lastUrlModule.current = urlModule;
      suppressUrlWrite.current = true;
      setActiveModule(isAnalysisModule(urlModule) ? urlModule : 'overview');
    }
  }, [searchParams, setActiveModule]);
  useEffect(() => {
    if (suppressUrlWrite.current) { suppressUrlWrite.current = false; return; }
    if (searchParams.get('module') === activeModule) return;
    const next = new URLSearchParams(searchParams);
    next.set('module', activeModule);
    lastUrlModule.current = activeModule;
    setSearchParams(next);
  }, [activeModule, searchParams, setSearchParams]);

  const handleClusterSelect = useCallback((id: number) => setSelectedCluster(id), [setSelectedCluster]);
  const handleSvSelect = useCallback((id: string) => setSelectedSvId(id), [setSelectedSvId]);
  const handleEventSelect = useCallback((id: string) => setSelectedEventId(id), [setSelectedEventId]);
  const handleGeneSelect = useCallback((gene: GeneAnnotation) => {
    setSelectedGene(gene.name);
    setSelectedRegion({ start: gene.start, end: gene.end, label: gene.name });
  }, [setSelectedGene, setSelectedRegion]);

  if (!isDemo && result.isLoading) return <ResultSkeleton />;
  if ((!isDemo && (result.isError || !result.data)) || !jobId) {
    return <main id="main-content" className="page-shell"><div role="alert" className="workspace-section border-coral/50 text-coral"><h1 className="text-lg font-semibold">Результат недоступен</h1><p className="mt-2 text-sm">{(result.error as Error | undefined)?.message ?? 'Не удалось загрузить результат.'}</p></div></main>;
  }

  const data = isDemo ? demoData : result.data!;
  return <main id="main-content" className="workspace-page">
    <div className="workspace-shell">
      <aside className="workspace-navigation"><ModuleNavigation active={activeModule} onSelect={setActiveModule} /></aside>
      <div className="workspace-main min-w-0">
        <SampleHeader data={data} jobId={jobId} isDemo={isDemo} />
        {activeModule === 'overview' ? <Overview
          data={data}
          activeLayers={activeLayers}
          selectedRegion={Boolean(selectedRegion)}
          onClusterSelect={handleClusterSelect}
          onSvSelect={handleSvSelect}
          onEventSelect={handleEventSelect}
          onGeneSelect={handleGeneSelect}
          onOpenModule={setActiveModule}
        /> : <AnalysisModulesPanel data={data} jobId={isDemo ? undefined : jobId} />}
      </div>
      <aside className="workspace-inspector"><SelectionInspector data={data} /></aside>
    </div>
  </main>;
}

function Overview({ data, activeLayers, selectedRegion, onClusterSelect, onSvSelect, onEventSelect, onGeneSelect, onOpenModule }: {
  data: MitoAnalysisData;
  activeLayers: ReturnType<typeof useMitoStore.getState>['activeLayers'];
  selectedRegion: boolean;
  onClusterSelect: (id: number) => void;
  onSvSelect: (id: string) => void;
  onEventSelect: (id: string) => void;
  onGeneSelect: (gene: GeneAnnotation) => void;
  onOpenModule: ReturnType<typeof useMitoStore.getState>['setActiveModule'];
}) {
  const coverage = summarizeCoverage(data);
  const variantCount = circularVariants(data).length;
  const moleculeCount = data.filter_stats.evidence_eligible_molecules ?? data.filter_stats.input_molecules;
  return <div className="grid min-w-0 gap-4">
    <section className="workspace-summary" aria-label="Сводка образца">
      <SummaryButton label={data.filter_stats.evidence_eligible_molecules === undefined ? 'Собранные молекулы' : 'Молекулы для evidence'} value={moleculeCount === undefined ? 'не рассчитано' : moleculeCount.toLocaleString()} detail="Отдельно от записей выравнивания" onClick={() => onOpenModule('evidence')} />
      <SummaryButton label="Средняя глубина" value={data.coverage.length ? coverage.meanDepth.toFixed(1) + '×' : 'не рассчитано'} detail="По референсным позициям" onClick={() => onOpenModule('coverage')} />
      <SummaryButton label="Малые варианты" value={variantCount.toLocaleString()} detail="Наблюдаемые события" onClick={() => onOpenModule('variants')} />
      <SummaryButton label="Перестройки" value={data.svs.length.toLocaleString()} detail="Структурные события" onClick={() => onOpenModule('rearrangements')} />
    </section>
    <figure className="workspace-section min-w-0">
      <div className="flex flex-wrap items-start justify-between gap-3">
        <div><h2 className="text-lg font-semibold">Карта митохондриального генома</h2><p className="mt-1 text-sm text-muted">Выберите ген или событие, чтобы увидеть данные справа.</p></div>
        <span className="text-xs text-muted">1–{data.metadata.reference_length.toLocaleString()} bp</span>
      </div>
      <GenomeToolbar genes={data.genes} />
      <CircosWrapper data={data} activeLayers={activeLayers} onClusterSelect={onClusterSelect} onSvSelect={onSvSelect} onEventSelect={onEventSelect} onGeneSelect={onGeneSelect} />
      <GenomeLegend maxDepth={coverage.maxDepth} />
      <figcaption className="mt-3 border-t border-line pt-3 text-xs leading-5 text-muted">Координаты относятся к кольцевому референсу. Доля события определяется по пригодным для оценки молекулам, если такие данные рассчитаны. Exact values remain available in the linked textual projections.</figcaption>
    </figure>
    {selectedRegion && <LinearGenomeView data={data} />}
  </div>;
}

function SummaryButton({ label, value, detail, onClick }: { label: string; value: string; detail: string; onClick: () => void }) {
  return <button type="button" onClick={onClick} className="workspace-summary-item group text-left"><span className="text-xs font-medium text-muted">{label}</span><span className="mt-2 block font-mono text-xl font-semibold tabular-nums text-text">{value}</span><span className="mt-1 block text-[11px] text-muted group-hover:text-aqua">{detail} →</span></button>;
}

function ResultSkeleton() {
  return <main id="main-content" className="workspace-page" aria-busy="true" aria-label="Загрузка результата"><div className="workspace-shell"><div className="workspace-navigation h-96 animate-pulse bg-panel2" /><div className="workspace-main grid gap-4"><div className="h-28 animate-pulse rounded-lg bg-panel2" /><div className="h-24 animate-pulse rounded-lg bg-panel2" /><div className="h-[520px] animate-pulse rounded-lg bg-panel2" /></div><div className="workspace-inspector h-72 animate-pulse bg-panel2" /></div></main>;
}
