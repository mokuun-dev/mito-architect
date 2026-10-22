import type { GeneAnnotation, MitoAnalysisData } from '@mito-architect/visualization-lib';
import { useQuery } from '@tanstack/react-query';
import { Orbit, Rows3 } from 'lucide-react';
import { useCallback, useEffect, useMemo, useRef } from 'react';
import { useParams, useSearchParams } from 'react-router-dom';
import AnalysisModulesPanel from '../components/AnalysisModulesPanel';
import CircosWrapper from '../components/CircosWrapper';
import { EventFiltersBar } from '../components/events/EventFiltersBar';
import { EventsTable } from '../components/events/EventsTable';
import { GenomeLegend } from '../components/genome/GenomeLegend';
import { GenomeView } from '../components/genome/GenomeView';
import { DetailDrawer } from '../components/inspector/DetailDrawer';
import { KpiStrip } from '../components/workspace/KpiStrip';
import { ModuleNavigation, parseAnalysisModule } from '../components/workspace/ModuleNavigation';
import { SampleHeader } from '../components/workspace/SampleHeader';
import { getResult } from '../lib/api';
import { demoData } from '../lib/demoData';
import { buildEventRows, filterEventRows, sortEventRows, type EventRow } from '../lib/eventRows';
import { useMitoStore, type OverviewView } from '../lib/store';

export default function ResultPage() {
  const { jobId } = useParams<{ jobId: string }>();
  const isDemo = jobId === 'demo';
  const activeModule = useMitoStore((state) => state.activeModule);
  const overviewView = useMitoStore((state) => state.overviewView);
  const filters = useMitoStore((state) => state.eventFilters);
  const selectedEventId = useMitoStore((state) => state.selectedEventId);
  const selectedSvId = useMitoStore((state) => state.selectedSvId);
  const selectedGene = useMitoStore((state) => state.selectedGene);
  const setActiveModule = useMitoStore((state) => state.setActiveModule);
  const setOverviewView = useMitoStore((state) => state.setOverviewView);
  const setData = useMitoStore((state) => state.setData);
  const setSelectedSvId = useMitoStore((state) => state.setSelectedSvId);
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);
  const setSelectedGene = useMitoStore((state) => state.setSelectedGene);
  const setSelectedRegion = useMitoStore((state) => state.setSelectedRegion);

  const result = useQuery({
    queryKey: ['result', jobId],
    queryFn: ({ signal }) => getResult(jobId!, signal),
    enabled: Boolean(jobId) && !isDemo
  });
  const data: MitoAnalysisData | undefined = isDemo ? demoData : result.data;

  useEffect(() => { setData(data); }, [data, setData]);
  useUrlState(Boolean(data));

  const rows = useMemo(() => (data ? buildEventRows(data) : []), [data]);
  const filtered = useMemo(() => filterEventRows(rows, filters), [filters, rows]);
  const sequence = useMemo(() => sortEventRows(filtered, 'position', false), [filtered]);

  const selectRow = useCallback((row: EventRow) => (row.kind === 'sv' ? setSelectedSvId(row.id) : setSelectedEventId(row.id)), [setSelectedEventId, setSelectedSvId]);
  const selectEvent = useCallback((id: string) => setSelectedEventId(id), [setSelectedEventId]);
  const selectSv = useCallback((id: string) => setSelectedSvId(id), [setSelectedSvId]);
  const selectGene = useCallback((gene: GeneAnnotation) => setSelectedGene(gene.name), [setSelectedGene]);
  const showInGenome = useCallback((start: number, end: number, label?: string) => {
    const length = data?.metadata.reference_length ?? 16569;
    const from = Math.max(1, Math.round(start));
    const to = Math.min(length, Math.round(end));
    setSelectedRegion(from <= 1 && to >= length ? undefined : { start: from, end: to, label });
    setOverviewView('genome');
  }, [data, setOverviewView, setSelectedRegion]);

  if (!isDemo && result.isLoading) return <ResultSkeleton />;
  if (!data || !jobId) {
    return (
      <main id="main-content" className="page-shell">
        <div role="alert" className="card border-coral/40 p-6">
          <h1 className="text-lg font-semibold">Result unavailable</h1>
          <p className="mt-2 text-sm text-muted">{(result.error as Error | undefined)?.message ?? 'The result could not be loaded.'}</p>
        </div>
      </main>
    );
  }

  const overview = activeModule === 'overview';
  const drawerOpen = overview && Boolean(selectedEventId || selectedSvId || selectedGene);
  const hasDenominators = rows.some((row) => row.callable !== undefined);
  const selectedId = selectedEventId ?? selectedSvId;

  return (
    <main id="main-content" className={'mx-auto max-w-[1600px] px-4 pb-16 transition-[padding] sm:px-6 ' + (drawerOpen ? '2xl:pr-[464px]' : '')}>
      <SampleHeader data={data} jobId={jobId} isDemo={isDemo} />
      <div className="mt-5">
        <KpiStrip data={data} rows={rows} onOpen={setActiveModule} />
      </div>
      <div className="sticky top-14 z-10 -mx-4 mt-4 bg-shell/95 px-4 backdrop-blur sm:-mx-6 sm:px-6">
        <ModuleNavigation active={activeModule} onSelect={setActiveModule} />
      </div>

      <div className="mt-5">
        {overview ? (
          <div className="grid gap-4">
            <div className="flex flex-col gap-3 lg:flex-row lg:items-center">
              <div className="segmented self-start" role="group" aria-label="Overview layout">
                <ViewButton view="table" current={overviewView} onSelect={setOverviewView} />
                <ViewButton view="genome" current={overviewView} onSelect={setOverviewView} />
              </div>
              <div className="min-w-0 flex-1">
                <EventFiltersBar rows={rows} shown={filtered.length} />
              </div>
            </div>

            {overviewView === 'table' ? (
              <div className="grid items-start gap-4 xl:grid-cols-[minmax(0,1fr)_300px]">
                <section className="card overflow-hidden" aria-label="Events">
                  <EventsTable rows={filtered} selectedId={selectedId} onSelect={selectRow} hasDenominators={hasDenominators} />
                  <p className="border-t border-line px-4 py-3 text-xs leading-5 text-muted">
                    {hasDenominators
                      ? 'HF is the heteroplasmic fraction: ALT molecules over callable molecules at the position, with a 95% confidence interval. QC facts are observational, not calibrated filters.'
                      : 'This result reports supporting records without callable denominators, so support is not a frequency. Re-run with schema 0.6 for heteroplasmic fractions.'}
                  </p>
                </section>
                <aside className="card xl:sticky xl:top-32" aria-labelledby="mini-map-title">
                  <div className="flex items-center justify-between gap-2 px-4 pt-3.5">
                    <h2 id="mini-map-title" className="text-[13px] font-semibold">Genome map</h2>
                    <button type="button" className="text-[13px] font-medium text-accent hover:underline" onClick={() => setOverviewView('genome')}>
                      Open genome view
                    </button>
                  </div>
                  <div className="mx-auto max-w-[280px] px-3">
                    <CircosWrapper
                      compact
                      data={data}
                      activeLayers={{ genes: true, coverage: true, clusters: false, svs: true, snps: true }}
                      selection={{ eventId: selectedEventId, svId: selectedSvId, gene: selectedGene }}
                      onSvSelect={selectSv}
                      onEventSelect={selectEvent}
                      onGeneSelect={selectGene}
                    />
                  </div>
                  <div className="border-t border-line px-4 py-3">
                    <GenomeLegend compact />
                  </div>
                </aside>
              </div>
            ) : (
              <GenomeView data={data} rows={filtered} onSelectEvent={selectEvent} onSelectSv={selectSv} onSelectGene={selectGene} />
            )}
          </div>
        ) : (
          <AnalysisModulesPanel data={data} jobId={isDemo ? undefined : jobId} />
        )}
      </div>

      {overview && (
        <DetailDrawer
          data={data}
          rows={rows}
          sequence={sequence}
          onSelectRow={selectRow}
          onSelectGene={selectGene}
          onShowInGenome={showInGenome}
        />
      )}
    </main>
  );
}

function ViewButton({ view, current, onSelect }: { view: OverviewView; current: OverviewView; onSelect: (view: OverviewView) => void }) {
  const Icon = view === 'table' ? Rows3 : Orbit;
  return (
    <button type="button" className="segmented-item" aria-pressed={current === view} onClick={() => onSelect(view)}>
      <Icon className="h-4 w-4" aria-hidden />
      {view === 'table' ? 'Table' : 'Genome'}
    </button>
  );
}

/**
 * Keeps `?module=`, `?view=`, and the selection (`?event=`, `?sv=`, `?gene=`)
 * in step with the store. URL changes (links, back/forward) update the store;
 * changes made in the UI write the URL, pushing history only when the
 * section or view changes.
 */
function useUrlState(dataLoaded: boolean) {
  const [searchParams, setSearchParams] = useSearchParams();
  const activeModule = useMitoStore((state) => state.activeModule);
  const overviewView = useMitoStore((state) => state.overviewView);
  const eventId = useMitoStore((state) => state.selectedEventId);
  const svId = useMitoStore((state) => state.selectedSvId);
  const gene = useMitoStore((state) => state.selectedGene);
  const setActiveModule = useMitoStore((state) => state.setActiveModule);
  const setOverviewView = useMitoStore((state) => state.setOverviewView);
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);
  const setSelectedSvId = useMitoStore((state) => state.setSelectedSvId);
  const setSelectedGene = useMitoStore((state) => state.setSelectedGene);

  const urlModule = parseAnalysisModule(searchParams.get('module')) ?? 'overview';
  const urlView: OverviewView = searchParams.get('view') === 'genome' ? 'genome' : 'table';
  const urlEvent = searchParams.get('event') ?? undefined;
  const urlSv = searchParams.get('sv') ?? undefined;
  const urlGene = searchParams.get('gene') ?? undefined;
  const placeKey = activeModule + '|' + overviewView;
  const storeKey = [placeKey, eventId, svId, gene].join('|');
  const urlKey = [urlModule + '|' + urlView, urlEvent, urlSv, urlGene].join('|');
  const previous = useRef({ storeKey, placeKey });

  useEffect(() => {
    setActiveModule(urlModule);
    setOverviewView(urlView);
  }, [setActiveModule, setOverviewView, urlModule, urlView]);

  // Loading a result clears the selection, so apply the URL selection after it.
  useEffect(() => {
    if (!dataLoaded) return;
    if (urlEvent) setSelectedEventId(urlEvent);
    else if (urlSv) setSelectedSvId(urlSv);
    else if (urlGene) setSelectedGene(urlGene);
  }, [dataLoaded, setSelectedEventId, setSelectedGene, setSelectedSvId, urlEvent, urlGene, urlSv]);

  useEffect(() => {
    const last = previous.current;
    previous.current = { storeKey, placeKey };
    if (last.storeKey === storeKey || storeKey === urlKey) return;
    const next = new URLSearchParams(searchParams);
    next.set('module', activeModule);
    const entries: Array<[string, string | undefined]> = [
      ['view', activeModule === 'overview' && overviewView === 'genome' ? 'genome' : undefined],
      ['event', eventId],
      ['sv', svId],
      ['gene', gene]
    ];
    for (const [name, value] of entries) {
      if (value) next.set(name, value);
      else next.delete(name);
    }
    setSearchParams(next, { replace: last.placeKey === placeKey });
    // Keyed on the store alone: only store-originated changes write the URL.
  }, [storeKey]);
}

function ResultSkeleton() {
  return (
    <main id="main-content" className="mx-auto max-w-[1600px] px-4 pt-6 sm:px-6" aria-busy="true" aria-label="Loading result">
      <div className="h-16 w-80 animate-pulse rounded-lg bg-panel2" />
      <div className="mt-5 h-20 animate-pulse rounded-xl bg-panel2" />
      <div className="mt-4 h-11 animate-pulse rounded-lg bg-panel2" />
      <div className="mt-5 h-[480px] animate-pulse rounded-xl bg-panel2" />
    </main>
  );
}
