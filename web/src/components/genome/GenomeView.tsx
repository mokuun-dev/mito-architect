import type { GeneAnnotation, LayerName, MitoAnalysisData } from '@mito-architect/visualization-lib';
import { ChevronLeft, ChevronRight, Layers3, Maximize2, Search, ZoomIn, ZoomOut } from 'lucide-react';
import { useMemo, useState } from 'react';
import type { EventRow } from '../../lib/eventRows';
import { useMitoStore } from '../../lib/store';
import CircosWrapper from '../CircosWrapper';
import { GenomeLegend } from './GenomeLegend';
import { LinearTracks } from './LinearTracks';

const LAYERS: Array<[LayerName, string]> = [
  ['genes', 'Genes'],
  ['coverage', 'Coverage'],
  ['snps', 'Variants'],
  ['svs', 'Structural events']
];

export function GenomeView({ data, rows, onSelectEvent, onSelectSv, onSelectGene }: {
  data: MitoAnalysisData;
  rows: EventRow[];
  onSelectEvent: (id: string) => void;
  onSelectSv: (id: string) => void;
  onSelectGene: (gene: GeneAnnotation) => void;
}) {
  const length = data.metadata.reference_length;
  const activeLayers = useMitoStore((state) => state.activeLayers);
  const setLayer = useMitoStore((state) => state.setLayer);
  const region = useMitoStore((state) => state.selectedRegion);
  const setRegion = useMitoStore((state) => state.setSelectedRegion);
  const eventId = useMitoStore((state) => state.selectedEventId);
  const svId = useMitoStore((state) => state.selectedSvId);
  const gene = useMitoStore((state) => state.selectedGene);
  const visible = region ?? { start: 1, end: length };
  const span = visible.end - visible.start + 1;

  const highlightMoleculeIds = useMemo(() => {
    const row = rows.find((item) => item.id === eventId || item.id === svId);
    return new Set(row?.variant?.supporting_molecule_ids ?? row?.variant?.supporting_reads ?? row?.sv?.supporting_reads ?? []);
  }, [eventId, rows, svId]);

  function show(start: number, end: number, label?: string) {
    const nextSpan = Math.min(length, Math.max(40, end - start + 1));
    const center = (start + end) / 2;
    let from = Math.round(center - nextSpan / 2);
    from = Math.max(1, Math.min(length - nextSpan + 1, from));
    const next = { start: from, end: from + nextSpan - 1, label };
    setRegion(next.start === 1 && next.end === length ? undefined : next);
  }

  return (
    <div className="grid gap-4 xl:grid-cols-[minmax(340px,420px)_minmax(0,1fr)]">
      <section className="card min-w-0" aria-labelledby="circular-map-title">
        <div className="card-header">
          <div>
            <h2 id="circular-map-title" className="card-title">Circular map</h2>
            <p className="card-subtitle">Origin at the top, clockwise. Ctrl + scroll zooms.</p>
          </div>
          <details className="relative">
            <summary className="btn-secondary list-none"><Layers3 className="h-4 w-4" aria-hidden />Layers</summary>
            <fieldset className="popover right-0 top-[calc(100%+6px)] w-56 py-1">
              <legend className="sr-only">Map layers</legend>
              {LAYERS.map(([name, label]) => (
                <label key={name} className="flex min-h-9 cursor-pointer items-center justify-between gap-4 px-3 text-sm hover:bg-panel2">
                  <span>{label}</span>
                  <input type="checkbox" checked={activeLayers[name]} onChange={(event) => setLayer(name, event.target.checked)} className="h-4 w-4 accent-[rgb(var(--color-accent))]" />
                </label>
              ))}
            </fieldset>
          </details>
        </div>
        <div className="p-4">
          <CircosWrapper
            data={data}
            activeLayers={activeLayers}
            selection={{ eventId, svId, gene }}
            onSvSelect={onSelectSv}
            onEventSelect={onSelectEvent}
            onGeneSelect={onSelectGene}
          />
        </div>
      </section>

      <section className="card min-w-0" aria-labelledby="linear-tracks-title">
        <div className="card-header">
          <div className="min-w-0">
            <h2 id="linear-tracks-title" className="card-title">Linear tracks</h2>
            <p className="card-subtitle tabular-nums">
              {region?.label ? <span className="font-medium text-text">{region.label} · </span> : null}
              {visible.start.toLocaleString('en-US')}–{visible.end.toLocaleString('en-US')} · {span.toLocaleString('en-US')} bp
            </p>
          </div>
          <div className="flex flex-wrap items-center gap-2">
            <GeneJump genes={data.genes} onJump={(target) => { onSelectGene(target); show(target.start - 60, target.end + 60, target.name); }} />
            <div className="flex items-center gap-1" role="group" aria-label="Region navigation">
              <button type="button" className="btn-icon" aria-label="Pan left" title="Pan left" disabled={visible.start <= 1} onClick={() => show(visible.start - span / 2, visible.end - span / 2)}><ChevronLeft className="h-4 w-4" aria-hidden /></button>
              <button type="button" className="btn-icon" aria-label="Pan right" title="Pan right" disabled={visible.end >= length} onClick={() => show(visible.start + span / 2, visible.end + span / 2)}><ChevronRight className="h-4 w-4" aria-hidden /></button>
              <button type="button" className="btn-icon" aria-label="Zoom in" title="Zoom in" disabled={span <= 40} onClick={() => show(visible.start + span / 4, visible.end - span / 4)}><ZoomIn className="h-4 w-4" aria-hidden /></button>
              <button type="button" className="btn-icon" aria-label="Zoom out" title="Zoom out" disabled={!region} onClick={() => show(visible.start - span / 2, visible.end + span / 2)}><ZoomOut className="h-4 w-4" aria-hidden /></button>
              <button type="button" className="btn-icon" aria-label="Show the whole reference" title="Whole reference" disabled={!region} onClick={() => setRegion(undefined)}><Maximize2 className="h-4 w-4" aria-hidden /></button>
            </div>
          </div>
        </div>
        <div className="px-4 pb-3 pt-4">
          <LinearTracks
            data={data}
            rows={rows}
            region={visible}
            onRegionChange={(next) => next && show(next.start, next.end)}
            selection={{ eventId, svId, gene }}
            highlightMoleculeIds={highlightMoleculeIds}
            onSelectEvent={onSelectEvent}
            onSelectSv={onSelectSv}
            onSelectGene={onSelectGene}
          />
        </div>
        <div className="grid gap-2 border-t border-line px-5 py-3">
          <GenomeLegend />
          <p className="text-xs leading-5 text-muted">
            Tracks share one reference axis. Variant height is the heteroplasmic fraction over callable molecules; hollow markers have
            no callable denominator, so their support is not a frequency. The ring shows the reference, not the shape of an individual molecule.
            Exact values remain available in the linked textual projections.
          </p>
        </div>
      </section>
    </div>
  );
}

function GeneJump({ genes, onJump }: { genes: GeneAnnotation[]; onJump: (gene: GeneAnnotation) => void }) {
  const [query, setQuery] = useState('');
  const normalized = query.trim().toLowerCase();
  const found = normalized ? genes.filter((gene) => gene.name.toLowerCase().includes(normalized)).slice(0, 8) : [];
  return (
    <div className="relative w-48">
      <label htmlFor="genome-gene-search" className="sr-only">Go to gene</label>
      <Search className="pointer-events-none absolute left-2.5 top-1/2 h-4 w-4 -translate-y-1/2 text-subtle" aria-hidden />
      <input
        id="genome-gene-search"
        value={query}
        onChange={(event) => setQuery(event.target.value)}
        onKeyDown={(event) => {
          if (event.key === 'Enter' && found[0]) {
            onJump(found[0]);
            setQuery('');
          }
          if (event.key === 'Escape') setQuery('');
        }}
        placeholder="Go to gene…"
        className="input pl-8"
        autoComplete="off"
      />
      {normalized && (
        <div className="popover left-0 top-[calc(100%+4px)] w-full overflow-hidden py-1" role="listbox" aria-label="Matching genes">
          {found.length ? found.map((gene) => (
            <button
              key={gene.name}
              type="button"
              role="option"
              aria-selected={false}
              onClick={() => { onJump(gene); setQuery(''); }}
              className="flex min-h-9 w-full items-center justify-between gap-2 px-3 text-left text-sm hover:bg-panel2"
            >
              <span className="font-medium">{gene.name}</span>
              <span className="font-mono text-[11px] text-muted">{gene.start}–{gene.end}</span>
            </button>
          )) : <p className="px-3 py-2 text-sm text-muted">No matching gene</p>}
        </div>
      )}
    </div>
  );
}
