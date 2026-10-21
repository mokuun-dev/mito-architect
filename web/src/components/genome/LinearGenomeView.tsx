import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { useMemo } from 'react';
import { circularVariants } from '@mito-architect/visualization-lib';
import { useMitoStore } from '../../lib/store';

export function LinearGenomeView({ data }: { data: MitoAnalysisData }) {
  const region = useMitoStore((state) => state.selectedRegion);
  const setRegion = useMitoStore((state) => state.setSelectedRegion);
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);
  const length = data.metadata.reference_length;
  const visible = region ?? { start: 1, end: length, label: 'Whole reference' };
  const start = Math.max(1, Math.min(visible.start, visible.end));
  const end = Math.min(length, Math.max(visible.start, visible.end));
  const span = Math.max(1, end - start + 1);
  const items = useMemo(() => circularVariants(data).filter((variant) => variant.position >= start && variant.position <= end), [data, end, start]);
  const coverage = useMemo(() => data.coverage.filter((bin) => bin.end >= start && bin.start <= end), [data.coverage, end, start]);
  const maxDepth = Math.max(1, ...coverage.map((bin) => bin.depth));
  const left = (position: number) => ((position - start) / span) * 100;

  return (
    <section className="overflow-hidden rounded-xl border border-line bg-panel shadow-tool" aria-labelledby="linear-genome-title">
      <div className="flex flex-wrap items-start justify-between gap-3 border-b border-line px-4 py-3">
        <div>
          <div className="section-kicker">Region detail</div>
          <h2 id="linear-genome-title" className="mt-1 text-base font-semibold">Linear evidence view</h2>
          <p className="mt-1 text-xs text-muted">{visible.label ?? 'Selected range'} · {start.toLocaleString()}–{end.toLocaleString()} bp</p>
        </div>
        <button type="button" onClick={() => setRegion(undefined)} disabled={!region} className="min-h-10 rounded-lg border border-line bg-panel2 px-3 text-sm text-muted hover:border-aqua hover:text-aqua disabled:opacity-50">
          Whole reference
        </button>
      </div>
      <div className="data-scroll p-4">
        <div className="relative min-w-[680px]">
          <div className="relative h-12 rounded-md border border-line bg-shell">
            {data.genes.filter((gene) => gene.end >= start && gene.start <= end).map((gene) => {
              const geneStart = Math.max(start, gene.start);
              const geneEnd = Math.min(end, gene.end);
              return <button key={gene.name} type="button" title={gene.name + ' ' + gene.start + '–' + gene.end} onClick={() => setRegion({ start: gene.start, end: gene.end, label: gene.name })}
                className="absolute top-2 h-8 overflow-hidden rounded bg-sky/70 px-1 text-left text-[10px] font-semibold text-shell hover:ring-2 hover:ring-aqua"
                style={{ left: left(geneStart) + '%', width: Math.max(0.8, ((geneEnd - geneStart + 1) / span) * 100) + '%' }}>{gene.name.replace('MT-', '')}</button>;
            })}
          </div>
          <div className="relative mt-3 h-24 overflow-hidden rounded-md border border-line bg-shell" aria-label="Coverage across selected region">
            {coverage.map((bin) => {
              const binStart = Math.max(start, bin.start);
              const binEnd = Math.min(end, bin.end);
              return <div key={bin.start + '-' + bin.end} title={bin.start + '–' + bin.end + ': ' + bin.depth + '×'} className={bin.depth === 0 ? 'absolute bottom-0 bg-coral/75' : 'absolute bottom-0 bg-aqua/70'}
                style={{ left: left(binStart) + '%', width: Math.max(0.15, ((binEnd - binStart + 1) / span) * 100) + '%', height: (bin.depth / maxDepth) * 100 + '%' }} />;
            })}
          </div>
          <div className="relative mt-3 h-8 border-t border-line">
            {items.map((variant) => <button key={variant.id} type="button" title={variant.id} onClick={() => setSelectedEventId(variant.id)}
              className="absolute top-2 h-3 w-3 -translate-x-1/2 rounded-full border border-shell bg-amber hover:ring-2 hover:ring-aqua"
              style={{ left: left(variant.position) + '%' }} />)}
          </div>
          <div className="flex justify-between font-mono text-xs text-muted"><span>{start.toLocaleString()}</span><span>{end.toLocaleString()}</span></div>
        </div>
      </div>
      <p className="border-t border-line px-4 py-3 text-xs leading-5 text-muted">Tracks share one coordinate axis: genes, read depth, then aggregated variants. A zero-height coverage bar is a measured zero in the delivered bins; unavailable positional callability is not drawn as zero.</p>
    </section>
  );
}
