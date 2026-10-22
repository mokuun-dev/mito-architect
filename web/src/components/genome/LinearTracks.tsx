import type { AlignmentFragment, CoverageBin, GeneAnnotation, MitoAnalysisData } from '@mito-architect/visualization-lib';
import { useLayoutEffect, useMemo, useRef, useState, type KeyboardEvent, type PointerEvent, type ReactNode } from 'react';
import type { EventRow } from '../../lib/eventRows';
import type { GenomeRegion } from '../../lib/store';
import { geneColor, variantShapePath } from './trackStyle';

export type TrackName = 'genes' | 'coverage' | 'variants' | 'svs' | 'molecules';

interface LinearTracksProps {
  data: MitoAnalysisData;
  rows: EventRow[];
  region: { start: number; end: number };
  onRegionChange?: (region?: GenomeRegion) => void;
  selection?: { eventId?: string; svId?: string; gene?: string };
  highlightMoleculeIds?: ReadonlySet<string>;
  onSelectEvent?: (id: string) => void;
  onSelectSv?: (id: string) => void;
  onSelectGene?: (gene: GeneAnnotation) => void;
  tracks?: TrackName[];
}

interface Layout {
  name: TrackName | 'axis';
  label: string;
  detail?: string;
  top: number;
  height: number;
}

const LABEL_WIDTH = 92;
const RIGHT_PAD = 12;
const TRACK_GAP = 12;
const MAX_SV_LANES = 6;
const MAX_MOLECULE_LANES = 24;

/**
 * Linear tracks over one reference region. Positions are 1-based bases; base
 * b spans the boundary interval [b - 1, b], so a region [s, e] maps the
 * boundaries s - 1 .. e onto the plot width.
 */
export function LinearTracks({
  data,
  rows,
  region,
  onRegionChange,
  selection = {},
  highlightMoleculeIds,
  onSelectEvent,
  onSelectSv,
  onSelectGene,
  tracks = ['genes', 'coverage', 'variants', 'svs', 'molecules']
}: LinearTracksProps) {
  const hostRef = useRef<HTMLDivElement>(null);
  const [width, setWidth] = useState(800);
  const [hover, setHover] = useState<{ x: number; y: number; content: ReactNode }>();
  const [crosshair, setCrosshair] = useState<number>();
  const [brush, setBrush] = useState<{ from: number; to: number }>();
  const referenceLength = data.metadata.reference_length;
  const start = Math.max(1, Math.min(region.start, region.end));
  const end = Math.min(referenceLength, Math.max(region.start, region.end));
  const span = end - start + 1;
  const plotWidth = Math.max(120, width - LABEL_WIDTH - RIGHT_PAD);
  const x = (boundary: number) => LABEL_WIDTH + ((boundary - (start - 1)) / span) * plotWidth;
  const baseAt = (px: number) => Math.min(end, Math.max(start, Math.floor(start - 1 + ((px - LABEL_WIDTH) / plotWidth) * span) + 1));

  // Measure before the first paint so tracks never render at a placeholder width.
  useLayoutEffect(() => {
    const host = hostRef.current;
    if (!host) return;
    setWidth(Math.round(host.getBoundingClientRect().width) || 800);
    const observer = new ResizeObserver(([entry]) => setWidth(Math.round(entry.contentRect.width) || 800));
    observer.observe(host);
    return () => observer.disconnect();
  }, []);

  const genes = useMemo(() => data.genes.filter((gene) => gene.end >= start && gene.start <= end), [data.genes, end, start]);
  const bins = useMemo(() => data.coverage.filter((bin) => bin.end >= start && bin.start <= end).sort((a, b) => a.start - b.start), [data.coverage, end, start]);
  const maxDepth = Math.max(1, ...bins.map((bin) => bin.depth));
  const variants = useMemo(() => rows.filter((row) => row.kind !== 'sv' && row.position >= start && row.position <= end), [end, rows, start]);
  const svLanes = useMemo(() => packStructural(rows.filter((row) => row.kind === 'sv'), start, end, referenceLength), [end, referenceLength, rows, start]);
  const moleculeLanes = useMemo(() => packMolecules(data.alignments ?? [], start, end), [data.alignments, end, start]);

  const showMolecules = tracks.includes('molecules') && (data.alignments?.length ?? 0) > 0;
  const layout: Layout[] = [];
  let cursor = 0;
  const add = (name: Layout['name'], label: string, height: number, detail?: string) => {
    layout.push({ name, label, detail, top: cursor, height });
    cursor += height + TRACK_GAP;
  };
  add('axis', '', 26);
  if (tracks.includes('genes')) add('genes', 'Genes', 26);
  if (tracks.includes('coverage')) add('coverage', 'Coverage', 64, `0–${maxDepth.toLocaleString('en-US')}×`);
  if (tracks.includes('variants')) add('variants', 'Variants', 76, 'HF 0–100%');
  if (tracks.includes('svs')) add('svs', 'Structural', Math.max(1, Math.min(MAX_SV_LANES, svLanes.lanes.length)) * 12 + 6);
  if (showMolecules) add('molecules', 'Molecules', Math.max(1, Math.min(MAX_MOLECULE_LANES, moleculeLanes.lanes.length)) * 7 + 6, `${moleculeLanes.total.toLocaleString('en-US')} in view`);
  const height = cursor - TRACK_GAP + 4;
  const track = (name: Layout['name']) => layout.find((item) => item.name === name);

  function showTooltip(event: PointerEvent<Element>, content: ReactNode) {
    const bounds = hostRef.current?.getBoundingClientRect();
    if (!bounds) return;
    setHover({ x: event.clientX - bounds.left, y: event.clientY - bounds.top, content });
  }

  function onPlotMove(event: PointerEvent<SVGSVGElement>) {
    const bounds = event.currentTarget.getBoundingClientRect();
    const px = event.clientX - bounds.left;
    if (brush) {
      setBrush({ ...brush, to: Math.min(LABEL_WIDTH + plotWidth, Math.max(LABEL_WIDTH, px)) });
      return;
    }
    if (px < LABEL_WIDTH || px > LABEL_WIDTH + plotWidth) {
      setCrosshair(undefined);
      return;
    }
    setCrosshair(px);
  }

  function finishBrush() {
    if (!brush) return;
    const from = Math.min(brush.from, brush.to);
    const to = Math.max(brush.from, brush.to);
    setBrush(undefined);
    if (to - from >= 4 && onRegionChange) {
      onRegionChange({ start: baseAt(from), end: baseAt(to) });
    }
  }

  const crosshairBase = crosshair === undefined ? undefined : baseAt(crosshair);
  const crosshairBin = crosshairBase === undefined ? undefined : bins.find((bin) => bin.start <= crosshairBase && bin.end >= crosshairBase);
  const crosshairGenes = crosshairBase === undefined ? [] : genes.filter((gene) => gene.start <= crosshairBase && gene.end >= crosshairBase);

  return (
    <div ref={hostRef} className="relative w-full select-none">
      <svg
        // Re-create the drawing when the width changes: Firefox can leave
        // stale geometry behind when many SVG attributes change at once.
        key={width}
        width={width}
        height={height}
        viewBox={`0 0 ${width} ${height}`}
        className="block max-w-full"
        role="group"
        aria-label={`Linear tracks for reference positions ${start.toLocaleString('en-US')} to ${end.toLocaleString('en-US')}`}
        onPointerMove={onPlotMove}
        onPointerUp={finishBrush}
        onPointerLeave={() => {
          setCrosshair(undefined);
          setHover(undefined);
          finishBrush();
        }}
      >
        {layout.filter((item) => item.name !== 'axis').map((item) => (
          <g key={item.name}>
            <text x={0} y={item.top + 12} className="fill-text text-[11px] font-semibold">{item.label}</text>
            {item.detail && <text x={0} y={item.top + 26} className="fill-muted text-[10px] tabular-nums">{item.detail}</text>}
            <line x1={LABEL_WIDTH} x2={LABEL_WIDTH + plotWidth} y1={item.top + item.height} y2={item.top + item.height} className="stroke-line" />
          </g>
        ))}

        <Axis
          top={track('axis')!.top}
          start={start}
          end={end}
          x={x}
          plotWidth={plotWidth}
          interactive={Boolean(onRegionChange)}
          onBrushStart={(px) => setBrush({ from: px, to: px })}
        />

        {track('genes') && genes.map((gene) => {
          const item = track('genes')!;
          const x0 = x(Math.max(gene.start, start) - 1);
          const x1 = x(Math.min(gene.end, end));
          const selected = selection.gene === gene.name;
          const label = gene.name.replace(/^MT-/, '');
          return (
            <g
              key={gene.name}
              role="button"
              tabIndex={onSelectGene ? 0 : -1}
              aria-label={`${gene.name}, ${gene.start}–${gene.end}, ${gene.biotype}, strand ${gene.strand}`}
              className={onSelectGene ? 'cursor-pointer' : undefined}
              onClick={() => onSelectGene?.(gene)}
              onKeyDown={(event) => activate(event, () => onSelectGene?.(gene))}
              onPointerMove={(event) => showTooltip(event, <><strong>{gene.name}</strong> · {gene.start.toLocaleString('en-US')}–{gene.end.toLocaleString('en-US')} · {biotypeLabel(gene.biotype)} · strand {gene.strand}</>)}
              onPointerLeave={() => setHover(undefined)}
            >
              <rect
                x={x0}
                y={item.top + 3}
                width={Math.max(1, x1 - x0)}
                height={item.height - 6}
                rx={3}
                fill={geneColor(gene.biotype)}
                fillOpacity={selected ? 0.32 : 0.16}
                stroke={selected ? 'rgb(var(--color-text))' : geneColor(gene.biotype)}
                strokeWidth={selected ? 1.5 : 1}
              />
              {x1 - x0 > label.length * 6.5 + 8 && (
                <text x={(x0 + x1) / 2} y={item.top + item.height / 2 + 4} textAnchor="middle" className="pointer-events-none fill-text text-[10.5px] font-semibold">
                  {label}
                </text>
              )}
            </g>
          );
        })}

        {track('coverage') && <CoverageTrack item={track('coverage')!} bins={bins} start={start} end={end} x={x} maxDepth={maxDepth} />}

        {track('variants') && (() => {
          const item = track('variants')!;
          const base = item.top + item.height;
          const usable = item.height - 12;
          return (
            <g>
              {[0.5, 1].map((fraction) => (
                <line key={fraction} x1={LABEL_WIDTH} x2={LABEL_WIDTH + plotWidth} y1={base - fraction * usable} y2={base - fraction * usable} className="stroke-line" strokeDasharray="2 4" />
              ))}
              <text x={LABEL_WIDTH - 6} y={base - usable + 3} textAnchor="end" className="fill-subtle text-[9px] tabular-nums">100%</text>
              <text x={LABEL_WIDTH - 6} y={base - usable / 2 + 3} textAnchor="end" className="fill-subtle text-[9px] tabular-nums">50%</text>
              {variants.map((row) => {
                const cx = x(row.position - 0.5);
                const measured = row.hf !== undefined;
                const cy = base - (measured ? row.hf! : 0.5) * usable;
                const selected = selection.eventId === row.id;
                return (
                  <g
                    key={row.id}
                    role="button"
                    tabIndex={onSelectEvent ? 0 : -1}
                    aria-label={`${row.typeLabel} ${row.position} ${row.change}${measured ? `, HF ${formatPercent(row.hf!)}` : ', no callable denominator'}`}
                    className="cursor-pointer"
                    onClick={() => onSelectEvent?.(row.id)}
                    onKeyDown={(event) => activate(event, () => onSelectEvent?.(row.id))}
                    onPointerMove={(event) => showTooltip(event, <VariantTooltip row={row} />)}
                    onPointerLeave={() => setHover(undefined)}
                  >
                    <line x1={cx} x2={cx} y1={base} y2={cy} stroke="rgb(var(--color-variant))" strokeOpacity={0.7} strokeWidth={1.5} />
                    <path
                      d={variantShapePath(row.type, 4.5)}
                      transform={`translate(${cx} ${cy})`}
                      fill={measured ? 'rgb(var(--color-variant))' : 'rgb(var(--color-panel))'}
                      stroke={measured ? 'rgb(var(--color-panel))' : 'rgb(var(--color-variant))'}
                      strokeWidth={1.5}
                    />
                    {selected && <circle cx={cx} cy={cy} r={8.5} fill="none" stroke="rgb(var(--color-text))" strokeWidth={2} />}
                    <circle cx={cx} cy={cy} r={9} fill="transparent" />
                  </g>
                );
              })}
            </g>
          );
        })()}

        {track('svs') && (() => {
          const item = track('svs')!;
          return (
            <g>
              {svLanes.lanes.slice(0, MAX_SV_LANES).map((lane, laneIndex) => lane.map(({ row, segments }) => {
                const selected = selection.svId === row.id;
                return (
                  <g
                    key={row.id}
                    role="button"
                    tabIndex={onSelectSv ? 0 : -1}
                    aria-label={`${row.typeLabel} ${row.change}, ${row.altSupport} supporting molecules`}
                    className="cursor-pointer"
                    onClick={() => onSelectSv?.(row.id)}
                    onKeyDown={(event) => activate(event, () => onSelectSv?.(row.id))}
                    onPointerMove={(event) => showTooltip(event, <><strong>{row.typeLabel}</strong> · {row.change}{row.length ? ` · ${row.length.toLocaleString('en-US')} bp` : ''} · {row.altSupport} supporting molecules</>)}
                    onPointerLeave={() => setHover(undefined)}
                  >
                    {segments.map(([from, to]) => (
                      <rect
                        key={from}
                        x={x(from - 1)}
                        y={item.top + 4 + laneIndex * 12}
                        width={Math.max(2, x(to) - x(from - 1))}
                        height={7}
                        rx={3.5}
                        fill="rgb(var(--color-sv))"
                        fillOpacity={selected ? 1 : 0.8}
                        stroke={selected ? 'rgb(var(--color-text))' : 'none'}
                        strokeWidth={1.5}
                      />
                    ))}
                  </g>
                );
              }))}
              {svLanes.lanes.length === 0 && <text x={LABEL_WIDTH + 4} y={item.top + 13} className="fill-subtle text-[11px]">No structural events in this region</text>}
              {svLanes.lanes.length > MAX_SV_LANES && (
                <text x={LABEL_WIDTH + plotWidth} y={item.top + item.height + 10} textAnchor="end" className="fill-muted text-[10px]">
                  {svLanes.lanes.slice(MAX_SV_LANES).reduce((sum, lane) => sum + lane.length, 0)} more not drawn
                </text>
              )}
            </g>
          );
        })()}

        {track('molecules') && (() => {
          const item = track('molecules')!;
          return (
            <g>
              {moleculeLanes.lanes.slice(0, MAX_MOLECULE_LANES).map((lane, laneIndex) => lane.map((molecule) => {
                const highlighted = highlightMoleculeIds?.has(molecule.id) ?? false;
                return (
                  <g
                    key={molecule.id}
                    onPointerMove={(event) => showTooltip(event, <><strong className="font-mono">{molecule.id}</strong> · {molecule.fragments.length} aligned fragment{molecule.fragments.length === 1 ? '' : 's'}{highlighted ? ' · supports the selected event' : ''}</>)}
                    onPointerLeave={() => setHover(undefined)}
                  >
                    {molecule.fragments.map((fragment) => (
                      <rect
                        key={fragment.id}
                        x={x(Math.max(start, fragment.start) - 1)}
                        y={item.top + 4 + laneIndex * 7}
                        width={Math.max(1, x(Math.min(end, fragment.end)) - x(Math.max(start, fragment.start) - 1))}
                        height={4}
                        rx={2}
                        fill={highlighted ? 'rgb(var(--color-variant))' : 'rgb(var(--color-coverage))'}
                        fillOpacity={highlighted ? 1 : fragment.supplementary ? 0.55 : 0.9}
                      />
                    ))}
                  </g>
                );
              }))}
              {moleculeLanes.lanes.length > MAX_MOLECULE_LANES && (
                <text x={LABEL_WIDTH + plotWidth} y={item.top + item.height + 10} textAnchor="end" className="fill-muted text-[10px]">
                  {moleculeLanes.lanes.slice(MAX_MOLECULE_LANES).reduce((sum, lane) => sum + lane.length, 0)} more molecules not drawn
                </text>
              )}
            </g>
          );
        })()}

        {crosshair !== undefined && !brush && (
          <line x1={crosshair} x2={crosshair} y1={0} y2={height} className="pointer-events-none stroke-text" strokeOpacity={0.35} strokeDasharray="3 3" />
        )}
        {brush && (
          <rect
            x={Math.min(brush.from, brush.to)}
            y={0}
            width={Math.abs(brush.to - brush.from)}
            height={height}
            className="pointer-events-none fill-accent stroke-accent"
            fillOpacity={0.1}
          />
        )}
      </svg>

      {(hover || crosshairBase !== undefined) && !brush && (
        <div
          className="pointer-events-none absolute z-10 max-w-[300px] rounded-lg border border-line bg-panel px-2.5 py-1.5 text-xs leading-5 text-text shadow-pop"
          style={tooltipPosition(hover ? hover : { x: crosshair ?? 0, y: 8 }, width)}
        >
          {hover ? hover.content : (
            <>
              <span className="font-mono font-semibold tabular-nums">{crosshairBase!.toLocaleString('en-US')}</span>
              {crosshairBin && <> · {crosshairBin.depth.toLocaleString('en-US')}× depth</>}
              {crosshairGenes.length > 0 && <> · {crosshairGenes.map((gene) => gene.name).join(', ')}</>}
            </>
          )}
        </div>
      )}
    </div>
  );
}

function Axis({ top, start, end, x, plotWidth, interactive, onBrushStart }: {
  top: number;
  start: number;
  end: number;
  x: (boundary: number) => number;
  plotWidth: number;
  interactive: boolean;
  onBrushStart: (px: number) => void;
}) {
  const ticks = niceTicks(start, end, Math.max(3, Math.floor(plotWidth / 96)));
  return (
    <g>
      <rect
        x={LABEL_WIDTH}
        y={top}
        width={plotWidth}
        height={26}
        rx={4}
        className={interactive ? 'cursor-col-resize fill-panel2' : 'fill-panel2'}
        onPointerDown={(event) => {
          if (!interactive) return;
          const bounds = (event.currentTarget.ownerSVGElement ?? event.currentTarget).getBoundingClientRect();
          onBrushStart(event.clientX - bounds.left);
        }}
      >
        {interactive && <title>Drag to zoom into a region</title>}
      </rect>
      {ticks.map((tick) => (
        <g key={tick} className="pointer-events-none">
          <line x1={x(tick - 0.5)} x2={x(tick - 0.5)} y1={top + 20} y2={top + 26} className="stroke-line-strong" />
          <text x={x(tick - 0.5)} y={top + 15} textAnchor="middle" className="fill-muted text-[10px] tabular-nums">
            {tick.toLocaleString('en-US')}
          </text>
        </g>
      ))}
      <text x={0} y={top + 16} className="fill-subtle text-[10px]">{interactive ? 'Drag to zoom' : 'Position'}</text>
    </g>
  );
}

function CoverageTrack({ item, bins, start, end, x, maxDepth }: {
  item: Layout;
  bins: CoverageBin[];
  start: number;
  end: number;
  x: (boundary: number) => number;
  maxDepth: number;
}) {
  const base = item.top + item.height;
  const y = (depth: number) => base - (depth / maxDepth) * (item.height - 4);
  let path = '';
  for (const bin of bins) {
    const from = x(Math.max(start, bin.start) - 1);
    const to = x(Math.min(end, bin.end));
    path += `M${from},${base}L${from},${y(bin.depth)}L${to},${y(bin.depth)}L${to},${base}Z`;
  }
  return (
    <g className="pointer-events-none">
      <path d={path} fill="rgb(var(--color-coverage))" fillOpacity={0.75} />
      {bins.filter((bin) => bin.depth === 0).map((bin) => (
        <rect
          key={bin.start}
          x={x(Math.max(start, bin.start) - 1)}
          y={base - 3}
          width={Math.max(1, x(Math.min(end, bin.end)) - x(Math.max(start, bin.start) - 1))}
          height={3}
          fill="rgb(var(--color-coral))"
        />
      ))}
      {bins.length === 0 && <text x={LABEL_WIDTH + 4} y={item.top + 14} className="fill-subtle text-[11px]">No coverage bins in this region</text>}
    </g>
  );
}

function VariantTooltip({ row }: { row: EventRow }) {
  return (
    <>
      <strong className="font-mono">{row.position.toLocaleString('en-US')} {row.change}</strong> · {row.typeLabel}
      {row.genes.length > 0 && <> · {row.genes.join(', ')}</>}
      <br />
      {row.hf !== undefined
        ? <>HF {formatPercent(row.hf)} ({row.altSupport}/{row.callable} molecules)</>
        : <>{row.altSupport} supporting records · no callable denominator</>}
    </>
  );
}

function packStructural(rows: EventRow[], start: number, end: number, length: number) {
  const lanes: Array<Array<{ row: EventRow; segments: Array<[number, number]> }>> = [];
  const occupied: Array<Array<[number, number]>> = [];
  const gap = Math.max(1, Math.round((end - start) / 200));
  for (const row of [...rows].sort((a, b) => a.position - b.position)) {
    // Events that wrap the origin (end < start) occupy two segments.
    const spans: Array<[number, number]> = row.end >= row.position ? [[row.position, row.end]] : [[row.position, length], [1, row.end]];
    const segments = spans
      .map(([from, to]) => [Math.max(from, start), Math.min(to, end)] as [number, number])
      .filter(([from, to]) => from <= to);
    if (segments.length === 0) continue;
    let lane = occupied.findIndex((intervals) =>
      segments.every(([from, to]) => intervals.every(([taken, until]) => to < taken || from > until)));
    if (lane === -1) {
      lane = lanes.length;
      lanes.push([]);
      occupied.push([]);
    }
    lanes[lane].push({ row, segments });
    occupied[lane].push(...segments.map(([from, to]) => [from - gap, to + gap] as [number, number]));
  }
  return { lanes };
}

interface MoleculeExtent {
  id: string;
  start: number;
  end: number;
  fragments: Array<{ id: string; start: number; end: number; supplementary: boolean }>;
}

function packMolecules(alignments: AlignmentFragment[], start: number, end: number) {
  const molecules = new Map<string, MoleculeExtent>();
  for (const fragment of alignments) {
    if (fragment.role === 'secondary' || fragment.role === 'unaligned_read' || fragment.reference_start <= 0) continue;
    const fragmentEnd = fragment.reference_start + referenceSpan(fragment.cigar) - 1;
    if (fragmentEnd < start || fragment.reference_start > end) continue;
    const molecule = molecules.get(fragment.molecule_id) ?? { id: fragment.molecule_id, start: Infinity, end: -Infinity, fragments: [] };
    molecule.start = Math.min(molecule.start, fragment.reference_start);
    molecule.end = Math.max(molecule.end, fragmentEnd);
    molecule.fragments.push({ id: fragment.id, start: fragment.reference_start, end: fragmentEnd, supplementary: fragment.role.includes('supplementary') });
    molecules.set(fragment.molecule_id, molecule);
  }
  const lanes: MoleculeExtent[][] = [];
  const laneEnds: number[] = [];
  for (const molecule of [...molecules.values()].sort((a, b) => a.start - b.start || a.id.localeCompare(b.id))) {
    let lane = laneEnds.findIndex((laneEnd) => laneEnd < molecule.start);
    if (lane === -1) {
      lane = lanes.length;
      lanes.push([]);
      laneEnds.push(0);
    }
    lanes[lane].push(molecule);
    laneEnds[lane] = molecule.end + Math.max(1, Math.round((end - start) / 200));
  }
  return { lanes, total: molecules.size };
}

/** Reference bases consumed by a CIGAR string (M, D, N, =, X). */
export function referenceSpan(cigar: string): number {
  let span = 0;
  for (const match of cigar.matchAll(/(\d+)([MIDNSHP=X])/g)) {
    if ('MDN=X'.includes(match[2])) span += Number(match[1]);
  }
  return Math.max(1, span);
}

export function niceTicks(start: number, end: number, count: number): number[] {
  const span = Math.max(1, end - start);
  const raw = span / Math.max(1, count);
  const magnitude = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((factor) => factor * magnitude).find((candidate) => candidate >= raw) ?? magnitude * 10;
  const ticks: number[] = [];
  for (let tick = Math.ceil(start / step) * step; tick <= end; tick += step) {
    if (tick >= 1) ticks.push(tick);
  }
  return ticks;
}

function tooltipPosition(point: { x: number; y: number }, width: number) {
  const left = point.x + 14;
  return left + 300 > width ? { right: Math.max(0, width - point.x + 14), top: point.y + 12 } : { left, top: point.y + 12 };
}

function activate(event: KeyboardEvent, action: () => void) {
  if (event.key === 'Enter' || event.key === ' ') {
    event.preventDefault();
    action();
  }
}

function biotypeLabel(biotype: string): string {
  return biotype === 'protein_coding' ? 'protein-coding' : biotype;
}

function formatPercent(value: number): string {
  return `${(value * 100).toFixed(value < 0.1 ? 2 : 1)}%`;
}
