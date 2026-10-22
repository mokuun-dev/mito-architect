import { ArrowDown, ArrowUp, ChevronsUpDown } from 'lucide-react';
import { useEffect, useMemo, useState } from 'react';
import { sortEventRows, type EventRow, type EventSortKey } from '../../lib/eventRows';
import { formatFraction } from '../../lib/format';

const PAGE_SIZE = 100;

export function EventsTable({ rows, selectedId, onSelect, hasDenominators }: {
  rows: EventRow[];
  selectedId?: string;
  onSelect: (row: EventRow) => void;
  hasDenominators: boolean;
}) {
  const [sortKey, setSortKey] = useState<EventSortKey>('position');
  const [descending, setDescending] = useState(false);
  const [page, setPage] = useState(0);
  const sorted = useMemo(() => sortEventRows(rows, sortKey, descending), [descending, rows, sortKey]);
  const pageCount = Math.max(1, Math.ceil(sorted.length / PAGE_SIZE));
  const safePage = Math.min(page, pageCount - 1);
  const visible = sorted.slice(safePage * PAGE_SIZE, (safePage + 1) * PAGE_SIZE);

  useEffect(() => setPage(0), [rows]);

  // Keep the selected row on screen when it is chosen elsewhere (map, drawer).
  useEffect(() => {
    if (!selectedId) return;
    const index = sorted.findIndex((row) => row.id === selectedId);
    if (index >= 0) setPage(Math.floor(index / PAGE_SIZE));
    window.requestAnimationFrame(() => document.getElementById(rowDomId(selectedId))?.scrollIntoView({ block: 'nearest' }));
  }, [selectedId, sorted]);

  function sortBy(key: EventSortKey) {
    if (key === sortKey) {
      setDescending((value) => !value);
    } else {
      setSortKey(key);
      setDescending(key === 'hf' || key === 'support' || key === 'callable');
    }
  }

  const header = (key: EventSortKey, label: string, align: 'left' | 'right' = 'left') => (
    <th scope="col" aria-sort={sortKey === key ? (descending ? 'descending' : 'ascending') : 'none'} className={align === 'right' ? 'text-right' : undefined}>
      <button type="button" onClick={() => sortBy(key)} className={'inline-flex items-center gap-1 uppercase hover:text-text ' + (sortKey === key ? 'text-text' : '')}>
        {label}
        {sortKey !== key ? <ChevronsUpDown className="h-3 w-3 opacity-60" aria-hidden /> : descending ? <ArrowDown className="h-3 w-3" aria-hidden /> : <ArrowUp className="h-3 w-3" aria-hidden />}
      </button>
    </th>
  );

  return (
    <div className="min-w-0">
      <div className="data-scroll">
        <table className="data-table min-w-[980px]">
          <thead>
            <tr>
              {header('position', 'Position')}
              <th scope="col">Change</th>
              {header('type', 'Type')}
              {header('gene', 'Gene')}
              <th scope="col">Consequence</th>
              {header('hf', 'HF')}
              {header('support', hasDenominators ? 'ALT / callable' : 'Support', 'right')}
              {hasDenominators && <th scope="col">95% CI</th>}
              {hasDenominators && <th scope="col">NUMT</th>}
              <th scope="col">QC</th>
              <th scope="col">Clinical</th>
            </tr>
          </thead>
          <tbody>
            {visible.length === 0 ? (
              <tr>
                <td colSpan={11} className="py-10 text-center text-muted">No events match the current filters.</td>
              </tr>
            ) : visible.map((row) => {
              const selected = row.id === selectedId;
              return (
                <tr
                  key={row.id}
                  id={rowDomId(row.id)}
                  onClick={() => onSelect(row)}
                  className={'cursor-pointer ' + (selected ? 'bg-accent/[0.07] shadow-[inset_3px_0_0_rgb(var(--color-accent))]' : 'hover:bg-panel2/70')}
                >
                  <td>
                    <button
                      type="button"
                      onClick={(event) => { event.stopPropagation(); onSelect(row); }}
                      aria-pressed={selected}
                      aria-label={`Inspect ${row.typeLabel} at ${row.position}`}
                      className="font-mono font-semibold tabular-nums text-accent hover:underline"
                    >
                      {row.position.toLocaleString('en-US')}
                    </button>
                  </td>
                  <td className="max-w-[160px] truncate font-mono text-[12px]" title={row.change}>{row.change}</td>
                  <td><TypeBadge row={row} /></td>
                  <td className="max-w-[150px] truncate" title={row.genes.join(', ')}>
                    {row.genes.length === 0 ? <span className="text-subtle">intergenic</span> : row.genes.length > 2 ? `${row.genes.slice(0, 2).join(', ')} +${row.genes.length - 2}` : row.genes.join(', ')}
                  </td>
                  <td className="max-w-[160px] truncate text-muted" title={row.consequence}>{row.consequence ?? '—'}</td>
                  <td><FrequencyCell row={row} /></td>
                  <td className="text-right font-mono text-[12px] tabular-nums">
                    {row.callable !== undefined ? <>{row.altSupport.toLocaleString('en-US')}<span className="text-subtle"> / {row.callable.toLocaleString('en-US')}</span></> : <>{row.altSupport.toLocaleString('en-US')}<span className="text-subtle"> {row.kind === 'sv' ? 'mol.' : 'rec.'}</span></>}
                  </td>
                  {hasDenominators && (
                    <td className="whitespace-nowrap font-mono text-[12px] tabular-nums text-muted">
                      {row.ciLow !== undefined && row.ciHigh !== undefined ? `${formatFraction(row.ciLow)}–${formatFraction(row.ciHigh)}` : '—'}
                    </td>
                  )}
                  {hasDenominators && <td><NumtBadge value={row.numt} /></td>}
                  <td>
                    {row.qcFlags.length === 0
                      ? <span className="text-subtle">—</span>
                      : <span className="badge-warning" title={row.qcFlags.join('\n')}>{row.qcFlags.length} flag{row.qcFlags.length === 1 ? '' : 's'}</span>}
                  </td>
                  <td className="max-w-[150px] truncate">
                    {row.clinical ? <span className="badge-neutral" title="Source assertion summary, not a clinical verdict">{row.clinical}</span> : <span className="text-subtle">—</span>}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
      {pageCount > 1 && (
        <div className="flex items-center justify-between gap-3 border-t border-line px-4 py-2.5 text-[13px] text-muted">
          <span className="tabular-nums">
            {(safePage * PAGE_SIZE + 1).toLocaleString('en-US')}–{Math.min(sorted.length, (safePage + 1) * PAGE_SIZE).toLocaleString('en-US')} of {sorted.length.toLocaleString('en-US')}
          </span>
          <div className="flex items-center gap-2">
            <button type="button" className="btn-secondary min-h-8" disabled={safePage === 0} onClick={() => setPage(safePage - 1)}>Previous</button>
            <span className="tabular-nums" aria-live="polite">{safePage + 1} / {pageCount}</span>
            <button type="button" className="btn-secondary min-h-8" disabled={safePage >= pageCount - 1} onClick={() => setPage(safePage + 1)}>Next</button>
          </div>
        </div>
      )}
    </div>
  );
}

export function TypeBadge({ row }: { row: EventRow }) {
  if (row.kind === 'sv') {
    return <span className="badge border-sv/30 bg-sv/10 text-text"><span className="h-1.5 w-1.5 rounded-full bg-sv" aria-hidden />{row.typeLabel}</span>;
  }
  return <span className="badge border-variant/30 bg-variant/10 text-text"><span className="h-1.5 w-1.5 rounded-full bg-variant" aria-hidden />{row.typeLabel}</span>;
}

function FrequencyCell({ row }: { row: EventRow }) {
  if (row.hf === undefined) {
    return <span className="text-subtle" title={row.kind === 'sv' ? 'Structural support has no negative-evidence denominator' : 'No callable denominator: support is not a frequency'}>—</span>;
  }
  const low = row.ciLow ?? row.hf;
  const high = row.ciHigh ?? row.hf;
  return (
    <span className="flex items-center gap-2" title={`HF ${formatFraction(row.hf)} (95% CI ${formatFraction(low)}–${formatFraction(high)})`}>
      <span className="relative h-1.5 w-14 shrink-0 overflow-hidden rounded-full bg-panel2" aria-hidden>
        <span className="absolute inset-y-0 rounded-full bg-variant/30" style={{ left: `${low * 100}%`, width: `${Math.max(0.5, (high - low) * 100)}%` }} />
        <span className="absolute inset-y-0 left-0 rounded-full bg-variant" style={{ width: `${Math.max(1, row.hf * 100)}%` }} />
      </span>
      <span className="font-mono text-[12px] tabular-nums">{formatFraction(row.hf)}</span>
    </span>
  );
}

function NumtBadge({ value }: { value?: string }) {
  if (!value) return <span className="text-subtle">—</span>;
  if (value === 'ASSESSABLE') return <span className="badge-success">Assessable</span>;
  if (value === 'NOT_ASSESSABLE') return <span className="badge-warning">Not assessable</span>;
  return <span className="badge-neutral">{value}</span>;
}

export function rowDomId(id: string): string {
  return 'event-row-' + id.replace(/[^a-zA-Z0-9_-]/g, '_');
}
