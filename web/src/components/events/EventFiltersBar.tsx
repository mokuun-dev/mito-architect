import { Search, X } from 'lucide-react';
import { useMemo } from 'react';
import { DEFAULT_EVENT_FILTERS, type EventKind, type EventRow } from '../../lib/eventRows';
import { useMitoStore } from '../../lib/store';

const KINDS: Array<[EventKind, string]> = [
  ['snv', 'SNV'],
  ['indel', 'Indel'],
  ['sv', 'Structural']
];

const HF_THRESHOLDS = [0, 0.01, 0.05, 0.1, 0.5, 0.9];

export function EventFiltersBar({ rows, shown }: { rows: EventRow[]; shown: number }) {
  const filters = useMitoStore((state) => state.eventFilters);
  const setFilters = useMitoStore((state) => state.setEventFilters);
  const reset = useMitoStore((state) => state.resetEventFilters);
  const genes = useMemo(() => [...new Set(rows.flatMap((row) => row.genes))].sort(), [rows]);
  const hasFrequencies = rows.some((row) => row.hf !== undefined);
  const hasFlags = rows.some((row) => row.qcFlags.length > 0);
  const active = JSON.stringify(filters) !== JSON.stringify(DEFAULT_EVENT_FILTERS);

  return (
    <div className="flex flex-wrap items-center gap-2" role="search" aria-label="Filter events">
      <div className="relative w-full sm:w-64">
        <Search className="pointer-events-none absolute left-2.5 top-1/2 h-4 w-4 -translate-y-1/2 text-subtle" aria-hidden />
        <label htmlFor="event-search" className="sr-only">Search events</label>
        <input
          id="event-search"
          value={filters.query}
          onChange={(event) => setFilters({ query: event.target.value })}
          placeholder="Position, allele, gene, QC flag…"
          className="input pl-8"
          autoComplete="off"
        />
      </div>
      <div className="flex flex-wrap items-center gap-1.5" role="group" aria-label="Event type">
        {KINDS.map(([kind, label]) => {
          const pressed = filters.kinds.includes(kind);
          return (
            <button
              key={kind}
              type="button"
              className="chip"
              aria-pressed={pressed}
              onClick={() => setFilters({ kinds: pressed ? filters.kinds.filter((item) => item !== kind) : [...filters.kinds, kind] })}
            >
              {label}
            </button>
          );
        })}
      </div>
      <label className="sr-only" htmlFor="event-gene">Gene</label>
      <select
        id="event-gene"
        value={filters.gene ?? ''}
        onChange={(event) => setFilters({ gene: event.target.value || undefined })}
        className="input w-auto pr-8"
      >
        <option value="">All genes</option>
        {genes.map((gene) => <option key={gene} value={gene}>{gene}</option>)}
      </select>
      {hasFrequencies && (
        <>
          <label className="sr-only" htmlFor="event-hf">Minimum heteroplasmic fraction</label>
          <select
            id="event-hf"
            value={filters.minHf}
            onChange={(event) => setFilters({ minHf: Number(event.target.value) })}
            className="input w-auto pr-8"
            title="Rows without a callable denominator are hidden when a minimum is set"
          >
            {HF_THRESHOLDS.map((threshold) => (
              <option key={threshold} value={threshold}>{threshold === 0 ? 'Any HF' : `HF ≥ ${threshold * 100}%`}</option>
            ))}
          </select>
        </>
      )}
      {hasFlags && (
        <button type="button" className="chip" aria-pressed={filters.hideFlagged} onClick={() => setFilters({ hideFlagged: !filters.hideFlagged })}>
          Hide QC-flagged
        </button>
      )}
      <span className="ml-auto text-[13px] tabular-nums text-muted" aria-live="polite">
        {shown.toLocaleString('en-US')} of {rows.length.toLocaleString('en-US')} events
      </span>
      {active && (
        <button type="button" className="btn-ghost min-h-8 px-2" onClick={reset}>
          <X className="h-4 w-4" aria-hidden />
          Clear filters
        </button>
      )}
    </div>
  );
}
