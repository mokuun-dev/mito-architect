import type { EvidenceObservation } from '@mito-architect/visualization-lib';
import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { useQuery } from '@tanstack/react-query';
import { useMemo } from 'react';
import { useState } from 'react';
import { FixedSizeList } from 'react-window';
import type { ListChildComponentProps } from 'react-window';
import type { EvidenceSearchFilters } from '../../lib/api';
import type { EvidenceSearchRow } from '../../lib/api';
import { searchEvidence } from '../../lib/api';
import { normalizeExactEvidenceFilter } from '../../lib/analysisNavigation';
import { EmptyRow } from './shared';

export const REMOTE_SEARCH_LIMIT = 100;

export function RemoteEvidenceSearch({
  jobId,
  moleculeIds,
  eventIds
}: {
  jobId: string;
  moleculeIds: string[];
  eventIds: string[];
}) {
  const moleculeSuggestions = useMemo(() => moleculeIds.slice(0, 500), [moleculeIds]);
  const eventSuggestions = useMemo(() => eventIds.slice(0, 500), [eventIds]);
  const [draft, setDraft] = useState<EvidenceSearchFilters>({});
  const [filters, setFilters] = useState<EvidenceSearchFilters>({});
  const [cursor, setCursor] = useState(0);
  const search = useQuery({
    queryKey: ['evidence-search', jobId, filters.moleculeId ?? '', filters.eventId ?? '', filters.state ?? '', cursor],
    queryFn: ({ signal }) => searchEvidence(jobId, filters, cursor, REMOTE_SEARCH_LIMIT, signal),
    staleTime: Number.POSITIVE_INFINITY,
    gcTime: 5 * 60 * 1000
  });
  const rows = search.data?.rows ?? [];
  const total = search.data?.total_matches ?? 0;
  const firstShown = rows.length === 0 ? 0 : cursor + 1;
  const lastShown = cursor + rows.length;
  const applyFilters = () => {
    setCursor(0);
    setFilters({
      moleculeId: normalizeExactEvidenceFilter(draft.moleculeId),
      eventId: normalizeExactEvidenceFilter(draft.eventId),
      state: draft.state || undefined
    });
  };
  return (
    <section className="overflow-hidden rounded-lg border border-line" aria-labelledby="global-evidence-search-title">
      <div className="border-b border-line bg-panel2 px-3 py-3">
        <h3 id="global-evidence-search-title" className="text-sm font-semibold">Global evidence search</h3>
        <p className="mt-0.5 text-xs text-muted">
          Exact filters are evaluated over every immutable page by the server index; empty filters browse all stored observations. Suggestions are bounded to 500 IDs, but any exact ID is accepted.
        </p>
        <form
          className="mt-3 grid gap-2 md:grid-cols-[1fr_1fr_0.7fr_auto_auto]"
          onSubmit={(event) => { event.preventDefault(); applyFilters(); }}
        >
          <label className="grid gap-1 text-xs text-muted">
            Molecule ID
            <input
              type="text"
              list="evidence-molecule-ids"
              value={draft.moleculeId ?? ''}
              onChange={(event) => setDraft((current) => ({ ...current, moleculeId: event.target.value }))}
              className="min-h-10 rounded-lg border border-line bg-shell px-3 py-2 text-sm text-text"
              placeholder="exact molecule ID"
            />
          </label>
          <datalist id="evidence-molecule-ids">
            {moleculeSuggestions.map((id) => <option key={id} value={id} />)}
          </datalist>
          <label className="grid gap-1 text-xs text-muted">
            Event ID
            <input
              type="text"
              list="evidence-event-ids"
              value={draft.eventId ?? ''}
              onChange={(event) => setDraft((current) => ({ ...current, eventId: event.target.value }))}
              className="min-h-10 rounded-lg border border-line bg-shell px-3 py-2 font-mono text-sm text-text"
              placeholder="exact event ID"
            />
          </label>
          <datalist id="evidence-event-ids">
            {eventSuggestions.map((id) => <option key={id} value={id} />)}
          </datalist>
          <label className="grid gap-1 text-xs text-muted">
            Stored state
            <select
              value={draft.state ?? ''}
              onChange={(event) => setDraft((current) => ({ ...current, state: event.target.value || undefined }))}
              className="min-h-10 rounded-lg border border-line bg-shell px-3 py-2 text-sm text-text"
            >
              <option value="">all states</option>
              {['ALTERNATE', 'REFERENCE', 'EVENT_ABSENT', 'LOW_QUALITY', 'CONFLICT'].map((state) => (
                <option key={state} value={state}>{state}</option>
              ))}
            </select>
          </label>
          <button type="submit" className="min-h-10 self-end rounded-lg border border-aqua bg-aqua/10 px-3 py-2 text-sm font-semibold text-aqua">
            Apply
          </button>
          <button
            type="button"
            className="min-h-10 self-end rounded-lg border border-line px-3 py-2 text-sm text-muted hover:text-text"
            onClick={() => { setDraft({}); setFilters({}); setCursor(0); }}
          >
            Clear
          </button>
        </form>
      </div>
      <div className="flex flex-col gap-2 border-b border-line px-3 py-2 text-xs text-muted sm:flex-row sm:items-center sm:justify-between" aria-live="polite">
        <span>{search.isFetching ? 'Searching…' : `${firstShown}–${lastShown} of ${total} exact matches`}</span>
        <span className="flex flex-wrap gap-2">
          <button
            type="button"
            className="min-h-9 rounded-lg border border-line px-2 py-1 disabled:opacity-40"
            disabled={cursor === 0 || search.isFetching}
            onClick={() => setCursor(Math.max(0, cursor - REMOTE_SEARCH_LIMIT))}
          >
            Previous matches
          </button>
          <button
            type="button"
            className="min-h-9 rounded-lg border border-line px-2 py-1 disabled:opacity-40"
            disabled={search.data?.next_cursor == null || search.isFetching}
            onClick={() => setCursor(search.data?.next_cursor ?? cursor)}
          >
            Next matches
          </button>
        </span>
      </div>
      {search.isLoading ? <EmptyRow text="Building the first bounded result page…" /> : search.isError ? (
        <div role="alert" className="px-3 py-4 text-sm text-coral">{(search.error as Error).message}</div>
      ) : rows.length === 0 ? <EmptyRow text="No stored observation matches the exact filters; this does not imply a reference state" /> : (
        <div className="data-scroll">
          <div className="min-w-[900px]">
            <div className="grid grid-cols-[0.55fr_1fr_1fr_0.7fr_0.9fr_1.1fr] gap-2 border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
              <div>Page:row</div><div>Molecule</div><div>Event</div><div>State</div><div>Alignment</div><div>Evidence</div>
            </div>
            <FixedSizeList
              height={Math.min(320, Math.max(42, rows.length * 42))}
              width="100%"
              itemCount={rows.length}
              itemSize={42}
              itemData={rows}
              itemKey={(index, items) => items[index].id}
            >
              {EvidenceSearchResultRow}
            </FixedSizeList>
          </div>
        </div>
      )}
    </section>
  );
}

export function EvidenceSearchResultRow({ index, style, data }: ListChildComponentProps<EvidenceSearchRow[]>) {
  const observation = data[index];
  return (
    <div style={style} className="grid grid-cols-[0.55fr_1fr_1fr_0.7fr_0.9fr_1.1fr] items-center gap-2 border-b border-line px-3 text-sm">
      <div className="tabular-nums text-muted">{observation.page_index}:{observation.row_index}</div>
      <div className="truncate font-medium" title={observation.molecule_id}>{observation.molecule_id}</div>
      <div className="truncate" title={observation.event_id}>{observation.event_id}</div>
      <div>{observation.state}</div>
      <div className="truncate" title={observation.alignment_id}>{observation.alignment_id}</div>
      <div className="truncate text-muted" title={observation.evidence_source}>{observation.evidence_source}</div>
    </div>
  );
}

export function RemoteEvidenceBrowser({
  page,
  pageIndex,
  pageCount,
  loading,
  error,
  onPageChange
}: {
  page?: NonNullable<MitoAnalysisData['observation_pages']>[number];
  pageIndex: number;
  pageCount: number;
  loading: boolean;
  error: Error | null;
  onPageChange: (page: number) => void;
}) {
  const rows = useMemo(() => page ? materializeObservationPage(page) : [], [page]);
  const lastPage = Math.max(0, pageCount - 1);
  return (
    <section className="overflow-hidden rounded-lg border border-line" aria-labelledby="remote-evidence-title">
      <div className="flex flex-wrap items-center justify-between gap-3 border-b border-line bg-panel2 px-3 py-2">
        <div>
          <h3 id="remote-evidence-title" className="text-sm font-semibold">Remote evidence pages</h3>
          <p className="mt-0.5 text-xs text-muted">Only the open immutable page is allocated and rendered in the browser.</p>
        </div>
        <div className="flex flex-wrap items-center gap-2 text-xs" aria-label="Evidence page navigation">
          <button type="button" className="min-h-9 rounded-lg border border-line px-2 py-1 disabled:opacity-40" aria-label="First evidence page" disabled={pageIndex === 0 || pageCount === 0} onClick={() => onPageChange(0)}>First</button>
          <button type="button" className="min-h-9 rounded-lg border border-line px-2 py-1 disabled:opacity-40" aria-label="Previous evidence page" disabled={pageIndex === 0 || pageCount === 0} onClick={() => onPageChange(pageIndex - 1)}>Prev</button>
          <label className="flex items-center gap-1 text-muted">
            Page
            <input
              className="min-h-9 w-16 rounded-lg border border-line bg-shell px-2 py-1 text-right text-text"
              type="number"
              min={1}
              max={Math.max(1, pageCount)}
              value={pageCount === 0 ? 0 : pageIndex + 1}
              disabled={pageCount === 0}
              aria-label="Evidence page number"
              onChange={(event) => {
                const requested = Number.parseInt(event.target.value, 10);
                if (Number.isFinite(requested)) onPageChange(Math.min(lastPage, Math.max(0, requested - 1)));
              }}
            />
            / {pageCount}
          </label>
          <button type="button" className="min-h-9 rounded-lg border border-line px-2 py-1 disabled:opacity-40" aria-label="Next evidence page" disabled={pageIndex >= lastPage || pageCount === 0} onClick={() => onPageChange(pageIndex + 1)}>Next</button>
          <button type="button" className="min-h-9 rounded-lg border border-line px-2 py-1 disabled:opacity-40" aria-label="Last evidence page" disabled={pageIndex >= lastPage || pageCount === 0} onClick={() => onPageChange(lastPage)}>Last</button>
        </div>
      </div>
      {loading ? <EmptyRow text="Loading immutable evidence page…" /> : error ? (
        <div role="alert" className="px-3 py-4 text-sm text-coral">{error.message}</div>
      ) : rows.length === 0 ? <EmptyRow text="Evidence page is empty" /> : (
        <div className="data-scroll">
          <div className="min-w-[780px]">
            <div className="grid grid-cols-[1fr_1fr_0.7fr_0.9fr_1.2fr] gap-2 border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
              <div>Molecule</div><div>Event</div><div>State</div><div>Alignment</div><div>Evidence</div>
            </div>
            <FixedSizeList
              height={Math.min(320, Math.max(42, rows.length * 42))}
              width="100%"
              itemCount={rows.length}
              itemSize={42}
              itemData={rows}
              itemKey={(index, items) => items[index].id}
            >
              {RemoteEvidenceRow}
            </FixedSizeList>
          </div>
        </div>
      )}
    </section>
  );
}

export function RemoteEvidenceRow({ index, style, data }: ListChildComponentProps<EvidenceObservation[]>) {
  const observation = data[index];
  return (
    <div style={style} className="grid grid-cols-[1fr_1fr_0.7fr_0.9fr_1.2fr] items-center gap-2 border-b border-line px-3 text-sm">
      <div className="truncate font-medium" title={observation.molecule_id}>{observation.molecule_id}</div>
      <div className="truncate" title={observation.event_id}>{observation.event_id}</div>
      <div>{observation.state}</div>
      <div className="truncate" title={observation.alignment_id}>{observation.alignment_id}</div>
      <div className="truncate text-muted" title={observation.evidence_source}>{observation.evidence_source}</div>
    </div>
  );
}

export function materializeObservationPage(
  page: NonNullable<MitoAnalysisData['observation_pages']>[number]
): EvidenceObservation[] {
  return Array.from({ length: page.count }, (_, index) => ({
    id: `observation:${page.offset + index}`,
    molecule_id: page.columns.molecule_id[index],
    event_id: page.columns.event_id[index],
    alignment_id: page.columns.alignment_id[index],
    state: page.columns.state[index],
    observed_allele: page.columns.observed_allele[index],
    base_quality: page.columns.base_quality[index],
    mapping_quality: page.columns.mapping_quality[index],
    strand: page.columns.strand[index],
    evidence_source: page.columns.evidence_source[index],
    read_position: page.columns.read_position[index]
  }));
}
