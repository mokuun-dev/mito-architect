import type { EvidenceObservation } from '@mito-architect/visualization-lib';
import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { useQuery } from '@tanstack/react-query';
import { useEffect } from 'react';
import { useMemo } from 'react';
import { useState } from 'react';
import { getEvidencePage } from '../../lib/api';
import { searchEvidence } from '../../lib/api';
import { useMitoStore } from '../../lib/store';
import { RemoteEvidenceSearch } from './RemoteEvidence';
import { RemoteEvidenceBrowser } from './RemoteEvidence';
import { Metric } from './shared';
import { EmptyRow } from './shared';
import { formatPercent } from './shared';

export function EvidencePanel({ data, jobId }: { data: MitoAnalysisData; jobId?: string }) {
  const callability = data.callability ?? [];
  const events = data.events ?? [];
  const phaseLinks = data.phase_links ?? [];
  const remoteMoleculeIds = useMemo(
    () => (data.molecules ?? []).slice(0, 500).map((molecule) => molecule.id),
    [data.molecules]
  );
  const remoteEventIds = useMemo(
    () => (data.events ?? []).slice(0, 500).map((event) => event.id),
    [data.events]
  );
  const moleculesById = new Map((data.molecules ?? []).map((molecule) => [molecule.id, molecule]));
  const selectedEventId = useMitoStore((state) => state.selectedEventId);
  const selectedMoleculeId = useMitoStore((state) => state.selectedMoleculeId);
  const selectedPhaseId = useMitoStore((state) => state.selectedPhaseId);
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);
  const setSelectedMoleculeId = useMitoStore((state) => state.setSelectedMoleculeId);
  const setSelectedPhaseId = useMitoStore((state) => state.setSelectedPhaseId);
  const remotePages = data.evidence_encoding?.observation_storage === 'remote_http_pages';
  const pageCount = data.evidence_encoding?.observation_page_count ?? 0;
  const [remotePageIndex, setRemotePageIndex] = useState(0);
  const remotePage = useQuery({
    queryKey: ['evidence-page', jobId, remotePageIndex],
    queryFn: ({ signal }) => getEvidencePage(jobId!, remotePageIndex, signal),
    enabled: remotePages && Boolean(jobId) && pageCount > 0,
    staleTime: Number.POSITIVE_INFINITY,
    gcTime: 5 * 60 * 1000
  });
  useEffect(() => {
    if (remotePageIndex >= pageCount) setRemotePageIndex(Math.max(0, pageCount - 1));
  }, [pageCount, remotePageIndex]);
  const selection = selectedPhaseId
    ? { kind: 'phase' as const, id: selectedPhaseId }
    : selectedEventId
      ? { kind: 'event' as const, id: selectedEventId }
      : selectedMoleculeId
        ? { kind: 'molecule' as const, id: selectedMoleculeId }
        : null;
  const selected = selection && (
    selection.kind === 'event'
      ? events.some((event) => event.id === selection.id)
      : selection.kind === 'molecule'
        ? callability.some((item) => item.molecule_id === selection.id)
        : phaseLinks.some((link) => link.id === selection.id)
  ) ? selection : null;
  const selectedRemoteObservations = useQuery({
    queryKey: ['evidence-selection', jobId, selected?.kind, selected?.id],
    queryFn: async ({ signal }) => {
      if (!selected || !jobId) return [];
      if (selected.kind === 'event') {
        return (await searchEvidence(jobId, { eventId: selected.id }, 0, 300, signal)).rows;
      }
      if (selected.kind === 'molecule') {
        return (await searchEvidence(jobId, { moleculeId: selected.id }, 0, 300, signal)).rows;
      }
      const phase = phaseLinks.find((link) => link.id === selected.id);
      if (!phase) return [];
      const responses = await Promise.all([
        searchEvidence(jobId, { eventId: phase.event_a_id }, 0, 300, signal),
        searchEvidence(jobId, { eventId: phase.event_b_id }, 0, 300, signal)
      ]);
      return Array.from(
        new Map(responses.flatMap((response) => response.rows).map((row) => [row.id, row])).values()
      )
        .sort((left, right) => left.page_index - right.page_index || left.row_index - right.row_index)
        .slice(0, 300);
    },
    enabled: remotePages && Boolean(jobId) && selected !== null,
    staleTime: Number.POSITIVE_INFINITY,
    gcTime: 5 * 60 * 1000
  });
  const selectedObservations = selected
    ? remotePages
      ? selectedRemoteObservations.data ?? []
      : matchingObservations(data, selected, 300)
    : [];
  const observationCount = data.evidence_encoding?.observation_count ?? data.observations?.length ?? 0;
  if (data.metadata.schema_version !== '0.6') {
    return (
      <div className="rounded-md border border-amber/50 bg-amber/10 p-4 text-sm text-amber">
        Molecule evidence is not present in this schema {data.metadata.schema_version ?? 'unknown'} result.
        Run analysis with the opt-in schema 0.6 evidence graph; a missing projection is not evidence that
        no molecular linkage exists.
      </div>
    );
  }

  return (
    <div className="grid gap-4">
      <div className="grid gap-3 sm:grid-cols-2 xl:grid-cols-4">
        <Metric label="Assembled molecules" value={(data.molecules?.length ?? 0).toString()} />
        <Metric label="Normalized events" value={events.length.toString()} />
        <Metric label="Phase links" value={phaseLinks.length.toString()} />
        <Metric label="Evidence rows" value={observationCount.toString()} />
      </div>

      {remotePages && (
        <>
          {jobId ? (
            <RemoteEvidenceSearch
              jobId={jobId}
              moleculeIds={remoteMoleculeIds}
              eventIds={remoteEventIds}
            />
          ) : (
            <div role="alert" className="rounded border border-amber/50 bg-amber/10 px-3 py-2 text-sm text-amber">
              Remote evidence search requires the originating job identifier.
            </div>
          )}
          <RemoteEvidenceBrowser
            page={remotePage.data}
            pageIndex={remotePageIndex}
            pageCount={pageCount}
            loading={remotePage.isLoading}
            error={remotePage.error as Error | null}
            onPageChange={setRemotePageIndex}
          />
        </>
      )}

      <div className="data-scroll rounded-lg border border-line">
        <div className="grid min-w-[820px] grid-cols-[1.05fr_0.7fr_0.65fr_1.1fr_1.4fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <div>Molecule</div><div>Status</div><div>Callable</div><div>Protocol decision</div><div>Alignment provenance</div>
        </div>
        {callability.length === 0 ? <EmptyRow text="Callability projection is empty" /> : callability.slice(0, 100).map((item) => {
          const molecule = moleculesById.get(item.molecule_id);
          const protocolDecision = [
            molecule?.identity_policy,
            ...(molecule?.protocol_flags ?? []),
            ...(molecule?.exclusion_reasons ?? []).map((reason) => `excluded:${reason}`)
          ].filter(Boolean).join(', ');
          return (
            <button
              key={item.molecule_id}
              type="button"
              aria-pressed={selected?.kind === 'molecule' && selected.id === item.molecule_id}
              onClick={() => {
                setSelectedEventId(undefined);
                setSelectedPhaseId(undefined);
                setSelectedMoleculeId(item.molecule_id);
              }}
              className="grid min-w-[820px] w-full grid-cols-[1.05fr_0.7fr_0.65fr_1.1fr_1.4fr] gap-2 border-b border-line px-3 py-2 text-left text-sm last:border-0 hover:bg-panel2 aria-pressed:bg-panel2"
            >
              <div className="truncate font-semibold" title={item.molecule_id}>{item.molecule_id}</div>
              <div className={molecule?.evidence_eligible && item.known ? 'text-aqua' : 'text-amber'}>{item.status}</div>
              <div>{item.known ? formatPercent(item.callable_fraction) : 'not assessable'}</div>
              <div className="truncate text-muted" title={protocolDecision || 'no protocol metadata'}>
                {protocolDecision || 'default identity'}
              </div>
              <div className="truncate text-muted" title={item.alignments.map((alignment) => `${alignment.alignment_id}:${alignment.status}`).join(', ')}>
                {item.alignments.map((alignment) => `${alignment.alignment_id} ${alignment.status}`).join(', ') || 'no eligible alignment'}
              </div>
            </button>
          );
        })}
      </div>

      <div className="data-scroll rounded-lg border border-line">
        <div className="grid min-w-[720px] grid-cols-[1.35fr_0.8fr_0.65fr_0.65fr_0.95fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <div>Event</div><div>Type</div><div>ALT</div><div>Callable</div><div>Assessability</div>
        </div>
        {events.length === 0 ? <EmptyRow text="No normalized events" /> : events.slice(0, 200).map((event) => (
          <button
            key={event.id}
            type="button"
            aria-pressed={selected?.kind === 'event' && selected.id === event.id}
            onClick={() => {
              setSelectedMoleculeId(undefined);
              setSelectedPhaseId(undefined);
              setSelectedEventId(event.id);
            }}
            className="grid min-w-[720px] w-full grid-cols-[1.35fr_0.8fr_0.65fr_0.65fr_0.95fr] gap-2 border-b border-line px-3 py-2 text-left text-sm last:border-0 hover:bg-panel2 aria-pressed:bg-panel2"
          >
            <div className="truncate font-semibold" title={event.id}>{eventLabel(event)}</div>
            <div>{event.type}</div>
            <div>{event.evidence_counts.alternate}</div>
            <div>{event.evidence_counts.callable}</div>
            <div
              className={event.assessability === 'REFERENCE_AND_ALTERNATE' ? 'text-aqua' : 'text-amber'}
              title={event.negative_evidence_rule}
            >
              {event.assessability === 'REFERENCE_AND_ALTERNATE' ? 'complete' : 'support-only'}
            </div>
          </button>
        ))}
      </div>

      <div className="data-scroll rounded-lg border border-line">
        <div className="grid min-w-[1040px] grid-cols-[1fr_1fr_0.48fr_0.48fr_0.48fr_0.48fr_0.65fr_0.8fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <div>Event A</div><div>Event B</div><div>Joint</div><div>A+/B−</div><div>A−/B+</div><div>Both ALT</div><div>Co-ALT CI</div><div>Inspect</div>
        </div>
        {phaseLinks.length === 0 ? <EmptyRow text="No phase links with co-observed alternate evidence" /> : phaseLinks.slice(0, 200).map((link) => (
          <div key={link.id} className="grid min-w-[1040px] grid-cols-[1fr_1fr_0.48fr_0.48fr_0.48fr_0.48fr_0.65fr_0.8fr] gap-2 border-b border-line px-3 py-2 text-sm last:border-0">
            <button type="button" className="truncate text-left hover:text-aqua" title={link.event_a_id} onClick={() => { setSelectedMoleculeId(undefined); setSelectedPhaseId(undefined); setSelectedEventId(link.event_a_id); }}>{link.event_a_id}</button>
            <button type="button" className="truncate text-left hover:text-aqua" title={link.event_b_id} onClick={() => { setSelectedMoleculeId(undefined); setSelectedPhaseId(undefined); setSelectedEventId(link.event_b_id); }}>{link.event_b_id}</button>
            <div>{link.jointly_callable}</div>
            <div>{link.a_alternate_b_absent}</div>
            <div>{link.a_absent_b_alternate}</div>
            <div title={phaseMoleculeTitle(link.supporting_molecule_indices, data)}>{link.both_alternate}</div>
            <div title={`uncertain n=${link.jointly_uncertain}; delta=${link.linkage_delta.toFixed(4)}`}>{formatPercent(link.co_alternate_ci95_low)}–{formatPercent(link.co_alternate_ci95_high)}</div>
            <button type="button" aria-pressed={selected?.kind === 'phase' && selected.id === link.id} className={link.assessability === 'COMPLETE_FOR_BOTH_EVENTS' ? 'text-left text-aqua hover:underline' : 'text-left text-amber hover:underline'} onClick={() => { setSelectedEventId(undefined); setSelectedMoleculeId(undefined); setSelectedPhaseId(link.id); }}>
              {link.assessability === 'COMPLETE_FOR_BOTH_EVENTS' ? 'complete pair' : 'conditioned'}
            </button>
          </div>
        ))}
      </div>
      <div className="data-scroll rounded-lg border border-line">
        <div className="flex items-center justify-between gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <span>{selected ? `${selected.kind}: ${selected.id}` : 'Evidence inspector'}</span>
          {selected && <button type="button" className="text-aqua hover:underline" onClick={() => { setSelectedEventId(undefined); setSelectedMoleculeId(undefined); setSelectedPhaseId(undefined); }}>clear</button>}
        </div>
        {!selected ? <EmptyRow text="Select a molecule, event, or phase endpoint to trace its evidence" /> : selectedRemoteObservations.isLoading ? (
          <EmptyRow text="Searching all immutable evidence pages…" />
        ) : selectedRemoteObservations.isError ? (
          <div role="alert" className="px-3 py-4 text-sm text-coral">
            {(selectedRemoteObservations.error as Error).message}
          </div>
        ) : selectedObservations.length === 0 ? (
          <EmptyRow text={remotePages ? 'No explicit observation in the global evidence index; sparse absence is not REF' : 'No explicit observation; sparse absence means NOT_CALLABLE'} />
        ) : (
          <>
            <div className="grid min-w-[760px] grid-cols-[1fr_1fr_0.65fr_0.8fr_1fr] gap-2 border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
              <div>Molecule</div><div>Event</div><div>State</div><div>Alignment</div><div>Evidence</div>
            </div>
            {selectedObservations.slice(0, 300).map((observation) => (
              <div key={observation.id} className="grid min-w-[760px] grid-cols-[1fr_1fr_0.65fr_0.8fr_1fr] gap-2 border-b border-line px-3 py-2 text-sm last:border-0">
                <button type="button" className="truncate text-left hover:text-aqua" title={observation.molecule_id} onClick={() => { setSelectedEventId(undefined); setSelectedPhaseId(undefined); setSelectedMoleculeId(observation.molecule_id); }}>{observation.molecule_id}</button>
                <button type="button" className="truncate text-left hover:text-aqua" title={observation.event_id} onClick={() => { setSelectedMoleculeId(undefined); setSelectedPhaseId(undefined); setSelectedEventId(observation.event_id); }}>{observation.event_id}</button>
                <div>{observation.state}</div>
                <div className="truncate" title={observation.alignment_id}>{observation.alignment_id}</div>
                <div className="truncate text-muted" title={observation.evidence_source}>{observation.evidence_source}</div>
              </div>
            ))}
          </>
        )}
      </div>
      {(callability.length > 100 || events.length > 200 || phaseLinks.length > 200 || observationCount > 300) && (
        <div className="rounded border border-amber/40 bg-amber/10 px-3 py-2 text-xs text-amber">
          Preview limits: {Math.min(callability.length, 100)}/{callability.length} molecules, {Math.min(events.length, 200)}/{events.length} events, {Math.min(phaseLinks.length, 200)}/{phaseLinks.length} phase links, and {Math.min(observationCount, 300)}/{observationCount} observations. Use exact molecule/event search to query the complete immutable evidence index.
        </div>
      )}
    </div>
  );
}

export function eventLabel(event: NonNullable<MitoAnalysisData['events']>[number]) {
  if (event.type === 'SNV' && event.start && event.ref && event.alt) {
    return `m.${event.start}${event.ref}>${event.alt}`;
  }
  if (event.start && event.end) {
    return `${event.start}-${event.end}`;
  }
  return event.id;
}

export function phaseMoleculeTitle(indices: number[], data: MitoAnalysisData): string {
  if (indices.length === 0) return 'no co-ALT support';
  const molecules = data.molecules ?? [];
  return indices.map((index) => molecules[index]?.id ?? `unresolved molecule index ${index}`).join(', ');
}

export function matchingObservations(
  data: MitoAnalysisData,
  selected: { kind: 'event' | 'molecule' | 'phase'; id: string },
  limit: number,
  pageOverride?: NonNullable<MitoAnalysisData['observation_pages']>
): EvidenceObservation[] {
  const matches: EvidenceObservation[] = [];
  const phase = selected.kind === 'phase'
    ? (data.phase_links ?? []).find((link) => link.id === selected.id)
    : undefined;
  const accept = (observation: EvidenceObservation) => {
    const matched = selected.kind === 'event'
      ? observation.event_id === selected.id
      : selected.kind === 'molecule'
        ? observation.molecule_id === selected.id
        : phase !== undefined && (observation.event_id === phase.event_a_id || observation.event_id === phase.event_b_id);
    if (matched && matches.length < limit) matches.push(observation);
  };
  if (data.observations) {
    for (const observation of data.observations) {
      accept(observation);
      if (matches.length === limit) return matches;
    }
    return matches;
  }
  for (const page of pageOverride ?? data.observation_pages ?? []) {
    for (let index = 0; index < page.count; index += 1) {
      accept({
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
      });
      if (matches.length === limit) return matches;
    }
  }
  return matches;
}
