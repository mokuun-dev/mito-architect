import type { AlignmentFragment } from '@mito-architect/visualization-lib';
import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { useState } from 'react';
import { useMitoStore } from '../../lib/store';
import { Metric } from './shared';
import { EmptyRow } from './shared';

export function RearrangementsPanel({ data }: { data: MitoAnalysisData }) {
  const paths = data.complex_events ?? [];
  const eventsById = new Map((data.events ?? []).map((event) => [event.id, event]));
  const [selectedPathId, setSelectedPathId] = useState<string | undefined>();
  const selectedPath = paths.find((path) => path.id === selectedPathId);
  const selectedEvent = selectedPath
    ? eventsById.get(selectedPath.event_id ?? selectedPath.id)
    : undefined;
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);
  const setSelectedMoleculeId = useMitoStore((state) => state.setSelectedMoleculeId);
  const setSelectedPhaseId = useMitoStore((state) => state.setSelectedPhaseId);
  const setSelectedSvId = useMitoStore((state) => state.setSelectedSvId);
  const supportingMoleculeIds = selectedPath
    ? selectedEvent?.supporting_molecule_ids ?? selectedPath.supporting_reads
    : [];
  const supportingMoleculeSet = new Set(supportingMoleculeIds);
  const supportingAlignments = (data.alignments ?? []).filter((alignment) =>
    supportingMoleculeSet.has(alignment.molecule_id)
  );
  const inspectEvent = (eventId: string) => {
    setSelectedMoleculeId(undefined);
    setSelectedPhaseId(undefined);
    setSelectedEventId(eventId);
  };
  const inspectMolecule = (moleculeId: string) => {
    setSelectedEventId(undefined);
    setSelectedPhaseId(undefined);
    setSelectedMoleculeId(moleculeId);
  };

  return (
    <div className="grid gap-4">
      <div className="grid gap-3 sm:grid-cols-3">
        <Metric label="Structural events" value={data.svs.length.toString()} />
        <Metric label="Multi-junction paths" value={paths.length.toString()} />
        <Metric
          label="Path-support molecules"
          value={new Set(paths.flatMap((path) => eventsById.get(path.event_id ?? path.id)?.supporting_molecule_ids ?? path.supporting_reads)).size.toString()}
        />
      </div>

      <section className="data-scroll rounded-lg border border-line" aria-labelledby="complex-paths-title">
        <div className="border-b border-line bg-panel2 px-3 py-2">
          <h3 id="complex-paths-title" className="text-sm font-semibold">Observed multi-junction paths</h3>
          <p className="mt-0.5 text-xs text-muted">
            Ordered split-alignment evidence only. Missing junctions and circular closure are never inferred.
          </p>
        </div>
        <div className="grid min-w-[720px] grid-cols-[1.6fr_0.55fr_0.55fr_0.7fr_0.75fr] gap-2 border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <div>Canonical path</div><div>Junctions</div><div>Segments</div><div>Molecules</div><div>Inspect</div>
        </div>
        {paths.length === 0 ? <EmptyRow text="No molecule contains two or more observed split-alignment junctions" /> : paths.map((path) => {
          const event = eventsById.get(path.event_id ?? path.id);
          const support = event?.supporting_molecule_ids ?? path.supporting_reads;
          const selected = selectedPath?.id === path.id;
          return (
            <button
              key={path.id}
              type="button"
              aria-pressed={selected}
              onClick={() => {
                setSelectedPathId(path.id);
                inspectEvent(path.event_id ?? path.id);
              }}
              className="grid min-w-[720px] w-full grid-cols-[1.6fr_0.55fr_0.55fr_0.7fr_0.75fr] gap-2 border-b border-line px-3 py-2 text-left text-sm last:border-0 hover:bg-panel2 aria-pressed:bg-panel2"
            >
              <span className="truncate font-mono text-xs" title={path.id}>{path.id}</span>
              <span>{path.junction_count}</span>
              <span>{path.segment_count}</span>
              <span>{support.length}</span>
              <span className="text-aqua">{selected ? 'selected' : 'trace path'}</span>
            </button>
          );
        })}
      </section>

      {selectedPath && (
        <section className="overflow-hidden rounded-md border border-aqua/50" aria-labelledby="selected-path-title">
          <div className="flex flex-wrap items-start justify-between gap-3 border-b border-line bg-panel2 px-3 py-3">
            <div className="min-w-0">
              <h3 id="selected-path-title" className="font-semibold">Selected rearrangement path</h3>
              <p className="mt-1 break-all font-mono text-xs text-aqua">{selectedPath.id}</p>
            </div>
            <button
              type="button"
              className="text-sm text-aqua hover:underline"
              onClick={() => { setSelectedPathId(undefined); setSelectedEventId(undefined); }}
            >
              Clear selection
            </button>
          </div>
          <div className="grid gap-4 p-3 lg:grid-cols-[minmax(0,1fr)_minmax(320px,0.8fr)]">
            <div>
              <div className="text-xs font-semibold uppercase tracking-normal text-muted">Ordered junction evidence</div>
              <div className="mt-2 flex items-center gap-1 overflow-x-auto rounded border border-line bg-shell p-2" aria-label="Observed path order">
                {selectedPath.junction_ids.map((junctionId, index) => <span key={`${junctionId}:diagram:${index}`} className="inline-flex shrink-0 items-center gap-1"><span className="rounded border border-magenta/50 bg-magenta/10 px-2 py-1 font-mono text-xs text-magenta">{junctionId}</span>{index + 1 < selectedPath.junction_ids.length && <span className="text-muted">→</span>}</span>)}
                <span className="ml-2 shrink-0 text-xs text-amber">open path; no closure inferred</span>
              </div>
              <ol className="mt-2 grid gap-2" aria-label="Junction traversal order">
                {selectedPath.junction_ids.map((junctionId, index) => {
                  const eventId = `sv:${junctionId}`;
                  const sv = data.svs.find((candidate) => candidate.id === junctionId);
                  return (
                    <li key={`${junctionId}:${index}`} className="grid grid-cols-[auto_1fr_auto] items-center gap-3 rounded border border-line bg-shell px-3 py-2">
                      <span className="grid h-7 w-7 place-items-center rounded-full border border-aqua/50 text-xs text-aqua">{index + 1}</span>
                      <span className="min-w-0">
                        <button
                          type="button"
                          className="block max-w-full truncate text-left font-mono text-xs hover:text-aqua"
                          title={junctionId}
                          onClick={() => {
                            setSelectedSvId(junctionId);
                            inspectEvent(sv?.event_id ?? eventId);
                          }}
                        >
                          {junctionId}
                        </button>
                        <span className="mt-0.5 block text-xs text-muted">
                          {sv ? `${sv.type} · ${sv.start}–${sv.end} · ${sv.length} bp` : 'component projection unavailable'}
                        </span>
                      </span>
                      <span className="font-mono text-xs text-amber" title="Observed strand transition">
                        {selectedPath.junction_orientations[index] ?? 'unknown'}
                      </span>
                    </li>
                  );
                })}
              </ol>
              <div className="mt-3 rounded border border-amber/50 bg-amber/10 px-3 py-2 text-xs text-amber">
                Canonicalization: {selectedPath.canonicalization}. This is an observed open traversal; rotation equivalence and molecular closure are not claimed.
              </div>
            </div>
            <div className="grid content-start gap-3">
              <TraceList
                label="Supporting molecules"
                values={supportingMoleculeIds}
                empty="No molecule reference was retained"
                onSelect={inspectMolecule}
              />
              <div className="rounded border border-line bg-shell p-3">
                <div className="text-xs font-semibold uppercase tracking-normal text-muted">Evidence contract</div>
                <dl className="mt-2 grid grid-cols-[auto_1fr] gap-x-3 gap-y-1 text-xs">
                  <dt className="text-muted">Event</dt><dd className="break-all font-mono">{selectedPath.event_id ?? selectedPath.id}</dd>
                  <dt className="text-muted">Assessability</dt><dd>{selectedEvent?.assessability ?? 'legacy projection'}</dd>
                  <dt className="text-muted">Negative evidence</dt><dd>{selectedEvent?.negative_evidence_rule ?? 'not encoded'}</dd>
                  <dt className="text-muted">Topology</dt><dd>observed partial path; closure not inferred</dd>
                </dl>
              </div>
            </div>
          </div>
          <AlignmentTraceTable alignments={supportingAlignments} />
        </section>
      )}

      <section className="data-scroll rounded-lg border border-line" aria-labelledby="simple-sv-title">
        <div className="border-b border-line bg-panel2 px-3 py-2">
          <h3 id="simple-sv-title" className="text-sm font-semibold">Canonical structural events</h3>
        </div>
        <div className="grid min-w-[820px] grid-cols-[1.2fr_0.7fr_0.7fr_0.6fr_0.65fr_1fr] gap-2 border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
          <div>Event</div><div>Type</div><div>Coordinates</div><div>Length</div><div>Molecules</div><div>Evidence / orientation</div>
        </div>
        {data.svs.length === 0 ? <EmptyRow text="No canonical structural event" /> : data.svs.map((sv) => {
          const event = eventsById.get(sv.event_id ?? `sv:${sv.id}`);
          return (
            <button
              key={sv.id}
              type="button"
              onClick={() => { setSelectedSvId(sv.id); inspectEvent(sv.event_id ?? `sv:${sv.id}`); }}
              className="grid min-w-[820px] w-full grid-cols-[1.2fr_0.7fr_0.7fr_0.6fr_0.65fr_1fr] gap-2 border-b border-line px-3 py-2 text-left text-sm last:border-0 hover:bg-panel2"
            >
              <span className="truncate font-mono text-xs" title={sv.id}>{sv.id}</span>
              <span>{sv.type}</span>
              <span>{sv.start}–{sv.end}</span>
              <span>{sv.length}</span>
              <span>{event?.supporting_molecule_ids.length ?? sv.supporting_reads.length}</span>
              <span className="truncate text-muted" title={[...(sv.evidence_sources ?? []), ...(sv.orientations ?? [])].join(', ')}>
                {(sv.evidence_sources ?? [sv.evidence_source ?? 'unknown']).join('+')} · {(sv.orientations ?? [sv.orientation ?? 'n/a']).join(', ')}
              </span>
            </button>
          );
        })}
      </section>
    </div>
  );
}

export function TraceList({
  label,
  values,
  empty,
  onSelect
}: {
  label: string;
  values: string[];
  empty: string;
  onSelect: (value: string) => void;
}) {
  return (
    <div className="rounded border border-line bg-shell p-3">
      <div className="text-xs font-semibold uppercase tracking-normal text-muted">{label}</div>
      <div className="mt-2 max-h-36 overflow-y-auto">
        {values.length === 0 ? <span className="text-xs text-muted">{empty}</span> : values.map((value) => (
          <button
            key={value}
            type="button"
            className="block w-full truncate rounded px-1 py-1 text-left font-mono text-xs hover:bg-panel2 hover:text-aqua"
            title={value}
            onClick={() => onSelect(value)}
          >
            {value}
          </button>
        ))}
      </div>
    </div>
  );
}

export function AlignmentTraceTable({ alignments }: { alignments: AlignmentFragment[] }) {
  return (
    <div className="data-scroll border-t border-line">
      <div className="px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">Source alignment fragments</div>
      <div className="grid min-w-[900px] grid-cols-[0.8fr_0.9fr_0.55fr_0.6fr_0.65fr_0.8fr_1.4fr] gap-2 border-t border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
        <div>Alignment</div><div>Molecule</div><div>Role</div><div>Strand</div><div>MAPQ</div><div>Start / CIGAR</div><div>SA provenance</div>
      </div>
      {alignments.length === 0 ? <EmptyRow text="Alignment-fragment provenance is unavailable in this result schema" /> : alignments.map((alignment) => (
        <div key={alignment.id} className="grid min-w-[900px] grid-cols-[0.8fr_0.9fr_0.55fr_0.6fr_0.65fr_0.8fr_1.4fr] gap-2 border-t border-line px-3 py-2 text-xs">
          <div className="truncate font-mono" title={alignment.id}>{alignment.id}</div>
          <div className="truncate" title={alignment.molecule_id}>{alignment.molecule_id}</div>
          <div>{alignment.role}</div>
          <div>{alignment.strand}</div>
          <div>{alignment.mapping_quality}</div>
          <div className="truncate" title={`${alignment.reference_name}:${alignment.reference_start} ${alignment.cigar}`}>
            {alignment.reference_start} · {alignment.cigar}
          </div>
          <div className="truncate font-mono text-muted" title={alignment.aux_tags?.SA ?? 'no SA tag'}>
            {alignment.aux_tags?.SA ?? 'no SA tag retained'}
          </div>
        </div>
      ))}
    </div>
  );
}
