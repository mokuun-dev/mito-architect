import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import type { MolecularArchitecture } from '@mito-architect/visualization-lib';
import { useMitoStore } from '../../lib/store';
import { Metric } from './shared';
import { EmptyRow } from './shared';
import { formatPercent } from './shared';

export function ArchitecturePanel({ data }: { data: MitoAnalysisData }) {
  const architectures = data.architectures ?? [];
  const inference = data.architecture_inference;
  const setSelectedMoleculeId = useMitoStore((state) => state.setSelectedMoleculeId);
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);
  const setActiveModule = useMitoStore((state) => state.setActiveModule);
  const inspectEvent = (eventId?: string) => {
    setSelectedMoleculeId(undefined);
    setSelectedEventId(eventId);
    setActiveModule('evidence');
  };
  const inspectMolecule = (moleculeId?: string) => {
    setSelectedEventId(undefined);
    setSelectedMoleculeId(moleculeId);
    setActiveModule('evidence');
  };

  if (data.metadata.schema_version !== '0.6' || !inference) {
    return (
      <div className="rounded-md border border-amber/50 bg-amber/10 p-4 text-sm text-amber">
        Architecture inference is unavailable in this result contract. This does not mean that the
        sample contains no distinct mtDNA molecule types.
      </div>
    );
  }

  return (
    <div className="grid gap-4">
      <div className="grid gap-3 sm:grid-cols-2 xl:grid-cols-5">
        <Metric label="Inference status" value={inference.status} />
        <Metric label="Candidates" value={inference.candidate_count.toString()} />
        <Metric label="Assigned" value={inference.assigned_molecules.toString()} />
        <Metric label="Ambiguous" value={inference.ambiguous_molecules.toString()} />
        <Metric label="Unassigned" value={inference.unassigned_molecules.toString()} />
      </div>

      <div className="rounded-md border border-sky/40 bg-sky/10 p-3 text-sm text-muted">
        These are candidate molecular architectures reconstructed from callable-aware evidence, not
        cell clones or clinical subpopulations. Missing molecule/event pairs remain NOT_CALLABLE.
        Stability is a deterministic full-dataset resampling estimate: each replicate reruns
        inference, then compares recovered signatures and abundance. It is not independent
        biological validation.
      </div>

      {inference.qc_flags.length > 0 && (
        <div className="flex flex-wrap gap-2">
          {inference.qc_flags.map((flag) => (
            <span key={flag} className="rounded border border-amber/50 bg-amber/10 px-2 py-1 text-xs text-amber">
              {flag}
            </span>
          ))}
        </div>
      )}

      {architectures.length === 0 ? (
        <EmptyRow text="No architecture reached the declared support and assignment gates." />
      ) : (
        <div className="grid gap-4">
          {architectures.map((architecture, index) => (
            <ArchitectureCard
              key={architecture.architecture_id}
              architecture={architecture}
              ordinal={index + 1}
              onMoleculeSelect={inspectMolecule}
              onEventSelect={inspectEvent}
            />
          ))}
        </div>
      )}
    </div>
  );
}

export function ArchitectureCard({
  architecture,
  ordinal,
  onMoleculeSelect,
  onEventSelect
}: {
  architecture: MolecularArchitecture;
  ordinal: number;
  onMoleculeSelect: (moleculeId?: string) => void;
  onEventSelect: (eventId?: string) => void;
}) {
  return (
    <article className="overflow-hidden rounded-xl border border-line bg-shell/30">
      <div className="flex flex-wrap items-start justify-between gap-3 border-b border-line bg-panel2 px-4 py-3">
        <div className="min-w-0">
          <div className="section-kicker">Candidate architecture MA-{ordinal.toString().padStart(2, '0')}</div>
          <div className="mt-1 break-all font-mono text-xs text-muted" title={architecture.architecture_id}>
            {architecture.architecture_id}
          </div>
        </div>
        <div className="text-right">
          <div className="font-mono text-xl font-semibold tabular-nums">
            {formatPercent(architecture.estimated_fraction)}
          </div>
          <div className="text-xs text-muted">
            95% CI {formatPercent(architecture.confidence_interval.low)}–{formatPercent(architecture.confidence_interval.high)}
          </div>
        </div>
      </div>

      <div className="grid gap-3 p-4 lg:grid-cols-[minmax(0,1.2fr)_minmax(320px,0.8fr)]">
        <div className="grid gap-3">
          <div>
            <div className="text-xs font-semibold uppercase tracking-normal text-muted">Defining signature</div>
            <div className="mt-2 flex flex-wrap gap-2">
              {architecture.defining_event_signature.length === 0 ? (
                <span className="text-sm text-muted">reference-like callable profile</span>
              ) : architecture.defining_event_signature.map((eventId) => (
                <button
                  key={eventId}
                  type="button"
                  onClick={() => onEventSelect(eventId)}
                  className="min-h-9 max-w-full truncate rounded-lg border border-aqua/40 bg-panel2 px-2 py-1 font-mono text-xs text-aqua hover:border-aqua"
                  title={`Inspect ${eventId} in Molecules & phase`}
                >
                  {eventId}
                </button>
              ))}
            </div>
          </div>

          <div>
            <div className="text-xs font-semibold uppercase tracking-normal text-muted">Optional events</div>
            <div className="mt-2 flex flex-wrap gap-2">
              {architecture.optional_event_ids.length === 0 ? (
                <span className="text-sm text-muted">none above the optional-profile threshold</span>
              ) : architecture.optional_event_ids.map((eventId) => (
                <button
                  key={eventId}
                  type="button"
                  onClick={() => onEventSelect(eventId)}
                  className="min-h-9 max-w-full truncate rounded-lg border border-line bg-panel2 px-2 py-1 font-mono text-xs text-muted hover:border-aqua hover:text-aqua"
                >
                  {eventId}
                </button>
              ))}
            </div>
          </div>

          <details className="rounded-lg border border-line bg-panel2/50 p-3 text-xs">
            <summary className="cursor-pointer font-semibold text-text">
              Callable event profile ({architecture.event_profile.length})
            </summary>
            <div className="data-scroll mt-3 max-h-64">
              <div className="grid min-w-[720px] grid-cols-[minmax(260px,1.5fr)_90px_90px_90px_110px_100px] gap-2 border-b border-line pb-2 font-semibold uppercase text-muted">
                <div>Event</div><div>State</div><div>ALT</div><div>Absent</div><div>Not callable</div><div>ALT fraction</div>
              </div>
              {architecture.event_profile.map((profile) => (
                <div key={profile.event_id} className="grid min-w-[720px] grid-cols-[minmax(260px,1.5fr)_90px_90px_90px_110px_100px] gap-2 border-b border-line/60 py-2 last:border-0">
                  <button type="button" onClick={() => onEventSelect(profile.event_id)} className="truncate text-left font-mono text-aqua" title={profile.event_id}>
                    {profile.event_id}
                  </button>
                  <div>{profile.consensus_state}</div>
                  <div>{profile.alternate}</div>
                  <div>{profile.absent}</div>
                  <div>{profile.not_callable}</div>
                  <div>{formatPercent(profile.alternate_fraction)}</div>
                </div>
              ))}
            </div>
          </details>
        </div>

        <div className="grid content-start gap-3 sm:grid-cols-2 lg:grid-cols-1 xl:grid-cols-2">
          <Metric label="Assigned molecules" value={architecture.molecule_count.toString()} />
          <Metric label="Median callable" value={formatPercent(architecture.median_callable_fraction)} />
          <Metric label="Mean assignment confidence" value={formatPercent(architecture.assignment_confidence.mean)} />
          <Metric label="Signature similarity" value={`${formatPercent(architecture.cluster_stability.value)} · ${architecture.cluster_stability.status}`} />
          <Metric label="Recovery" value={formatPercent(architecture.cluster_stability.recovery_rate)} />
          <Metric label="Abundance SD" value={formatPercent(architecture.cluster_stability.abundance_standard_deviation)} />
          <div className="rounded-md border border-line bg-panel p-3 sm:col-span-2 lg:col-span-1 xl:col-span-2">
            <div className="text-xs font-semibold uppercase tracking-normal text-muted">Representative molecules</div>
            <div className="mt-2 flex flex-wrap gap-2">
              {architecture.representative_molecule_ids.map((moleculeId) => (
                <button
                  key={moleculeId}
                  type="button"
                  onClick={() => onMoleculeSelect(moleculeId)}
                  className="min-h-9 max-w-full truncate rounded-lg border border-line px-2 py-1 font-mono text-xs hover:border-aqua hover:text-aqua"
                  title={`Inspect ${moleculeId} in Molecules & phase`}
                >
                  {moleculeId}
                </button>
              ))}
            </div>
          </div>
          <div className="rounded-md border border-line bg-panel p-3 text-xs text-muted sm:col-span-2 lg:col-span-1 xl:col-span-2">
            {architecture.ambiguous_molecule_ids.length} ambiguous molecules near this candidate · haplogroup:{' '}
            {architecture.haplogroup_evidence.best} ({architecture.haplogroup_evidence.status}, {architecture.haplogroup_evidence.quality.toFixed(1)}%)
          </div>
        </div>
      </div>
    </article>
  );
}
