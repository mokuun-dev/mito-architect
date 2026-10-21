import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { EmptyRow } from './shared';

export function HaplogroupPanel({ data }: { data: MitoAnalysisData }) {
  return (
    <div className="data-scroll rounded-lg border border-line">
      <div className="grid min-w-[760px] grid-cols-[0.7fr_0.8fr_0.55fr_0.8fr_1.4fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
        <div>Cluster</div>
        <div>Haplogroup</div>
        <div>Method score</div>
        <div>Molecules</div>
        <div>Alternatives / molecule signature</div>
      </div>
      {data.clusters.length === 0 ? (
        <EmptyRow text="No clusters" />
      ) : (
        data.clusters.map((cluster) => (
          <div
            key={cluster.id}
            className="grid min-w-[760px] grid-cols-[0.7fr_0.8fr_0.55fr_0.8fr_1.4fr] gap-2 border-b border-line px-3 py-2 text-sm last:border-0"
          >
            <div className="font-semibold">{cluster.label}</div>
            <div className={cluster.haplogroup === 'unassigned' ? 'text-muted' : 'text-aqua'}>
              {cluster.haplogroup ?? 'unassigned'}
              {cluster.haplogroup_assignment?.contamination_warning && (
                <span className="ml-2 text-amber" title="Competing macrohaplogroup assignments">
                  mixed?
                </span>
              )}
            </div>
            <div title="Ranking score emitted by the haplogroup method; it is not a probability or a confidence interval.">
              {cluster.haplogroup_assignment ? cluster.haplogroup_assignment.quality.toFixed(1) : '—'}
            </div>
            <div>{cluster.size}</div>
            <div className="truncate text-muted">
              {cluster.haplogroup_assignment?.candidates.slice(1).map((candidate) =>
                `${candidate.name} ${candidate.score.toFixed(1)}%`
              ).join(', ') || 'no ranked alternatives'}
              {' | '}
              {cluster.sv_signature.length === 0 ? 'no SVs' : cluster.sv_signature.map((sv) => `${sv.sv_id} n=${sv.support}`).join(', ')}
              {(cluster.complex_event_signature?.length ?? 0) > 0
                ? ` | complex paths: ${cluster.complex_event_signature?.map((event) => `n=${event.support}`).join(', ')}`
                : ''}
              {cluster.haplogroup_assignment
                ? ` | markers: ${cluster.haplogroup_assignment.observed_markers?.length ?? 0}, callable bp: ${cluster.haplogroup_assignment.callable_ranges?.reduce((total, range) => total + range.end - range.start + 1, 0) ?? 'unknown'}`
                : ''}
            </div>
          </div>
        ))
      )}
    </div>
  );
}
