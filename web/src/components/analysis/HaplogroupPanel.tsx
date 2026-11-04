import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { EmptyRow } from './shared';

export function HaplogroupPanel({ data }: { data: MitoAnalysisData }) {
  const calls = data.precision_calls;
  const lineages = calls?.haplotypes ?? [];
  return (
    <div className="grid gap-5">
      {lineages.length > 0 && (
        <section aria-labelledby="phased-lineages">
          <h3 id="phased-lineages" className="text-sm font-semibold">Phased lineages</h3>
          <p className="mt-1 max-w-[80ch] text-[13px] leading-5 text-muted">
            Reads are grouped by the alleles they carry at heteroplasmic sites, treating uncovered sites as unknown.
            {lineages.length > 1
              ? ' A second lineage can reflect heteroplasmy, a mixed sample, or contamination; the reads do not say which.'
              : ' Molecules here form a single lineage.'}
          </p>
          <ul className="mt-3 grid gap-2">
            {lineages.map((lineage) => (
              <li key={lineage.id} className="rounded-lg border border-line px-4 py-3">
                <div className="flex flex-wrap items-baseline justify-between gap-x-4 gap-y-1">
                  <span className="text-[15px] font-semibold text-accent">{lineage.haplogroup ?? 'unassigned'}</span>
                  <span className="text-[13px] tabular-nums text-muted">
                    {(lineage.fraction * 100).toFixed(1)}% of informative molecules · {lineage.molecules.toLocaleString('en-US')} placed
                    {lineage.haplogroup_quality !== undefined && ` · score ${lineage.haplogroup_quality.toFixed(1)}`}
                  </span>
                </div>
                <div className="mt-2 h-1.5 overflow-hidden rounded-full bg-panel2" aria-hidden>
                  <div className="h-full rounded-full bg-accent" style={{ width: `${Math.max(1, lineage.fraction * 100)}%` }} />
                </div>
                {lineage.defining_variants.length > 0 && (
                  <p className="mt-2 text-xs text-muted">
                    {lineage.defining_variants.length} variants distinguish it from the largest lineage
                  </p>
                )}
              </li>
            ))}
          </ul>
        </section>
      )}
      <HaplogroupClusters data={data} />
    </div>
  );
}

function HaplogroupClusters({ data }: { data: MitoAnalysisData }) {
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
            <div className={cluster.haplogroup === 'unassigned' ? 'text-muted' : 'text-accent'}>
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
