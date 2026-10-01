import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { Metric } from './shared';
import { summarizeCoverage } from './summaries';
import { formatPercent } from './shared';
import { coverageDistribution } from '../../lib/coverageModel';

export function CoveragePanel({
  data,
  metrics
}: {
  data: MitoAnalysisData;
  metrics: ReturnType<typeof summarizeCoverage>;
}) {
  const distribution = coverageDistribution(data.coverage);
  const maxDistribution = Math.max(1, ...distribution.map((item) => item.count));
  return (
    <div className="grid gap-4 lg:grid-cols-[320px_minmax(0,1fr)]">
      <div className="grid gap-3 sm:grid-cols-3 lg:grid-cols-1">
        <Metric label="Mean depth" value={`${metrics.meanDepth.toFixed(2)}×`} />
        <Metric label={`${metrics.pctLabel === 'sites' ? 'Sites' : 'Bins'} >20×`} value={formatPercent(metrics.pctGt20 / 100)} />
        <Metric label="Max depth" value={metrics.maxDepth.toString()} />
        <Metric label="Zero / low bins" value={`${metrics.zeroBins} / ${metrics.lowBins}`} />
      </div>
      <div className="grid gap-3 rounded-md border border-line bg-panel2 p-3">
        <div>
          <div className="text-xs font-semibold uppercase tracking-normal text-muted">Depth distribution</div>
          <p className="mt-1 text-xs text-muted">X: depth category; Y: number of reference positions represented by bins.</p>
        </div>
        <div className="grid grid-cols-5 items-end gap-2" aria-label="Coverage depth distribution">
          {distribution.map((item) => (
            <div key={item.label} className="grid gap-1 text-center">
              <div className="flex h-24 items-end rounded-sm bg-shell px-1"><div className="w-full rounded-t-sm bg-sky" style={{ height: `${(item.count / maxDistribution) * 100}%` }} title={`${item.label}: ${item.count} positions`} /></div>
              <div className="text-[11px] text-muted">{item.label}</div>
              <div className="font-mono text-[11px] text-text">{item.count}</div>
            </div>
          ))}
        </div>
        <p className="text-xs leading-5 text-muted">The full positional profile is shown in the linear evidence view. Bins are storage intervals, not individual base calls.</p>
      </div>
    </div>
  );
}
