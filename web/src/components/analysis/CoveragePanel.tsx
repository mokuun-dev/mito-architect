import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { Maximize2 } from 'lucide-react';
import { useState } from 'react';
import { coverageDistribution } from '../../lib/coverageModel';
import type { GenomeRegion } from '../../lib/store';
import { LinearTracks } from '../genome/LinearTracks';
import { formatPercent, Metric } from './shared';
import { summarizeCoverage } from './summaries';

export function CoveragePanel({
  data,
  metrics
}: {
  data: MitoAnalysisData;
  metrics: ReturnType<typeof summarizeCoverage>;
}) {
  const [region, setRegion] = useState<GenomeRegion>();
  const length = data.metadata.reference_length;
  const visible = region ?? { start: 1, end: length };
  const distribution = coverageDistribution(data.coverage);
  const maxDistribution = Math.max(1, ...distribution.map((item) => item.count));
  return (
    <div className="grid gap-5">
      <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-4">
        <Metric label="Mean depth" value={`${metrics.meanDepth.toFixed(2)}×`} />
        <Metric label={`${metrics.pctLabel === 'sites' ? 'Sites' : 'Bins'} >20×`} value={formatPercent(metrics.pctGt20 / 100)} />
        <Metric label="Max depth" value={`${metrics.maxDepth.toLocaleString('en-US')}×`} />
        <Metric label="Zero / low (<20×) bins" value={`${metrics.zeroBins} / ${metrics.lowBins}`} />
      </div>

      <section className="rounded-lg border border-line" aria-labelledby="coverage-profile-title">
        <div className="flex flex-wrap items-center justify-between gap-2 border-b border-line px-4 py-3">
          <div>
            <h3 id="coverage-profile-title" className="text-[13px] font-semibold">Depth profile</h3>
            <p className="text-xs text-muted tabular-nums">
              {visible.start.toLocaleString('en-US')}–{visible.end.toLocaleString('en-US')} bp · drag on the axis to zoom
            </p>
          </div>
          <button type="button" className="btn-secondary min-h-8" disabled={!region} onClick={() => setRegion(undefined)}>
            <Maximize2 className="h-4 w-4" aria-hidden />
            Whole reference
          </button>
        </div>
        <div className="px-4 py-3">
          <LinearTracks
            data={data}
            rows={[]}
            region={visible}
            tracks={['genes', 'coverage']}
            onRegionChange={(next) => setRegion(next && next.end - next.start >= 20 ? next : undefined)}
          />
        </div>
        <p className="border-t border-line px-4 py-2.5 text-xs leading-5 text-muted">
          Bins are storage intervals, not individual base calls. A red mark is a measured zero; positions without reported
          callability are not drawn as zero.
        </p>
      </section>

      <section className="rounded-lg border border-line p-4" aria-labelledby="coverage-distribution-title">
        <h3 id="coverage-distribution-title" className="text-[13px] font-semibold">Depth distribution</h3>
        <p className="mt-0.5 text-xs text-muted">Reference positions in each depth band.</p>
        <div className="mt-4 grid grid-cols-5 items-end gap-3" role="img" aria-label={'Coverage depth distribution: ' + distribution.map((item) => `${item.label} ${item.count} positions`).join(', ')}>
          {distribution.map((item) => (
            <div key={item.label} className="grid gap-1.5 text-center" title={`${item.label}: ${item.count.toLocaleString('en-US')} positions`}>
              <div className="flex h-28 items-end">
                <div className="w-full rounded-t-[4px] bg-coverage" style={{ height: `${Math.max(item.count ? 2 : 0, (item.count / maxDistribution) * 100)}%` }} />
              </div>
              <div className="font-mono text-xs font-medium tabular-nums">{item.count.toLocaleString('en-US')}</div>
              <div className="text-[11px] text-muted">{item.label}</div>
            </div>
          ))}
        </div>
      </section>
    </div>
  );
}
