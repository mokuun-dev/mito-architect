import type { MitoAnalysisData, PrecisionCalls, PrecisionDiagnostics } from '@mito-architect/visualization-lib';
import { useQuery } from '@tanstack/react-query';
import { useMemo, useState } from 'react';
import { useSearchParams } from 'react-router-dom';
import { fetchDemoIndex, resolveDemoSample, type TruthScore } from '../../lib/demoPanel';
import { filterLabel } from '../../lib/precisionProjection';
import { EmptyRow, Metric } from './shared';

const REFERENCE_LENGTH_FALLBACK = 16569;

export function CallQualityPanel({ data }: { data: MitoAnalysisData }) {
  const calls = data.precision_calls;
  if (!calls?.assessed || !calls.diagnostics) {
    return (
      <p className="text-[13px] leading-6 text-muted">
        The statistical caller did not run on this result{calls?.reason ? `: ${calls.reason}` : ''}. It needs aligned reads and an engine with this feature.
      </p>
    );
  }
  const diagnostics = calls.diagnostics;
  const length = data.metadata.reference_length ?? REFERENCE_LENGTH_FALLBACK;
  const callable = (calls.callable_ranges ?? []).reduce((total, [start, end]) => total + end - start + 1, 0);
  return (
    <div className="grid gap-6">
      <p className="max-w-[85ch] text-[13px] leading-6 text-muted">
        Every mismatch and gap starts as a hypothesis. A call has to exceed an error background measured in this sample, on both strands, with the false discovery rate held at 1% over all hypotheses tested. Thresholds are not tuned per sample, so noisier data simply supports fewer low-level calls.
      </p>
      <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-4">
        <Metric label="Molecules used" value={diagnostics.molecules_used.toLocaleString('en-US')} />
        <Metric label="Excluded as divergent" value={diagnostics.molecules_excluded_divergent.toLocaleString('en-US')} />
        <Metric label="Median depth" value={`${Math.round(diagnostics.median_depth).toLocaleString('en-US')}×`} />
        <Metric label="Callable positions" value={`${callable.toLocaleString('en-US')} of ${length.toLocaleString('en-US')}`} />
        <Metric label="SNV hypotheses tested" value={diagnostics.snv_hypotheses.toLocaleString('en-US')} />
        <Metric label="Indel candidates" value={diagnostics.indel_candidates.toLocaleString('en-US')} />
        <Metric label="Substitution rate per base" value={`${(diagnostics.mean_substitution_rate * 100).toFixed(3)}%`} />
        <Metric label="Indel noise, 99th percentile" value={`${(diagnostics.indel_noise_p99 * 100).toFixed(2)}%`} />
      </div>
      {callable < length && <CallableNote calls={calls} length={length} />}
      <TruthCheck />
      <ErrorProfile diagnostics={diagnostics} />
      <SetAside calls={calls} />
    </div>
  );
}

function CallableNote({ calls, length }: { calls: PrecisionCalls; length: number }) {
  const gaps: Array<[number, number]> = [];
  let cursor = 1;
  for (const [start, end] of calls.callable_ranges ?? []) {
    if (start > cursor) gaps.push([cursor, start - 1]);
    cursor = end + 1;
  }
  if (cursor <= length) gaps.push([cursor, length]);
  return (
    <p className="notice-warning" role="status">
      <span>
        <strong className="font-semibold">Not callable: </strong>
        {gaps.slice(0, 6).map(([start, end]) => `${start.toLocaleString('en-US')}–${end.toLocaleString('en-US')}`).join(', ')}
        {gaps.length > 6 ? ` and ${gaps.length - 6} more` : ''}. These positions have too little depth to assess, which says nothing about whether a variant is there.
      </span>
    </p>
  );
}

function TruthCheck() {
  const [searchParams] = useSearchParams();
  const index = useQuery({ queryKey: ['demo-index'], queryFn: ({ signal }) => fetchDemoIndex(signal), staleTime: Infinity });
  const isDemo = window.location.pathname.endsWith('/result/demo');
  const entry = isDemo && index.data ? resolveDemoSample(index.data, searchParams.get('sample')) : undefined;
  const truth = entry?.truth_check;
  if (!truth) return null;
  const rows: Array<[string, TruthScore]> = [
    ['SNVs', truth.snv],
    ['SNVs outside long C-tracts', truth.snv_easy],
    ['Small indels', truth.indel],
    ['Small indels outside long C-tracts', truth.indel_easy],
    ['Large deletions', truth.large_deletion]
  ];
  return (
    <section aria-labelledby="truth-check-title">
      <h3 id="truth-check-title" className="text-sm font-semibold">Check against the known truth</h3>
      <p className="mt-1 max-w-[85ch] text-[13px] leading-5 text-muted">
        These reads were simulated, so the true variants are known. This measures the caller against a simplified error model and does not establish accuracy on real sequencers. In long C-tracts one molecule can be written several equally valid ways, so a correct call may not match the truth's spelling.
      </p>
      <div className="data-scroll mt-3 rounded-lg border border-line">
        <table className="w-full min-w-[620px] text-left text-[13px]">
          <thead className="border-b border-line bg-panel2 text-xs font-semibold text-muted">
            <tr>
              <th scope="col" className="px-3 py-2 font-semibold">Class</th>
              <th scope="col" className="px-3 py-2 text-right font-semibold">Found</th>
              <th scope="col" className="px-3 py-2 text-right font-semibold">Missed</th>
              <th scope="col" className="px-3 py-2 text-right font-semibold">False</th>
              <th scope="col" className="px-3 py-2 text-right font-semibold">Mean level error</th>
            </tr>
          </thead>
          <tbody className="tabular-nums">
            {rows.map(([label, score]) => (
              <tr key={label} className="border-b border-line last:border-0">
                <th scope="row" className="px-3 py-2 font-medium">{label}</th>
                <td className="px-3 py-2 text-right">{score.tp}</td>
                <td className="px-3 py-2 text-right">{score.fn}</td>
                <td className="px-3 py-2 text-right">{score.fp}</td>
                <td className="px-3 py-2 text-right">{score.fraction_mae === null ? '—' : `${(score.fraction_mae * 100).toFixed(1)} pp`}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </section>
  );
}

const SERIES = [
  { key: 'homopolymer_deletion_rate', label: 'Deletion', color: 'rgb(var(--color-gene-coding))' },
  { key: 'homopolymer_insertion_rate', label: 'Insertion', color: 'rgb(var(--color-gene-trna))' }
] as const;

/** Rates at the top of the estimate's range come from a trend, not observation. */
const EXTRAPOLATED_RATE = 0.5;

function ErrorProfile({ diagnostics }: { diagnostics: PrecisionDiagnostics }) {
  const [hover, setHover] = useState<{ run: number; label: string; rate: number; x: number; y: number }>();
  const series = useMemo(
    () => SERIES.map((item) => ({
      ...item,
      points: diagnostics[item.key].filter((point) => point.rate < EXTRAPOLATED_RATE && point.rate > 0)
    })).filter((item) => item.points.length > 0),
    [diagnostics]
  );
  const allPoints = series.flatMap((item) => item.points);
  if (allPoints.length < 2) return null;

  const width = 620;
  const height = 240;
  const margin = { top: 16, right: 72, bottom: 38, left: 54 };
  const runs = allPoints.map((point) => point.run);
  const minRun = Math.min(...runs);
  const maxRun = Math.max(...runs);
  const decadeLow = Math.floor(Math.log10(Math.min(...allPoints.map((point) => point.rate))));
  const decadeHigh = Math.ceil(Math.log10(Math.max(...allPoints.map((point) => point.rate))));
  const x = (run: number) => margin.left + ((run - minRun) / Math.max(1, maxRun - minRun)) * (width - margin.left - margin.right);
  const y = (rate: number) =>
    margin.top + ((decadeHigh - Math.log10(rate)) / Math.max(1, decadeHigh - decadeLow)) * (height - margin.top - margin.bottom);
  const decades: number[] = [];
  for (let decade = decadeLow; decade <= decadeHigh; decade += 1) decades.push(decade);
  const ticks: number[] = [];
  for (let run = minRun; run <= maxRun; run += 1) ticks.push(run);

  return (
    <section aria-labelledby="error-profile-title">
      <h3 id="error-profile-title" className="text-sm font-semibold">Length errors measured in this sample</h3>
      <p className="mt-1 max-w-[85ch] text-[13px] leading-5 text-muted">
        Share of molecules with a one-base length error, by homopolymer length. The caller uses these rates as the background for indels and as gap costs when it re-scores reads.
      </p>
      <div className="mt-2 flex items-center gap-4 text-xs text-muted" aria-hidden>
        {series.map((item) => (
          <span key={item.key} className="inline-flex items-center gap-1.5">
            <span className="h-2.5 w-2.5 rounded-full" style={{ background: item.color }} />
            {item.label}
          </span>
        ))}
      </div>
      <div className="relative mt-1 max-w-[620px]">
        <svg viewBox={`0 0 ${width} ${height}`} className="h-auto w-full" role="img" aria-label="Length error rate rises with homopolymer length; exact values are in the table below.">
          {decades.map((decade) => (
            <g key={decade}>
              <line x1={margin.left} x2={width - margin.right} y1={y(10 ** decade)} y2={y(10 ** decade)} stroke="rgb(var(--color-line))" strokeWidth={1} />
              <text x={margin.left - 8} y={y(10 ** decade) + 4} textAnchor="end" className="fill-muted text-[11px]">
                {decade >= -2 ? `${10 ** (decade + 2)}%` : `${(10 ** (decade + 2)).toFixed(-decade - 2)}%`}
              </text>
            </g>
          ))}
          {ticks.map((run) => (
            <text key={run} x={x(run)} y={height - margin.bottom + 18} textAnchor="middle" className="fill-muted text-[11px]">{run}</text>
          ))}
          <text x={(margin.left + width - margin.right) / 2} y={height - 4} textAnchor="middle" className="fill-muted text-[11px]">Homopolymer length, bases</text>
          {series.map((item) => {
            const last = item.points[item.points.length - 1];
            return (
              <g key={item.key}>
                <polyline
                  fill="none"
                  stroke={item.color}
                  strokeWidth={2}
                  strokeLinejoin="round"
                  points={item.points.map((point) => `${x(point.run)},${y(point.rate)}`).join(' ')}
                />
                {item.points.map((point) => (
                  <g key={point.run}>
                    <circle cx={x(point.run)} cy={y(point.rate)} r={4} fill={item.color} stroke="rgb(var(--color-panel))" strokeWidth={2} />
                    <circle
                      cx={x(point.run)}
                      cy={y(point.rate)}
                      r={14}
                      fill="transparent"
                      onMouseEnter={() => setHover({ run: point.run, label: item.label, rate: point.rate, x: x(point.run), y: y(point.rate) })}
                      onMouseLeave={() => setHover(undefined)}
                    />
                  </g>
                ))}
                <text x={x(last.run) + 8} y={y(last.rate) + 4} className="fill-text text-[11px] font-medium">{item.label}</text>
              </g>
            );
          })}
        </svg>
        {hover && (
          <div
            className="pointer-events-none absolute z-10 rounded-md border border-line bg-panel px-2 py-1 text-xs shadow-pop"
            style={{ left: `${(hover.x / width) * 100}%`, top: `${(hover.y / height) * 100}%`, transform: 'translate(-50%, -130%)' }}
          >
            <span className="font-medium">{hover.label}</span>, run of {hover.run}: {(hover.rate * 100).toFixed(2)}%
          </div>
        )}
      </div>
      <details className="mt-2 text-[13px]">
        <summary className="cursor-pointer font-medium text-accent">Values as a table</summary>
        <div className="data-scroll mt-2 rounded-lg border border-line">
          <table className="w-full min-w-[320px] text-left">
            <thead className="border-b border-line bg-panel2 text-xs text-muted">
              <tr>
                <th scope="col" className="px-3 py-1.5 font-semibold">Run length</th>
                {series.map((item) => <th key={item.key} scope="col" className="px-3 py-1.5 text-right font-semibold">{item.label}</th>)}
              </tr>
            </thead>
            <tbody className="tabular-nums">
              {ticks.map((run) => (
                <tr key={run} className="border-b border-line last:border-0">
                  <th scope="row" className="px-3 py-1.5 font-medium">{run}</th>
                  {series.map((item) => {
                    const point = item.points.find((candidate) => candidate.run === run);
                    return <td key={item.key} className="px-3 py-1.5 text-right">{point ? `${(point.rate * 100).toFixed(2)}%` : '—'}</td>;
                  })}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </details>
    </section>
  );
}

function SetAside({ calls }: { calls: PrecisionCalls }) {
  const rows = useMemo(() => {
    const snvs = (calls.snvs ?? []).filter((call) => call.status !== 'PASS').map((call) => ({
      id: call.id,
      label: `${call.position.toLocaleString('en-US')} ${call.ref}>${call.alt}`,
      depth: `${call.alt_depth.toLocaleString('en-US')} of ${call.depth.toLocaleString('en-US')}`,
      filters: call.filters,
      p: call.p_value
    }));
    const indels = (calls.indels ?? []).filter((call) => call.status !== 'PASS').map((call) => ({
      id: call.id,
      label: call.hgvs,
      depth: `${call.raw_alt_depth.toLocaleString('en-US')} of ${call.depth.toLocaleString('en-US')}`,
      filters: call.filters,
      p: call.p_value
    }));
    const deletions = (calls.deletions ?? []).filter((call) => call.status !== 'PASS').map((call) => ({
      id: call.id,
      label: `del ${call.start.toLocaleString('en-US')}–${call.end.toLocaleString('en-US')}`,
      depth: `${call.junction_molecules.toLocaleString('en-US')} junctions`,
      filters: call.filters,
      p: call.p_value
    }));
    return [...snvs, ...indels, ...deletions].sort((a, b) => a.p - b.p).slice(0, 25);
  }, [calls]);
  return (
    <section aria-labelledby="set-aside-title">
      <h3 id="set-aside-title" className="text-sm font-semibold">Strongest candidates that were set aside</h3>
      <p className="mt-1 max-w-[85ch] text-[13px] leading-5 text-muted">
        A rejected candidate may still be real, but the reads here do not separate it from error. Review these before treating an expected variant as absent.
      </p>
      <div className="data-scroll mt-2 rounded-lg border border-line">
        {rows.length === 0 ? <EmptyRow text="No candidate came close enough to report." /> : (
          <table className="w-full min-w-[620px] text-left text-[13px]">
            <thead className="border-b border-line bg-panel2 text-xs text-muted">
              <tr>
                <th scope="col" className="px-3 py-1.5 font-semibold">Candidate</th>
                <th scope="col" className="px-3 py-1.5 font-semibold">Molecules</th>
                <th scope="col" className="px-3 py-1.5 font-semibold">Why it was set aside</th>
              </tr>
            </thead>
            <tbody>
              {rows.map((row) => (
                <tr key={row.id} className="border-b border-line last:border-0">
                  <th scope="row" className="whitespace-nowrap px-3 py-1.5 font-mono font-medium">{row.label}</th>
                  <td className="whitespace-nowrap px-3 py-1.5 tabular-nums">{row.depth}</td>
                  <td className="px-3 py-1.5 text-muted">{row.filters.map(filterLabel).join('; ')}</td>
                </tr>
              ))}
            </tbody>
          </table>
        )}
      </div>
    </section>
  );
}
