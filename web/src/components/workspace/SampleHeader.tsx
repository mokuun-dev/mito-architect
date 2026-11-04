import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { AlertTriangle, ChevronDown, Info } from 'lucide-react';
import { useSearchParams } from 'react-router-dom';
import ExportButtons from '../ExportButtons';
import { platformLabel, type DemoIndex, type DemoSampleEntry } from '../../lib/demoPanel';

export function SampleHeader({ data, jobId, isDemo, demoIndex, demoEntry }: {
  data: MitoAnalysisData;
  jobId: string;
  isDemo: boolean;
  demoIndex?: DemoIndex;
  demoEntry?: DemoSampleEntry;
}) {
  const numt = data.filter_stats.numt_assessment;
  const metadata = data.metadata;
  return (
    <header className="grid gap-4 pt-6">
      <div className="flex flex-wrap items-start justify-between gap-4">
        <div className="min-w-0">
          <div className="flex flex-wrap items-center gap-2">
            <span className="section-kicker">Analysis result</span>
            {isDemo && <span className="badge-warning">{demoEntry?.simulated === false ? 'Demo · real public reads' : 'Demo · simulated reads'}</span>}
            <span className="badge-neutral font-mono">schema {metadata.schema_version ?? '—'}</span>
          </div>
          <h1 className="mt-1.5 break-words text-2xl font-semibold tracking-tight">{demoEntry?.title ?? metadata.sample}</h1>
          <p className="mt-1 text-[13px] text-muted">
            <span className="font-mono">{metadata.reference_accession ?? 'Custom reference'}</span>
            {' · '}{metadata.reference_length.toLocaleString('en-US')} bp{' · '}engine {metadata.engine_version}
          </p>
        </div>
        <div className="flex flex-wrap items-center gap-2">
          <details className="relative">
            <summary className="btn-secondary list-none">
              <Info className="h-4 w-4" aria-hidden />
              Run details
              <ChevronDown className="h-4 w-4 text-muted" aria-hidden />
            </summary>
            <div className="popover right-0 top-[calc(100%+6px)] w-[min(380px,calc(100vw-32px))] p-4">
              <dl className="grid grid-cols-[auto_minmax(0,1fr)] gap-x-4 gap-y-2 text-[13px]">
                <dt className="text-muted">Engine</dt><dd className="font-mono">{metadata.engine_version}</dd>
                <dt className="text-muted">Schema</dt><dd className="font-mono">{metadata.schema_version ?? '—'}</dd>
                <dt className="text-muted">Input</dt><dd className="break-all font-mono">{metadata.input_path ?? '—'}</dd>
                <dt className="text-muted">Reference</dt><dd className="break-all font-mono">{metadata.reference_path ?? '—'}</dd>
                {metadata.calling_parameters && (
                  <>
                    <dt className="text-muted">Min MAPQ</dt><dd className="font-mono">{metadata.calling_parameters.min_mapping_quality}</dd>
                    <dt className="text-muted">Min base Q</dt><dd className="font-mono">{metadata.calling_parameters.min_base_quality}</dd>
                  </>
                )}
              </dl>
              <p className="mt-3 border-t border-line pt-3 text-xs leading-5 text-muted">
                Validating the result format does not confirm the biological correctness of its conclusions.
              </p>
            </div>
          </details>
          <ExportButtons jobId={jobId} jsonData={data} htmlAvailable={!isDemo} />
        </div>
      </div>
      {isDemo && demoEntry && demoIndex && <DemoNotice entry={demoEntry} index={demoIndex} />}
      {numt && !numt.specificity_assessable && (
        <div className="notice-warning" role="status">
          <AlertTriangle className="mt-0.5 h-4 w-4 shrink-0 text-amber" aria-hidden />
          <span>
            <strong className="font-semibold">NUMT specificity was not assessed.</strong>{' '}
            Without a verified competitive nuclear-plus-mitochondrial alignment, low-fraction variants may come from nuclear copies.
            {numt.provenance_reason && <span className="text-muted"> ({numt.provenance_reason})</span>}
          </span>
        </div>
      )}
    </header>
  );
}

function DemoNotice({ entry, index }: { entry: DemoSampleEntry; index: DemoIndex }) {
  const [searchParams, setSearchParams] = useSearchParams();
  const selectSample = (id: string) => {
    // Selections and filters belong to the previous sample.
    const next = new URLSearchParams(searchParams);
    next.set('sample', id);
    for (const key of ['event', 'sv', 'gene']) next.delete(key);
    setSearchParams(next);
  };
  return (
    <section className="grid gap-3 rounded-lg border border-line bg-panel px-4 py-3" aria-label="Demo sample">
      <div className="flex flex-wrap items-center gap-x-3 gap-y-2">
        <label htmlFor="demo-sample" className="text-[13px] font-medium">Sample</label>
        <select
          id="demo-sample"
          className="h-9 max-w-full rounded-lg border border-line bg-panel2 px-2 text-[13px]"
          value={entry.id}
          onChange={(event) => selectSample(event.target.value)}
        >
          {index.samples.map((sample) => <option key={sample.id} value={sample.id}>{sample.title}</option>)}
        </select>
        <span className="text-[13px] text-muted">{platformLabel(entry.platform)} · {Math.round(entry.median_depth).toLocaleString('en-US')}× median depth</span>
      </div>
      <p className="max-w-[90ch] text-[13px] leading-5 text-muted">{entry.summary}</p>
      <p className="max-w-[90ch] text-xs leading-5 text-muted">
        {entry.simulated
          ? `Background: public GenBank genome${entry.accessions.length > 1 ? 's' : ''} ${entry.accessions.join(', ')}. Reads are simulated, and any clinical variant or large deletion is added in silico, so it is not a finding about that person.`
          : <>Reads from the public run <a className="text-accent hover:underline" href={entry.source} target="_blank" rel="noreferrer">{entry.accessions.join(', ')}</a>. No truth set exists for this sample.</>}
      </p>
    </section>
  );
}
