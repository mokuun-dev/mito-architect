import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { AlertTriangle, ChevronDown, Info } from 'lucide-react';
import ExportButtons from '../ExportButtons';

export function SampleHeader({ data, jobId, isDemo }: { data: MitoAnalysisData; jobId: string; isDemo: boolean }) {
  const numt = data.filter_stats.numt_assessment;
  const metadata = data.metadata;
  return (
    <header className="grid gap-4 pt-6">
      <div className="flex flex-wrap items-start justify-between gap-4">
        <div className="min-w-0">
          <div className="flex flex-wrap items-center gap-2">
            <span className="section-kicker">Analysis result</span>
            {isDemo && <span className="badge-warning">Demo data</span>}
            <span className="badge-neutral font-mono">schema {metadata.schema_version ?? '—'}</span>
          </div>
          <h1 className="mt-1.5 break-words text-2xl font-semibold tracking-tight">{metadata.sample}</h1>
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
