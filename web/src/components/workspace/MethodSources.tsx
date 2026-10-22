import type { MitoAnalysisData } from '@mito-architect/visualization-lib';

export function MethodSources({ data }: { data: MitoAnalysisData }) {
  const parameters = data.metadata.calling_parameters;
  return (
    <div className="grid gap-6">
      <p className="max-w-[72ch] text-[13px] leading-6 text-muted">
        Versions and files that make this analysis reproducible. A cited resource does not independently confirm every conclusion drawn from it.
      </p>
      <dl className="grid gap-x-6 gap-y-2 text-[13px] sm:grid-cols-[180px_minmax(0,1fr)]">
        <dt className="text-muted">Engine</dt><dd className="font-mono">{data.metadata.engine_version}</dd>
        <dt className="text-muted">Result schema</dt><dd className="font-mono">{data.metadata.schema_version ?? '—'}</dd>
        <dt className="text-muted">Reference</dt><dd className="break-all font-mono">{data.metadata.reference_accession ?? data.metadata.reference_path ?? '—'}</dd>
        <dt className="text-muted">Input</dt><dd className="break-all font-mono">{data.metadata.input_path ?? '—'}</dd>
        {parameters && (
          <>
            <dt className="text-muted">Mapping quality ≥</dt><dd className="font-mono">{parameters.min_mapping_quality}</dd>
            <dt className="text-muted">Base quality ≥</dt><dd className="font-mono">{parameters.min_base_quality}</dd>
            {parameters.numt_threshold !== undefined && <><dt className="text-muted">NUMT threshold</dt><dd className="font-mono">{parameters.numt_threshold}</dd></>}
          </>
        )}
      </dl>
      {(data.metadata.algorithm_notes ?? []).length > 0 && (
        <section>
          <h3 className="section-kicker">Algorithm notes</h3>
          <ul className="mt-2 list-disc space-y-1 pl-5 text-[13px] text-muted">
            {data.metadata.algorithm_notes?.map((note) => <li key={note}>{note}</li>)}
          </ul>
        </section>
      )}
      <section>
        <h3 className="section-kicker">Versioned resources</h3>
        {!data.metadata.resources?.length ? <p className="mt-2 text-[13px] text-muted">No resource information was reported.</p> : (
          <ul className="mt-2 divide-y divide-line rounded-lg border border-line">
            {data.metadata.resources.map((resource) => (
              <li key={resource.name + ':' + resource.version} className="grid gap-1 px-4 py-3 text-[13px] sm:grid-cols-[minmax(0,1fr)_auto] sm:items-start">
                <div className="min-w-0">
                  <div className="font-semibold">{resource.name} <span className="font-mono font-normal text-muted">{resource.version}</span></div>
                  <div className="mt-1 break-all font-mono text-xs text-muted">SHA-256 {resource.sha256}</div>
                  <div className="mt-1 text-xs text-muted">Retrieved {resource.retrieved || '—'} · License {resource.license || '—'}</div>
                </div>
                {/^https:\/\//.test(resource.source) && (
                  <a href={resource.source} target="_blank" rel="noopener noreferrer" className="font-medium text-accent underline-offset-2 hover:underline">Source ↗</a>
                )}
              </li>
            ))}
          </ul>
        )}
      </section>
    </div>
  );
}
