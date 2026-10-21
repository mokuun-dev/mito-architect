import type { VariantSummary } from './summaries';
import { EmptyRow } from './shared';
import { formatPercent } from './shared';

export function ClinicalPanel({ variants }: { variants: VariantSummary[] }) {
  return (
    <div className="data-scroll rounded-lg border border-line">
      <div className="grid min-w-[860px] grid-cols-[1fr_0.55fr_0.65fr_1fr_0.75fr_0.6fr] gap-2 border-b border-line bg-panel2 px-3 py-2 text-xs font-semibold uppercase tracking-normal text-muted">
        <div>Variant</div>
        <div>Support</div>
        <div>Frequency</div>
        <div>Clinical summary</div>
        <div>Conflict</div>
        <div>Records</div>
      </div>
      {variants.length === 0 ? (
        <EmptyRow text="No variants" />
      ) : (
        variants.map((variant) => {
          const annotation = variant.annotation;
          const assertions = annotation?.assertions ?? [];
          const conflicting = annotation?.conflict_status === 'conflicting';
          return (
            <div key={variant.key} className="border-b border-line last:border-0">
              <div className="grid min-w-[860px] grid-cols-[1fr_0.55fr_0.65fr_1fr_0.75fr_0.6fr] gap-2 px-3 py-2 text-sm">
                <div className="font-semibold">{variant.label}</div>
                <div>{variant.support}</div>
                <div>{formatPercent(variant.frequency)}</div>
                <div className={annotation ? 'text-amber' : 'text-muted'}>
                  {annotation?.consensus_significance ?? annotation?.pathogenicity ?? 'not annotated'}
                </div>
                <div className={conflicting ? 'font-semibold text-coral' : annotation ? 'text-text' : 'text-muted'}>
                  {annotation?.conflict_status ?? 'not assessed'}
                </div>
                <div>{assertions.length || '—'}</div>
              </div>
              {annotation && (
                <details className="border-t border-line/70 bg-shell/40 px-3 py-2 text-xs text-muted">
                  <summary className="cursor-pointer font-semibold text-text">
                    Source assertions and provenance
                  </summary>
                  <div className="mt-2 grid gap-2">
                    {assertions.length === 0 ? (
                      <div>Legacy summary only; assertion-level provenance is unavailable.</div>
                    ) : (
                      assertions.map((assertion, index) => (
                        <div
                          key={`${assertion.source}:${assertion.assertion_id}:${index}`}
                          className="rounded border border-line bg-panel2 p-2"
                        >
                          <div className="flex flex-wrap items-center gap-x-2 gap-y-1 text-text">
                            <span className="font-semibold">{assertion.source}</span>
                            <span>{assertion.clinical_significance || 'significance not provided'}</span>
                            <span className="text-muted">({assertion.normalized_significance})</span>
                          </div>
                          <div className="mt-1">
                            {assertion.disease || 'disease not provided'}
                            {assertion.review_status ? ` | review: ${assertion.review_status}` : ''}
                            {assertion.assertion_date ? ` | asserted: ${assertion.assertion_date}` : ''}
                            {assertion.allele_id ? ` | allele ${assertion.allele_id}` : ''}
                          </div>
                          <div className="mt-1 break-words">
                            snapshot {assertion.resource_version || 'not provided'} / {assertion.retrieved_at || 'retrieval unknown'}
                            {assertion.references.length ? ` | ${assertion.references.join(', ')}` : ''}
                            {assertion.source_url && (
                              <>
                                {' | '}
                                <a
                                  className="text-aqua underline underline-offset-2"
                                  href={assertion.source_url}
                                  target="_blank"
                                  rel="noreferrer"
                                >
                                  source record
                                </a>
                              </>
                            )}
                          </div>
                        </div>
                      ))
                    )}
                  </div>
                </details>
              )}
            </div>
          );
        })
      )}
    </div>
  );
}
