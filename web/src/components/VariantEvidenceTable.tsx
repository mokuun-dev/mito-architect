import type { AggregateVariant, MitoAnalysisData } from '@mito-architect/visualization-lib';
import { circularVariants, legacyEventId } from '@mito-architect/visualization-lib';
import { ArrowDownUp, Search } from 'lucide-react';
import { useMemo, useState } from 'react';
import { useMitoStore } from '../lib/store';

interface VariantEvidenceTableProps {
  data: MitoAnalysisData;
}

type SortKey = 'position' | 'heteroplasmy' | 'alt_depth' | 'callable_depth';

export default function VariantEvidenceTable({ data }: VariantEvidenceTableProps) {
  const [sortKey, setSortKey] = useState<SortKey>('position');
  const [query, setQuery] = useState('');
  const selectedEventId = useMitoStore((state) => state.selectedEventId);
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);
  const setSelectedMoleculeId = useMitoStore((state) => state.setSelectedMoleculeId);
  const setSelectedPhaseId = useMitoStore((state) => state.setSelectedPhaseId);
  const variants = data.variants ?? [];

  const rows = useMemo(() => {
    const normalizedQuery = query.trim().toLowerCase();
    return [...variants]
      .filter((variant) => {
        if (!normalizedQuery) return true;
        return [
          variant.event_id,
          variant.type,
          variant.gene,
          variant.consequence,
          variant.ref,
          variant.alt,
          ...(variant.qc_flags ?? [])
        ]
          .filter(Boolean)
          .some((value) => String(value).toLowerCase().includes(normalizedQuery));
      })
      .sort((a, b) => compareVariants(a, b, sortKey));
  }, [query, sortKey, variants]);

  const selected = variants.find((variant) => (variant.event_id ?? legacyEventId(variant)) === selectedEventId);

  if (!data.variants) return <LegacyVariantTable data={data} selectedEventId={selectedEventId} onSelect={setSelectedEventId} />;

  return (
    <section className="overflow-hidden rounded-xl border border-line bg-panel shadow-tool" aria-labelledby="variant-evidence-title">
      <div className="flex flex-col gap-3 border-b border-line px-4 py-3 md:flex-row md:items-center md:justify-between">
        <div>
          <h2 id="variant-evidence-title" className="text-sm font-semibold">
            Unified variant evidence
          </h2>
          <p className="mt-1 text-xs text-muted">
            SNV and small-indel counts are derived from the same molecule/event observations.
          </p>
        </div>
        <label className="flex min-h-11 min-w-0 items-center gap-2 rounded-lg border border-line bg-panel2 px-3 py-2 text-sm focus-within:border-aqua focus-within:shadow-focus md:w-80">
          <Search className="h-4 w-4 shrink-0 text-muted" aria-hidden />
          <span className="sr-only">Filter variants</span>
          <input
            value={query}
            onChange={(event) => setQuery(event.target.value)}
            placeholder="Event, allele, gene, QC flag"
            className="min-w-0 flex-1 bg-transparent outline-none placeholder:text-muted"
          />
        </label>
      </div>

      <div className="data-scroll">
        <table className="w-full min-w-[1040px] border-collapse text-sm">
          <thead className="bg-panel2 text-xs uppercase tracking-normal text-muted">
            <tr>
              <SortableHeader label="Event" sortKey={sortKey} />
              <SortableHeader label="Position" value="position" sortKey={sortKey} onSort={setSortKey} />
              <th className="px-3 py-2 text-left font-semibold">Alleles</th>
              <SortableHeader label="HF" value="heteroplasmy" sortKey={sortKey} onSort={setSortKey} />
              <SortableHeader label="ALT" value="alt_depth" sortKey={sortKey} onSort={setSortKey} />
              <SortableHeader label="Callable" value="callable_depth" sortKey={sortKey} onSort={setSortKey} />
              <th className="px-3 py-2 text-left font-semibold">95% CI</th>
              <th className="px-3 py-2 text-left font-semibold">NUMT</th>
              <th className="px-3 py-2 text-left font-semibold">QC facts</th>
            </tr>
          </thead>
          <tbody>
            {rows.length === 0 ? (
              <tr>
                <td colSpan={9} className="px-4 py-6 text-center text-muted">
                  No variants match the current filter.
                </td>
              </tr>
            ) : (
              rows.map((variant) => {
                const selectedRow = (variant.event_id ?? legacyEventId(variant)) === selectedEventId;
                return (
                  <tr
                    key={variant.event_id ?? `${variant.position}:${variant.ref}:${variant.alt}`}
                    className={selectedRow ? 'border-b border-line bg-aqua/10' : 'border-b border-line hover:bg-panel2/70'}
                  >
                    <td className="max-w-64 px-3 py-2">
                      <button
                        type="button"
                        onClick={() => {
                          setSelectedMoleculeId(undefined);
                          setSelectedPhaseId(undefined);
                          setSelectedEventId(variant.event_id ?? legacyEventId(variant));
                        }}
                        aria-pressed={selectedRow}
                        className="block w-full truncate text-left font-semibold text-aqua hover:underline focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-aqua"
                        title={variant.event_id ?? 'Legacy SNP projection'}
                      >
                        {variant.type ?? 'SNV'} · {variant.event_id ?? 'legacy event'}
                      </button>
                    </td>
                    <td className="px-3 py-2 tabular-nums">{variant.position}</td>
                    <td className="px-3 py-2 font-mono text-xs">{variant.ref} → {variant.alt}</td>
                    <td className="px-3 py-2 tabular-nums">{formatPercent(variant.heteroplasmy)}</td>
                    <td className="px-3 py-2 tabular-nums">{variant.alt_depth}</td>
                    <td className="px-3 py-2 tabular-nums">{variant.callable_depth}</td>
                    <td className="px-3 py-2 tabular-nums text-muted">
                      {formatPercent(variant.ci95_low)}–{formatPercent(variant.ci95_high)}
                    </td>
                    <td className={variant.numt_assessability === 'ASSESSABLE' ? 'px-3 py-2 text-aqua' : 'px-3 py-2 text-amber'}>
                      {variant.numt_assessability ?? 'legacy'}
                    </td>
                    <td className="px-3 py-2">
                      <div className="flex max-w-80 flex-wrap gap-1">
                        {(variant.qc_flags ?? []).length === 0 ? (
                          <span className="text-muted">none recorded</span>
                        ) : (
                          variant.qc_flags?.map((flag) => (
                            <span key={flag} className="rounded border border-line bg-panel2 px-1.5 py-0.5 text-[11px] text-muted">
                              {flag}
                            </span>
                          ))
                        )}
                      </div>
                    </td>
                  </tr>
                );
              })
            )}
          </tbody>
        </table>
      </div>

      <div className="border-t border-line px-4 py-3 text-xs text-muted" aria-live="polite">
        {rows.length} of {variants.length} variants shown. QC facts are observational; calibrated hard filters have not been applied.
      </div>

      {selected && (
        <div className="grid gap-4 border-t border-line bg-panel2/40 p-4 lg:grid-cols-[minmax(0,1fr)_minmax(280px,0.55fr)]">
          <div className="grid gap-3 sm:grid-cols-2 xl:grid-cols-4">
            <EvidenceMetric label="Reference molecules" value={selected.ref_depth.toString()} />
            <EvidenceMetric label="Other allele" value={selected.other_depth.toString()} />
            <EvidenceMetric label="Low quality" value={(selected.low_quality_depth ?? 0).toString()} />
            <EvidenceMetric label="Conflicts" value={(selected.conflict_depth ?? 0).toString()} />
            <EvidenceMetric label="ALT strand F/R" value={selected.strand_support ? `${selected.strand_support.alt_forward}/${selected.strand_support.alt_reverse}` : 'not recorded'} />
            <EvidenceMetric label="ALT mean MAPQ" value={selected.mapping_quality?.alternate.mean?.toFixed(1) ?? 'not estimable'} />
            <EvidenceMetric label="Homopolymer run" value={selected.homopolymer_context?.run_length.toString() ?? 'not recorded'} />
            <EvidenceMetric label="VCF" value={selected.vcf_representable === false ? 'not lossless' : 'representable'} />
          </div>
          <div className="min-w-0 rounded-md border border-line bg-shell p-3">
            <div className="text-xs font-semibold uppercase tracking-normal text-muted">Supporting molecules</div>
            <div className="mt-2 max-h-44 overflow-y-auto font-mono text-xs">
              {(selected.supporting_molecule_ids ?? selected.supporting_reads).length === 0 ? (
                <div className="text-muted">No alternate-supporting molecule IDs.</div>
              ) : (
                (selected.supporting_molecule_ids ?? selected.supporting_reads).map((moleculeId) => (
                  <button
                    key={moleculeId}
                    type="button"
                    onClick={() => {
                      setSelectedEventId(undefined);
                      setSelectedPhaseId(undefined);
                      setSelectedMoleculeId(moleculeId);
                    }}
                    className="block w-full truncate rounded px-1 py-1 text-left hover:bg-panel hover:text-aqua focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-aqua"
                    title={`Inspect molecule ${moleculeId}`}
                  >
                    {moleculeId}
                  </button>
                ))
              )}
            </div>
          </div>
        </div>
      )}
    </section>
  );
}

function LegacyVariantTable({ data, selectedEventId, onSelect }: { data: MitoAnalysisData; selectedEventId?: string; onSelect: (id?: string) => void }) {
  const [query, setQuery] = useState('');
  const [page, setPage] = useState(0);
  const items = useMemo(() => circularVariants(data).filter((item) => !query || `${item.id} ${item.position}`.toLowerCase().includes(query.toLowerCase())).sort((a, b) => a.position - b.position), [data, query]);
  const pageSize = 100;
  const pageCount = Math.max(1, Math.ceil(items.length / pageSize));
  const safePage = Math.min(page, pageCount - 1);
  const visible = items.slice(safePage * pageSize, (safePage + 1) * pageSize);
  return <section className="overflow-hidden rounded-lg border border-line bg-panel" aria-label="Наблюдаемые варианты">
    <div className="flex flex-wrap items-center justify-between gap-3 border-b border-line p-4"><div><h3 className="text-sm font-semibold">Наблюдаемые варианты</h3><p className="mt-1 text-xs text-muted">Схема {data.metadata.schema_version ?? 'не указана'}: пригодный знаменатель для события не представлен. Число поддерживающих записей не является частотой.</p></div><label className="text-xs text-muted">Поиск <input value={query} onChange={(event) => { setQuery(event.target.value); setPage(0); }} className="ml-2 rounded border border-line bg-panel px-2 py-1 text-sm text-text" /></label></div>
    <div className="data-scroll"><table className="w-full min-w-[500px] border-collapse text-sm"><thead className="bg-panel2 text-left text-xs text-muted"><tr><th className="px-3 py-2">Координата</th><th className="px-3 py-2">Замена</th><th className="px-3 py-2">Тип</th><th className="px-3 py-2">Поддержка</th></tr></thead><tbody>{visible.map((item) => <tr key={item.id} className={selectedEventId === item.id ? 'border-t border-line bg-aqua/10' : 'border-t border-line hover:bg-panel2'}><td className="px-3 py-2"><button type="button" onClick={() => onSelect(item.id)} className="font-mono font-semibold text-aqua hover:underline">{item.position}</button></td><td className="px-3 py-2 font-mono">{item.ref} → {item.alt}</td><td className="px-3 py-2">{item.type}</td><td className="px-3 py-2 tabular-nums">{item.support} записей</td></tr>)}</tbody></table></div>
    <div className="flex flex-wrap items-center justify-between gap-2 border-t border-line px-4 py-2 text-xs text-muted"><span>Показано {items.length ? safePage * pageSize + 1 : 0}–{Math.min((safePage + 1) * pageSize, items.length)} из {items.length} найденных; всего {circularVariants(data).length}.</span><div className="flex items-center gap-2"><button type="button" disabled={safePage === 0} onClick={() => setPage(safePage - 1)} className="rounded border border-line px-2 py-1 text-text disabled:opacity-40">Назад</button><span aria-live="polite">{safePage + 1}/{pageCount}</span><button type="button" disabled={safePage >= pageCount - 1} onClick={() => setPage(safePage + 1)} className="rounded border border-line px-2 py-1 text-text disabled:opacity-40">Далее</button></div></div>
  </section>;
}

function SortableHeader({
  label,
  value,
  sortKey,
  onSort
}: {
  label: string;
  value?: SortKey;
  sortKey: SortKey;
  onSort?: (value: SortKey) => void;
}) {
  if (!value || !onSort) {
    return <th className="px-3 py-2 text-left font-semibold">{label}</th>;
  }
  return (
    <th className="px-3 py-2 text-left font-semibold">
      <button
        type="button"
        onClick={() => onSort(value)}
        aria-pressed={sortKey === value}
        className={sortKey === value ? 'inline-flex items-center gap-1 text-aqua' : 'inline-flex items-center gap-1 hover:text-aqua'}
      >
        {label}
        <ArrowDownUp className="h-3 w-3" aria-hidden />
      </button>
    </th>
  );
}

function EvidenceMetric({ label, value }: { label: string; value: string }) {
  return (
    <div className="rounded-md border border-line bg-panel p-3">
      <div className="text-xs uppercase tracking-normal text-muted">{label}</div>
      <div className="mt-1 break-words font-semibold tabular-nums">{value}</div>
    </div>
  );
}

function compareVariants(a: AggregateVariant, b: AggregateVariant, key: SortKey): number {
  const difference = Number(a[key]) - Number(b[key]);
  return difference || (a.event_id ?? '').localeCompare(b.event_id ?? '');
}

function formatPercent(value: number): string {
  return `${(value * 100).toFixed(2)}%`;
}
