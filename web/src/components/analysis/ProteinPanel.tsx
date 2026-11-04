import { Search } from 'lucide-react';
import { useEffect, useMemo, useState } from 'react';
import { COMPLEX_STRUCTURES } from '../../lib/proteinStructureSources';
import { useMitoStore } from '../../lib/store';
import NglProteinViewer from '.././NglProteinViewer';
import type { ProteinViewerVariant } from '.././NglProteinViewer';
import type { VariantSummary } from './summaries';
import { Metric } from './shared';
import { EmptyRow } from './shared';
import { formatPercent } from './shared';

/** One protein change: a variant can hold two where genes overlap. */
export interface ProteinRow {
  key: string;
  variant: VariantSummary;
  gene?: string;
  residue?: string;
  consequence?: string;
  grantham?: number;
  structure?: VariantSummary['structure'];
}

const NON_PROTEIN_CONSEQUENCES = new Set(['tRNA_variant', 'rRNA_variant', 'control_region_variant', 'intergenic_variant']);
const SILENT_CONSEQUENCES = new Set(['synonymous_variant', 'start_retained_variant', 'stop_retained_variant']);

export function proteinRows(variants: VariantSummary[]): ProteinRow[] {
  const rows: ProteinRow[] = [];
  for (const variant of variants) {
    if (variant.effects) {
      for (const effect of variant.effects) {
        if (effect.biotype !== 'protein_coding') continue;
        rows.push({
          key: `${variant.key}|${effect.gene ?? ''}`,
          variant,
          gene: effect.gene,
          residue: effect.hgvs_p,
          consequence: effect.consequence,
          grantham: effect.grantham,
          structure: effect.structure
        });
      }
    } else if ((variant.protein || variant.consequence) && !NON_PROTEIN_CONSEQUENCES.has(variant.consequence ?? '')) {
      rows.push({
        key: variant.key,
        variant,
        gene: variant.protein ?? variant.gene,
        residue: variant.residue,
        consequence: variant.consequence,
        structure: variant.structure
      });
    }
  }
  return rows;
}

function toViewerVariant(row: ProteinRow): ProteinViewerVariant {
  return {
    key: row.key,
    label: row.variant.label,
    gene: row.gene,
    protein: row.gene,
    residue: row.residue,
    frequency: row.variant.frequency,
    structure: row.structure
  };
}

function consequenceLabel(consequence: string | undefined): string {
  return consequence ? consequence.replace(/_variant$/, '').replace(/_/g, ' ') : 'consequence not annotated';
}

export function ProteinPanel({ variants }: { variants: VariantSummary[] }) {
  const allRows = useMemo(() => proteinRows(variants), [variants]);
  const [selectedKey, setSelectedKey] = useState<string>();
  const [complex, setComplex] = useState<string>('all');
  const [query, setQuery] = useState('');
  const [includeSilent, setIncludeSilent] = useState(false);
  const selectedEventId = useMitoStore((state) => state.selectedEventId);
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);

  const rows = useMemo(() => {
    const needle = query.trim().toLowerCase();
    return allRows.filter((row) => {
      if (!includeSilent && SILENT_CONSEQUENCES.has(row.consequence ?? '')) return false;
      if (complex !== 'all' && row.structure?.structure_id?.toUpperCase() !== complex) return false;
      if (!needle) return true;
      return [row.gene, row.residue, row.variant.label, row.consequence].some((value) => value?.toLowerCase().includes(needle));
    });
  }, [allRows, complex, includeSilent, query]);

  const selected = rows.find((row) => row.key === selectedKey) ?? rows[0];
  const mapped = allRows.filter((row) => row.structure?.structure_id && row.structure.residue_index !== undefined);
  const viewerVariants = useMemo(() => allRows.filter((row) => !SILENT_CONSEQUENCES.has(row.consequence ?? '') || row.key === selected?.key).map(toViewerVariant), [allRows, selected?.key]);

  useEffect(() => {
    const match = selectedEventId ? allRows.find((row) => row.variant.key === selectedEventId) : undefined;
    if (match && !match.key.startsWith(`${selectedKey?.split('|')[0] ?? ''}|`)) {
      setSelectedKey(match.key);
    }
    // Only an outside selection moves this list; its own clicks already did.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [allRows, selectedEventId]);

  const choose = (key: string) => {
    setSelectedKey(key);
    const row = allRows.find((item) => item.key === key);
    if (row) setSelectedEventId(row.variant.key);
  };

  return (
    <div className="grid min-w-0 gap-4 lg:grid-cols-[300px_minmax(0,1fr)] xl:grid-cols-[340px_minmax(0,1fr)]">
      <div className="grid content-start gap-3">
        <div className="grid grid-cols-2 gap-3">
          <Metric label="Protein changes" value={allRows.length.toString()} />
          <Metric label="On a verified structure" value={mapped.length.toString()} />
        </div>
        <div className="overflow-hidden rounded-lg border border-line bg-panel2">
          <div className="grid gap-2 border-b border-line p-2.5">
            <div className="relative">
              <Search className="pointer-events-none absolute left-2.5 top-1/2 h-4 w-4 -translate-y-1/2 text-subtle" aria-hidden />
              <input
                value={query}
                onChange={(event) => setQuery(event.target.value)}
                placeholder="Gene or residue"
                aria-label="Filter protein changes by gene or residue"
                className="h-9 w-full rounded-lg border border-line bg-panel pl-8 pr-2 text-[13px]"
              />
            </div>
            <div className="grid gap-2">
              <select
                value={complex}
                onChange={(event) => setComplex(event.target.value)}
                aria-label="Filter by respiratory complex"
                className="h-9 w-full min-w-0 rounded-lg border border-line bg-panel px-2 text-[13px]"
              >
                <option value="all">All complexes</option>
                {COMPLEX_STRUCTURES.map((item) => <option key={item.pdbId} value={item.pdbId}>{item.complex}</option>)}
              </select>
              <label className="inline-flex min-h-6 cursor-pointer items-center gap-1.5 text-xs text-muted">
                <input type="checkbox" checked={includeSilent} onChange={(event) => setIncludeSilent(event.target.checked)} />
                Include synonymous changes
              </label>
            </div>
          </div>
          <div className="max-h-[560px] overflow-auto scrollbar-thin">
            {rows.length === 0 ? (
              <EmptyRow text={allRows.length === 0 ? 'No protein-changing variant in this sample' : 'No protein change matches these filters'} />
            ) : (
              rows.map((row) => (
                <button
                  key={row.key}
                  type="button"
                  onClick={() => choose(row.key)}
                  aria-pressed={selected?.key === row.key}
                  className={[
                    'block min-h-11 w-full border-b px-3 py-2 text-left text-[13px] last:border-0',
                    selected?.key === row.key ? 'border-accent bg-accent/10' : 'border-line hover:bg-panel'
                  ].join(' ')}
                >
                  <span className="flex items-baseline justify-between gap-2">
                    <span className="truncate font-semibold">{row.gene} <span className="font-mono font-normal">{row.residue}</span></span>
                    <span className="shrink-0 tabular-nums text-muted">{formatPercent(row.variant.frequency)}</span>
                  </span>
                  <span className="mt-0.5 block truncate text-xs text-muted">
                    {row.variant.label} · {consequenceLabel(row.consequence)}
                    {row.grantham !== undefined ? ` · Grantham ${row.grantham}` : ''}
                    {row.structure ? '' : ' · no model'}
                  </span>
                </button>
              ))
            )}
          </div>
        </div>
      </div>
      <div className="min-w-0 rounded-lg border border-line bg-panel2 p-3 sm:p-4">
        {selected ? (
          <div className="grid gap-4">
            <NglProteinViewer
              variant={toViewerVariant(selected)}
              companions={viewerVariants}
              onSelectVariant={choose}
            />
            <dl className="grid gap-3 sm:grid-cols-2 xl:grid-cols-4">
              <Metric label="Variant" value={selected.variant.label} />
              <Metric label="Residue" value={selected.residue ?? 'not annotated'} />
              <Metric label="Complex" value={selected.structure?.complex ?? 'no verified model'} />
              <Metric label="Frequency" value={formatPercent(selected.variant.frequency)} />
            </dl>
            <MoleculeDistribution variant={selected.variant} />
          </div>
        ) : (
          <div className="text-sm text-muted">
            {allRows.length === 0
              ? 'Protein changes appear here once a coding variant passes calling. Each is shown on an experimental human structure of its respiratory complex.'
              : 'Adjust the filters to select a protein change.'}
          </div>
        )}
      </div>
    </div>
  );
}

export function MoleculeDistribution({ variant }: { variant: VariantSummary }) {
  const clusters = new Map<number, number>();
  for (const molecule of variant.molecules) {
    clusters.set(molecule.clusterId, (clusters.get(molecule.clusterId) ?? 0) + 1);
  }
  const groups = [...clusters.entries()].sort(([left], [right]) => left - right);
  const maximum = Math.max(1, ...groups.map(([, count]) => count));

  return (
    <div className="rounded-md border border-line bg-panel p-3">
      <div className="flex flex-wrap items-center justify-between gap-2">
        <div>
          <div className="text-xs font-semibold uppercase tracking-normal text-muted">Molecule-level support</div>
          <div className="mt-1 text-sm text-text">
            {variant.molecules.length} observed molecules across {groups.length} read cluster{groups.length === 1 ? '' : 's'}
          </div>
        </div>
        <div className="text-xs text-muted">
          {variant.callableDepth == null
            ? 'Observed fraction; locus-callable depth unavailable'
            : `HF ${formatPercent(variant.frequency)} | DP ${variant.callableDepth} | 95% CI ${variant.ci95Low === undefined || variant.ci95High === undefined ? 'not recorded' : `${formatPercent(variant.ci95Low)}-${formatPercent(variant.ci95High)}`}`}
        </div>
      </div>
      {groups.length > 0 && (
        <div className="mt-3 grid gap-2 sm:grid-cols-2 xl:grid-cols-3">
          {groups.map(([clusterId, count]) => (
            <div key={clusterId} className="rounded border border-line bg-panel2 px-2.5 py-2 text-xs">
              <div className="flex items-center justify-between gap-2">
                <span className="font-semibold text-text">{clusterId < 0 ? 'Outliers' : `Cluster ${clusterId + 1}`}</span>
                <span className="text-accent">n={count}</span>
              </div>
              <div className="mt-2 h-1.5 overflow-hidden rounded-full bg-shell">
                <div className="h-full rounded-full bg-accent" style={{ width: `${(count / maximum) * 100}%` }} />
              </div>
            </div>
          ))}
        </div>
      )}
      <details className="mt-3 text-xs text-muted">
        <summary className="cursor-pointer select-none font-semibold text-text">Supporting molecule IDs</summary>
        <div className="mt-2 max-h-28 overflow-auto rounded border border-line bg-shell p-2 font-mono">
          {variant.molecules.map((molecule) => (
            <div key={`${molecule.clusterId}:${molecule.id}`} className="truncate" title={molecule.id}>
              {molecule.id}
            </div>
          ))}
        </div>
      </details>
      {(variant.strandSupport || variant.alleleQuality || variant.readPosition) && (
        <div className="mt-3 grid gap-2 sm:grid-cols-2 xl:grid-cols-4">
          <Metric
            label="Alt strand F/R"
            value={variant.strandSupport ? `${variant.strandSupport.alt_forward}/${variant.strandSupport.alt_reverse}` : '—'}
          />
          <Metric
            label="Strand delta"
            value={variant.strandBiasDelta == null ? 'not estimable' : variant.strandBiasDelta.toFixed(3)}
          />
          <Metric
            label="Alt mean Phred"
            value={variant.alleleQuality?.alternate.mean_phred?.toFixed(1) ?? 'not estimable'}
          />
          <Metric
            label="Read-position delta"
            value={variant.readPosition?.bias_delta == null ? 'not estimable' : variant.readPosition.bias_delta.toFixed(3)}
          />
        </div>
      )}
    </div>
  );
}
