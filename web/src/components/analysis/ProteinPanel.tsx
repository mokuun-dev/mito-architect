import { useEffect } from 'react';
import { useState } from 'react';
import { useMitoStore } from '../../lib/store';
import NglProteinViewer from '.././NglProteinViewer';
import type { ProteinViewerVariant } from '.././NglProteinViewer';
import type { VariantSummary } from './summaries';
import { Metric } from './shared';
import { EmptyRow } from './shared';
import { formatPercent } from './shared';

export function ProteinPanel({ variants }: { variants: VariantSummary[] }) {
  const proteinVariants = variants.filter((variant) => variant.protein || variant.consequence);
  const structurallyMapped = proteinVariants.filter((variant) => variant.structure?.structure_id && variant.structure.residue_index !== undefined);
  const [selectedKey, setSelectedKey] = useState<string>();
  const selectedEventId = useMitoStore((state) => state.selectedEventId);
  const setSelectedEventId = useMitoStore((state) => state.setSelectedEventId);
  const selected = proteinVariants.find((variant) => variant.key === selectedKey) ?? proteinVariants[0];

  useEffect(() => {
    if (selectedEventId && proteinVariants.some((variant) => variant.key === selectedEventId)) {
      setSelectedKey(selectedEventId);
    } else if (!selectedKey && proteinVariants[0]) {
      setSelectedKey(proteinVariants[0].key);
    } else if (selectedKey && !proteinVariants.some((variant) => variant.key === selectedKey)) {
      setSelectedKey(proteinVariants[0]?.key);
    }
  }, [proteinVariants, selectedEventId, selectedKey]);

  return (
    <div className="grid min-w-0 gap-4 lg:grid-cols-[280px_minmax(0,1fr)] xl:grid-cols-[320px_minmax(0,1fr)]">
      <div className="grid content-start gap-3">
        <Metric label="Protein consequences" value={proteinVariants.length.toString()} />
        <Metric label="Verified structure coordinates" value={structurallyMapped.length.toString()} />
        <Metric label="Viewer" value={selected?.structure?.structure_id ?? 'no curated structure'} />
        <Metric label="Example" value={selected ? selected.residue ?? selected.label : 'no protein variant'} />
        <div className="overflow-hidden rounded-lg border border-line bg-panel2">
          <div className="border-b border-line px-3 py-2 text-xs font-semibold uppercase tracking-[0.1em] text-muted">
            Protein variants
          </div>
          <div className="max-h-52 overflow-auto scrollbar-thin">
            {proteinVariants.length === 0 ? (
              <EmptyRow text="No protein-mapped variants" />
            ) : (
              proteinVariants.map((variant) => (
                <button
                  key={variant.key}
                  type="button"
                  onClick={() => {
                    setSelectedKey(variant.key);
                    setSelectedEventId(variant.key);
                  }}
                  className={[
                    'min-h-11 w-full border-b px-3 py-2 text-left text-sm last:border-0',
                    selected?.key === variant.key ? 'border-aqua bg-aqua/10' : 'border-line hover:bg-panel'
                  ].join(' ')}
                >
                  <div className="font-semibold">{variant.protein ?? variant.gene}</div>
                  <div className="mt-1 truncate text-muted">{variant.residue ?? variant.label}</div>
                </button>
              ))
            )}
          </div>
        </div>
      </div>
      <div className="min-w-0 rounded-lg border border-line bg-panel2 p-3 sm:p-4">
        {selected ? (
          <div className="grid gap-4">
            <NglProteinViewer variant={toProteinViewerVariant(selected)} />
            <dl className="grid gap-3 sm:grid-cols-2 xl:grid-cols-4">
              <Metric label="Variant" value={selected.label} />
              <Metric label="Residue" value={selected.residue ?? 'mapped'} />
              <Metric label="Structure" value={selected.structure?.structure_id ?? 'local model'} />
              <Metric label="Frequency" value={formatPercent(selected.frequency)} />
            </dl>
            <MoleculeDistribution variant={selected} />
          </div>
        ) : (
          <div className="text-sm text-muted">
            Non-synonymous consequence calls and local structure mappings are required before residue highlights are populated.
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
                <span className="text-aqua">n={count}</span>
              </div>
              <div className="mt-2 h-1.5 overflow-hidden rounded-full bg-shell">
                <div className="h-full rounded-full bg-aqua" style={{ width: `${(count / maximum) * 100}%` }} />
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

export function toProteinViewerVariant(variant: VariantSummary): ProteinViewerVariant {
  return {
    key: variant.key,
    label: variant.label,
    gene: variant.gene,
    protein: variant.protein,
    residue: variant.residue,
    frequency: variant.frequency,
    structure: variant.structure
  };
}
