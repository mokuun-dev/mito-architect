import type { MitoAnalysisData, ReadFeature } from '@mito-architect/visualization-lib';
import { ArrowDownUp } from 'lucide-react';
import { useMemo, useState } from 'react';
import { FixedSizeList, ListChildComponentProps } from 'react-window';
import { useMitoStore } from '../lib/store';

interface VariantTableProps {
  data: MitoAnalysisData;
}

type SortKey = 'id' | 'cluster_id' | 'mean_quality' | 'length';
const TABLE_GRID =
  'grid-cols-[minmax(260px,1.7fr)_120px_120px_120px_minmax(240px,1.6fr)]';

export default function VariantTable({ data }: VariantTableProps) {
  const [sortKey, setSortKey] = useState<SortKey>('cluster_id');
  const selectedCluster = useMitoStore((state) => state.selectedCluster);
  const minQuality = useMitoStore((state) => state.minQuality);

  const rows = useMemo(() => {
    return [...data.reads]
      .filter((read) => !read.filtered_numt)
      .filter((read) => selectedCluster === undefined || read.cluster_id === selectedCluster)
      .filter((read) => read.mean_quality >= minQuality)
      .sort((a, b) => compareReads(a, b, sortKey));
  }, [data.reads, minQuality, selectedCluster, sortKey]);

  return (
    <section className="overflow-hidden rounded-xl border border-line bg-panel" aria-labelledby="read-variant-title">
      <div className="flex flex-wrap items-center justify-between gap-3 border-b border-line px-4 py-3">
        <div>
          <h2 id="read-variant-title" className="text-sm font-semibold">Reads and variants</h2>
          <p className="mt-1 text-xs text-muted">Virtualized alignment-record projection for rapid inspection.</p>
        </div>
        <div className="flex items-center gap-2 text-xs text-muted">
          <ArrowDownUp className="h-4 w-4" aria-hidden />
          {rows.length} rows
        </div>
      </div>
      <div className="data-scroll">
        <div className="min-w-[900px]" role="table" aria-rowcount={rows.length + 1}>
          <div role="row" className={`grid ${TABLE_GRID} gap-2 border-b border-line bg-panel2 px-4 py-2 text-xs font-semibold uppercase tracking-[0.08em] text-muted`}>
            <Header align="left" label="Read" value="id" sortKey={sortKey} setSortKey={setSortKey} />
            <Header align="center" label="Cluster" value="cluster_id" sortKey={sortKey} setSortKey={setSortKey} />
            <Header align="right" label="Q mean" value="mean_quality" sortKey={sortKey} setSortKey={setSortKey} />
            <Header align="right" label="Length" value="length" sortKey={sortKey} setSortKey={setSortKey} />
            <div role="columnheader" className="text-left">SNPs / SVs</div>
          </div>
          <FixedSizeList height={340} width="100%" itemCount={rows.length} itemSize={44} itemData={rows}>
            {Row}
          </FixedSizeList>
        </div>
      </div>
    </section>
  );
}

function Header({
  label,
  value,
  sortKey,
  setSortKey,
  align
}: {
  label: string;
  value: SortKey;
  sortKey: SortKey;
  setSortKey: (value: SortKey) => void;
  align: 'left' | 'center' | 'right';
}) {
  return (
    <div role="columnheader">
      <button
        type="button"
        onClick={() => setSortKey(value)}
        aria-pressed={sortKey === value}
        className={[
          'min-h-8 w-full rounded-sm',
          align === 'left' ? 'text-left' : align === 'center' ? 'text-center' : 'text-right',
          sortKey === value ? 'text-aqua' : 'text-muted hover:text-text'
        ].join(' ')}
      >
        {label}
      </button>
    </div>
  );
}

function Row({ index, style, data }: ListChildComponentProps<ReadFeature[]>) {
  const read = data[index];
  return (
    <div
      style={style}
      role="row"
      aria-rowindex={index + 2}
      className={`grid ${TABLE_GRID} items-center gap-2 border-b border-line px-4 text-sm hover:bg-panel2/70`}
    >
      <div role="cell" className="truncate font-mono text-xs font-medium">{read.id}</div>
      <div role="cell" className="text-center">{read.cluster_id < 0 ? 'Outlier' : `C${read.cluster_id + 1}`}</div>
      <div role="cell" className="text-right font-mono tabular-nums">{read.mean_quality.toFixed(1)}</div>
      <div role="cell" className="text-right font-mono tabular-nums">{read.length}</div>
      <div role="cell" className="truncate text-muted">
        {read.snps.length} SNPs
        {read.sv_ids.length > 0 ? ` | ${read.sv_ids.join(', ')}` : ''}
      </div>
    </div>
  );
}

function compareReads(a: ReadFeature, b: ReadFeature, sortKey: SortKey): number {
  if (sortKey === 'id') return a.id.localeCompare(b.id);
  return Number(a[sortKey]) - Number(b[sortKey]);
}
