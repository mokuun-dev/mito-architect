import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import type { EventRow } from '../../lib/eventRows';
import { summarizeCoverage } from '../../lib/coverageModel';
import type { AnalysisModule } from '../../lib/store';

/** Sample-level headline numbers; each tile opens the module that explains it. */
export function KpiStrip({ data, rows, onOpen }: { data: MitoAnalysisData; rows: EventRow[]; onOpen: (module: AnalysisModule) => void }) {
  const coverage = summarizeCoverage(data);
  const molecules = data.filter_stats.evidence_eligible_molecules ?? data.filter_stats.input_molecules;
  const variants = rows.filter((row) => row.kind !== 'sv');
  const measured = variants.filter((row) => row.hf !== undefined);
  const heteroplasmic = measured.filter((row) => row.hf! < 0.95).length;
  const svCount = rows.length - variants.length;
  const haplogroup = [...data.clusters]
    .filter((cluster) => cluster.haplogroup && cluster.haplogroup !== 'unassigned')
    .sort((a, b) => b.size - a.size)[0];

  return (
    <section className="grid grid-cols-2 gap-px overflow-hidden rounded-xl border border-line bg-line shadow-card lg:grid-cols-5" aria-label="Sample summary">
      <Tile
        label={data.filter_stats.evidence_eligible_molecules === undefined ? 'Molecules' : 'Evidence molecules'}
        value={molecules === undefined ? '—' : molecules.toLocaleString('en-US')}
        detail={`${data.filter_stats.numt_filtered_reads.toLocaleString('en-US')} NUMT-filtered reads`}
        onClick={() => onOpen('evidence')}
      />
      <Tile
        label="Mean depth"
        value={data.coverage.length ? `${coverage.meanDepth.toFixed(1)}×` : '—'}
        detail={`${coverage.pctGt20.toFixed(1)}% of ${coverage.pctLabel} >20×`}
        onClick={() => onOpen('coverage')}
      />
      <Tile
        label="Small variants"
        value={variants.length.toLocaleString('en-US')}
        detail={measured.length ? `${heteroplasmic} with HF below 95%` : 'No callable denominators'}
        onClick={() => onOpen('overview')}
      />
      <Tile
        label="Structural events"
        value={svCount.toLocaleString('en-US')}
        detail={`${(data.complex_events?.length ?? 0).toLocaleString('en-US')} complex paths`}
        onClick={() => onOpen('rearrangements')}
      />
      <Tile
        wide
        label="Haplogroup"
        value={haplogroup?.haplogroup ?? '—'}
        detail={haplogroup ? `Largest of ${data.clusters.length} cluster${data.clusters.length === 1 ? '' : 's'}${haplogroup.haplogroup_assignment ? ` · score ${haplogroup.haplogroup_assignment.quality.toFixed(1)}` : ''}` : 'Not assigned'}
        onClick={() => onOpen('haplogroups')}
      />
    </section>
  );
}

function Tile({ label, value, detail, onClick, wide = false }: { label: string; value: string; detail: string; onClick: () => void; wide?: boolean }) {
  return (
    <button
      type="button"
      onClick={onClick}
      className={'group min-w-0 bg-panel px-4 py-3.5 text-left hover:bg-panel2 ' + (wide ? 'col-span-2 lg:col-span-1' : '')}
    >
      <span className="block text-xs font-medium text-muted">{label}</span>
      <span className="mt-1 block truncate text-xl font-semibold tracking-tight">{value}</span>
      <span className="mt-0.5 block truncate text-xs text-muted group-hover:text-text">{detail}</span>
    </button>
  );
}
