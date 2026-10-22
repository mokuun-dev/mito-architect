import type { GeneAnnotation, MitoAnalysisData } from '@mito-architect/visualization-lib';
import { ArrowRight, ChevronDown, ChevronUp, X } from 'lucide-react';
import { useEffect, useMemo, type ReactNode } from 'react';
import type { EventRow } from '../../lib/eventRows';
import { formatCount, formatFraction } from '../../lib/format';
import { useMitoStore, type AnalysisModule } from '../../lib/store';
import { TypeBadge } from '../events/EventsTable';
import { LinearTracks } from '../genome/LinearTracks';
import { geneColor } from '../genome/trackStyle';

interface DetailDrawerProps {
  data: MitoAnalysisData;
  rows: EventRow[];
  /** Rows in table order, for previous/next stepping. */
  sequence: EventRow[];
  onSelectRow: (row: EventRow) => void;
  onSelectGene: (gene: GeneAnnotation) => void;
  onShowInGenome: (start: number, end: number, label?: string) => void;
}

/** Non-modal slide-over with the evidence behind the selected event or gene. */
export function DetailDrawer({ data, rows, sequence, onSelectRow, onSelectGene, onShowInGenome }: DetailDrawerProps) {
  const eventId = useMitoStore((state) => state.selectedEventId);
  const svId = useMitoStore((state) => state.selectedSvId);
  const geneName = useMitoStore((state) => state.selectedGene);
  const clear = useMitoStore((state) => state.clearInspection);
  const row = rows.find((item) => item.id === (eventId ?? svId));
  const gene = geneName ? data.genes.find((item) => item.name === geneName) : undefined;
  const open = Boolean(row || gene);
  const index = row ? sequence.findIndex((item) => item.id === row.id) : -1;

  useEffect(() => {
    if (!open) return;
    const onKey = (event: KeyboardEvent) => {
      if (event.key === 'Escape') clear();
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [clear, open]);

  return (
    <aside
      aria-labelledby={open ? 'detail-drawer-title' : undefined}
      aria-label={open ? undefined : 'Details'}
      className={[
        'fixed bottom-0 right-0 top-14 z-20 flex w-full max-w-[440px] flex-col border-l border-line bg-panel transition-transform duration-200',
        open ? 'translate-x-0 shadow-pop' : 'invisible translate-x-full'
      ].join(' ')}
    >
      {open && (
        <>
          <header className="flex items-start justify-between gap-3 border-b border-line px-5 py-4">
            <div className="min-w-0">
              <div className="flex flex-wrap items-center gap-2">
                {row ? <TypeBadge row={row} /> : <span className="badge-neutral">Gene</span>}
                {row && row.genes.length > 0 && <span className="truncate text-[13px] text-muted">{row.genes.slice(0, 3).join(', ')}{row.genes.length > 3 ? ` +${row.genes.length - 3}` : ''}</span>}
              </div>
              <h2 id="detail-drawer-title" className="mt-1.5 break-all font-mono text-lg font-semibold tracking-tight">
                {row ? (row.kind === 'sv' ? row.change : `${row.position.toLocaleString('en-US')} ${row.change}`) : gene!.name}
              </h2>
            </div>
            <div className="flex shrink-0 items-center gap-1">
              {row && (
                <>
                  <button type="button" className="btn-icon h-8 w-8" aria-label="Previous event" title="Previous event" disabled={index <= 0} onClick={() => onSelectRow(sequence[index - 1])}><ChevronUp className="h-4 w-4" aria-hidden /></button>
                  <button type="button" className="btn-icon h-8 w-8" aria-label="Next event" title="Next event" disabled={index < 0 || index >= sequence.length - 1} onClick={() => onSelectRow(sequence[index + 1])}><ChevronDown className="h-4 w-4" aria-hidden /></button>
                </>
              )}
              <button type="button" className="btn-icon h-8 w-8" aria-label="Close details" title="Close (Esc)" onClick={clear}><X className="h-4 w-4" aria-hidden /></button>
            </div>
          </header>
          <div className="flex-1 overflow-y-auto scrollbar-thin">
            {row && <Locus data={data} rows={rows} row={row} onSelectRow={onSelectRow} />}
            {row?.kind === 'sv' ? <StructuralDetails row={row} onShowInGenome={onShowInGenome} />
              : row ? <VariantDetails data={data} row={row} onShowInGenome={onShowInGenome} />
              : gene ? <GeneDetails gene={gene} rows={rows} onSelectRow={onSelectRow} onSelectGene={onSelectGene} onShowInGenome={onShowInGenome} />
              : null}
          </div>
        </>
      )}
    </aside>
  );
}

function Locus({ data, rows, row, onSelectRow }: { data: MitoAnalysisData; rows: EventRow[]; row: EventRow; onSelectRow: (row: EventRow) => void }) {
  const length = data.metadata.reference_length;
  const region = useMemo(() => {
    if (row.kind === 'sv' && row.end >= row.position) {
      const pad = Math.max(50, Math.round((row.end - row.position) * 0.08));
      return { start: Math.max(1, row.position - pad), end: Math.min(length, row.end + pad) };
    }
    if (row.kind === 'sv') return { start: 1, end: length };
    return { start: Math.max(1, row.position - 120), end: Math.min(length, row.position + 120) };
  }, [length, row]);
  const byId = useMemo(() => new Map(rows.map((item) => [item.id, item])), [rows]);
  return (
    <section className="border-b border-line px-5 py-4" aria-label="Locus">
      <LinearTracks
        data={data}
        rows={rows}
        region={region}
        tracks={['genes', 'coverage', row.kind === 'sv' ? 'svs' : 'variants']}
        selection={{ eventId: row.kind === 'sv' ? undefined : row.id, svId: row.kind === 'sv' ? row.id : undefined }}
        onSelectEvent={(id) => { const next = byId.get(id); if (next) onSelectRow(next); }}
        onSelectSv={(id) => { const next = byId.get(id); if (next) onSelectRow(next); }}
      />
    </section>
  );
}

function VariantDetails({ data, row, onShowInGenome }: { data: MitoAnalysisData; row: EventRow; onShowInGenome: DetailDrawerProps['onShowInGenome'] }) {
  const variant = row.variant;
  const openModule = useOpenModule();
  const setSelectedMoleculeId = useMitoStore((state) => state.setSelectedMoleculeId);
  const molecules = variant ? variant.supporting_molecule_ids ?? variant.supporting_reads : [];
  const schema06 = data.metadata.schema_version === '0.6';
  return (
    <>
      <Section title="Frequency">
        {row.hf !== undefined ? (
          <div>
            <div className="flex items-baseline gap-3">
              <span className="text-3xl font-semibold tracking-tight">{formatFraction(row.hf)}</span>
              {row.ciLow !== undefined && row.ciHigh !== undefined && (
                <span className="text-[13px] text-muted">95% CI {formatFraction(row.ciLow)}–{formatFraction(row.ciHigh)}</span>
              )}
            </div>
            <p className="mt-1 text-[13px] text-muted">
              {formatCount(row.altSupport)} ALT of {formatCount(row.callable)} callable molecules
            </p>
          </div>
        ) : (
          <p className="notice-warning">
            No callable denominator was reported for this position. {formatCount(row.altSupport)} supporting records are not a frequency.
          </p>
        )}
      </Section>

      {variant && (
        <Section title="Evidence">
          <dl className="grid grid-cols-2 gap-x-4 gap-y-2.5 text-[13px]">
            <Fact label="REF molecules" value={formatCount(variant.ref_depth)} />
            <Fact label="Other allele" value={formatCount(variant.other_depth)} />
            <Fact label="Low quality" value={formatCount(variant.low_quality_depth)} />
            <Fact label="Conflicting" value={formatCount(variant.conflict_depth)} />
            <Fact label="ALT strand F / R" value={variant.strand_support ? `${variant.strand_support.alt_forward} / ${variant.strand_support.alt_reverse}` : '—'} />
            <Fact label="Strand bias Δ" value={variant.strand_bias_delta === undefined || variant.strand_bias_delta === null ? '—' : variant.strand_bias_delta.toFixed(2)} />
            <Fact label="ALT mean MAPQ" value={variant.mapping_quality?.alternate.mean?.toFixed(1) ?? '—'} />
            <Fact label="ALT mean base Q" value={variant.allele_quality?.alternate.mean_phred?.toFixed(1) ?? '—'} />
            <Fact label="Read-position bias Δ" value={variant.read_position?.bias_delta?.toFixed(2) ?? '—'} />
            <Fact label="Homopolymer run" value={variant.homopolymer_context ? String(variant.homopolymer_context.run_length) : '—'} />
            <Fact label="NUMT" value={variant.numt_assessability === 'ASSESSABLE' ? 'Assessable' : variant.numt_assessability === 'NOT_ASSESSABLE' ? 'Not assessable' : variant.numt_assessability ?? '—'} />
            <Fact label="VCF" value={variant.vcf_representable === false ? 'Not lossless' : 'Representable'} />
          </dl>
          {row.qcFlags.length > 0 && (
            <div className="mt-3 flex flex-wrap gap-1.5" aria-label="QC flags">
              {row.qcFlags.map((flag) => <span key={flag} className="badge-warning font-mono">{flag}</span>)}
            </div>
          )}
          <p className="mt-3 text-xs leading-5 text-muted">QC facts are observational; calibrated hard filters have not been applied.</p>
        </Section>
      )}

      {(row.consequence || variant?.protein || variant?.annotation) && (
        <Section title="Annotation">
          {(row.consequence || variant?.protein) && (
            <dl className="mb-2.5 grid gap-2.5 text-[13px]">
              {row.consequence && <Fact label="Consequence" value={row.consequence} />}
              {variant?.protein && <Fact label="Protein change" value={`${variant.protein}${variant.residue ? ` ${variant.residue}` : ''}`} />}
            </dl>
          )}
          {variant?.annotation && <ClinicalSummary annotation={variant.annotation} />}
        </Section>
      )}

      {molecules.length > 0 && (
        <Section title={`Supporting molecules (${molecules.length.toLocaleString('en-US')})`}>
          <ul className="max-h-48 overflow-y-auto rounded-lg border border-line font-mono text-xs scrollbar-thin">
            {molecules.slice(0, 200).map((id) => (
              <li key={id} className="border-b border-line last:border-0">
                <button
                  type="button"
                  className="block w-full truncate px-3 py-1.5 text-left hover:bg-panel2 hover:text-accent disabled:hover:bg-transparent disabled:hover:text-text"
                  disabled={!schema06}
                  title={schema06 ? `Inspect molecule ${id}` : 'Molecule records are available in schema 0.6 results'}
                  onClick={() => { setSelectedMoleculeId(id); openModule('evidence'); }}
                >
                  {id}
                </button>
              </li>
            ))}
          </ul>
          {molecules.length > 200 && <p className="mt-2 text-xs text-muted">The first 200 are listed.</p>}
        </Section>
      )}

      <Actions>
        <Action onClick={() => onShowInGenome(row.position - 150, row.position + 150, `${row.position} ${row.change}`)}>Show in genome view</Action>
        {schema06 && <Action onClick={() => openModule('evidence')}>Open molecules & phasing</Action>}
        {variant?.annotation && <Action onClick={() => openModule('clinical')}>Open annotations</Action>}
      </Actions>
    </>
  );
}

function StructuralDetails({ row, onShowInGenome }: { row: EventRow; onShowInGenome: DetailDrawerProps['onShowInGenome'] }) {
  const sv = row.sv!;
  const openModule = useOpenModule();
  return (
    <>
      <Section title="Event">
        <dl className="grid grid-cols-2 gap-x-4 gap-y-2.5 text-[13px]">
          <Fact label="Span" value={row.change} />
          <Fact label="Length" value={sv.length ? `${sv.length.toLocaleString('en-US')} bp` : '—'} />
          <Fact label="Supporting molecules" value={formatCount(sv.supporting_reads.length)} />
          <Fact label="Known event" value={sv.known_event ? 'Yes' : 'No'} />
          <Fact label="Evidence" value={(sv.evidence_sources ?? (sv.evidence_source ? [sv.evidence_source] : [])).join(', ') || '—'} />
          <Fact label="Orientation" value={(sv.orientations ?? (sv.orientation ? [sv.orientation] : [])).join(', ') || '—'} />
          {sv.segment_count !== undefined && <Fact label="Segments" value={String(sv.segment_count)} />}
        </dl>
        <p className="mt-3 text-xs leading-5 text-muted">
          Structural support is not converted into a frequency without a matching denominator of negative observations.
        </p>
      </Section>
      {row.genes.length > 0 && (
        <Section title={`Affected genes (${row.genes.length})`}>
          <p className="text-[13px] leading-6">{row.genes.join(', ')}</p>
        </Section>
      )}
      {sv.annotation && (
        <Section title="Annotation">
          <ClinicalSummary annotation={sv.annotation} />
        </Section>
      )}
      <Actions>
        <Action onClick={() => (row.end >= row.position ? onShowInGenome(row.position, row.end, row.typeLabel) : onShowInGenome(1, Number.MAX_SAFE_INTEGER))}>Show in genome view</Action>
        <Action onClick={() => openModule('rearrangements')}>Open rearrangements</Action>
      </Actions>
    </>
  );
}

function GeneDetails({ gene, rows, onSelectRow, onSelectGene, onShowInGenome }: {
  gene: GeneAnnotation;
  rows: EventRow[];
  onSelectRow: (row: EventRow) => void;
  onSelectGene: (gene: GeneAnnotation) => void;
  onShowInGenome: DetailDrawerProps['onShowInGenome'];
}) {
  const setFilters = useMitoStore((state) => state.setEventFilters);
  const setView = useMitoStore((state) => state.setOverviewView);
  const events = rows.filter((row) => row.genes.includes(gene.name));
  return (
    <>
      <Section title="Gene">
        <dl className="grid grid-cols-2 gap-x-4 gap-y-2.5 text-[13px]">
          <Fact label="Range" value={`${gene.start.toLocaleString('en-US')}–${gene.end.toLocaleString('en-US')}`} />
          <Fact label="Length" value={`${(gene.end - gene.start + 1).toLocaleString('en-US')} bp`} />
          <Fact label="Strand" value={gene.strand} />
          <div>
            <dt className="text-muted">Biotype</dt>
            <dd className="mt-0.5 flex items-center gap-1.5 font-medium">
              <span className="h-2.5 w-2.5 rounded-sm" style={{ background: geneColor(gene.biotype) }} aria-hidden />
              {gene.biotype === 'protein_coding' ? 'Protein-coding' : gene.biotype}
            </dd>
          </div>
        </dl>
      </Section>
      <Section title={`Events in gene (${events.length})`}>
        {events.length === 0 ? <p className="text-[13px] text-muted">No events overlap this gene.</p> : (
          <ul className="grid gap-1">
            {events.slice(0, 12).map((row) => (
              <li key={row.id}>
                <button type="button" onClick={() => onSelectRow(row)} className="flex w-full items-center justify-between gap-3 rounded-lg px-2 py-1.5 text-left text-[13px] hover:bg-panel2">
                  <span className="font-mono">{row.kind === 'sv' ? row.change : `${row.position} ${row.change}`}</span>
                  <span className="text-muted">{row.hf !== undefined ? formatFraction(row.hf) : row.typeLabel}</span>
                </button>
              </li>
            ))}
          </ul>
        )}
      </Section>
      <Actions>
        <Action onClick={() => { setFilters({ gene: gene.name }); setView('table'); onSelectGene(gene); }}>{`Filter the table to ${gene.name}`}</Action>
        <Action onClick={() => onShowInGenome(gene.start - 60, gene.end + 60, gene.name)}>Show in genome view</Action>
      </Actions>
    </>
  );
}

function ClinicalSummary({ annotation }: { annotation: NonNullable<EventRow['variant']>['annotation'] }) {
  if (!annotation) return null;
  const sources = annotation.sources ?? (annotation.source ? [annotation.source] : []);
  return (
    <>
      <dl className="grid gap-2.5 text-[13px]">
        <Fact label="Clinical significance" value={annotation.consensus_significance ?? annotation.pathogenicity ?? '—'} />
        {annotation.conflict_status && <Fact label="Assertion agreement" value={annotation.conflict_status.replace(/_/g, ' ')} />}
        {annotation.phenotype && <Fact label="Phenotype" value={annotation.phenotype} />}
        {sources.length > 0 && <Fact label="Sources" value={sources.join(', ')} />}
      </dl>
      <p className="mt-2.5 text-xs leading-5 text-muted">Summarised source assertions, not a clinical verdict.</p>
    </>
  );
}

function useOpenModule() {
  const setActiveModule = useMitoStore((state) => state.setActiveModule);
  return (module: AnalysisModule) => {
    setActiveModule(module);
    window.scrollTo({ top: 0 });
  };
}

function Section({ title, children }: { title: string; children: ReactNode }) {
  return (
    <section className="border-b border-line px-5 py-4">
      <h3 className="section-kicker mb-3">{title}</h3>
      {children}
    </section>
  );
}

function Fact({ label, value }: { label: string; value: string }) {
  return (
    <div className="min-w-0">
      <dt className="text-muted">{label}</dt>
      <dd className="mt-0.5 break-words font-medium tabular-nums">{value}</dd>
    </div>
  );
}

function Actions({ children }: { children: ReactNode }) {
  return <div className="grid gap-1 px-3 py-3">{children}</div>;
}

function Action({ children, onClick }: { children: string; onClick: () => void }) {
  return (
    <button type="button" onClick={onClick} className="flex min-h-9 items-center justify-between rounded-lg px-2 text-left text-[13px] font-medium text-accent hover:bg-panel2">
      {children}
      <ArrowRight className="h-4 w-4" aria-hidden />
    </button>
  );
}
