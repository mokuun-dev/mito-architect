import type {
  AggregateVariant,
  GeneAnnotation,
  MitoAnalysisData,
  ReadFeature,
  StructuralVariant
} from '@mito-architect/visualization-lib';

export type EventKind = 'snv' | 'indel' | 'sv';
export type EventSortKey = 'position' | 'type' | 'gene' | 'hf' | 'support' | 'callable';

/** One row of the unified event table: a small variant or a structural event. */
export interface EventRow {
  /** Selection identifier: the event ID of a variant or the ID of an SV. */
  id: string;
  kind: EventKind;
  type: string;
  typeLabel: string;
  position: number;
  end: number;
  length?: number;
  ref?: string;
  alt?: string;
  change: string;
  genes: string[];
  consequence?: string;
  /** Only present when a callable denominator was reported for the event. */
  hf?: number;
  ciLow?: number;
  ciHigh?: number;
  altSupport: number;
  callable?: number;
  numt?: string;
  qcFlags: string[];
  clinical?: string;
  variant?: AggregateVariant;
  sv?: StructuralVariant;
}

export interface EventFilters {
  query: string;
  kinds: EventKind[];
  gene?: string;
  minHf: number;
  hideFlagged: boolean;
}

export const DEFAULT_EVENT_FILTERS: EventFilters = {
  query: '',
  kinds: [],
  gene: undefined,
  minHf: 0,
  hideFlagged: false
};

export function buildEventRows(data: Pick<MitoAnalysisData, 'variants' | 'reads' | 'svs' | 'genes' | 'metadata'>): EventRow[] {
  const length = data.metadata.reference_length;
  const variants = data.variants
    ? data.variants.map((variant) => aggregateRow(variant, data.genes))
    : legacyRows(data.reads, data.genes);
  const svs = data.svs.map((sv) => structuralRow(sv, data.genes, length));
  return [...variants, ...svs].sort((a, b) => a.position - b.position || a.id.localeCompare(b.id));
}

export function filterEventRows(rows: EventRow[], filters: EventFilters): EventRow[] {
  const query = filters.query.trim().toLowerCase();
  return rows.filter((row) => {
    if (filters.kinds.length > 0 && !filters.kinds.includes(row.kind)) return false;
    if (filters.gene && !row.genes.includes(filters.gene)) return false;
    if (filters.minHf > 0 && (row.hf === undefined || row.hf < filters.minHf)) return false;
    if (filters.hideFlagged && row.qcFlags.length > 0) return false;
    if (!query) return true;
    return [row.id, row.change, row.typeLabel, String(row.position), row.consequence, row.clinical, ...row.genes, ...row.qcFlags]
      .some((value) => value !== undefined && value.toLowerCase().includes(query));
  });
}

/** Sorts by the key; rows without a value sort last in both directions. */
export function sortEventRows(rows: EventRow[], key: EventSortKey, descending: boolean): EventRow[] {
  const value = (row: EventRow): number | string | undefined => {
    switch (key) {
      case 'position': return row.position;
      case 'type': return row.typeLabel;
      case 'gene': return row.genes[0];
      case 'hf': return row.hf;
      case 'support': return row.altSupport;
      case 'callable': return row.callable;
    }
  };
  return [...rows].sort((a, b) => {
    const left = value(a);
    const right = value(b);
    if (left === undefined || right === undefined) {
      if (left === right) return a.position - b.position;
      return left === undefined ? 1 : -1;
    }
    const order = typeof left === 'number' && typeof right === 'number'
      ? left - right
      : String(left).localeCompare(String(right));
    return (descending ? -order : order) || a.position - b.position || a.id.localeCompare(b.id);
  });
}

export function genesAt(genes: GeneAnnotation[], start: number, end: number, length: number): string[] {
  // Circular events may wrap the origin (end < start).
  const spans = end >= start ? [[start, end]] : [[start, length], [1, end]];
  return genes
    .filter((gene) => spans.some(([from, to]) => gene.start <= to && gene.end >= from))
    .map((gene) => gene.name);
}

function aggregateRow(variant: AggregateVariant, genes: GeneAnnotation[]): EventRow {
  const type = variant.type ?? 'SNV';
  const end = variant.end ?? variant.position;
  const hasDenominator = variant.callable_depth > 0;
  return {
    id: variant.event_id ?? `snp:${variant.position}:${variant.ref}:${variant.alt}`,
    kind: type === 'SNV' ? 'snv' : 'indel',
    type,
    typeLabel: variantTypeLabel(type),
    position: variant.position,
    end,
    length: variant.length,
    ref: variant.ref,
    alt: variant.alt,
    change: `${variant.ref}>${variant.alt}`,
    genes: variant.gene ? [variant.gene] : genesAt(genes, variant.position, end, Number.MAX_SAFE_INTEGER),
    consequence: variant.consequence,
    hf: hasDenominator ? variant.heteroplasmy : undefined,
    ciLow: hasDenominator ? variant.ci95_low : undefined,
    ciHigh: hasDenominator ? variant.ci95_high : undefined,
    altSupport: variant.alt_depth,
    callable: variant.callable_depth,
    numt: variant.numt_assessability,
    qcFlags: variant.qc_flags ?? [],
    clinical: clinicalLabel(variant.annotation),
    variant
  };
}

function legacyRows(reads: ReadFeature[], genes: GeneAnnotation[]): EventRow[] {
  const rows = new Map<string, EventRow>();
  for (const read of reads) {
    if (read.filtered_numt) continue;
    for (const snp of read.snps) {
      const id = `snp:${snp.position}:${snp.ref}:${snp.alt}`;
      const existing = rows.get(id);
      if (existing) {
        existing.altSupport += 1;
        continue;
      }
      rows.set(id, {
        id,
        kind: 'snv',
        type: 'SNV',
        typeLabel: 'SNV',
        position: snp.position,
        end: snp.position,
        ref: snp.ref,
        alt: snp.alt,
        change: `${snp.ref}>${snp.alt}`,
        genes: snp.gene ? [snp.gene] : genesAt(genes, snp.position, snp.position, Number.MAX_SAFE_INTEGER),
        consequence: snp.consequence,
        altSupport: 1,
        qcFlags: [],
        clinical: clinicalLabel(snp.annotation)
      });
    }
  }
  return [...rows.values()];
}

function structuralRow(sv: StructuralVariant, genes: GeneAnnotation[], length: number): EventRow {
  return {
    id: sv.id,
    kind: 'sv',
    type: sv.type,
    typeLabel: structuralTypeLabel(sv.type),
    position: sv.start,
    end: sv.end,
    length: sv.length,
    change: `${sv.start.toLocaleString('en-US')}–${sv.end.toLocaleString('en-US')}`,
    genes: genesAt(genes, sv.start, sv.end, length),
    altSupport: sv.supporting_reads.length,
    qcFlags: [],
    clinical: clinicalLabel(sv.annotation),
    sv
  };
}

function variantTypeLabel(type: string): string {
  if (type === 'SNV') return 'SNV';
  if (type === 'SMALL_INSERTION') return 'Insertion';
  if (type === 'SMALL_DELETION') return 'Deletion';
  return titleCase(type.replace(/_/g, ' ').toLowerCase());
}

const SV_ABBREVIATIONS: Record<string, string> = { del: 'deletion', ins: 'insertion', dup: 'duplication', inv: 'inversion' };

/** "SV deletion", "SV soft clip right": distinct from small indels of the same kind. */
function structuralTypeLabel(type: string): string {
  const name = SV_ABBREVIATIONS[type.toLowerCase()] ?? type.toLowerCase().replace(/^sv_/, '').replace(/_/g, ' ');
  return `SV ${name}`;
}

function titleCase(value: string): string {
  return value.charAt(0).toUpperCase() + value.slice(1);
}

function clinicalLabel(annotation: AggregateVariant['annotation']): string | undefined {
  return annotation?.consensus_significance ?? annotation?.pathogenicity ?? undefined;
}
