import type {
  AggregateVariant,
  MitoAnalysisData,
  PrecisionDeletion,
  PrecisionEffect,
  PrecisionIndel,
  PrecisionSnv,
  StructuralVariant
} from '@mito-architect/visualization-lib';

/** Effect shown for an event overall: a mapped protein change, if there is one. */
function primaryEffect(effects: PrecisionEffect[]): PrecisionEffect | undefined {
  return (
    effects.find((effect) => effect.structure && effect.consequence !== 'synonymous_variant') ??
    effects.find((effect) => effect.biotype === 'protein_coding') ??
    effects[0]
  );
}

function effectFields(effects: PrecisionEffect[]) {
  const effect = primaryEffect(effects);
  const coding = effect?.biotype === 'protein_coding';
  return {
    effects,
    gene: effect?.gene,
    consequence: effect?.consequence,
    protein: coding ? effect?.gene : undefined,
    residue: coding ? effect?.hgvs_p : undefined,
    structure: effect?.structure ? { ...effect.structure, residue_label: effect.hgvs_p } : undefined
  };
}

function snvVariant(call: PrecisionSnv, original: AggregateVariant | undefined): AggregateVariant {
  const reference = original?.ref_depth ?? Math.max(0, call.depth - call.alt_depth);
  return {
    ...original,
    type: 'SNV',
    position: call.position,
    ref: call.ref,
    alt: call.alt,
    alt_depth: call.alt_depth,
    ref_depth: reference,
    other_depth: original?.other_depth ?? Math.max(0, call.depth - call.alt_depth - reference),
    callable_depth: call.depth,
    heteroplasmy: call.fraction,
    ci95_low: call.ci95_low,
    ci95_high: call.ci95_high,
    supporting_reads: original?.supporting_reads ?? [],
    strand_support: {
      alt_forward: call.strands.alt_forward,
      alt_reverse: call.strands.alt_reverse,
      ref_forward: original?.strand_support?.ref_forward ?? 0,
      ref_reverse: original?.strand_support?.ref_reverse ?? 0,
      other_forward: original?.strand_support?.other_forward ?? 0,
      other_reverse: original?.strand_support?.other_reverse ?? 0
    },
    homopolymer_context: { reference_base: call.ref, run_length: call.repeat.homopolymer_run },
    filter_status: call.status,
    qc_flags: call.flags,
    precision: call,
    ...effectFields(call.effects),
    annotation: original?.annotation
  };
}

function indelVariant(call: PrecisionIndel): AggregateVariant {
  const depth = Math.round(call.effective_depth);
  const alternate = Math.round(call.effective_alt_depth);
  return {
    event_id: call.id,
    type: call.type === 'deletion' ? 'SMALL_DELETION' : 'SMALL_INSERTION',
    position: call.vcf_position,
    start: call.position,
    end: call.type === 'deletion' ? call.position + call.length - 1 : call.position,
    length: call.length,
    ref: call.vcf_ref,
    alt: call.vcf_alt,
    hgvs: call.hgvs,
    vcf_position: call.vcf_position,
    alt_depth: alternate,
    ref_depth: Math.max(0, depth - alternate),
    other_depth: 0,
    callable_depth: depth,
    heteroplasmy: call.fraction,
    ci95_low: call.ci95_low,
    ci95_high: call.ci95_high,
    supporting_reads: call.supporting_molecules,
    homopolymer_context: { reference_base: null, run_length: call.repeat.homopolymer_run },
    filter_status: call.status,
    qc_flags: call.flags,
    precision: call,
    ...effectFields(call.effects)
  };
}

function deletionVariant(call: PrecisionDeletion): StructuralVariant {
  return {
    id: `deletion:${call.start}-${call.end}`,
    type: 'deletion',
    start: call.start,
    end: call.end,
    length: call.length,
    known_event: call.known_common_deletion,
    supporting_reads: call.supporting_molecules,
    evidence_sources: ['consensus_junction'],
    precision: call
  };
}

/** Types of raw event the deletion caller does not re-call: rearrangements. */
const RETAINED_SV_TYPES = new Set(['inversion', 'duplication', 'ambiguous_adjacency', 'insertion']);

/**
 * The result as the interface should present it: statistically supported calls
 * in the variant and structural projections. A result without the statistical
 * caller, or one where it could not run, is returned unchanged.
 */
export function projectPrecisionCalls(data: MitoAnalysisData): MitoAnalysisData {
  const calls = data.precision_calls;
  if (!calls?.assessed) return data;
  const originals = new Map((data.variants ?? []).map((variant) => [`${variant.position}:${variant.ref}:${variant.alt}`, variant]));
  const snvs = (calls.snvs ?? [])
    .filter((call) => call.status === 'PASS')
    .map((call) => snvVariant(call, originals.get(`${call.position}:${call.ref}:${call.alt}`)));
  const indels = (calls.indels ?? []).filter((call) => call.status === 'PASS').map(indelVariant);
  const deletions = (calls.deletions ?? []).filter((call) => call.status === 'PASS').map(deletionVariant);
  const retained = data.svs.filter(
    (sv) => RETAINED_SV_TYPES.has(sv.type) && sv.supporting_reads.length >= 3 && (sv.type !== 'insertion' || sv.length >= 50)
  );
  return {
    ...data,
    variants: [...snvs, ...indels].sort((a, b) => a.position - b.position),
    svs: [...deletions, ...retained]
  };
}

/** Raw events the statistical caller set aside; zero if nothing is hidden. */
export function unsupportedEventCount(data: MitoAnalysisData, projected: MitoAnalysisData): number {
  return Math.max(0, (data.variants?.length ?? 0) - (projected.variants?.filter((variant) => variant.type === 'SNV').length ?? 0));
}

const FILTER_LABELS: Record<string, string> = {
  not_significant: 'Not significant after multiple-testing correction',
  below_error_floor: 'Within twice the error background',
  strand_bias: 'Strong imbalance between strands',
  one_strand_support: 'No support on one strand',
  low_alternate_base_quality: 'Alternate bases have low quality',
  read_end_bias: 'Alternate bases cluster at read ends',
  gap_adjacent: 'Alternate bases sit mostly next to gaps, as an aligned length error does',
  repeat_alignment_artifact: 'Within the slip rate of a long homopolymer, where slipped reads align as mismatches',
  low_support: 'Too few molecules',
  low_lod: 'Allele likelihood too weak',
  realignment_only_support: 'Few reads show the gap in their alignment',
  single_strand_support: 'All junction molecules on one strand',
  linearization_origin: 'At the rCRS origin, where reads are split'
};

export function filterLabel(filter: string): string {
  return FILTER_LABELS[filter] ?? filter.replace(/_/g, ' ');
}
