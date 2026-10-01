import type { AggregateVariant, MitoAnalysisData, ReadFeature, SnpCall } from '../types';

export interface CircularVariant {
  id: string;
  position: number;
  ref: string;
  alt: string;
  type: string;
  support: number;
  callableDepth?: number;
  frequency?: number;
}

/** Canonical aggregate variants, with a lossless legacy read-level fallback. */
export function circularVariants(data: MitoAnalysisData): CircularVariant[] {
  if (data.variants?.length) {
    return data.variants.map((variant) => fromAggregateVariant(variant));
  }
  return legacyVariants(data.reads);
}

export function legacyEventId(snp: SnpCall): string {
  return `snp:${snp.position}:${snp.ref}:${snp.alt}`;
}

function fromAggregateVariant(variant: AggregateVariant): CircularVariant {
  return {
    id: variant.event_id ?? legacyEventId(variant),
    position: variant.position,
    ref: variant.ref,
    alt: variant.alt,
    type: variant.type ?? 'SNV',
    support: variant.alt_depth,
    callableDepth: variant.callable_depth,
    frequency: variant.heteroplasmy
  };
}

function legacyVariants(reads: ReadFeature[]): CircularVariant[] {
  const variants = new Map<string, CircularVariant>();
  for (const read of reads) {
    if (read.filtered_numt) continue;
    for (const snp of read.snps) {
      const id = legacyEventId(snp);
      const previous = variants.get(id);
      if (previous) previous.support += 1;
      else variants.set(id, { id, position: snp.position, ref: snp.ref, alt: snp.alt, type: 'SNV', support: 1 });
    }
  }
  const denominator = Math.max(1, reads.filter((read) => !read.filtered_numt).length);
  return [...variants.values()].map((variant) => ({ ...variant, frequency: variant.support / denominator }));
}
