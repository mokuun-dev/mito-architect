import type { AggregateVariant } from '@mito-architect/visualization-lib';
import type { ClinicalAnnotation } from '@mito-architect/visualization-lib';
import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import type { ProteinStructureMapping } from '@mito-architect/visualization-lib';
import type { ReadFeature } from '@mito-architect/visualization-lib';
export { summarizeCoverage } from '../../lib/coverageModel';

export interface VariantSummary {
  key: string;
  label: string;
  type: string;
  position: number;
  support: number;
  frequency: number;
  gene?: string;
  consequence?: string;
  protein?: string;
  residue?: string;
  annotation?: ClinicalAnnotation;
  structure?: ProteinStructureMapping;
  molecules: Array<{ id: string; clusterId: number }>;
  callableDepth?: number;
  ci95Low?: number;
  ci95High?: number;
  strandSupport?: AggregateVariant['strand_support'];
  strandBiasDelta?: number | null;
  alleleQuality?: AggregateVariant['allele_quality'];
  readPosition?: AggregateVariant['read_position'];
}

export function summarizeVariants(reads: ReadFeature[], data: MitoAnalysisData): VariantSummary[] {
  const denominator = Math.max(1, reads.length);
  const snps = new Map<string, VariantSummary>();
  const readClusterById = new Map(reads.map((read) => [read.id, read.cluster_id]));

  if (data.variants) {
    for (const snp of data.variants) {
      const key = snp.event_id ?? `snp:${snp.position}:${snp.ref}:${snp.alt}`;
      const type = snp.type ?? 'SNV';
      snps.set(key, {
        key,
        label: type === 'SNV'
          ? `${snp.position} ${snp.ref}>${snp.alt}`
          : `${type.replace('SMALL_', '').toLowerCase()} ${snp.position} ${snp.ref}>${snp.alt}`,
        type,
        position: snp.position,
        support: snp.alt_depth,
        frequency: snp.heteroplasmy,
        gene: snp.gene,
        consequence: snp.consequence,
        protein: snp.protein,
        residue: snp.residue,
        annotation: snp.annotation,
        structure: snp.structure,
        molecules: (snp.supporting_molecule_ids ?? snp.supporting_reads).map((id) => ({
          id,
          clusterId: readClusterById.get(id) ?? -1
        })),
        callableDepth: snp.callable_depth,
        ci95Low: snp.ci95_low,
        ci95High: snp.ci95_high,
        strandSupport: snp.strand_support,
        strandBiasDelta: snp.strand_bias_delta,
        alleleQuality: snp.allele_quality,
        readPosition: snp.read_position
      });
    }
  } else {
    for (const read of reads) {
      for (const snp of read.snps) {
      const key = `snp:${snp.position}:${snp.ref}:${snp.alt}`;
      const existing = snps.get(key);
      if (existing) {
        existing.support += 1;
        existing.frequency = existing.support / denominator;
        existing.molecules.push({ id: read.id, clusterId: read.cluster_id });
      } else {
        snps.set(key, {
          key,
          label: `${snp.position} ${snp.ref}>${snp.alt}`,
          type: 'SNP',
          position: snp.position,
          support: 1,
          frequency: 1 / denominator,
          gene: snp.gene,
          consequence: snp.consequence,
          protein: snp.protein,
          residue: snp.residue,
          annotation: snp.annotation,
          structure: snp.structure,
          molecules: [{ id: read.id, clusterId: read.cluster_id }]
        });
      }
    }
    }
  }

  const svs = data.svs.map((sv) => ({
    key: `sv:${sv.id}`,
    label: `${sv.type} ${sv.start}-${sv.end}`,
    type: sv.type,
    position: sv.start,
    support: sv.supporting_reads.length,
    // A structural-event negative-evidence denominator is not guaranteed by legacy payloads.
    frequency: Number.NaN,
    annotation: sv.annotation,
    molecules: sv.supporting_reads.map((id) => ({
      id,
      clusterId: readClusterById.get(id) ?? -1
    }))
  }));

  return [...snps.values(), ...svs].sort((a, b) => a.position - b.position || a.label.localeCompare(b.label));
}
