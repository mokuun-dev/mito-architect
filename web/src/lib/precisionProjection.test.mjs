import assert from 'node:assert/strict';
import test from 'node:test';

import { filterLabel, projectPrecisionCalls, unsupportedEventCount } from './precisionProjection.ts';

const strands = { alt_forward: 5, alt_reverse: 5, depth_forward: 50, depth_reverse: 50 };
const repeat = { homopolymer_run: 1, unit: 1, span: 1 };
const atp6 = {
  gene: 'MT-ATP6',
  biotype: 'protein_coding',
  consequence: 'missense_variant',
  hgvs_p: 'p.Leu156Arg',
  residue: 156,
  structure: { structure_id: '8H9S', chain: 'N', residue_index: 156, complex: 'Complex V ATP synthase' }
};

function snv(position, status, effects = []) {
  return {
    id: `snv:${position}:T:G`, position, ref: 'T', alt: 'G', status, filters: [], flags: [], fraction: 0.1,
    ci95_low: 0.05, ci95_high: 0.2, depth: 100, alt_depth: 10, strands, p_value: 1e-9, q_value: 1e-7,
    background_rate: 0.001, alt_mean_quality: 20, repeat, effects
  };
}

function result(precision) {
  return {
    metadata: { reference_length: 16569 },
    filter_stats: {},
    genes: [],
    coverage: [],
    clusters: [],
    reads: [],
    svs: [
      { id: 'deletion:10-90', type: 'deletion', start: 10, end: 90, length: 81, known_event: false, supporting_reads: ['a'] },
      { id: 'soft_clip_left:5+30', type: 'soft_clip_left', start: 5, end: 5, length: 30, known_event: false, supporting_reads: ['a', 'b', 'c'] },
      { id: 'inversion:100-400', type: 'inversion', start: 100, end: 400, length: 301, known_event: false, supporting_reads: ['a', 'b', 'c'] },
      { id: 'inversion:500-800', type: 'inversion', start: 500, end: 800, length: 301, known_event: false, supporting_reads: ['a'] }
    ],
    variants: [
      { position: 8993, ref: 'T', alt: 'G', alt_depth: 10, ref_depth: 90, other_depth: 0, callable_depth: 100, heteroplasmy: 0.1, ci95_low: 0, ci95_high: 1,
        supporting_reads: ['read-1'], annotation: { pathogenicity: 'pathogenic' } },
      { position: 9000, ref: 'T', alt: 'G', alt_depth: 1, ref_depth: 99, other_depth: 0, callable_depth: 100, heteroplasmy: 0.01, ci95_low: 0, ci95_high: 1, supporting_reads: [] }
    ],
    precision_calls: precision
  };
}

test('a result without the statistical caller is returned unchanged', () => {
  const legacy = result(undefined);
  assert.equal(projectPrecisionCalls(legacy), legacy);
  const unassessed = result({ schema_version: '1.0', assessed: false });
  assert.equal(projectPrecisionCalls(unassessed), unassessed);
});

test('only supported calls reach the variant projection, with their evidence', () => {
  const projected = projectPrecisionCalls(result({
    schema_version: '1.0', assessed: true,
    snvs: [snv(8993, 'PASS', [atp6]), snv(9000, 'FILTERED')],
    indels: [],
    deletions: []
  }));
  assert.equal(projected.variants.length, 1);
  const [variant] = projected.variants;
  assert.equal(variant.heteroplasmy, 0.1);
  assert.equal(variant.residue, 'p.Leu156Arg');
  assert.equal(variant.structure.chain, 'N');
  assert.deepEqual(variant.supporting_reads, ['read-1']);
  assert.equal(variant.annotation.pathogenicity, 'pathogenic');
  assert.equal(unsupportedEventCount(result(undefined), projected), 1);
});

test('non-coding effects carry no protein change', () => {
  const projected = projectPrecisionCalls(result({
    schema_version: '1.0', assessed: true,
    snvs: [snv(8993, 'PASS', [{ gene: 'MT-TL1', biotype: 'tRNA', consequence: 'tRNA_variant' }])],
    indels: [], deletions: []
  }));
  assert.equal(projected.variants[0].gene, 'MT-TL1');
  assert.equal(projected.variants[0].protein, undefined);
  assert.equal(projected.variants[0].residue, undefined);
});

test('called deletions replace raw gaps; well-supported rearrangements remain', () => {
  const deletion = {
    id: 'pdeletion:8483-13459', start: 8483, end: 13459, left_aligned_start: 8470, left_aligned_end: 13446, length: 4977,
    microhomology: 'ACCTCCCTCACCA', status: 'PASS', filters: [], flags: [], fraction: 0.2, ci95_low: 0.15, ci95_high: 0.25,
    junction_molecules: 20, junction_forward: 10, junction_reverse: 10, wild_type_left: 80, wild_type_right: 80,
    breakpoint_variants: 3, depth_ratio: 0.8, p_value: 1e-20, q_value: 1e-20, known_common_deletion: true,
    impact: { genes_lost: [], genes_partial: [], trnas_affected: 5, heavy_strand_origin_lost: false, light_strand_origin_lost: false, major_arc: true, contact_zone: true },
    supporting_molecule_count: 20, supporting_molecules: ['m1']
  };
  const projected = projectPrecisionCalls(result({ schema_version: '1.0', assessed: true, snvs: [], indels: [], deletions: [deletion] }));
  assert.deepEqual(projected.svs.map((sv) => sv.id), ['deletion:8483-13459', 'inversion:100-400']);
  assert.equal(projected.svs[0].precision.fraction, 0.2);
  assert.equal(projected.svs[0].known_event, true);
});

test('filter reasons read as sentences, and unknown ones stay legible', () => {
  assert.equal(filterLabel('strand_bias'), 'Strong imbalance between strands');
  assert.equal(filterLabel('new_filter_name'), 'new filter name');
});
