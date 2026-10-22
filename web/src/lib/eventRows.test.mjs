import assert from 'node:assert/strict';
import test from 'node:test';

import {
  buildEventRows,
  DEFAULT_EVENT_FILTERS,
  filterEventRows,
  genesAt,
  sortEventRows
} from './eventRows.ts';

const genes = [
  { name: 'MT-TL1', start: 3230, end: 3304, strand: '+', biotype: 'tRNA' },
  { name: 'MT-ND1', start: 3307, end: 4262, strand: '+', biotype: 'protein_coding' },
  { name: 'MT-CYB', start: 14747, end: 15887, strand: '+', biotype: 'protein_coding' }
];
const metadata = { reference_length: 16569 };

function aggregate(position, ref, alt, overrides = {}) {
  return {
    event_id: `snv:${position}:${ref}:${alt}`,
    type: 'SNV',
    position,
    ref,
    alt,
    alt_depth: 10,
    ref_depth: 90,
    other_depth: 0,
    callable_depth: 100,
    heteroplasmy: 0.1,
    ci95_low: 0.05,
    ci95_high: 0.17,
    supporting_reads: [],
    ...overrides
  };
}

test('aggregate variants keep frequencies only with a callable denominator', () => {
  const rows = buildEventRows({
    metadata,
    genes,
    reads: [],
    svs: [],
    variants: [aggregate(3243, 'A', 'G'), aggregate(3500, 'C', 'T', { callable_depth: 0, heteroplasmy: 0 })]
  });
  assert.equal(rows.length, 2);
  assert.equal(rows[0].hf, 0.1);
  assert.deepEqual(rows[0].genes, ['MT-TL1']);
  assert.equal(rows[1].hf, undefined);
  assert.equal(rows[1].ciLow, undefined);
});

test('legacy reads aggregate support without inventing a frequency', () => {
  const read = (id, filtered) => ({ id, filtered_numt: filtered, cluster_id: 0, snps: [{ position: 3243, ref: 'A', alt: 'G' }] });
  const rows = buildEventRows({ metadata, genes, svs: [], reads: [read('a', false), read('b', false), read('c', true)] });
  assert.equal(rows.length, 1);
  assert.equal(rows[0].id, 'snp:3243:A:G');
  assert.equal(rows[0].altSupport, 2);
  assert.equal(rows[0].hf, undefined);
});

test('structural events list every gene they overlap, across the origin', () => {
  assert.deepEqual(genesAt(genes, 15000, 3250, 16569), ['MT-TL1', 'MT-CYB']);
  const rows = buildEventRows({
    metadata,
    genes,
    reads: [],
    variants: [],
    svs: [{ id: 'sv1', type: 'deletion', start: 3300, end: 15000, length: 11700, known_event: false, supporting_reads: ['m1', 'm2'] }]
  });
  assert.equal(rows[0].kind, 'sv');
  assert.equal(rows[0].altSupport, 2);
  assert.deepEqual(rows[0].genes, ['MT-TL1', 'MT-ND1', 'MT-CYB']);
});

test('filters combine kind, gene, frequency and QC constraints', () => {
  const rows = buildEventRows({
    metadata,
    genes,
    reads: [],
    variants: [
      aggregate(3243, 'A', 'G', { heteroplasmy: 0.4 }),
      aggregate(3400, 'C', 'T', { heteroplasmy: 0.02, qc_flags: ['LOW_ALT_SUPPORT'] }),
      aggregate(3500, 'C', 'CA', { type: 'SMALL_INSERTION', callable_depth: 0 })
    ],
    svs: [{ id: 'sv1', type: 'deletion', start: 3300, end: 15000, length: 11700, known_event: false, supporting_reads: [] }]
  });
  const ids = (filters) => filterEventRows(rows, { ...DEFAULT_EVENT_FILTERS, ...filters }).map((row) => row.id);
  assert.deepEqual(ids({ kinds: ['sv'] }), ['sv1']);
  assert.deepEqual(ids({ kinds: ['indel'] }), ['snv:3500:C:CA']);
  assert.deepEqual(ids({ gene: 'MT-TL1' }), ['snv:3243:A:G', 'sv1']);
  assert.deepEqual(ids({ minHf: 0.05 }), ['snv:3243:A:G']);
  assert.deepEqual(ids({ hideFlagged: true, kinds: ['snv'] }), ['snv:3243:A:G']);
  assert.deepEqual(ids({ query: 'low_alt' }), ['snv:3400:C:T']);
});

test('sorting keeps rows without a value last in both directions', () => {
  const rows = buildEventRows({
    metadata,
    genes,
    reads: [],
    svs: [],
    variants: [
      aggregate(100, 'A', 'G', { heteroplasmy: 0.2 }),
      aggregate(200, 'A', 'G', { callable_depth: 0 }),
      aggregate(300, 'A', 'G', { heteroplasmy: 0.9 })
    ]
  });
  assert.deepEqual(sortEventRows(rows, 'hf', true).map((row) => row.position), [300, 100, 200]);
  assert.deepEqual(sortEventRows(rows, 'hf', false).map((row) => row.position), [100, 300, 200]);
});
