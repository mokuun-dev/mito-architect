import assert from 'node:assert/strict';
import test from 'node:test';
import { circularVariants } from './dataAdapter.ts';

test('circular map prefers normalized aggregate variants over legacy read SNPs', () => {
  const variants = circularVariants({
    variants: [{ event_id: 'snv:73:A:G', position: 73, ref: 'A', alt: 'G', alt_depth: 4, callable_depth: 10, heteroplasmy: 0.4 }],
    reads: [{ filtered_numt: false, snps: [{ position: 73, ref: 'A', alt: 'G' }, { position: 263, ref: 'A', alt: 'G' }] }]
  });
  assert.deepEqual(variants, [{ id: 'snv:73:A:G', position: 73, ref: 'A', alt: 'G', type: 'SNV', support: 4, callableDepth: 10, frequency: 0.4 }]);
});

test('legacy aggregation excludes NUMT-filtered records', () => {
  const variants = circularVariants({
    reads: [
      { filtered_numt: false, snps: [{ position: 73, ref: 'A', alt: 'G' }] },
      { filtered_numt: true, snps: [{ position: 73, ref: 'A', alt: 'G' }] }
    ]
  });
  assert.equal(variants[0].support, 1);
  assert.equal(variants[0].frequency, 1);
});
