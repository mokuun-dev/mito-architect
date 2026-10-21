import assert from 'node:assert/strict';
import test from 'node:test';
import { coverageDistribution, summarizeCoverage } from './coverageModel.ts';

test('coverage fallback weights bins by represented reference bases', () => {
  const summary = summarizeCoverage({
    coverage: [
      { start: 1, end: 1, depth: 100 },
      { start: 2, end: 101, depth: 0 }
    ]
  });
  assert.equal(summary.meanDepth, 100 / 101);
  assert.equal(summary.pctGt20, 1 / 101 * 100);
  assert.equal(summary.pctLabel, 'sites');
});

test('coverage distribution retains zero-depth intervals as zero', () => {
  assert.deepEqual(coverageDistribution([
    { start: 1, end: 10, depth: 0 },
    { start: 11, end: 15, depth: 25 }
  ]), [
    { label: '0×', count: 10 },
    { label: '1–9×', count: 0 },
    { label: '10–19×', count: 0 },
    { label: '20–49×', count: 5 },
    { label: '50×+', count: 0 }
  ]);
});
