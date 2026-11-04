import assert from 'node:assert/strict';
import test from 'node:test';

import {
  normalizePdbId,
  rcsbCandidates
} from './proteinStructureSources.ts';

test('RCSB sources are deterministic and local-first', () => {
  const candidates = rcsbCandidates('8h9s', 'N');
  assert.equal(candidates.length, 3);
  assert.equal(candidates[0].url, '/structures/8h9s.bcif');
  assert.equal(candidates[0].chain, 'N');
  assert.equal(candidates[0].confidenceSemantics, 'bfactor');
  assert.equal(candidates[1].url, 'https://models.rcsb.org/8h9s.bcif');
  assert.equal(candidates[2].url, 'https://files.rcsb.org/download/8H9S.cif');
});

test('structure identifiers are validated and normalized', () => {
  assert.equal(normalizePdbId(' 9i4i '), '9I4I');
  assert.equal(normalizePdbId('Q00360'), undefined);
  assert.equal(normalizePdbId('AF-P00846-F1'), undefined);
});

test('chain IDs are refused rather than replaced', async () => {
  const { validChainId } = await import('./proteinStructureSources.ts');
  assert.equal(validChainId('r'), 'r');
  assert.equal(validChainId('AA'), 'AA');
  assert.equal(validChainId('N or 1'), undefined);
  assert.equal(validChainId(undefined), undefined);
});

test('the structure table matches the verified catalog', async () => {
  const { readFile } = await import('node:fs/promises');
  const { COMPLEX_STRUCTURES } = await import('./proteinStructureSources.ts');
  const catalog = await readFile(new URL('../../../core/data/mt_structure_catalog.tsv', import.meta.url), 'utf8');
  const rows = catalog.trim().split('\n').slice(1).map((line) => line.split('\t'));
  assert.equal(rows.length, 13);
  for (const [gene, , pdbId, chain] of rows) {
    const entry = COMPLEX_STRUCTURES.find((item) => item.pdbId === pdbId);
    assert.ok(entry?.subunits.some((subunit) => subunit.gene === gene && subunit.chain === chain), gene);
  }
});
