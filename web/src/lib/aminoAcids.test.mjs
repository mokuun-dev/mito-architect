import assert from 'node:assert/strict';
import test from 'node:test';

import { aminoAcid, describeSubstitution, granthamDistance, parseProteinChange, residueMatches } from './aminoAcids.ts';

test('protein changes emitted by the engine parse', () => {
  const missense = parseProteinChange('p.Leu156Arg');
  assert.equal(missense.kind, 'missense');
  assert.equal(missense.position, 156);
  assert.equal(missense.reference.three, 'Leu');
  assert.equal(missense.alternate.three, 'Arg');
  assert.equal(parseProteinChange('p.Leu156=').kind, 'synonymous');
  assert.equal(parseProteinChange('p.Trp8Ter').kind, 'nonsense');
  assert.equal(parseProteinChange('p.Leu156fs').kind, 'frameshift');
  assert.equal(parseProteinChange('not a change'), undefined);
});

test('Grantham distances match the published values', () => {
  assert.equal(granthamDistance(aminoAcid('L'), aminoAcid('R')), 102);
  assert.equal(granthamDistance(aminoAcid('R'), aminoAcid('H')), 29);
  assert.equal(describeSubstitution(aminoAcid('L'), aminoAcid('R')).granthamClass, 'moderately radical');
});

test('residue identity check accepts only the expected residue', () => {
  assert.equal(residueMatches(aminoAcid('Leu'), 'LEU'), true);
  assert.equal(residueMatches(aminoAcid('Leu'), 'ARG'), false);
  assert.equal(residueMatches(aminoAcid('Met'), 'FME'), true);
  assert.equal(residueMatches(undefined, 'LEU'), false);
});
