import assert from 'node:assert/strict';
import test from 'node:test';
import { inclusiveIntervalToSpans, positionToAngle } from './coordinates.ts';

test('last inclusive base reaches the circular endpoint', () => {
  assert.equal(positionToAngle(16569, 16569), Math.PI * 2);
  assert.deepEqual(inclusiveIntervalToSpans(16569, 16569, 16569), [{ start: 16568, end: 16569 }]);
});

test('declared wrap interval is split at the reference boundary', () => {
  assert.deepEqual(inclusiveIntervalToSpans(16560, 10, 16569, true), [
    { start: 16559, end: 16569 },
    { start: 0, end: 10 }
  ]);
});
