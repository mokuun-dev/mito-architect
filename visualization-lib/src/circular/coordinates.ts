/**
 * Coordinate helpers for a circular reference.
 *
 * Public result payloads use inclusive, one-based positions. Geometry is
 * calculated from boundaries (0..referenceLength), so the final base reaches
 * the end of the circle instead of being clamped to its own start.
 */
export interface CircularInterval {
  start: number;
  end: number;
}

export function clampBase(position: number, referenceLength: number): number {
  return Math.max(1, Math.min(referenceLength, Math.round(position)));
}

export function baseBoundary(position: number, referenceLength: number): number {
  return Math.max(0, Math.min(referenceLength, position));
}

/** Converts a one-based inclusive interval into one or two zero-based spans. */
export function inclusiveIntervalToSpans(
  start: number,
  end: number,
  referenceLength: number,
  wrapsReference = false
): CircularInterval[] {
  const first = clampBase(start, referenceLength) - 1;
  const last = clampBase(end, referenceLength);
  if (wrapsReference && first >= last) {
    return [{ start: first, end: referenceLength }, { start: 0, end: last }];
  }
  return [{ start: Math.min(first, last - 1), end: Math.max(first + 1, last) }];
}

export function positionToAngle(position: number, referenceLength: number): number {
  return (baseBoundary(position, referenceLength) / Math.max(1, referenceLength)) * Math.PI * 2;
}

export function intervalMidpoint(start: number, end: number, referenceLength: number): number {
  const spans = inclusiveIntervalToSpans(start, end, referenceLength);
  if (spans.length === 1) return (spans[0].start + spans[0].end) / 2;
  const length = spans.reduce((sum, span) => sum + span.end - span.start, 0);
  return (spans[0].start + length / 2) % referenceLength;
}
