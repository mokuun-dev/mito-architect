/** Formats a 0–1 fraction as a percentage with precision suited to its size. */
export function formatFraction(value: number): string {
  if (!Number.isFinite(value)) return '—';
  const percent = value * 100;
  if (percent !== 0 && Math.abs(percent) < 1) return `${percent.toFixed(2)}%`;
  return `${percent.toFixed(1)}%`;
}

export function formatCount(value: number | undefined): string {
  return value === undefined ? '—' : value.toLocaleString('en-US');
}
