

export function Metric({ label, value }: { label: string; value: string }) {
  return (
    <div className="min-w-0 rounded-lg border border-line bg-panel2 p-3">
      <div className="text-xs font-semibold uppercase tracking-[0.08em] text-muted">{label}</div>
      <div className="mt-1 break-words font-mono text-lg font-semibold leading-snug [overflow-wrap:anywhere]">
        {value}
      </div>
    </div>
  );
}

export function EmptyRow({ text }: { text: string }) {
  return <div className="px-3 py-3 text-sm text-muted">{text}</div>;
}

export function formatPercent(value: number): string {
  return Number.isFinite(value) ? `${(value * 100).toFixed(1)}%` : 'не оценивается';
}
