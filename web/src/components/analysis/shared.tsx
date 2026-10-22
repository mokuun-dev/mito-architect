

export function Metric({ label, value }: { label: string; value: string }) {
  return (
    <div className="min-w-0 rounded-lg border border-line bg-panel px-4 py-3">
      <div className="text-xs font-medium text-muted">{label}</div>
      <div className="mt-1 break-words text-lg font-semibold leading-snug tracking-tight [overflow-wrap:anywhere]">
        {value}
      </div>
    </div>
  );
}

export function EmptyRow({ text }: { text: string }) {
  return <div className="px-4 py-6 text-center text-[13px] text-muted">{text}</div>;
}

export function formatPercent(value: number): string {
  return Number.isFinite(value) ? `${(value * 100).toFixed(1)}%` : 'not estimable';
}
