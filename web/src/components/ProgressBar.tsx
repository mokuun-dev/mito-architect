import type { JobStatusValue } from '../lib/api';

interface ProgressBarProps {
  status?: JobStatusValue;
  progress?: number;
  error?: string | null;
}

export default function ProgressBar({ status = 'queued', progress = 0, error }: ProgressBarProps) {
  const stages = ['queued', 'processing', 'done'] as const;
  const currentStage = status === 'done' ? 2 : status === 'processing' ? 1 : status === 'queued' ? 0 : -1;
  const boundedProgress = Math.max(0, Math.min(100, progress));
  return (
    <section className="lab-panel p-5 sm:p-6" aria-labelledby="analysis-progress-title" aria-live="polite">
      <div className="flex items-start justify-between gap-4">
        <div>
          <div className="section-kicker">Native pipeline</div>
          <h2 id="analysis-progress-title" className="mt-1 text-lg font-semibold">Analysis progress</h2>
          <p role={error ? 'alert' : 'status'} className={error ? 'mt-1 text-sm text-coral' : 'mt-1 text-sm text-muted'}>
            {error ?? statusLabel(status)}
          </p>
        </div>
        <div className="font-mono text-2xl font-semibold tabular-nums">{Math.round(boundedProgress)}%</div>
      </div>
      <div
        className="mt-5 h-2 overflow-hidden rounded-full bg-panel2"
        role="progressbar"
        aria-label="Analysis completion"
        aria-valuemin={0}
        aria-valuemax={100}
        aria-valuenow={Math.round(boundedProgress)}
      >
        <div
          className="h-full rounded-full bg-aqua transition-all duration-300"
          style={{ width: `${boundedProgress}%` }}
        />
      </div>
      <div className="mt-5 grid grid-cols-3 gap-2">
        {stages.map((stage, index) => (
          <div
            key={stage}
            className={[
              'min-h-11 rounded-lg border px-3 py-2 text-center font-mono text-[11px] font-semibold uppercase tracking-[0.08em]',
              index <= currentStage
                ? 'border-aqua/60 bg-aqua/10 text-aqua'
                : 'border-line text-muted'
            ].join(' ')}
          >
            {stage}
          </div>
        ))}
      </div>
    </section>
  );
}

function statusLabel(status: JobStatusValue): string {
  switch (status) {
    case 'queued':
      return 'Queued for background processing.';
    case 'processing':
      return 'Extracting features, calling SVs, and clustering reads.';
    case 'done':
      return 'Result is ready.';
    case 'error':
      return 'Analysis failed.';
    case 'cancelled':
      return 'Analysis was cancelled.';
  }
}
