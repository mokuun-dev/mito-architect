import { Check } from 'lucide-react';
import type { JobStatusValue } from '../lib/api';

interface ProgressBarProps {
  status?: JobStatusValue;
  progress?: number;
  error?: string | null;
}

const STAGES = [
  { id: 'queued', label: 'Queued' },
  { id: 'processing', label: 'Analysing' },
  { id: 'done', label: 'Result ready' }
] as const;

export default function ProgressBar({ status = 'queued', progress = 0, error }: ProgressBarProps) {
  const currentStage = status === 'done' ? 2 : status === 'processing' ? 1 : status === 'queued' ? 0 : -1;
  const boundedProgress = Math.max(0, Math.min(100, progress));
  const failed = status === 'error' || status === 'cancelled' || Boolean(error);
  return (
    <section className="card p-5 sm:p-6" aria-labelledby="analysis-progress-title" aria-live="polite">
      <div className="flex items-start justify-between gap-4">
        <div>
          <h2 id="analysis-progress-title" className="text-[15px] font-semibold">Analysis progress</h2>
          <p role={error ? 'alert' : 'status'} className={error ? 'mt-1 text-sm text-coral' : 'mt-1 text-sm text-muted'}>
            {error ?? statusLabel(status)}
          </p>
        </div>
        <div className="font-mono text-2xl font-semibold tabular-nums">{Math.round(boundedProgress)}%</div>
      </div>
      <div
        className="mt-5 h-1.5 overflow-hidden rounded-full bg-panel2"
        role="progressbar"
        aria-label="Analysis completion"
        aria-valuemin={0}
        aria-valuemax={100}
        aria-valuenow={Math.round(boundedProgress)}
      >
        <div
          className={'h-full rounded-full transition-all duration-300 ' + (failed ? 'bg-coral' : 'bg-accent')}
          style={{ width: `${boundedProgress}%` }}
        />
      </div>
      <ol className="mt-5 grid grid-cols-3 gap-3">
        {STAGES.map((stage, index) => {
          const complete = index < currentStage || status === 'done';
          const active = index === currentStage && status !== 'done';
          return (
            <li key={stage.id} className="flex items-center gap-2 text-[13px]">
              <span
                className={[
                  'grid h-5 w-5 shrink-0 place-items-center rounded-full border text-[10px] font-semibold',
                  complete ? 'border-accent bg-accent text-onaccent' : active ? 'border-accent text-accent' : 'border-line text-subtle'
                ].join(' ')}
              >
                {complete ? <Check className="h-3 w-3" aria-hidden /> : index + 1}
              </span>
              <span className={complete || active ? 'font-medium text-text' : 'text-muted'}>{stage.label}</span>
            </li>
          );
        })}
      </ol>
    </section>
  );
}

function statusLabel(status: JobStatusValue): string {
  switch (status) {
    case 'queued':
      return 'Waiting for an analysis slot.';
    case 'processing':
      return 'Reading alignments, calling variants and structural events, and assembling molecules.';
    case 'done':
      return 'The result is ready.';
    case 'error':
      return 'The analysis failed.';
    case 'cancelled':
      return 'The analysis was cancelled.';
  }
}
