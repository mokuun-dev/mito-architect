import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { AlertTriangle, ChevronDown } from 'lucide-react';
import ExportButtons from '../ExportButtons';

export function SampleHeader({ data, jobId, isDemo }: { data: MitoAnalysisData; jobId: string; isDemo: boolean }) {
  const numt = data.filter_stats.numt_assessment;
  return <header className="workspace-sample-header">
    <div className="min-w-0">
      <div className="flex flex-wrap items-center gap-2 text-xs text-muted">
        <span>Результат анализа</span>
        {isDemo && <span className="rounded border border-amber/40 bg-amber/10 px-2 py-0.5 font-semibold text-amber">DEMO</span>}
      </div>
      <h1 className="mt-1 break-words text-2xl font-semibold leading-8 tracking-tight">{data.metadata.sample}</h1>
      <p className="mt-1 font-mono text-xs text-muted">{data.metadata.reference_accession ?? 'Собственный референс'} · {data.metadata.reference_length.toLocaleString()} bp</p>
    </div>
    <div className="flex flex-wrap items-center gap-2">
      <details className="relative">
        <summary className="workspace-header-action list-none">Паспорт результата <ChevronDown className="h-4 w-4" aria-hidden /></summary>
        <div className="workspace-header-popover">
          <dl className="grid grid-cols-[auto_1fr] gap-x-4 gap-y-2 text-xs">
            <dt className="text-muted">Движок</dt><dd className="font-mono">{data.metadata.engine_version}</dd>
            <dt className="text-muted">Схема</dt><dd className="font-mono">{data.metadata.schema_version ?? 'не указана'}</dd>
            <dt className="text-muted">Вход</dt><dd className="break-all font-mono">{data.metadata.input_path ?? 'не указан'}</dd>
            <dt className="text-muted">Референс</dt><dd className="break-all font-mono">{data.metadata.reference_path ?? 'не указан'}</dd>
          </dl>
          <p className="mt-3 border-t border-line pt-2 text-xs text-muted">Проверка формата результата не подтверждает биологическую достоверность выводов.</p>
        </div>
      </details>
      <ExportButtons jobId={jobId} jsonData={data} htmlAvailable={!isDemo} />
    </div>
    {numt && !numt.specificity_assessable && <div className="workspace-data-notice" role="status">
      <AlertTriangle className="h-4 w-4 shrink-0" aria-hidden />
      <span>Специфичность NUMT не оценена. Низкочастотные варианты требуют осторожной интерпретации.</span>
    </div>}
  </header>;
}
