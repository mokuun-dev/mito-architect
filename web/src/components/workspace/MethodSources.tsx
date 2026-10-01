import type { MitoAnalysisData } from '@mito-architect/visualization-lib';

export function MethodSources({ data }: { data: MitoAnalysisData }) {
  return <section className="workspace-section" aria-labelledby="method-heading">
    <h2 id="method-heading" className="text-lg font-semibold">Метод и источники</h2>
    <p className="mt-2 max-w-[70ch] text-sm text-muted">Версии и файлы позволяют воспроизвести анализ. Наличие ссылки на ресурс не означает независимого подтверждения каждого вывода.</p>
    <dl className="mt-5 grid gap-2 text-sm sm:grid-cols-[180px_1fr]">
      <dt className="text-muted">Движок</dt><dd className="font-mono">{data.metadata.engine_version}</dd>
      <dt className="text-muted">Схема результата</dt><dd className="font-mono">{data.metadata.schema_version ?? 'не указана'}</dd>
      <dt className="text-muted">Референс</dt><dd className="font-mono">{data.metadata.reference_accession ?? data.metadata.reference_path ?? 'не указан'}</dd>
      <dt className="text-muted">Тип входа</dt><dd className="break-all font-mono">{data.metadata.input_path ?? 'не указан'}</dd>
    </dl>
    {(data.metadata.algorithm_notes ?? []).length > 0 && <div className="mt-5 border-t border-line pt-4"><h3 className="text-sm font-semibold">Примечания алгоритма</h3><ul className="mt-2 list-disc space-y-1 pl-5 text-sm text-muted">{data.metadata.algorithm_notes?.map((note) => <li key={note}>{note}</li>)}</ul></div>}
    <div className="mt-5 border-t border-line pt-4"><h3 className="text-sm font-semibold">Версионированные ресурсы</h3>
      {!data.metadata.resources?.length ? <p className="mt-2 text-sm text-muted">Сведения о ресурсах не представлены.</p> : <ul className="mt-3 divide-y divide-line border-y border-line">{data.metadata.resources.map((resource) => <li key={resource.name + ':' + resource.version} className="grid gap-1 py-3 text-sm sm:grid-cols-[minmax(0,1fr)_auto] sm:items-start">
        <div><div className="font-semibold">{resource.name} <span className="font-mono font-normal text-muted">{resource.version}</span></div><div className="mt-1 break-all font-mono text-xs text-muted">SHA-256: {resource.sha256}</div><div className="mt-1 text-xs text-muted">Получен: {resource.retrieved || 'не указано'} · Лицензия: {resource.license || 'не указана'}</div></div>
        {/^https:\/\//.test(resource.source) && <a href={resource.source} target="_blank" rel="noopener noreferrer" className="font-medium text-aqua underline-offset-2 hover:underline">Открыть источник ↗</a>}
      </li>)}</ul>}</div>
  </section>;
}
