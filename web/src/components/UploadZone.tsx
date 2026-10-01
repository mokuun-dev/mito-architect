import { FileCheck2, FileUp, HardDriveUpload, LockKeyhole, Upload } from 'lucide-react';
import { DragEvent, useId, useRef, useState } from 'react';

const ACCEPTED = ['.fastq', '.fq', '.sam', '.bam', '.cram'];

interface UploadZoneProps {
  onFile: (file: File, evidenceGraph: boolean) => void;
  defaultEvidenceGraph?: boolean;
}

export default function UploadZone({ onFile, defaultEvidenceGraph = false }: UploadZoneProps) {
  const inputRef = useRef<HTMLInputElement>(null);
  const headingId = useId();
  const helpId = useId();
  const errorId = useId();
  const [dragging, setDragging] = useState(false);
  const [error, setError] = useState<string>();
  const [evidenceGraph, setEvidenceGraph] = useState(defaultEvidenceGraph);

  function handleFile(file?: File) {
    if (!file) return;
    const lower = file.name.toLowerCase();
    if (!ACCEPTED.some((extension) => lower.endsWith(extension))) {
      setError('Expected FASTQ, SAM, BAM, or CRAM input.');
      return;
    }
    setError(undefined);
    onFile(file, evidenceGraph);
  }

  function onDrop(event: DragEvent<HTMLElement>) {
    event.preventDefault();
    setDragging(false);
    handleFile(event.dataTransfer.files.item(0) ?? undefined);
  }

  return (
    <section
      aria-labelledby={headingId}
      aria-describedby={[helpId, error ? errorId : undefined].filter(Boolean).join(' ')}
      onDragOver={(event) => {
        event.preventDefault();
        setDragging(true);
      }}
      onDragLeave={() => setDragging(false)}
      onDrop={onDrop}
      className={[
        'glass-panel relative grid place-items-center overflow-hidden rounded-xl border border-dashed px-5 py-6 text-center shadow-tool sm:px-8 sm:py-7',
        dragging ? 'border-aqua bg-aqua/10 shadow-focus' : 'border-line hover:border-aqua/50'
      ].join(' ')}
    >
      <div className="pointer-events-none absolute inset-x-[12%] top-0 h-px bg-gradient-to-r from-transparent via-aqua/70 to-transparent" />
      <div className="relative max-w-xl">
        <div className="mx-auto mb-3 grid h-12 w-12 place-items-center rounded-xl border border-aqua/40 bg-aqua/10 shadow-focus">
          {dragging ? <HardDriveUpload className="h-8 w-8 text-aqua" aria-hidden /> : <FileUp className="h-8 w-8 text-aqua" aria-hidden />}
        </div>
        <div className="section-kicker">New analysis</div>
        <h2 id={headingId} className="mt-1 text-2xl font-semibold tracking-tight">Load long-read evidence</h2>
        <p id={helpId} className="mt-2 text-sm leading-6 text-muted">
          Выберите один файл выравнивания или FASTQ. Для покрытия, вариантов и перестроек
          нужны координаты из SAM, BAM или CRAM; FASTQ сам по себе их не содержит.
        </p>
        <fieldset className="mt-4 grid gap-2 text-left">
          <legend className="mb-2 text-sm font-semibold">Режим анализа</legend>
          <label className="flex cursor-pointer gap-3 rounded-lg border border-line bg-panel px-3 py-3 text-sm">
            <input type="radio" name="analysis-schema" checked={!evidenceGraph} onChange={() => setEvidenceGraph(false)} className="mt-1 accent-aqua" />
            <span><strong className="block">Базовый · схема 0.5</strong><span className="mt-1 block text-xs leading-5 text-muted">Покрытие, наблюдаемые варианты и перестройки. Группы по признакам не являются подтверждёнными молекулярными архитектурами.</span></span>
          </label>
          <label className="flex cursor-pointer gap-3 rounded-lg border border-line bg-panel px-3 py-3 text-sm">
            <input type="radio" name="analysis-schema" checked={evidenceGraph} onChange={() => setEvidenceGraph(true)} className="mt-1 accent-aqua" />
            <span><strong className="block">Расширенный · схема 0.6</strong><span className="mt-1 block text-xs leading-5 text-muted">Граф доказательств, связи событий и кандидаты архитектур. На шумных длинных чтениях возможен отказ по лимиту ресурсов.</span></span>
          </label>
        </fieldset>
        <div className="mt-4 flex flex-wrap justify-center gap-3">
          <button
            type="button"
            onClick={() => inputRef.current?.click()}
            className="inline-flex min-h-11 items-center gap-2 rounded-lg border border-aqua bg-aqua px-5 py-3 text-sm font-semibold text-shell shadow-focus hover:bg-aqua/90"
          >
            <Upload className="h-4 w-4" aria-hidden />
            Choose file
          </button>
        </div>
        <div className="mt-3 flex flex-wrap justify-center gap-2" aria-label="Accepted file extensions">
          {ACCEPTED.map((extension) => (
            <span key={extension} className="status-pill font-mono uppercase">
              {extension.slice(1)}
            </span>
          ))}
        </div>
        <div className="mt-3 grid gap-2 text-left text-xs text-muted sm:grid-cols-2">
          <span className="flex items-center gap-2"><FileCheck2 className="h-4 w-4 shrink-0 text-leaf" aria-hidden />Content and extension validation</span>
          <span className="flex items-center gap-2"><LockKeyhole className="h-4 w-4 shrink-0 text-sky" aria-hidden />Per-job isolated processing</span>
        </div>
        {error && <p id={errorId} role="alert" className="mt-4 rounded-lg border border-coral/40 bg-coral/10 px-3 py-2 text-sm text-coral">{error}</p>}
        <input
          ref={inputRef}
          type="file"
          className="hidden"
          accept={ACCEPTED.join(',')}
          onChange={(event) => handleFile(event.target.files?.item(0) ?? undefined)}
        />
      </div>
    </section>
  );
}
