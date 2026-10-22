import { FileUp, Upload } from 'lucide-react';
import { DragEvent, useId, useRef, useState } from 'react';

const ACCEPTED = ['.bam', '.cram', '.sam', '.fastq', '.fq'];

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
      setError('Expected a BAM, CRAM, SAM, or FASTQ file.');
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
      className="card p-5 sm:p-6"
    >
      <h2 id={headingId} className="text-lg font-semibold tracking-tight">New analysis</h2>
      <p id={helpId} className="mt-1 text-[13px] leading-5 text-muted">
        Coverage, variants, and rearrangements need reference coordinates, so upload an aligned SAM, BAM, or CRAM file.
        FASTQ reads are accepted but carry no coordinates.
      </p>

      <fieldset className="mt-5">
        <legend className="text-[13px] font-semibold">Analysis mode</legend>
        <div className="mt-2 grid gap-3 sm:grid-cols-2">
          <ModeOption
            checked={!evidenceGraph}
            onSelect={() => setEvidenceGraph(false)}
            title="Standard"
            schema="Schema 0.5"
            detail="Coverage, observed variants, rearrangements, and haplogroups."
          />
          <ModeOption
            checked={evidenceGraph}
            onSelect={() => setEvidenceGraph(true)}
            title="Evidence graph"
            schema="Schema 0.6"
            detail="Adds molecule links, phasing, and candidate architectures. Noisy long reads may hit the resource limit."
          />
        </div>
      </fieldset>

      <div
        onDragOver={(event) => {
          event.preventDefault();
          setDragging(true);
        }}
        onDragLeave={() => setDragging(false)}
        onDrop={onDrop}
        className={[
          'mt-5 grid place-items-center rounded-xl border-2 border-dashed px-6 py-10 text-center',
          dragging ? 'border-accent bg-accent/5' : 'border-line bg-panel2/50'
        ].join(' ')}
      >
        <span className="grid h-11 w-11 place-items-center rounded-full bg-accent/10 text-accent">
          <FileUp className="h-5 w-5" aria-hidden />
        </span>
        <p className="mt-3 text-sm font-medium">Drop a file here, or</p>
        <button type="button" onClick={() => inputRef.current?.click()} className="btn-primary mt-3">
          <Upload className="h-4 w-4" aria-hidden />
          Choose file
        </button>
        <p className="mt-3 font-mono text-[11px] uppercase tracking-wide text-subtle" aria-label="Accepted file extensions">
          {ACCEPTED.map((extension) => extension.slice(1)).join(' · ')}
        </p>
      </div>

      {error && (
        <p id={errorId} role="alert" className="notice mt-4 border-coral/30 bg-coral/10 text-coral">
          {error}
        </p>
      )}
      <input
        ref={inputRef}
        type="file"
        className="hidden"
        accept={ACCEPTED.join(',')}
        onChange={(event) => handleFile(event.target.files?.item(0) ?? undefined)}
      />
    </section>
  );
}

function ModeOption({ checked, onSelect, title, schema, detail }: { checked: boolean; onSelect: () => void; title: string; schema: string; detail: string }) {
  return (
    <label
      className={[
        'flex cursor-pointer gap-3 rounded-lg border p-3.5 text-sm',
        checked ? 'border-accent bg-accent/5 shadow-focus' : 'border-line hover:border-line-strong'
      ].join(' ')}
    >
      <input type="radio" name="analysis-schema" checked={checked} onChange={onSelect} className="mt-0.5 h-4 w-4 accent-[rgb(var(--color-accent))]" />
      <span>
        <span className="flex items-center gap-2 font-semibold">
          {title}
          <span className="badge-neutral font-mono">{schema}</span>
        </span>
        <span className="mt-1 block text-[13px] leading-5 text-muted">{detail}</span>
      </span>
    </label>
  );
}
