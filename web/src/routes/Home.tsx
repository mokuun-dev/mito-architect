import {
  Activity,
  ArrowRight,
  DatabaseZap,
  FileStack,
  GitBranch,
  Layers3,
  ShieldCheck
} from 'lucide-react';
import { useNavigate } from 'react-router-dom';
import UploadZone from '../components/UploadZone';
import { demoData } from '../lib/demoData';
import { useMitoStore } from '../lib/store';

export default function Home() {
  const navigate = useNavigate();
  const setSelectedFile = useMitoStore((state) => state.setSelectedFile);
  const setSelectedEvidenceGraph = useMitoStore((state) => state.setSelectedEvidenceGraph);
  const setData = useMitoStore((state) => state.setData);
  const setJobId = useMitoStore((state) => state.setJobId);

  return (
    <main id="main-content" className="page-shell surface-grid">
      <section className="grid gap-5 xl:grid-cols-[minmax(0,1fr)_380px]">
        <div className="grid content-start gap-5">
          <section className="lab-panel overflow-hidden p-5">
            <div className="section-kicker">Long-read mitochondrial genomics</div>
            <h1 className="mt-2 max-w-4xl text-3xl font-semibold tracking-tight sm:text-4xl">
              Исследуйте митохондриальный геном от варианта до молекулы.
            </h1>
            <p className="mt-3 max-w-3xl text-sm leading-6 text-muted sm:text-base sm:leading-7">
              В базовом режиме изучайте покрытие и наблюдаемые события. Расширенный режим
              связывает события с молекулами и строит кандидаты архитектур с явными ограничениями.
            </p>
            <div className="mt-3 flex flex-wrap gap-2" aria-label="Core analysis projections">
              <span className="status-pill"><Activity className="h-3.5 w-3.5 text-aqua" aria-hidden />Variant-centric</span>
              <span className="status-pill"><GitBranch className="h-3.5 w-3.5 text-sky" aria-hidden />Molecule-centric</span>
              <span className="status-pill"><Layers3 className="h-3.5 w-3.5 text-magenta" aria-hidden />Architecture-centric</span>
            </div>
          </section>

          <UploadZone
            onFile={(file, evidenceGraph) => {
              setSelectedEvidenceGraph(evidenceGraph);
              setSelectedFile(file);
              navigate('/upload');
            }}
          />
        </div>

        <aside className="grid content-start gap-4" aria-label="Analysis overview">
          <button
            type="button"
            onClick={() => {
              setData(demoData);
              setJobId('demo');
              setSelectedFile(undefined);
              navigate('/result/demo');
            }}
            className="glass-panel group rounded-xl border border-line p-4 text-left shadow-tool hover:border-aqua hover:shadow-focus"
          >
            <span className="flex items-center justify-between gap-3">
              <span>
                <span className="section-kicker">Annotated sandbox</span>
                <span className="mt-1 block text-base font-semibold">Open the traceable demo</span>
                <span className="mt-1 block text-xs leading-5 text-muted">Explore rCRS variants, SV evidence, phasing, clinical assertions, and protein mapping without uploading data.</span>
              </span>
              <span className="grid h-11 w-11 shrink-0 place-items-center rounded-lg border border-aqua/40 bg-aqua/10 text-aqua transition-transform group-hover:translate-x-0.5">
                <ArrowRight className="h-5 w-5" aria-hidden />
              </span>
            </span>
          </button>

          <div className="glass-panel rounded-xl border border-line p-4 shadow-tool">
            <div className="section-kicker mb-3">Evidence model</div>
            <div className="grid gap-2" aria-label="Evidence processing stages">
              <EvidenceStep index="01" label="Alignment evidence" detail="Primary + supplementary records" />
              <EvidenceStep index="02" label="Callable observations" detail="ALT, REF, absent, ambiguous" />
              <EvidenceStep index="03" label="Linked projections" detail="Variants, molecules, architectures" />
            </div>
          </div>

          <div className="glass-panel rounded-xl border border-line p-4 shadow-tool">
            <div className="section-kicker mb-3">Analysis capabilities</div>
            <div className="grid gap-2">
              <SignalRow icon={<FileStack className="h-4 w-4" />} label="FASTQ / SAM / BAM / CRAM" status="native" />
              <SignalRow icon={<ShieldCheck className="h-4 w-4" />} label="NUMT assessability" status="explicit" />
              <SignalRow icon={<Activity className="h-4 w-4" />} label="SNV + indel + complex SV" status="linked" />
              <SignalRow icon={<DatabaseZap className="h-4 w-4" />} label="Versioned annotations" status="provenance" />
            </div>
          </div>
        </aside>
      </section>
    </main>
  );
}

function SignalRow({ icon, label, status }: { icon: JSX.Element; label: string; status: string }) {
  return (
    <div className="flex min-h-11 items-center justify-between gap-3 rounded-lg border border-line bg-panel2 px-3 py-2">
      <span className="flex items-center gap-2 text-sm">
        <span className="text-aqua">{icon}</span>
        {label}
      </span>
      <span className="font-mono text-[11px] text-muted">{status}</span>
    </div>
  );
}

function EvidenceStep({ index, label, detail }: { index: string; label: string; detail: string }) {
  return (
    <div className="grid grid-cols-[36px_minmax(0,1fr)] gap-3 rounded-lg border border-line bg-panel2 p-3">
      <span className="grid h-9 w-9 place-items-center rounded-md border border-aqua/30 bg-aqua/10 font-mono text-xs font-semibold text-aqua">
        {index}
      </span>
      <span className="min-w-0">
        <span className="block text-sm font-semibold">{label}</span>
        <span className="mt-0.5 block text-xs leading-5 text-muted">{detail}</span>
      </span>
    </div>
  );
}
