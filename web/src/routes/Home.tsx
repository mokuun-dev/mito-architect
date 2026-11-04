import { useQuery } from '@tanstack/react-query';
import { ArrowRight, FileStack, GitBranch, Layers3, ScanLine, ShieldCheck } from 'lucide-react';
import type { ReactNode } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import UploadZone from '../components/UploadZone';
import { fetchDemoIndex, platformLabel } from '../lib/demoPanel';
import { useMitoStore } from '../lib/store';

export default function Home() {
  const navigate = useNavigate();
  const setSelectedFile = useMitoStore((state) => state.setSelectedFile);
  const setSelectedEvidenceGraph = useMitoStore((state) => state.setSelectedEvidenceGraph);
  const panel = useQuery({ queryKey: ['demo-index'], queryFn: ({ signal }) => fetchDemoIndex(signal), staleTime: Infinity });

  return (
    <main id="main-content" className="page-shell">
      <section className="max-w-3xl">
        <p className="section-kicker text-accent">Long-read mitochondrial genomics</p>
        <h1 className="mt-2 text-3xl font-semibold tracking-tight sm:text-[40px] sm:leading-[1.15]">
          From aligned reads to variants, molecules, and rearrangements.
        </h1>
        <p className="mt-3 text-[15px] leading-7 text-muted">
          Upload a long-read alignment of human mtDNA. Mito-Architect reports coverage, heteroplasmic
          variants with callable denominators, structural events, and haplogroups, and links each
          finding back to the molecules that support it.
        </p>
      </section>

      <div className="mt-8 grid items-start gap-6 lg:grid-cols-[minmax(0,1fr)_360px]">
        <UploadZone
          onFile={(file, evidenceGraph) => {
            setSelectedEvidenceGraph(evidenceGraph);
            setSelectedFile(file);
            navigate('/upload');
          }}
        />

        <aside className="grid content-start gap-4" aria-label="Getting started">
          <section className="card p-5" aria-labelledby="how-it-works">
            <h2 id="how-it-works" className="text-[15px] font-semibold">How it works</h2>
            <ol className="mt-4 grid gap-4">
              <Step index={1} title="Read alignment evidence" detail="Primary and supplementary records are assembled into molecules." />
              <Step index={2} title="Classify observations" detail="Each position is ALT, REF, absent, or ambiguous on every molecule." />
              <Step index={3} title="Project the results" detail="Variants, molecules, and architectures stay linked to their evidence." />
            </ol>
          </section>

          <section className="card p-5" aria-labelledby="capabilities">
            <h2 id="capabilities" className="text-[15px] font-semibold">Capabilities</h2>
            <ul className="mt-3 grid gap-2.5 text-[13px]">
              <Capability icon={<FileStack className="h-4 w-4" />}>SAM, BAM, CRAM, and FASTQ input</Capability>
              <Capability icon={<ScanLine className="h-4 w-4" />}>SNVs, small indels, and complex SVs</Capability>
              <Capability icon={<GitBranch className="h-4 w-4" />}>Molecule-level phasing</Capability>
              <Capability icon={<Layers3 className="h-4 w-4" />}>Candidate molecular architectures</Capability>
              <Capability icon={<ShieldCheck className="h-4 w-4" />}>Explicit NUMT assessability</Capability>
            </ul>
          </section>
        </aside>
      </div>

      <section className="mt-10" aria-labelledby="demo-panel-title">
        <div className="max-w-3xl">
          <h2 id="demo-panel-title" className="text-xl font-semibold tracking-tight">Explore example samples</h2>
          <p className="mt-1.5 text-[13px] leading-6 text-muted">
            Actual engine output for mitochondrial genomes of different people from public records, on simulated long reads, plus one real nanopore run.
            Clinical variants and large deletions are added in silico; they are not findings about those people.
          </p>
        </div>
        {panel.isError ? (
          <p className="notice-warning mt-4" role="status">{(panel.error as Error).message}</p>
        ) : (
          <ul className="mt-4 grid gap-3 sm:grid-cols-2 xl:grid-cols-3">
            {panel.data?.samples.map((sample) => (
              <li key={sample.id}>
                <Link to={`/result/demo?sample=${encodeURIComponent(sample.id)}`} className="card group flex h-full flex-col gap-2 p-4 hover:border-accent/50">
                  <span className="flex items-start justify-between gap-3">
                    <span className="text-[14px] font-semibold">{sample.title}</span>
                    <ArrowRight className="mt-0.5 h-4 w-4 shrink-0 text-accent transition-transform group-hover:translate-x-0.5" aria-hidden />
                  </span>
                  <span className="text-[13px] leading-5 text-muted">{sample.summary}</span>
                  <span className="mt-auto flex flex-wrap items-center gap-1.5 pt-1">
                    <span className={sample.simulated ? 'badge-neutral' : 'badge-warning'}>{sample.simulated ? platformLabel(sample.platform) : 'Real reads'}</span>
                    {sample.haplogroup && <span className="badge-neutral font-mono">{sample.haplogroup}</span>}
                  </span>
                </Link>
              </li>
            ))}
          </ul>
        )}
      </section>
    </main>
  );
}

function Step({ index, title, detail }: { index: number; title: string; detail: string }) {
  return (
    <li className="grid grid-cols-[28px_minmax(0,1fr)] gap-3">
      <span className="grid h-7 w-7 place-items-center rounded-full border border-line bg-panel2 font-mono text-xs font-semibold text-muted">{index}</span>
      <span>
        <span className="block text-[13px] font-semibold">{title}</span>
        <span className="mt-0.5 block text-[13px] leading-5 text-muted">{detail}</span>
      </span>
    </li>
  );
}

function Capability({ icon, children }: { icon: ReactNode; children: string }) {
  return (
    <li className="flex items-center gap-2.5">
      <span className="text-accent" aria-hidden>{icon}</span>
      {children}
    </li>
  );
}
