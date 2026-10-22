import { useMutation, useQuery } from '@tanstack/react-query';
import { ArrowLeft, FileText } from 'lucide-react';
import { useEffect } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import ProgressBar from '../components/ProgressBar';
import UploadZone from '../components/UploadZone';
import { getStatus, uploadFileWithOptions } from '../lib/api';
import { useMitoStore } from '../lib/store';

export default function UploadPage() {
  const navigate = useNavigate();
  const selectedFile = useMitoStore((state) => state.selectedFile);
  const selectedEvidenceGraph = useMitoStore((state) => state.selectedEvidenceGraph);
  const setSelectedFile = useMitoStore((state) => state.setSelectedFile);
  const setSelectedEvidenceGraph = useMitoStore((state) => state.setSelectedEvidenceGraph);
  const jobId = useMitoStore((state) => state.jobId);
  const setJobId = useMitoStore((state) => state.setJobId);

  const upload = useMutation({
    mutationFn: ({ file, evidenceGraph }: { file: File; evidenceGraph: boolean }) => uploadFileWithOptions(file, { evidenceGraph }),
    onSuccess: (response) => setJobId(response.job_id)
  });

  useEffect(() => {
    if (selectedFile && !jobId && !upload.isPending && !upload.isSuccess && !upload.isError) {
      upload.mutate({ file: selectedFile, evidenceGraph: selectedEvidenceGraph });
    }
  }, [jobId, selectedEvidenceGraph, selectedFile, upload]);

  const status = useQuery({
    queryKey: ['status', jobId],
    queryFn: ({ signal }) => getStatus(jobId!, signal),
    enabled: Boolean(jobId),
    refetchInterval: (query) =>
      query.state.data && ['done', 'error', 'cancelled'].includes(query.state.data.status)
        ? false
        : 1000
  });

  useEffect(() => {
    if (status.data?.status === 'done') {
      navigate(`/result/${status.data.job_id}`);
    }
  }, [navigate, status.data]);

  if (!selectedFile && !jobId) {
    return (
      <main id="main-content" className="page-shell grid max-w-[960px] gap-6">
        <div>
          <BackLink>Back</BackLink>
          <h1 className="mt-3 text-2xl font-semibold tracking-tight">Start a new analysis</h1>
          <p className="mt-1 text-sm text-muted">Select one long-read input for this job.</p>
        </div>
        <UploadZone
          defaultEvidenceGraph={selectedEvidenceGraph}
          onFile={(file, evidenceGraph) => {
            setSelectedEvidenceGraph(evidenceGraph);
            setSelectedFile(file);
            setJobId(undefined);
          }}
        />
      </main>
    );
  }

  const schema = upload.data?.result_schema ?? (selectedEvidenceGraph ? '0.6' : '0.5');
  return (
    <main id="main-content" className="page-shell grid max-w-[960px] gap-6">
      <BackLink>New analysis</BackLink>
      <section className="card flex flex-wrap items-center justify-between gap-4 p-5">
        <div className="flex min-w-0 items-center gap-3">
          <span className="grid h-10 w-10 shrink-0 place-items-center rounded-lg bg-panel2 text-muted">
            <FileText className="h-5 w-5" aria-hidden />
          </span>
          <div className="min-w-0">
            <h1 className="truncate text-lg font-semibold tracking-tight">{selectedFile?.name ?? `Job ${jobId}`}</h1>
            <p className="mt-0.5 text-[13px] text-muted">
              {selectedFile ? formatBytes(selectedFile.size) + ' · ' : ''}
              {schema === '0.6' ? 'Evidence graph' : 'Standard'} analysis · schema {schema}
            </p>
          </div>
        </div>
        {jobId && <span className="badge-neutral font-mono">{jobId}</span>}
      </section>
      <ProgressBar
        status={status.data?.status ?? (upload.isError ? 'error' : upload.isPending ? 'queued' : 'processing')}
        progress={status.data?.progress ?? (upload.isPending ? 5 : 0)}
        error={(upload.error as Error | undefined)?.message ?? status.data?.error?.message}
      />
    </main>
  );
}

function BackLink({ children }: { children: string }) {
  return (
    <Link to="/" className="inline-flex w-fit items-center gap-1.5 rounded-lg text-sm font-medium text-muted hover:text-text">
      <ArrowLeft className="h-4 w-4" aria-hidden />
      {children}
    </Link>
  );
}

function formatBytes(bytes: number): string {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KiB`;
  if (bytes < 1024 * 1024 * 1024) return `${(bytes / (1024 * 1024)).toFixed(1)} MiB`;
  return `${(bytes / (1024 * 1024 * 1024)).toFixed(2)} GiB`;
}
