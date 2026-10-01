import { useMutation, useQuery } from '@tanstack/react-query';
import { ArrowLeft, FileUp } from 'lucide-react';
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
      <main id="main-content" className="page-shell grid max-w-[1100px] gap-5">
        <div>
          <Link to="/" className="inline-flex min-h-11 items-center gap-2 rounded-md text-sm text-muted hover:text-aqua">
            <ArrowLeft className="h-4 w-4" aria-hidden />
            Back to workspace
          </Link>
          <h1 className="mt-3 text-3xl font-semibold tracking-tight">Start a new analysis</h1>
          <p className="mt-2 text-sm text-muted">Select one supported long-read input for this job.</p>
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

  return (
    <main id="main-content" className="page-shell grid max-w-[1100px] gap-5">
      <Link to="/" className="inline-flex min-h-11 w-fit items-center gap-2 rounded-md text-sm text-muted hover:text-aqua">
        <ArrowLeft className="h-4 w-4" aria-hidden />
        New analysis
      </Link>
      <section className="lab-panel p-5 sm:p-6">
        <div className="flex flex-col gap-4 sm:flex-row sm:items-center sm:justify-between">
          <div className="min-w-0">
            <div className="section-kicker">Active job</div>
            <h1 className="mt-2 truncate text-2xl font-semibold tracking-tight">{selectedFile?.name ?? `Job ${jobId}`}</h1>
            <p className="mt-2 text-sm text-muted">
              Режим {upload.data?.result_schema ?? (selectedEvidenceGraph ? '0.6' : '0.5')} · проверка входа и воспроизводимый результат.
            </p>
          </div>
          <span className="status-pill shrink-0">
            <FileUp className="h-4 w-4 text-aqua" aria-hidden />
            {selectedFile ? formatBytes(selectedFile.size) : 'uploaded'}
          </span>
        </div>
      </section>
      <ProgressBar
        status={status.data?.status ?? (upload.isError ? 'error' : upload.isPending ? 'queued' : 'processing')}
        progress={status.data?.progress ?? (upload.isPending ? 5 : 0)}
        error={(upload.error as Error | undefined)?.message ?? status.data?.error?.message}
      />
    </main>
  );
}

function formatBytes(bytes: number): string {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KiB`;
  if (bytes < 1024 * 1024 * 1024) return `${(bytes / (1024 * 1024)).toFixed(1)} MiB`;
  return `${(bytes / (1024 * 1024 * 1024)).toFixed(2)} GiB`;
}
