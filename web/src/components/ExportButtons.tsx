import type { MitoAnalysisData } from '@mito-architect/visualization-lib';
import { ChevronDown, Download, FileCode2, FileJson, Table } from 'lucide-react';
import type { ReactNode } from 'react';
import { downloadReportUrl, downloadResultJsonUrl } from '../lib/api';

interface ExportButtonsProps {
  jobId: string;
  jsonData: MitoAnalysisData;
  htmlAvailable?: boolean;
}

export default function ExportButtons({ jobId, jsonData, htmlAvailable = true }: ExportButtonsProps) {
  function downloadJson() {
    downloadBlob(JSON.stringify(jsonData, null, 2), 'application/json', `mito-${jobId}.json`);
  }

  function downloadVariantTsv() {
    const header = [
      'event_id', 'type', 'position', 'ref', 'alt', 'alt_molecules', 'ref_molecules',
      'other_molecules', 'callable_molecules', 'heteroplasmy', 'ci95_low', 'ci95_high',
      'numt_assessability', 'qc_flags', 'supporting_molecule_ids'
    ];
    const rows = (jsonData.variants ?? []).map((variant) => [
      variant.event_id ?? `legacy:${variant.position}:${variant.ref}:${variant.alt}`,
      variant.type ?? 'SNV',
      variant.position,
      variant.ref,
      variant.alt,
      variant.alt_depth,
      variant.ref_depth,
      variant.other_depth,
      variant.callable_depth,
      variant.heteroplasmy.toFixed(8),
      variant.ci95_low.toFixed(8),
      variant.ci95_high.toFixed(8),
      variant.numt_assessability ?? '',
      (variant.qc_flags ?? []).join(';'),
      (variant.supporting_molecule_ids ?? variant.supporting_reads).join(';')
    ].map(tsvField).join('\t'));
    downloadBlob(
      `${header.join('\t')}\n${rows.join('\n')}\n`,
      'text/tab-separated-values;charset=utf-8',
      `mito-${jobId}-variants.tsv`
    );
  }

  return (
    <details className="relative">
      <summary className="btn-primary list-none">
        <Download className="h-4 w-4" aria-hidden />
        Export
        <ChevronDown className="h-4 w-4 opacity-80" aria-hidden />
      </summary>
      <div className="popover right-0 top-[calc(100%+6px)] w-64 py-1" role="group" aria-label="Export analysis result">
        {htmlAvailable ? (
          <MenuLink href={downloadResultJsonUrl(jobId)} download={`mito-${jobId}.json`} icon={<FileJson className="h-4 w-4" />} detail="Complete result document">Full JSON</MenuLink>
        ) : (
          <MenuButton onClick={downloadJson} icon={<FileJson className="h-4 w-4" />} detail="Complete result document">JSON</MenuButton>
        )}
        <MenuButton
          onClick={downloadVariantTsv}
          icon={<Table className="h-4 w-4" />}
          detail={jsonData.variants ? 'Variants with denominators' : 'Requires aggregate variants'}
          disabled={!jsonData.variants}
        >
          Variant TSV
        </MenuButton>
        {htmlAvailable && (
          <MenuLink href={downloadReportUrl(jobId)} icon={<FileCode2 className="h-4 w-4" />} detail="Self-contained report">HTML report</MenuLink>
        )}
      </div>
    </details>
  );
}

const ITEM = 'flex w-full items-start gap-3 px-3 py-2 text-left hover:bg-panel2 disabled:opacity-50 disabled:hover:bg-transparent';

function MenuLink({ href, download, icon, detail, children }: { href: string; download?: string; icon: ReactNode; detail: string; children: string }) {
  return (
    <a href={href} download={download} className={ITEM}>
      <span className="mt-0.5 text-muted" aria-hidden>{icon}</span>
      <span><span className="block text-sm font-medium">{children}</span><span className="block text-xs text-muted">{detail}</span></span>
    </a>
  );
}

function MenuButton({ onClick, icon, detail, children, disabled }: { onClick: () => void; icon: ReactNode; detail: string; children: string; disabled?: boolean }) {
  return (
    <button type="button" onClick={onClick} disabled={disabled} className={ITEM}>
      <span className="mt-0.5 text-muted" aria-hidden>{icon}</span>
      <span><span className="block text-sm font-medium">{children}</span><span className="block text-xs text-muted">{detail}</span></span>
    </button>
  );
}

function tsvField(value: string | number): string {
  return String(value).replace(/[\t\r\n]+/g, ' ').trim();
}

function downloadBlob(contents: string, type: string, fileName: string) {
  const blob = new Blob([contents], { type });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = fileName;
  link.click();
  URL.revokeObjectURL(url);
}
