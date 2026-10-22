import type { ReactNode } from 'react';
import { geneColor, variantShapePath } from './trackStyle';

/** Legend for the circular map and the linear tracks; every color is named. */
export function GenomeLegend({ compact = false }: { compact?: boolean }) {
  return (
    <ul className={'flex flex-wrap items-center gap-x-4 gap-y-1.5 text-[11px] text-muted' + (compact ? '' : ' sm:text-xs')} aria-label="Map legend">
      <LegendItem label="Protein-coding"><Swatch color={geneColor('protein_coding')} /></LegendItem>
      <LegendItem label="tRNA"><Swatch color={geneColor('tRNA')} /></LegendItem>
      <LegendItem label="rRNA"><Swatch color={geneColor('rRNA')} /></LegendItem>
      <LegendItem label="Coverage"><Swatch color="rgb(var(--color-coverage))" /></LegendItem>
      {!compact && <LegendItem label="No coverage"><span className="inline-block h-[3px] w-3 bg-coral" /></LegendItem>}
      <LegendItem label={compact ? 'Variant' : 'SNV · insertion · deletion'}>
        <svg width={compact ? 12 : 40} height={12} aria-hidden>
          <path d={variantShapePath('SNV', 4)} transform="translate(6 6)" fill="rgb(var(--color-variant))" />
          {!compact && <path d={variantShapePath('SMALL_INSERTION', 4)} transform="translate(20 6.5)" fill="rgb(var(--color-variant))" />}
          {!compact && <path d={variantShapePath('SMALL_DELETION', 4)} transform="translate(34 6)" fill="rgb(var(--color-variant))" />}
        </svg>
      </LegendItem>
      {!compact && (
        <LegendItem label="No callable denominator">
          <svg width={12} height={12} aria-hidden><path d={variantShapePath('SNV', 4)} transform="translate(6 6)" fill="rgb(var(--color-panel))" stroke="rgb(var(--color-variant))" strokeWidth={1.5} /></svg>
        </LegendItem>
      )}
      <LegendItem label="Structural event"><span className="inline-block h-1.5 w-4 rounded-full bg-sv" /></LegendItem>
    </ul>
  );
}

function LegendItem({ label, children }: { label: string; children: ReactNode }) {
  return (
    <li className="inline-flex items-center gap-1.5">
      {children}
      <span>{label}</span>
    </li>
  );
}

function Swatch({ color }: { color: string }) {
  return <span className="inline-block h-2.5 w-2.5 rounded-sm" style={{ background: color }} />;
}
