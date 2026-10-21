import { useMitoStore } from '../../lib/store';

export function GenomeLegend({ maxDepth }: { maxDepth: number }) {
  const layers = useMitoStore((state) => state.activeLayers);
  return <div className="workspace-genome-legend" aria-label="Условные обозначения карты">
    {layers.genes && <><span><i className="legend-swatch bg-sky" />Белок-кодирующие</span><span><i className="legend-swatch bg-[#a78bfa]" />tRNA</span><span><i className="legend-swatch bg-aqua" />rRNA</span></>}
    {layers.coverage && <span><i className="legend-swatch bg-aqua/40" />Глубина 0–{maxDepth}×</span>}
    {layers.snps && <span><i className="legend-dot bg-amber" />Малый вариант</span>}
    {layers.svs && <span><i className="legend-swatch bg-magenta" />Перестройка</span>}
    <span className="text-muted">Кольцо показывает референс, а не форму отдельной молекулы</span>
  </div>;
}
