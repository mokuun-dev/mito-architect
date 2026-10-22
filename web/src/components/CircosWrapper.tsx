import { MitoCircos } from '@mito-architect/visualization-lib';
import type { GeneAnnotation, LayerName, MitoAnalysisData, MitoCircosSelection } from '@mito-architect/visualization-lib';
import { Download, Minus, Plus, RotateCcw } from 'lucide-react';
import { useEffect, useRef } from 'react';

interface CircosWrapperProps {
  data: MitoAnalysisData;
  activeLayers: Record<LayerName, boolean>;
  selection: MitoCircosSelection;
  compact?: boolean;
  onSvSelect: (svId: string) => void;
  onEventSelect: (eventId: string) => void;
  onGeneSelect: (gene: GeneAnnotation) => void;
}

export default function CircosWrapper({
  data,
  activeLayers,
  selection,
  compact = false,
  onSvSelect,
  onEventSelect,
  onGeneSelect
}: CircosWrapperProps) {
  const hostRef = useRef<HTMLDivElement>(null);
  const instanceRef = useRef<MitoCircos>();
  const { eventId, svId, gene } = selection;

  useEffect(() => {
    const host = hostRef.current;
    if (!host) return;
    const instance = new MitoCircos(host, data, { activeLayers, showControls: false, compact }).render();
    instanceRef.current = instance;

    const handleSv = (event: Event) => onSvSelect((event as CustomEvent<{ id: string }>).detail.id);
    const handleEvent = (event: Event) => onEventSelect((event as CustomEvent<{ id: string }>).detail.id);
    const handleGene = (event: Event) => onGeneSelect((event as CustomEvent<GeneAnnotation>).detail);
    const resizeObserver = new ResizeObserver(() => instance.resize());
    resizeObserver.observe(host);
    host.addEventListener('mito:sv-select', handleSv);
    host.addEventListener('mito:event-select', handleEvent);
    host.addEventListener('mito:gene-select', handleGene);
    return () => {
      host.removeEventListener('mito:sv-select', handleSv);
      host.removeEventListener('mito:event-select', handleEvent);
      host.removeEventListener('mito:gene-select', handleGene);
      resizeObserver.disconnect();
      instance.destroy();
      instanceRef.current = undefined;
    };
    // Layers and selection are applied incrementally by the effects below.
  }, [data, compact, onEventSelect, onGeneSelect, onSvSelect]);

  useEffect(() => {
    instanceRef.current?.setLayers(activeLayers);
  }, [activeLayers]);

  useEffect(() => {
    instanceRef.current?.setSelection({ eventId, svId, gene });
  }, [data, compact, eventId, svId, gene]);

  function downloadSvg() {
    const svg = instanceRef.current?.exportSVG();
    if (!svg) return;
    const url = URL.createObjectURL(new Blob([svg], { type: 'image/svg+xml;charset=utf-8' }));
    const link = document.createElement('a');
    link.href = url;
    link.download = 'mito-' + data.metadata.sample.replace(/[^a-zA-Z0-9_-]/g, '_') + '-map.svg';
    link.click();
    window.setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  return (
    <div className="relative min-w-0" role="group" aria-label="Circular map of the mitochondrial genome">
      <div ref={hostRef} className="w-full" />
      {!compact && (
        <div className="absolute right-0 top-0 flex items-center gap-1" role="group" aria-label="Map zoom and export">
          <button type="button" onClick={() => instanceRef.current?.zoomBy(1.25)} aria-label="Zoom in" title="Zoom in" className="btn-icon h-8 w-8"><Plus className="h-4 w-4" aria-hidden /></button>
          <button type="button" onClick={() => instanceRef.current?.zoomBy(0.8)} aria-label="Zoom out" title="Zoom out" className="btn-icon h-8 w-8"><Minus className="h-4 w-4" aria-hidden /></button>
          <button type="button" onClick={() => instanceRef.current?.resetZoom()} aria-label="Reset zoom" title="Reset zoom" className="btn-icon h-8 w-8"><RotateCcw className="h-4 w-4" aria-hidden /></button>
          <button type="button" onClick={downloadSvg} aria-label="Download map as SVG" title="Download SVG" className="btn-icon h-8 w-8"><Download className="h-4 w-4" aria-hidden /></button>
        </div>
      )}
    </div>
  );
}
