import { MitoCircos } from '@mito-architect/visualization-lib';
import type { LayerName, MitoAnalysisData } from '@mito-architect/visualization-lib';
import { useEffect, useRef } from 'react';
import type { GeneAnnotation } from '@mito-architect/visualization-lib';
import { Download, Minus, Plus, RotateCcw } from 'lucide-react';

interface CircosWrapperProps {
  data: MitoAnalysisData;
  activeLayers: Record<LayerName, boolean>;
  onClusterSelect: (clusterId: number) => void;
  onSvSelect: (svId: string) => void;
  onEventSelect: (eventId: string) => void;
  onGeneSelect: (gene: GeneAnnotation) => void;
}

export default function CircosWrapper({
  data,
  activeLayers,
  onClusterSelect,
  onSvSelect,
  onEventSelect,
  onGeneSelect
}: CircosWrapperProps) {
  const hostRef = useRef<HTMLDivElement>(null);
  const instanceRef = useRef<MitoCircos>();

  useEffect(() => {
    const host = hostRef.current;
    if (!host) return;
    const instance = new MitoCircos(host, data, {
      activeLayers,
      showControls: false
    }).render();
    instanceRef.current = instance;

    const handleCluster = (event: Event) => {
      const detail = (event as CustomEvent<{ id: number }>).detail;
      onClusterSelect(detail.id);
    };
    const handleSv = (event: Event) => {
      const detail = (event as CustomEvent<{ id: string }>).detail;
      onSvSelect(detail.id);
    };
    const handleEvent = (event: Event) => {
      const detail = (event as CustomEvent<{ id: string }>).detail;
      onEventSelect(detail.id);
    };
    const handleGene = (event: Event) => onGeneSelect((event as CustomEvent<GeneAnnotation>).detail);
    const resizeObserver = new ResizeObserver(() => instance.resize());
    resizeObserver.observe(host);
    host.addEventListener('mito:cluster-select', handleCluster);
    host.addEventListener('mito:sv-select', handleSv);
    host.addEventListener('mito:event-select', handleEvent);
    host.addEventListener('mito:gene-select', handleGene);
    return () => {
      host.removeEventListener('mito:cluster-select', handleCluster);
      host.removeEventListener('mito:sv-select', handleSv);
      host.removeEventListener('mito:event-select', handleEvent);
      host.removeEventListener('mito:gene-select', handleGene);
      resizeObserver.disconnect();
      instance.destroy();
    };
  }, [data, onClusterSelect, onEventSelect, onGeneSelect, onSvSelect]);

  useEffect(() => {
    instanceRef.current?.setLayers(activeLayers);
  }, [activeLayers]);

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
    <div className="min-w-0" role="group" aria-label="Круговая карта митохондриального генома">
      <div ref={hostRef} className="min-h-[380px] w-full" />
      <div className="mt-2 flex flex-wrap items-center justify-between gap-2">
        <p className="text-xs text-muted">Начало сверху · координаты по часовой стрелке · Ctrl + колесо для масштаба</p>
        <div className="flex items-center gap-1" role="group" aria-label="Масштаб и экспорт карты">
          <button type="button" onClick={() => instanceRef.current?.zoomBy(1.25)} aria-label="Увеличить карту" className="workspace-icon-button"><Plus className="h-4 w-4" aria-hidden /></button>
          <button type="button" onClick={() => instanceRef.current?.zoomBy(0.8)} aria-label="Уменьшить карту" className="workspace-icon-button"><Minus className="h-4 w-4" aria-hidden /></button>
          <button type="button" onClick={() => instanceRef.current?.resetZoom()} aria-label="Сбросить масштаб карты" className="workspace-icon-button"><RotateCcw className="h-4 w-4" aria-hidden /></button>
          <button type="button" onClick={downloadSvg} aria-label="Скачать карту SVG" className="workspace-icon-button"><Download className="h-4 w-4" aria-hidden /></button>
        </div>
      </div>
    </div>
  );
}
