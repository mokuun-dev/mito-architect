import type { GeneAnnotation, LayerName } from '@mito-architect/visualization-lib';
import { Layers3, Search } from 'lucide-react';
import { useMemo, useState } from 'react';
import { useMitoStore } from '../../lib/store';

const LAYERS: Array<[LayerName, string]> = [
  ['genes', 'Гены'], ['coverage', 'Покрытие'], ['snps', 'Варианты'], ['svs', 'Перестройки']
];

export function GenomeToolbar({ genes }: { genes: GeneAnnotation[] }) {
  const [query, setQuery] = useState('');
  const activeLayers = useMitoStore((state) => state.activeLayers);
  const setLayer = useMitoStore((state) => state.setLayer);
  const setSelectedGene = useMitoStore((state) => state.setSelectedGene);
  const setSelectedRegion = useMitoStore((state) => state.setSelectedRegion);
  const found = useMemo(() => query.trim() ? genes.filter((gene) => gene.name.toLowerCase().includes(query.trim().toLowerCase())).slice(0, 8) : [], [query, genes]);

  return <div className="workspace-genome-toolbar">
    <div className="relative min-w-0 flex-1">
      <label htmlFor="genome-gene-search" className="sr-only">Найти ген</label>
      <Search className="pointer-events-none absolute left-3 top-1/2 h-4 w-4 -translate-y-1/2 text-muted" aria-hidden />
      <input id="genome-gene-search" value={query} onChange={(event) => setQuery(event.target.value)} placeholder="Найти ген, например MT-ND5" className="workspace-search-input" />
      {query.trim() && <div className="workspace-search-results" role="status">
        {found.length ? found.map((gene) => <button key={gene.name} type="button" onClick={() => { setSelectedGene(gene.name); setSelectedRegion({ start: gene.start, end: gene.end, label: gene.name }); setQuery(''); }} className="flex min-h-10 w-full items-center justify-between gap-2 border-b border-line px-3 text-left text-sm last:border-0 hover:bg-panel2">
          <span className="font-medium">{gene.name}</span><span className="font-mono text-xs text-muted">{gene.start}–{gene.end}</span>
        </button>) : <p className="px-3 py-2 text-sm text-muted">Ген не найден</p>}
      </div>}
    </div>
    <details className="relative shrink-0"><summary className="workspace-header-action list-none"><Layers3 className="h-4 w-4" aria-hidden /> Слои</summary>
      <fieldset className="workspace-layer-menu"><legend className="sr-only">Слои карты</legend>{LAYERS.map(([name, label]) => <label key={name} className="flex min-h-10 cursor-pointer items-center justify-between gap-4 px-3 text-sm hover:bg-panel2"><span>{label}</span><input type="checkbox" checked={activeLayers[name]} onChange={(event) => setLayer(name, event.target.checked)} className="h-4 w-4 accent-aqua" /></label>)}</fieldset>
    </details>
  </div>;
}
