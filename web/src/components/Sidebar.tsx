import type { ClusterSummary, GeneAnnotation, LayerName } from '@mito-architect/visualization-lib';
import { Layers, Search, SlidersHorizontal } from 'lucide-react';
import { useMemo } from 'react';
import { useMitoStore } from '../lib/store';

const LAYERS: Array<[LayerName, string]> = [
  ['genes', 'Genes'],
  ['coverage', 'Coverage'],
  ['clusters', 'Clusters'],
  ['svs', 'SVs'],
  ['snps', 'SNPs']
];

interface SidebarProps {
  genes: GeneAnnotation[];
  clusters: ClusterSummary[];
}

export default function Sidebar({ genes, clusters }: SidebarProps) {
  const activeLayers = useMitoStore((state) => state.activeLayers);
  const setLayer = useMitoStore((state) => state.setLayer);
  const geneQuery = useMitoStore((state) => state.geneQuery);
  const setGeneQuery = useMitoStore((state) => state.setGeneQuery);
  const minQuality = useMitoStore((state) => state.minQuality);
  const setMinQuality = useMitoStore((state) => state.setMinQuality);
  const selectedCluster = useMitoStore((state) => state.selectedCluster);
  const setSelectedCluster = useMitoStore((state) => state.setSelectedCluster);
  const setSelectedGene = useMitoStore((state) => state.setSelectedGene);
  const setSelectedRegion = useMitoStore((state) => state.setSelectedRegion);

  const filteredGenes = useMemo(() => {
    const query = geneQuery.trim().toLowerCase();
    return query ? genes.filter((gene) => gene.name.toLowerCase().includes(query)) : genes;
  }, [geneQuery, genes]);

  return (
    <aside className="grid min-w-0 content-start gap-4 xl:sticky xl:top-20" aria-label="Genome view controls">
      <fieldset className="glass-panel rounded-xl border border-line p-4 shadow-tool">
        <legend className="sr-only">Genome map layers</legend>
        <h2 className="mb-3 flex items-center gap-2 text-sm font-semibold">
          <Layers className="h-4 w-4 text-aqua" aria-hidden />
          Map layers
        </h2>
        <div className="grid grid-cols-2 gap-2">
          {LAYERS.map(([layer, label]) => (
            <label
              key={layer}
              className={[
                'flex min-h-11 cursor-pointer items-center justify-between gap-3 rounded-lg border px-3 py-2 text-sm',
                activeLayers[layer] ? 'border-aqua/70 bg-aqua/10 text-text' : 'border-line bg-panel2 text-muted'
              ].join(' ')}
            >
              <span>{label}</span>
              <input
                type="checkbox"
                checked={activeLayers[layer]}
                onChange={(event) => setLayer(layer, event.currentTarget.checked)}
                className="h-4 w-4 accent-aqua"
              />
            </label>
          ))}
        </div>
      </fieldset>

      <section className="glass-panel rounded-xl border border-line p-4 shadow-tool">
        <h2 className="mb-3 flex items-center gap-2 text-sm font-semibold">
          <Search className="h-4 w-4 text-aqua" aria-hidden />
          Gene Search
        </h2>
        <label htmlFor="gene-search" className="sr-only">Filter mitochondrial genes</label>
        <input
          id="gene-search"
          value={geneQuery}
          onChange={(event) => setGeneQuery(event.target.value)}
          placeholder="MT-ND5"
          className="min-h-11 w-full rounded-lg border border-line bg-panel2 px-3 py-2 text-sm placeholder:text-muted focus:border-aqua focus:shadow-focus"
        />
        <div className="mt-3 max-h-40 overflow-auto pr-1 text-sm scrollbar-thin" aria-live="polite">
          {filteredGenes.map((gene) => (
            <button key={gene.name} type="button" onClick={() => { setSelectedGene(gene.name); setSelectedRegion({ start: gene.start, end: gene.end, label: gene.name }); }} className="flex w-full justify-between border-b border-line py-2 text-left last:border-0 hover:text-aqua">
              <span className="font-medium">{gene.name}</span>
              <span className="text-muted">{gene.start}-{gene.end}</span>
            </button>
          ))}
        </div>
      </section>

      <section className="glass-panel rounded-xl border border-line p-4 shadow-tool">
        <h2 className="mb-3 flex items-center gap-2 text-sm font-semibold">
          <SlidersHorizontal className="h-4 w-4 text-aqua" aria-hidden />
          Evidence filter
        </h2>
        <div className="flex items-center justify-between gap-3">
          <label htmlFor="minimum-quality" className="text-xs font-semibold uppercase tracking-[0.1em] text-muted">Minimum mean quality</label>
          <output htmlFor="minimum-quality" className="font-mono text-sm font-semibold text-aqua">Q{minQuality}</output>
        </div>
        <input
          id="minimum-quality"
          type="range"
          min={0}
          max={45}
          step={1}
          value={minQuality}
          onChange={(event) => setMinQuality(Number(event.currentTarget.value))}
          aria-valuetext={`Minimum mean quality ${minQuality}`}
          className="mt-3 min-h-8 w-full accent-aqua"
        />
      </section>

      <section className="glass-panel rounded-xl border border-line p-4 shadow-tool">
        <h2 className="mb-3 text-sm font-semibold">Read clusters</h2>
        <div className="grid gap-2">
          <button
            type="button"
            onClick={() => setSelectedCluster(undefined)}
            className={[
              'min-h-11 rounded-lg border px-3 py-2 text-left text-sm',
              selectedCluster === undefined ? 'border-aqua bg-aqua/10 text-aqua shadow-focus' : 'border-line bg-panel2 text-muted'
            ].join(' ')}
            aria-pressed={selectedCluster === undefined}
          >
            All clusters
          </button>
          {clusters.map((cluster) => (
            <button
              key={cluster.id}
              type="button"
              onClick={() => setSelectedCluster(cluster.id)}
              className={[
                'min-h-11 rounded-lg border px-3 py-2 text-left text-sm',
                selectedCluster === cluster.id ? 'border-aqua bg-aqua/10 text-aqua shadow-focus' : 'border-line bg-panel2 text-muted hover:border-aqua/50'
              ].join(' ')}
              aria-pressed={selectedCluster === cluster.id}
            >
              <span className="flex items-center justify-between gap-3">
                <span>
                  <span className="font-semibold">{cluster.label}</span>
                  <span className="ml-2 text-muted">{cluster.size} records</span>
                </span>
                {cluster.haplogroup ? (
                  <span className="text-xs text-muted">{cluster.haplogroup}</span>
                ) : null}
              </span>
            </button>
          ))}
        </div>
      </section>
    </aside>
  );
}
