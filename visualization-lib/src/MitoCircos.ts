import * as d3 from 'd3';
import { circularVariants, type CircularVariant } from './circular/dataAdapter';
import { intervalMidpoint, positionToAngle } from './circular/coordinates';
import type { GeneAnnotation, LayerName, LayerState, MitoAnalysisData, MitoCircosOptions, StructuralVariant } from './types';

const DEFAULT_LAYERS: LayerState = { genes: true, coverage: true, clusters: false, svs: true, snps: true };
const GENE_COLORS: Record<string, string> = { protein_coding: 'var(--mito-gene-coding)', tRNA: 'var(--mito-gene-trna)', rRNA: 'var(--mito-gene-rrna)' };
const THEME_TOKENS = ['shell', 'panel', 'panel-2', 'line', 'text', 'muted', 'accent', 'gene-coding', 'gene-trna', 'gene-rrna', 'variant', 'sv', 'coverage'];

export interface MitoCircosSelection {
  eventId?: string;
  svId?: string;
  gene?: string;
}

/** Interactive circular reference overview; precise values remain in linked tables. */
export class MitoCircos {
  private readonly container: HTMLElement;
  private data: MitoAnalysisData;
  private readonly options: MitoCircosOptions;
  private layers: LayerState;
  private svg?: d3.Selection<SVGSVGElement, unknown, null, undefined>;
  private plot?: d3.Selection<SVGGElement, unknown, null, undefined>;
  private drawing?: d3.Selection<SVGGElement, unknown, null, undefined>;
  private tooltip?: HTMLDivElement;
  private zoomBehavior?: d3.ZoomBehavior<SVGSVGElement, unknown>;
  private transform = d3.zoomIdentity;
  private selection: MitoCircosSelection = {};
  private width = 760;
  private height = 620;
  private radius = 250;

  constructor(container: HTMLElement, data: MitoAnalysisData, options: MitoCircosOptions = {}) {
    this.container = container;
    this.data = data;
    this.options = options;
    this.layers = { ...DEFAULT_LAYERS, ...options.activeLayers };
  }

  render(): this {
    if (this.svg) return this.redraw();
    this.measure();
    const root = d3.select(this.container).append('div').attr('class', 'mito-circos mito-circos-' + (this.options.theme ?? 'dark') + (this.options.compact ? ' mito-circos-compact' : ''));
    const frame = root.append('div').attr('class', 'mito-circos-frame');
    this.svg = frame.append('svg').attr('viewBox', '0 0 ' + this.width + ' ' + this.height).attr('role', 'img')
      .attr('aria-label', 'Circular mtDNA evidence overview for ' + this.data.metadata.sample).attr('class', 'mito-circos-svg');
    this.tooltip = frame.append('div').attr('class', 'mito-circos-tooltip').node() ?? undefined;
    this.plot = this.svg.append('g');
    if (!this.options.compact) {
      this.zoomBehavior = d3.zoom<SVGSVGElement, unknown>().scaleExtent([0.7, 6])
        .filter((event) => event.type !== 'wheel' || event.ctrlKey)
        .on('zoom', (event) => {
          this.transform = event.transform;
          this.positionPlot();
        });
      this.svg.call(this.zoomBehavior);
    }
    this.injectStyles();
    return this.redraw();
  }

  updateData(data: MitoAnalysisData): this {
    this.data = data;
    this.transform = d3.zoomIdentity;
    return this.svg ? this.redraw() : this.render();
  }

  resize(): this {
    if (!this.svg) return this;
    const previousWidth = this.width;
    const previousHeight = this.height;
    this.measure();
    if (this.width === previousWidth && this.height === previousHeight) return this;
    this.svg.attr('viewBox', '0 0 ' + this.width + ' ' + this.height);
    return this.redraw();
  }

  destroy(): void {
    this.container.innerHTML = '';
    this.svg = undefined;
    this.plot = undefined;
    this.drawing = undefined;
    this.tooltip = undefined;
  }

  exportSVG(): string {
    const node = this.svg?.node();
    if (!node) return '';
    const clone = node.cloneNode(true) as SVGSVGElement;
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    clone.setAttribute('class', 'mito-circos mito-circos-svg');
    const rootStyle = getComputedStyle(document.documentElement);
    for (const name of THEME_TOKENS) {
      clone.style.setProperty('--color-' + name, rootStyle.getPropertyValue('--color-' + name));
    }
    const styles = document.getElementById('mito-circos-styles')?.textContent;
    if (styles) {
      const embedded = document.createElementNS('http://www.w3.org/2000/svg', 'style');
      embedded.textContent = styles;
      clone.insertBefore(embedded, clone.firstChild);
    }
    const title = document.createElementNS('http://www.w3.org/2000/svg', 'title');
    title.textContent = 'mtDNA circular reference map: ' + this.data.metadata.sample;
    clone.insertBefore(title, clone.firstChild);
    const description = document.createElementNS('http://www.w3.org/2000/svg', 'desc');
    description.textContent = 'Reference ' + (this.data.metadata.reference_accession ?? 'custom') + ', length ' + this.data.metadata.reference_length + ' bp. Layers: ' + this.getActiveLayers().join(', ') + '. Coordinates are a reference projection; exact values are in the analysis result.';
    clone.insertBefore(description, title.nextSibling);
    return new XMLSerializer().serializeToString(clone);
  }

  getActiveLayers(): LayerName[] {
    return (Object.keys(this.layers) as LayerName[]).filter((layer) => this.layers[layer]);
  }

  zoomBy(factor: number): this {
    if (this.svg && this.zoomBehavior) this.svg.call(this.zoomBehavior.scaleBy, factor);
    return this;
  }

  resetZoom(): this {
    if (this.svg && this.zoomBehavior) this.svg.call(this.zoomBehavior.transform, d3.zoomIdentity);
    return this;
  }

  setLayer(layer: LayerName, active: boolean): this {
    this.layers[layer] = active;
    return this.redraw();
  }

  setLayers(layers: Partial<LayerState>): this {
    this.layers = { ...this.layers, ...layers };
    return this.redraw();
  }

  /** Highlights one variant, structural event, or gene. */
  setSelection(selection: MitoCircosSelection): this {
    this.selection = { ...selection };
    return this.redraw();
  }

  private redraw(): this {
    if (!this.plot) return this;
    this.drawing?.remove();
    this.drawing = this.plot.append('g').attr('class', 'mito-circos-drawing');
    this.positionPlot();
    this.drawBackbone();
    if (this.layers.coverage) this.drawCoverage();
    if (this.layers.genes) this.drawGenes();
    if (this.layers.svs) this.drawStructuralVariants();
    if (this.layers.snps) this.drawVariants();
    // Cluster membership has no genomic coordinate; it belongs in the architecture panel.
    this.drawSelection();
    this.drawCenterLabel();
    return this;
  }

  private positionPlot(): void {
    this.plot?.attr('transform', 'translate(' + this.width / 2 + ', ' + this.height / 2 + ') ' + this.transform.toString());
  }

  private measure(): void {
    const bounds = this.container.getBoundingClientRect();
    const compact = Boolean(this.options.compact);
    this.width = Math.max(compact ? 220 : 320, this.options.width ?? (bounds.width || 760));
    this.height = Math.max(compact ? 220 : 420, this.options.height ?? (compact ? this.width : Math.min(760, this.width * 0.82)));
    this.radius = Math.max(compact ? 70 : 130, Math.min(this.width, this.height) / 2 - (compact ? 30 : 64));
  }

  private drawBackbone(): void {
    this.drawing?.append('circle').attr('r', this.radius + 17).attr('fill', 'none').attr('stroke', 'var(--mito-border)').attr('stroke-width', 1.5);
    const divisions = this.options.compact ? 4 : 12;
    const step = Math.max(1000, Math.ceil(this.data.metadata.reference_length / divisions / 1000) * 1000);
    const ticks = d3.range(0, this.data.metadata.reference_length, step);
    const track = this.drawing?.append('g');
    track?.selectAll('line').data(ticks).join('line')
      .attr('x1', (d) => this.point(d, this.radius + 15)[0]).attr('y1', (d) => this.point(d, this.radius + 15)[1])
      .attr('x2', (d) => this.point(d, this.radius + 23)[0]).attr('y2', (d) => this.point(d, this.radius + 23)[1])
      .attr('stroke', 'var(--mito-muted)').attr('stroke-opacity', 0.7);
    track?.selectAll('text').data(ticks).join('text').attr('class', 'mito-circos-tick-label')
      .attr('x', (d) => this.point(d, this.radius + (this.options.compact ? 31 : 37))[0]).attr('y', (d) => this.point(d, this.radius + (this.options.compact ? 31 : 37))[1])
      .attr('text-anchor', 'middle').attr('dominant-baseline', 'middle').text((d) => d === 0 ? '1' : Math.round(d / 1000) + 'k');
  }

  private drawCoverage(): void {
    const maxDepth = Math.max(1, d3.max(this.data.coverage, (d) => d.depth) ?? 1);
    const base = this.radius - this.ring(34);
    const arc = d3.arc<{ start: number; end: number; depth: number }>().innerRadius(base)
      .outerRadius((d) => base + Math.max(0, d.depth) / maxDepth * this.ring(38))
      .startAngle((d) => this.boundaryAngle(d.start - 1)).endAngle((d) => this.boundaryAngle(d.end));
    this.drawing?.append('g').selectAll('path').data(this.data.coverage).join('path')
      .attr('d', (d) => arc(d) ?? '').attr('fill', 'var(--mito-coverage)').attr('fill-opacity', 0.7).attr('stroke', 'var(--mito-coverage)').attr('stroke-opacity', 0.9).attr('stroke-width', 0.5)
      .on('pointermove', (event, d) => this.showTooltip(event, 'Coverage ' + d.start + '–' + d.end + ': ' + d.depth + '×')).on('pointerleave', () => this.hideTooltip());
  }

  private drawGenes(): void {
    const arc = d3.arc<GeneAnnotation>().innerRadius(this.radius - this.ring(78)).outerRadius(this.radius - this.ring(54)).padAngle(0.004).cornerRadius(this.options.compact ? 2 : 4)
      .startAngle((d) => this.boundaryAngle(d.start - 1)).endAngle((d) => this.boundaryAngle(d.end));
    const track = this.drawing?.append('g');
    track?.selectAll('path').data(this.data.genes).join('path').attr('d', (d) => arc(d) ?? '')
      .attr('fill', (d) => GENE_COLORS[d.biotype] ?? 'var(--mito-muted)').attr('fill-opacity', 0.85).attr('stroke', 'var(--mito-bg)').attr('stroke-width', 1).attr('tabindex', 0).attr('role', 'button')
      .attr('aria-label', (d) => d.name + ', ' + d.start + ' to ' + d.end + ', strand ' + d.strand)
      .on('pointermove', (event, d) => this.showTooltip(event, d.name + ' · ' + d.start + '–' + d.end + ' · strand ' + d.strand)).on('pointerleave', () => this.hideTooltip())
      .on('click', (_, d) => this.dispatch('mito:gene-select', d)).on('keydown', (event, d) => this.activate(event, () => this.dispatch('mito:gene-select', d)));
    const labelCount = this.options.compact ? 0 : Math.min(12, Math.max(6, Math.floor(this.radius / 20)));
    const labels = [...this.data.genes].sort((a, b) => b.end - b.start - (a.end - a.start)).slice(0, labelCount);
    track?.selectAll('text').data(labels).join('text').attr('class', 'mito-circos-label')
      .attr('x', (d) => this.point(intervalMidpoint(d.start, d.end, this.data.metadata.reference_length), this.radius - this.ring(42))[0])
      .attr('y', (d) => this.point(intervalMidpoint(d.start, d.end, this.data.metadata.reference_length), this.radius - this.ring(42))[1])
      .attr('text-anchor', 'middle').attr('dominant-baseline', 'middle').text((d) => d.name.replace('MT-', ''));
  }

  private drawStructuralVariants(): void {
    this.drawing?.append('g').selectAll('path').data(this.data.svs).join('path').attr('d', (d) => this.chordPath(d)).attr('fill', 'none')
      .attr('stroke', 'var(--mito-sv)').attr('stroke-width', (d) => Math.min(this.options.compact ? 5 : 10, 1.5 + Math.sqrt(d.supporting_reads.length) * 1.5)).attr('stroke-opacity', 0.8).attr('stroke-linecap', 'round').attr('tabindex', 0).attr('role', 'button')
      .attr('aria-label', (d) => d.type + ' ' + d.start + ' to ' + d.end + '; ' + d.supporting_reads.length + ' supporting molecules')
      .on('pointermove', (event, d) => this.showTooltip(event, d.type + ' ' + d.start + '–' + d.end + '; ' + d.supporting_reads.length + ' supporting molecules')).on('pointerleave', () => this.hideTooltip())
      .on('click', (_, d) => this.dispatch('mito:sv-select', d)).on('keydown', (event, d) => this.activate(event, () => this.dispatch('mito:sv-select', d)));
  }

  private drawVariants(): void {
    const variants = circularVariants(this.data);
    this.drawing?.append('g').selectAll('path').data(variants).join('path')
      .attr('transform', (d) => { const point = this.point(d.position - 1, this.variantRadius()); return 'translate(' + point[0] + ',' + point[1] + ')'; })
      .attr('d', (d) => d3.symbol().type(d.type.includes('INS') ? d3.symbolTriangle : d.type.includes('DEL') ? d3.symbolDiamond : d3.symbolCircle).size((d.type === 'SNV' ? 45 : 64) * (this.options.compact ? 0.6 : 1))())
      .attr('fill', 'var(--mito-variant)').attr('stroke', 'var(--mito-bg)').attr('stroke-width', 1).attr('tabindex', 0).attr('role', 'button')
      .attr('aria-label', (d) => d.id + '; ' + d.position + ' ' + d.ref + ' to ' + d.alt + '; support ' + d.support)
      .on('pointermove', (event, d) => this.showTooltip(event, d.position + ' ' + d.ref + '>' + d.alt + ' · ALT ' + d.support + (d.callableDepth !== undefined ? ' / callable ' + d.callableDepth : ''))).on('pointerleave', () => this.hideTooltip())
      .on('click', (_, d) => this.dispatch('mito:event-select', d)).on('keydown', (event, d) => this.activate(event, () => this.dispatch('mito:event-select', d)));
  }

  private drawSelection(): void {
    const { eventId, svId, gene } = this.selection;
    const length = this.data.metadata.reference_length;
    let position: number | undefined;
    const group = this.drawing?.append('g').attr('class', 'mito-circos-selection').attr('pointer-events', 'none');
    if (eventId && this.layers.snps) {
      const variant = circularVariants(this.data).find((item) => item.id === eventId);
      if (variant) {
        position = variant.position;
        const [x, y] = this.point(variant.position - 1, this.variantRadius());
        group?.append('circle').attr('cx', x).attr('cy', y).attr('r', this.options.compact ? 6 : 9)
          .attr('fill', 'none').attr('stroke', 'var(--mito-text)').attr('stroke-width', 2);
      }
    }
    if (svId && this.layers.svs) {
      const sv = this.data.svs.find((item) => item.id === svId);
      if (sv) {
        position = sv.start;
        group?.append('path').attr('d', this.chordPath(sv)).attr('fill', 'none').attr('stroke', 'var(--mito-text)')
          .attr('stroke-width', this.options.compact ? 2 : 3).attr('stroke-linecap', 'round');
      }
    }
    if (gene && this.layers.genes) {
      const match = this.data.genes.find((item) => item.name === gene);
      if (match) {
        position = intervalMidpoint(match.start, match.end, length);
        const arc = d3.arc<GeneAnnotation>().innerRadius(this.radius - this.ring(80)).outerRadius(this.radius - this.ring(52)).cornerRadius(4)
          .startAngle((d) => this.boundaryAngle(d.start - 1)).endAngle((d) => this.boundaryAngle(d.end));
        group?.append('path').attr('d', arc(match) ?? '').attr('fill', 'none').attr('stroke', 'var(--mito-text)').attr('stroke-width', 2);
      }
    }
    if (position !== undefined) {
      // A pointer outside the tick ring shows where the selection sits.
      const tip = this.point(position - 1, this.radius + 13);
      const left = this.point(position - 1 - length * 0.008, this.radius + 22);
      const right = this.point(position - 1 + length * 0.008, this.radius + 22);
      group?.append('path').attr('d', 'M' + tip.join(',') + 'L' + left.join(',') + 'L' + right.join(',') + 'Z').attr('fill', 'var(--mito-text)');
    }
  }

  private drawCenterLabel(): void {
    if (this.options.compact) return;
    const group = this.drawing?.append('g');
    group?.append('text').attr('text-anchor', 'middle').attr('y', -12).attr('class', 'mito-circos-title').text('mtDNA');
    group?.append('text').attr('text-anchor', 'middle').attr('y', 13).attr('class', 'mito-circos-subtitle').text(this.data.metadata.reference_length.toLocaleString() + ' bp');
  }

  private boundaryAngle(boundary: number): number { return positionToAngle(boundary, this.data.metadata.reference_length); }
  private point(boundary: number, radius: number): [number, number] { const angle = this.boundaryAngle(boundary) - Math.PI / 2; return [Math.cos(angle) * radius, Math.sin(angle) * radius]; }
  /** Scales a radial offset designed for a 250px radius to the current radius. */
  private ring(offset: number): number { return offset * Math.min(1, this.radius / 250); }
  private variantRadius(): number { return this.radius - this.ring(116); }
  private chordPath(sv: StructuralVariant): string { const start = this.point(sv.start - 1, this.radius - this.ring(132)); const end = this.point(sv.end - 1, this.radius - this.ring(132)); return 'M ' + start[0].toFixed(2) + ' ' + start[1].toFixed(2) + ' Q 0 0 ' + end[0].toFixed(2) + ' ' + end[1].toFixed(2); }
  private activate(event: KeyboardEvent, action: () => void): void { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); action(); } }
  private showTooltip(event: PointerEvent, text: string): void { if (!this.tooltip) return; const host = this.container.getBoundingClientRect(); this.tooltip.textContent = text; this.tooltip.style.opacity = '1'; this.tooltip.style.left = event.clientX - host.left + 14 + 'px'; this.tooltip.style.top = event.clientY - host.top + 14 + 'px'; }
  private hideTooltip(): void { if (this.tooltip) this.tooltip.style.opacity = '0'; }
  private dispatch(name: string, detail: GeneAnnotation | StructuralVariant | CircularVariant): void { this.container.dispatchEvent(new CustomEvent(name, { detail, bubbles: true })); }

  private injectStyles(): void {
    if (document.getElementById('mito-circos-styles')) return;
    const style = document.createElement('style');
    style.id = 'mito-circos-styles';
    style.textContent = `
      .mito-circos { --mito-bg: rgb(var(--color-panel)); --mito-panel: rgb(var(--color-panel)); --mito-border: rgb(var(--color-line)); --mito-text: rgb(var(--color-text)); --mito-muted: rgb(var(--color-muted)); --mito-coverage: rgb(var(--color-coverage)); --mito-variant: rgb(var(--color-variant)); --mito-sv: rgb(var(--color-sv)); --mito-gene-coding: rgb(var(--color-gene-coding)); --mito-gene-trna: rgb(var(--color-gene-trna)); --mito-gene-rrna: rgb(var(--color-gene-rrna)); width: 100%; color: var(--mito-text); }
      .mito-circos-frame { position: relative; overflow: hidden; background: var(--mito-panel); }
      .mito-circos-svg { display: block; width: 100%; height: auto; cursor: grab; }
      .mito-circos-svg:active { cursor: grabbing; }
      .mito-circos-compact .mito-circos-svg { cursor: default; }
      .mito-circos-svg [role=button]:focus-visible { outline: 2px solid var(--mito-variant); }
      .mito-circos-label { fill: var(--mito-muted); font: 600 10px/1 'Inter Variable',ui-sans-serif,system-ui,sans-serif; }
      .mito-circos-tick-label { fill: var(--mito-muted); font: 500 9px/1 'Inter Variable',ui-sans-serif,system-ui,sans-serif; }
      .mito-circos-title { fill: var(--mito-text); font: 650 18px/1 'Inter Variable',ui-sans-serif,system-ui,sans-serif; }
      .mito-circos-subtitle { fill: var(--mito-muted); font: 500 12px/1 'Inter Variable',ui-sans-serif,system-ui,sans-serif; }
      .mito-circos-tooltip { position: absolute; pointer-events: none; opacity: 0; z-index: 4; max-width: 280px; border: 1px solid var(--mito-border); border-radius: 6px; background: var(--mito-panel); color: var(--mito-text); padding: 7px 9px; box-shadow: 0 10px 24px rgb(0 0 0 / .16); font: 500 12px/1.35 'Inter Variable',ui-sans-serif,system-ui,sans-serif; transition: opacity 120ms ease; }
    `;
    document.head.appendChild(style);
  }
}
