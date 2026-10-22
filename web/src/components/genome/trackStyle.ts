/** Shared encodings for the genome map, linear tracks, legends, and tables. */

export function geneColor(biotype: string): string {
  if (biotype === 'protein_coding') return 'rgb(var(--color-gene-coding))';
  if (biotype === 'tRNA') return 'rgb(var(--color-gene-trna))';
  if (biotype === 'rRNA') return 'rgb(var(--color-gene-rrna))';
  return 'rgb(var(--color-muted))';
}

/** Marker shape by variant class: circle (SNV), triangle (insertion), diamond (deletion). */
export function variantShapePath(type: string, radius: number): string {
  const r = radius;
  if (type.includes('INS')) {
    return `M0,${-r * 1.15}L${r * 1.05},${r * 0.7}L${-r * 1.05},${r * 0.7}Z`;
  }
  if (type.includes('DEL')) {
    return `M0,${-r * 1.2}L${r * 1.2},0L0,${r * 1.2}L${-r * 1.2},0Z`;
  }
  return `M${-r},0A${r},${r} 0 1 0 ${r},0A${r},${r} 0 1 0 ${-r},0Z`;
}
