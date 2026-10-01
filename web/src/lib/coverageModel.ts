import type { CoverageBin, MitoAnalysisData } from '@mito-architect/visualization-lib';

export interface CoverageSummary {
  meanDepth: number;
  pctGt20: number;
  maxDepth: number;
  pctLabel: 'sites' | 'bins';
  zeroBins: number;
  lowBins: number;
}

export function summarizeCoverage(data: MitoAnalysisData): CoverageSummary {
  if (data.coverage_metrics) {
    const hasSites = data.coverage_metrics.pct_sites_gt20x !== undefined;
    return {
      meanDepth: data.coverage_metrics.mean_depth,
      pctGt20: hasSites ? data.coverage_metrics.pct_sites_gt20x! : data.coverage_metrics.pct_bins_gt20x ?? 0,
      maxDepth: data.coverage_metrics.max_depth,
      pctLabel: hasSites ? 'sites' : 'bins',
      zeroBins: data.coverage.filter((bin) => bin.depth === 0).length,
      lowBins: data.coverage.filter((bin) => bin.depth > 0 && bin.depth < 20).length
    };
  }
  const bins = data.coverage;
  const bases = bins.reduce((sum, bin) => sum + binLength(bin), 0);
  const weightedDepth = bins.reduce((sum, bin) => sum + bin.depth * binLength(bin), 0);
  const maxDepth = Math.max(0, ...bins.map((bin) => bin.depth));
  return {
    meanDepth: bases === 0 ? 0 : weightedDepth / bases,
    pctGt20: bases === 0 ? 0 : bins.filter((bin) => bin.depth > 20).reduce((sum, bin) => sum + binLength(bin), 0) * 100 / bases,
    maxDepth,
    pctLabel: 'sites',
    zeroBins: bins.filter((bin) => bin.depth === 0).length,
    lowBins: bins.filter((bin) => bin.depth > 0 && bin.depth < 20).length
  };
}

export function binLength(bin: CoverageBin): number {
  return Math.max(0, bin.end - bin.start + 1);
}

export function coverageDistribution(bins: CoverageBin[]): Array<{ label: string; count: number }> {
  const groups = [
    ['0×', (depth: number) => depth === 0],
    ['1–9×', (depth: number) => depth >= 1 && depth < 10],
    ['10–19×', (depth: number) => depth >= 10 && depth < 20],
    ['20–49×', (depth: number) => depth >= 20 && depth < 50],
    ['50×+', (depth: number) => depth >= 50]
  ] as const;
  return groups.map(([label, matches]) => ({ label, count: bins.filter((bin) => matches(bin.depth)).reduce((sum, bin) => sum + binLength(bin), 0) }));
}
