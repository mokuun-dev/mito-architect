# Benchmarks

Measured on 2026-11-04 with `scripts/mtsim/benchmark.py`. Each sample mixes two
public genomes (AF346963 and AF347015). The minor genome carries eight small indels
at its level, and a third component of the same size carries one large deletion.

These reads are simulated. The numbers measure the caller against
`scripts/mtsim/mtsim.py` and catch regressions; they do not establish accuracy on a
real sequencer. Real nanopore indel noise was heavier-tailed than this simulator
produces, so expect lower indel precision on real reads than shown here.

## Reading the table

- **Minor** is the level of the minor genome, and of the deleted molecules.
- **SNV** and **Indel** exclude long C-tracts (for example 300-317 and 16180-16195),
  where one molecule has several equally valid spellings.
- **Raw F1** applies a fixed threshold to the raw per-position table, the only
  option before this caller. The raw table has no indels. It is omitted at 4000x,
  where it is too large to keep.
- **Deletion F1** requires the exact length, so a call one base off counts as missed.
- **Deletion level error** is the mean absolute difference from the true share of
  molecules. It grows with the level because the simulator draws reads per molecule,
  which favours shorter circles. Real libraries differ by protocol.

## Results

| Platform | Minor | Depth | NUMT-like | SNV F1 | SNV P/R | Indel F1 | Indel P/R | Deletion F1 | Deletion level error | Raw F1 (1%, 3 reads) | Raw F1 (2%, 5 reads) | Time, s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| hifi | 0.3 | 1000 | 0 | 1.000 | 1.000/1.000 | 1.000 | 1.000/1.000 | 1.00 | 0.0818 | 1.000 | 1.000 | 0.94 |
| hifi | 0.1 | 1000 | 0 | 1.000 | 1.000/1.000 | 1.000 | 1.000/1.000 | 1.00 | 0.0373 | 1.000 | 1.000 | 0.93 |
| hifi | 0.05 | 1000 | 0 | 1.000 | 1.000/1.000 | 1.000 | 1.000/1.000 | 1.00 | 0.026 | 1.000 | 1.000 | 0.92 |
| hifi | 0.02 | 1000 | 0 | 1.000 | 1.000/1.000 | 1.000 | 1.000/1.000 | 1.00 | 0.018 | 1.000 | 0.959 | 0.82 |
| hifi | 0.01 | 1000 | 0 | 0.990 | 1.000/0.980 | 1.000 | 1.000/1.000 | 1.00 | 0.002 | 0.949 | 0.725 | 1.32 |
| hifi | 0.02 | 4000 | 0 | 1.000 | 1.000/1.000 | 1.000 | 1.000/1.000 | 1.00 | 0.0029 | - | - | 4.97 |
| hifi | 0.01 | 4000 | 0 | 1.000 | 1.000/1.000 | 1.000 | 1.000/1.000 | 1.00 | 0.0019 | - | - | 5.2 |
| hifi | 0.005 | 4000 | 0 | 0.990 | 1.000/0.980 | 1.000 | 1.000/1.000 | 1.00 | 0.001 | - | - | 4.38 |
| hifi | 0.1 | 100 | 0 | 1.000 | 1.000/1.000 | 1.000 | 1.000/1.000 | 1.00 | 0.0378 | 0.915 | 0.891 | 0.25 |
| hifi | 0.05 | 1000 | 0.02 | 1.000 | 1.000/1.000 | 1.000 | 1.000/1.000 | 1.00 | 0.0084 | 0.245 | 0.259 | 1.45 |
| ont_r10 | 0.3 | 1000 | 0 | 0.990 | 0.981/1.000 | 1.000 | 1.000/1.000 | 1.00 | 0.1052 | 0.689 | 1.000 | 1.08 |
| ont_r10 | 0.1 | 1000 | 0 | 1.000 | 1.000/1.000 | 1.000 | 1.000/1.000 | 1.00 | 0.0546 | 0.734 | 1.000 | 1.17 |
| ont_r10 | 0.05 | 1000 | 0 | 0.990 | 1.000/0.980 | 1.000 | 1.000/1.000 | 1.00 | 0.0373 | 0.667 | 1.000 | 0.98 |
| ont_r10 | 0.02 | 1000 | 0 | 0.990 | 1.000/0.980 | 1.000 | 1.000/1.000 | 1.00 | 0.0077 | 0.689 | 0.756 | 1.42 |
| ont_r10 | 0.01 | 1000 | 0 | 0.990 | 1.000/0.980 | 0.857 | 1.000/0.750 | 1.00 | 0.0087 | 0.632 | 0.762 | 0.95 |
| ont_r10 | 0.02 | 4000 | 0 | 0.990 | 1.000/0.980 | 0.941 | 0.889/1.000 | 1.00 | 0.0015 | - | - | 3.78 |
| ont_r10 | 0.01 | 4000 | 0 | 0.980 | 1.000/0.961 | 0.875 | 0.875/0.875 | 1.00 | 0.0022 | - | - | 5.26 |
| ont_r10 | 0.005 | 4000 | 0 | 0.841 | 1.000/0.726 | 0.769 | 1.000/0.625 | 1.00 | 0.0002 | - | - | 5.72 |
| ont_r10 | 0.1 | 100 | 0 | 0.990 | 1.000/0.980 | 1.000 | 1.000/1.000 | 1.00 | 0.0629 | 0.879 | 1.000 | 0.3 |
| ont_r10 | 0.05 | 1000 | 0.02 | 0.990 | 1.000/0.980 | 1.000 | 1.000/1.000 | 1.00 | 0.0042 | 0.217 | 0.265 | 1.52 |
| ont_r9 | 0.3 | 1000 | 0 | 0.981 | 0.962/1.000 | 0.842 | 0.727/1.000 | 0.00 | - | 0.105 | 0.545 | 1.8 |
| ont_r9 | 0.1 | 1000 | 0 | 0.971 | 0.962/0.980 | 1.000 | 1.000/1.000 | 1.00 | 0.0466 | 0.111 | 0.551 | 2.12 |
| ont_r9 | 0.05 | 1000 | 0 | 0.980 | 0.980/0.980 | 0.889 | 0.800/1.000 | 1.00 | 0.0136 | 0.113 | 0.545 | 1.63 |
| ont_r9 | 0.02 | 1000 | 0 | 0.866 | 0.913/0.824 | 0.750 | 0.750/0.750 | 1.00 | 0.0141 | 0.110 | 0.471 | 1.89 |
| ont_r9 | 0.01 | 1000 | 0 | 0.742 | 0.868/0.647 | 0.615 | 0.800/0.500 | 1.00 | 0.0079 | 0.106 | 0.337 | 1.53 |
| ont_r9 | 0.02 | 4000 | 0 | 0.893 | 0.820/0.980 | 0.500 | 0.417/0.625 | 1.00 | 0.0061 | - | - | 5.05 |
| ont_r9 | 0.01 | 4000 | 0 | 0.800 | 0.778/0.824 | 0.706 | 0.667/0.750 | 1.00 | 0.0 | - | - | 6.76 |
| ont_r9 | 0.005 | 4000 | 0 | 0.674 | 0.829/0.569 | 0.333 | 0.500/0.250 | 1.00 | 0.001 | - | - | 4.54 |
| ont_r9 | 0.1 | 100 | 0 | 0.970 | 0.980/0.961 | 0.667 | 1.000/0.500 | 1.00 | 0.0714 | 0.562 | 0.951 | 0.36 |
| ont_r9 | 0.05 | 1000 | 0.02 | 0.971 | 0.962/0.980 | 0.889 | 0.800/1.000 | 1.00 | 0.005 | 0.085 | 0.202 | 2.16 |

## Real reads

The first 3,000 reads of the public GM12878 nanopore run SRR18110025 (median depth
1,547x). No truth set exists for this sample, so this shows behaviour and not accuracy.

| Measure | Raw projections | Statistical calls |
| --- | --- | --- |
| SNV records | 26,255 | 18 |
| Small indels | not reported | 0 |
| Read clusters | 2,463 | 1 phased lineage |
| Haplogroup | a different label per cluster | H13a1a1a, score 99.2 |
| Result size | 111 MB | 2.8 MB |
| Run time, 8 threads | 5.6 s | 3.3 s |
| Peak memory | 1,070 MB | 688 MB |

The haplogroup is not confirmed by an independent source. The low-level SNVs at
this depth (2645G>A, 14535C>A, 14918G>A) and 16023G>A at 55% have no independent
confirmation either, and may be systematic errors.
