# Statistical variant calling

The engine reports a `precision_calls` section next to the raw per-read
projections. A raw mismatch or gap is only a hypothesis; a call has to stand out
from an error background measured in the same sample, with the false discovery
rate controlled at 1% over all hypotheses. These are research calls, not
clinical findings.

## Methods

**SNVs.** Each position is counted once per molecule, across primary and
supplementary alignments. Errors follow a beta-binomial model per read-oriented
trinucleotide context and homopolymer class (deepSNV, Shearwater). Each strand is
tested separately and combined with Fisher's method, and a call must also show
support on the weaker strand, because sequencer errors depend on the strand read.
Benjamini-Hochberg adjustment counts every hypothesis, including untested ones.
Calls are also rejected when the alternate bases have low quality, cluster at read
ends, or sit mostly next to gaps. An aligner often writes a length error as a gap
plus a mismatch beside it, and beside a long homopolymer as a pair of mismatches
with no gap, as at rCRS 309-310, so low levels there are not reported.

**Small indels.** Placements are normalized on the circular reference (Tan et
al.). Error rates per repeat length come from the sample's own reads. Plausible
candidates are rescored with a banded pair HMM against the reference and competing
alleles at the locus, and an EM fit gives the allele fraction and LOD (the Mutect2
model). The significance test uses molecules that show the gap in their alignment,
the same count the error model is fitted on; only events of 7 bp or more use rescored
counts, because aligners fragment those. A call must also exceed the sample's empirical noise tail, a Pareto law fit
by a right-censored Hill estimator, because motif-specific errors recur at the same
site and a sampling model cannot describe them. The caller rejects short events
that most reads do not show in their alignment, since those fit noisy bases better
than they fit the allele. The 99th percentile of indel noise is reported as
`indel_noise_p99`, so a reader can judge any low-level call.

**Large deletions.** Junctions from CIGAR gaps and split alignments are corrected
for query overlap and moved to their leftmost placement. They are clustered in
(placement, length), since microhomology moves a placement without changing the
length. Breakpoints are refined by maximum likelihood, and the level is junction
molecules over junction plus mean wild-type molecules spanning the breakpoints
(as in mitoSAlt and LostArc). Annotation covers genes lost, junction fusion and
frame, replication origins, and the 6-9 / 13-16 kb contact zone (Shamanskiy et al.).

**Haplotypes.** A Bernoulli mixture over supported heteroplasmic SNVs phases reads,
treating uncovered sites as missing. BIC selects the number of components, and a
new component must differ at three or more sites. A second component can reflect
heteroplasmy, a mixed sample or contamination; the data do not say which.

**NUMT-like reads.** Molecules far above the sample's median divergence are
excluded. This does not establish NUMT specificity, which needs competitive
nuclear-plus-mitochondrial alignment.

## What the calls cannot show

- Low-level calls have no independent confirmation here. On real data, compare them
  with another platform or a replicate before relying on them.
- Sequence in long C-tracts has several equally valid spellings, so a level there is
  less reliable, and calls carry a flag saying so.
- Deletion levels count molecules. Long-range PCR favours shorter molecules, so an
  amplified sample can read high.
- A read name no longer marks a read as NUMT except with `--allow-development-tags`.

## Options

`--compact-min-molecules N` (CLI) and `MITO_COMPACT_MIN_MOLECULES` (server, default
200) keep only supported events in the read, cluster and event projections at or
above N molecules, and use the phased haplotypes as read clusters.

## Testing

`scripts/mtsim/` simulates reads from real genomes with a known truth set, and
`scripts/mtsim/benchmark.py` scores the caller. The simulator's error model is
simplified: on real nanopore reads, indel noise was heavier-tailed than it
simulated, which is why the empirical tail exists. Simulated scores show relative
behaviour and regressions and do not establish accuracy on real sequencers.

## Sources

Retrieved from PubMed.

- Gerstung et al. 2012, deepSNV. [doi:10.1038/ncomms1814](https://doi.org/10.1038/ncomms1814)
- Gerstung et al. 2014, Shearwater. [doi:10.1093/bioinformatics/btt750](https://doi.org/10.1093/bioinformatics/btt750)
- Tan et al. 2015, variant normalization. [doi:10.1093/bioinformatics/btv112](https://doi.org/10.1093/bioinformatics/btv112)
- Behera et al. 2025, DRAGEN. [doi:10.1038/s41587-024-02382-1](https://doi.org/10.1038/s41587-024-02382-1)
- Edge and Bansal 2019, Longshot. [doi:10.1038/s41467-019-12493-y](https://doi.org/10.1038/s41467-019-12493-y)
- Lüth et al. 2022, low-frequency mtDNA calling on nanopore. [doi:10.3389/fgene.2022.887644](https://doi.org/10.3389/fgene.2022.887644)
- Lüth et al. 2026, HiFi mtDNA heteroplasmy. [doi:10.1186/s12864-026-13426-y](https://doi.org/10.1186/s12864-026-13426-y)
- Weissensteiner et al. 2024, mtDNA-Server 2. [doi:10.1093/nar/gkae296](https://doi.org/10.1093/nar/gkae296)
- Basu et al. 2020, mitoSAlt. [doi:10.1371/journal.pgen.1009242](https://doi.org/10.1371/journal.pgen.1009242)
- Lujan et al. 2020, LostArc. [doi:10.1186/s13059-020-02138-5](https://doi.org/10.1186/s13059-020-02138-5)
- Fearn et al. 2026, NanoDel. [doi:10.1093/bioinformatics/btag684](https://doi.org/10.1093/bioinformatics/btag684)
- Shamanskiy et al. 2023, contact zone. [doi:10.1186/s12915-023-01606-1](https://doi.org/10.1186/s12915-023-01606-1)
- Zakarian et al. 2026, MitoScope (preprint). [doi:10.64898/2026.04.13.718259](https://doi.org/10.64898/2026.04.13.718259)
- Schönherr et al. 2023, Haplogrep 3. [doi:10.1093/nar/gkad284](https://doi.org/10.1093/nar/gkad284)
- van Oven and Kayser 2009, PhyloTree. [doi:10.1002/humu.20921](https://doi.org/10.1002/humu.20921)
- Ingman et al. 2000, public genomes. [doi:10.1038/35047064](https://doi.org/10.1038/35047064)
- Grantham 1974, amino acid distance. [doi:10.1126/science.185.4154.862](https://doi.org/10.1126/science.185.4154.862)
