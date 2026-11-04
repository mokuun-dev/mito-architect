#pragma once
#include "common.hpp"
#include "reads.hpp"
#include "resources.hpp"
#include "variants.hpp"

// Statistical mtDNA variant calling on top of the molecule assembly.
//
// Every aligned mismatch and indel is only a hypothesis here. A call has to
// stand out from a sequencing error background that is estimated from the
// sample itself, with control of the false discovery rate over all tested
// hypotheses. Method references are collected in docs/variant-calling.md.
namespace mito::detail {

// ---------------------------------------------------------------- statistics

/** Thread-safe natural log of the gamma function (std::lgamma writes signgam).
 */
[[nodiscard]] double log_gamma(double value) noexcept;

/** P(K >= k) for K ~ BetaBinomial(n, mean, rho), with rho the intraclass
 * correlation. rho == 0 is the binomial limit. */
[[nodiscard]] double beta_binomial_survival(std::size_t k, std::size_t n,
                                            double mean, double rho) noexcept;

/** Fisher's method for two p-values (chi-squared with 4 degrees of freedom). */
[[nodiscard]] double fisher_combined(double first, double second) noexcept;

/** Two-sided Fisher exact test for the 2x2 table [[a, b], [c, d]]. */
[[nodiscard]] double fisher_exact_two_sided(std::size_t a, std::size_t b,
                                            std::size_t c,
                                            std::size_t d) noexcept;

/** Benjamini-Hochberg adjusted p-values, in input order. total_tests counts
 * hypotheses that were not evaluated because their p-value is one. */
[[nodiscard]] std::vector<double>
benjamini_hochberg(const std::vector<double> &p_values,
                   std::size_t total_tests);

struct SiteCount {
  std::uint32_t k = 0;
  std::uint32_t n = 0;
};

struct BetaBinomialFit {
  double mean = 0.0;
  double rho = 0.0;
  std::size_t sites = 0;
};

/** Robust moment fit (Kleinman 1973). The largest trim_fraction of sites by
 * rate are removed first, so true variants do not inflate the error model. */
[[nodiscard]] BetaBinomialFit
fit_beta_binomial(std::vector<SiteCount> sites, double trim_fraction,
                  double mean_floor, double rho_floor);

/** Upper tail of a rate distribution across a genome, as a Pareto law above a
 * high quantile. Most sites only carry noise, so this tail is the noise itself,
 * including the recurrent motif errors that a sampling model cannot describe. */
struct EmpiricalTail {
  double threshold = 0.0;
  double alpha = 8.0;
  double tail_fraction = 0.01;
};

/** Hill (1975) estimator over the top tail_fraction of rates, where zero_sites
 * counts sites with no observation. The very largest rates, where true variants
 * would be, are treated as right-censored. */
[[nodiscard]] EmpiricalTail fit_empirical_tail(std::vector<double> nonzero_rates,
                                               std::size_t zero_sites,
                                               double minimum_threshold);

/** P(a noise site reaches this rate) under the fitted tail. */
[[nodiscard]] double empirical_tail_probability(const EmpiricalTail &tail,
                                                double rate) noexcept;

/** Wilson 95% interval using the overdispersion-adjusted effective depth. */
[[nodiscard]] std::pair<double, double>
overdispersed_wilson(double k, double n, double rho) noexcept;

/** Two-sided rank-sum (Mann-Whitney) p-value from two histograms over the same
 * ordered bins. mean_shift receives mean_bin(first) - mean_bin(second). */
[[nodiscard]] double rank_sum_p_value(const std::uint32_t *first,
                                      const std::uint32_t *second,
                                      std::size_t bins,
                                      double *mean_shift) noexcept;

// ----------------------------------------------------------- sequence context

/** Tandem-repeat context of each reference position (1-based, circular). */
struct RepeatContext {
  std::vector<std::uint16_t> homopolymer_run;
  std::vector<std::uint8_t> unit;
  std::vector<std::uint16_t> span;
};

[[nodiscard]] RepeatContext build_repeat_context(const std::string &sequence);

/** A circular indel with its full interval of equivalent placements.
 *
 * left and right are the outermost equivalent positions. For a deletion they
 * are the first deleted base; for an insertion, the base preceding it. */
struct NormalizedCircularIndel {
  bool deletion = false;
  std::size_t left = 0;
  std::size_t right = 0;
  std::string left_allele;
  std::string right_allele;
};

[[nodiscard]] NormalizedCircularIndel
normalize_circular_indel(const std::string &reference, std::size_t position,
                         std::size_t deleted_length,
                         std::string_view inserted);

// ------------------------------------------------------------------- pair HMM

/** Context-dependent gap-open probabilities by homopolymer run length. */
struct PairHmmParameters {
  static constexpr std::size_t kMaximumRun = 16U;
  std::array<double, kMaximumRun + 1U> deletion_open{};
  std::array<double, kMaximumRun + 1U> insertion_open{};
  double gap_extend = 0.25;
  double minimum_error = 1e-4;
};

struct PairHmmScratch {
  std::vector<double> match;
  std::vector<double> insertion;
  std::vector<double> deletion;
  std::vector<double> next_match;
  std::vector<double> next_insertion;
  std::vector<double> next_deletion;
  std::vector<std::uint8_t> hap_run;
};

/** Banded forward log10 likelihood of a read window given a haplotype. */
[[nodiscard]] double
pair_hmm_log10_likelihood(std::string_view read, std::string_view qualities,
                          std::string_view haplotype,
                          const PairHmmParameters &parameters,
                          std::size_t band_half_width,
                          PairHmmScratch &scratch);

struct AlleleFractionFit {
  double fraction = 0.0;
  /** Sum of per-read alternate-allele posteriors at the fitted fraction. */
  double effective_alternate = 0.0;
  /** log10 likelihood ratio of the fitted fraction against fraction = 0. */
  double lod = 0.0;
};

/** EM fit of the mixture f * P(read | alt) + (1 - f) * P(read | ref) from
 * per-read log10 likelihoods (the Mutect2 somatic genotyping model). */
[[nodiscard]] AlleleFractionFit
fit_allele_fraction(const std::vector<double> &log10_reference,
                    const std::vector<double> &log10_alternate,
                    const std::vector<double> &weights);

// ----------------------------------------------------------------- annotation

struct StructureMapping {
  std::string pdb_id;
  std::string chain;
  std::size_t residue = 0;
  std::string complex_name;
};

struct VariantEffect {
  std::string gene;
  std::string biotype;
  std::string consequence;
  std::string hgvs_protein;
  std::string codon_change;
  std::size_t amino_acid_position = 0;
  char reference_amino_acid = '\0';
  char alternate_amino_acid = '\0';
  std::optional<int> grantham;
  std::optional<StructureMapping> structure;
};

/** Functional effect of a substitution in every overlapping rCRS feature. */
[[nodiscard]] std::vector<VariantEffect>
annotate_substitution(const std::string &reference, std::size_t position,
                      char alternate);

/** Functional effect of a small indel, placed by the HGVS 3' rule. */
[[nodiscard]] std::vector<VariantEffect>
annotate_small_indel(const std::string &reference, std::size_t position,
                     std::size_t reference_length_removed,
                     std::size_t inserted_length, bool insertion);

struct LargeDeletionImpact {
  std::vector<std::string> genes_lost;
  std::vector<std::string> genes_partial;
  std::size_t transfer_rnas_affected = 0;
  /** A partial protein-coding gene whose reading frame shifts at the junction,
   * or the 5' and 3' partial genes when the junction fuses two. */
  std::string fusion;
  bool heavy_strand_origin_lost = false;
  bool light_strand_origin_lost = false;
  bool major_arc = false;
  /** One breakpoint in 6-9 kb and the other in 13-16 kb (Shamanskiy 2023). */
  bool contact_zone = false;
};

/** Impact of removing the circular interval [start, end]. */
[[nodiscard]] LargeDeletionImpact
annotate_large_deletion(std::size_t start, std::size_t end,
                        std::size_t reference_length);

/** Reference amino acid sequence from translating one gene of rCRS. */
[[nodiscard]] std::string translate_rcrs_gene(const std::string &reference,
                                              std::string_view gene);

/** Physicochemical distance between two amino acids (Grantham 1974). */
[[nodiscard]] std::optional<int> grantham_distance(char first,
                                                   char second) noexcept;

/** Three-letter code for a one-letter amino acid, Ter for '*'. */
[[nodiscard]] std::string_view amino_acid_three_letter(char code) noexcept;

// ------------------------------------------------------------------ the caller

struct PrecisionRepeat {
  std::size_t homopolymer_run = 1;
  std::size_t unit = 1;
  std::size_t span = 1;
};

struct PrecisionStrandSupport {
  std::size_t alternate_forward = 0;
  std::size_t alternate_reverse = 0;
  std::size_t depth_forward = 0;
  std::size_t depth_reverse = 0;
};

struct PrecisionSnv {
  std::size_t position = 0;
  char reference = 'N';
  char alternate = 'N';
  bool pass = false;
  std::vector<std::string> filters;
  std::vector<std::string> flags;
  double fraction = 0.0;
  double ci95_low = 0.0;
  double ci95_high = 0.0;
  std::size_t depth = 0;
  std::size_t alternate_depth = 0;
  PrecisionStrandSupport strands;
  double p_value = 1.0;
  double q_value = 1.0;
  double background_rate = 0.0;
  double alternate_mean_quality = 0.0;
  PrecisionRepeat repeat;
  std::vector<VariantEffect> effects;
};

struct PrecisionIndel {
  std::string id;
  bool deletion = false;
  std::size_t left = 0;
  std::size_t right = 0;
  std::size_t length = 0;
  std::string allele;
  std::size_t vcf_position = 0;
  std::string vcf_reference;
  std::string vcf_alternate;
  std::string hgvs;
  bool pass = false;
  std::vector<std::string> filters;
  std::vector<std::string> flags;
  double fraction = 0.0;
  double ci95_low = 0.0;
  double ci95_high = 0.0;
  std::size_t depth = 0;
  std::size_t raw_alternate_depth = 0;
  double effective_alternate_depth = 0.0;
  double effective_depth = 0.0;
  PrecisionStrandSupport strands;
  double lod = 0.0;
  double p_value = 1.0;
  double q_value = 1.0;
  double background_rate = 0.0;
  bool rescored = false;
  PrecisionRepeat repeat;
  std::vector<VariantEffect> effects;
  std::vector<std::string> supporting_molecules;
};

struct PrecisionDeletion {
  std::string id;
  /** First and last deleted base under the HGVS 3' rule. */
  std::size_t start = 0;
  std::size_t end = 0;
  /** The same event at its leftmost equivalent placement. */
  std::size_t left_start = 0;
  std::size_t left_end = 0;
  std::size_t length = 0;
  std::string microhomology;
  bool pass = false;
  std::vector<std::string> filters;
  std::vector<std::string> flags;
  /** Molecules with the junction over junction plus mean wild-type spanning
   * molecules at the two breakpoints; nullopt when no breakpoint is spanned. */
  std::optional<double> fraction;
  double ci95_low = 0.0;
  double ci95_high = 0.0;
  std::size_t junction_molecules = 0;
  std::size_t junction_forward = 0;
  std::size_t junction_reverse = 0;
  std::size_t wild_type_left = 0;
  std::size_t wild_type_right = 0;
  /** Observed placements merged into this event. */
  std::size_t breakpoint_variants = 0;
  /** Median depth inside the deleted arc over median depth outside it. */
  std::optional<double> depth_ratio;
  double p_value = 1.0;
  double q_value = 1.0;
  bool known_common_deletion = false;
  LargeDeletionImpact impact;
  std::vector<std::string> supporting_molecules;
};

/** A read-phased molecular subpopulation from the mixture model. It may reflect
 * heteroplasmy, a mixed sample or contamination; the data do not say which. */
struct PrecisionHaplotype {
  int id = 0;
  double fraction = 0.0;
  std::size_t molecules = 0;
  /** Supported variants enriched in this component against the largest one. */
  std::vector<std::string> defining_variants;
  /** Supported majority markers of this component, in PhyloTree notation. */
  ReadFeature feature;
  std::optional<ClusterHaplogroupAssignment> haplogroup;
};

struct PrecisionDiagnostics {
  std::size_t molecules_used = 0;
  std::size_t molecules_excluded_divergent = 0;
  std::size_t reference_positions_covered = 0;
  double median_depth = 0.0;
  double mean_substitution_rate = 0.0;
  double median_molecule_divergence = 0.0;
  double divergence_limit = 0.0;
  std::size_t snv_hypotheses = 0;
  std::size_t indel_candidates = 0;
  std::size_t indels_rescored = 0;
  /** 99th percentile of one-base indel rates outside repeats in this sample. */
  double indel_noise_p99 = 0.0;
  /** Rate of one-base length errors by run length, including extrapolated runs.
   * Internal to calling; reports use the measured rates below. */
  std::array<double, 13> slip_rate_by_run{};
  std::size_t junction_molecules = 0;
  std::size_t singleton_junctions = 0;
  std::vector<std::pair<std::size_t, double>> homopolymer_deletion_rate;
  std::vector<std::pair<std::size_t, double>> homopolymer_insertion_rate;
};

struct PrecisionCallResult {
  bool assessed = false;
  std::vector<PrecisionSnv> snvs;
  std::vector<PrecisionIndel> indels;
  std::vector<PrecisionDeletion> deletions;
  PrecisionDiagnostics diagnostics;
  /** Sample consensus: major allele by position, rCRS where nothing is seen. */
  std::string consensus;
  /** Molecules spanning each position, on both strands. */
  std::vector<std::uint32_t> depth_by_position;
  /** Positions with enough depth to assess; everything else is NOT_CALLABLE. */
  std::vector<std::pair<std::size_t, std::size_t>> callable_ranges;
  /** Haplogroup of the major molecules, from supported variants only. */
  std::optional<ClusterHaplogroupAssignment> consensus_haplogroup;
  std::vector<PrecisionHaplotype> haplotypes;
  /** Haplotype of each assembled molecule; -1 when its reads do not cover enough
   * informative sites to place it. */
  std::vector<int> molecule_haplotype;
};

[[nodiscard]] PrecisionCallResult
call_precision_variants(const MoleculeAssemblyResult &assembly,
                        const std::vector<ReadFeature> &features,
                        const std::string &reference,
                        const AnalysisConfig &config);

/** A one-molecule stand-in for the major haplotype, for haplogroup assignment.
 * Its markers are the supported majority variants, in PhyloTree notation. */
[[nodiscard]] ReadFeature
consensus_haplogroup_feature(const PrecisionCallResult &calls,
                             const std::string &reference,
                             const MoleculeAssemblyResult &assembly);

/** The precision_calls JSON object; deterministic for a given result. */
[[nodiscard]] std::string
render_precision_calls_json(const PrecisionCallResult &result,
                            const AnalysisConfig &config);

/** Coverage over all alignments of each molecule. Reads crossing the origin are
 * split, so counting primary alignments alone understates depth. */
[[nodiscard]] CoverageResult
coverage_from_molecule_depth(const std::vector<std::uint32_t> &depth, std::size_t bin_count);

/** Replaces token clusters with the phased haplotype of each molecule. */
void assign_haplotype_clusters(const PrecisionCallResult &calls,
                               std::vector<ReadFeature> &features);

/** Keeps supported events in the read and event projections. Used for noisy
 * long-read samples, where raw per-read errors swamp both. */
void apply_compact_projection(const PrecisionCallResult &calls,
                              const std::string &reference,
                              std::vector<ReadFeature> &features,
                              std::map<std::string, SvCall> &svs,
                              std::map<std::string, ComplexSvCall> &complex_events);

} // namespace mito::detail
