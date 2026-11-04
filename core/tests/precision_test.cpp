// Unit tests for the statistical calling methods, with expected values worked
// out independently of this code (exact sums, published examples, or rCRS facts).
#include "detail/pipeline.hpp"
#include "detail/precision.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

namespace {

using namespace mito;
using namespace mito::detail;

int failures = 0;

void require(const bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAILED: " << message << "\n";
    ++failures;
  }
}

bool near(const double value, const double expected, const double tolerance) {
  return std::abs(value - expected) <= tolerance;
}

void test_log_gamma() {
  require(near(log_gamma(1.0), 0.0, 1e-12), "log_gamma(1)");
  require(near(log_gamma(5.0), std::log(24.0), 1e-12), "log_gamma(5)");
  require(near(log_gamma(0.5), 0.5 * std::log(3.14159265358979323846), 1e-12), "log_gamma(1/2)");
  require(near(log_gamma(1001.0), 5912.128178488163, 1e-8), "log_gamma(1001)");
}

void test_beta_binomial() {
  // Exact binomial sum: P(X >= 3) for n = 10, p = 0.1 is 0.0701908.
  require(near(beta_binomial_survival(3U, 10U, 0.1, 0.0), 0.0701908, 1e-6), "binomial limit");
  require(beta_binomial_survival(0U, 10U, 0.1, 0.1) == 1.0, "k = 0 is certain");
  require(beta_binomial_survival(11U, 10U, 0.1, 0.1) == 0.0, "k > n is impossible");
  // With rho = 0.5 and mean 0.5, Beta(1, 1) makes K uniform on 0..n.
  require(near(beta_binomial_survival(3U, 9U, 0.5, 1.0 / 3.0), 0.7, 1e-9), "uniform beta-binomial");
  // Overdispersion raises the tail of a rare event.
  require(beta_binomial_survival(20U, 500U, 0.01, 0.01) >
              beta_binomial_survival(20U, 500U, 0.01, 0.0),
          "overdispersion widens the tail");
  // A probability, even in extreme tails.
  const double tiny = beta_binomial_survival(5000U, 10000U, 1e-4, 1e-4);
  require(tiny >= 0.0 && tiny < 1e-100, "deep tail stays a probability");
}

void test_multiple_testing() {
  // Hand-computed Benjamini-Hochberg example.
  const auto adjusted = benjamini_hochberg({0.01, 0.04, 0.03, 0.2}, 4U);
  require(near(adjusted[0], 0.04, 1e-12), "BH rank 1");
  require(near(adjusted[1], 0.16 / 3.0, 1e-12), "BH rank 3");
  require(near(adjusted[2], 0.16 / 3.0, 1e-12), "BH keeps adjusted values monotone");
  require(near(adjusted[3], 0.2, 1e-12), "BH rank 4");
  // Hypotheses that were never evaluated still count against the discoveries.
  const auto wider = benjamini_hochberg({0.001}, 1000U);
  require(near(wider[0], 1.0, 1e-12), "BH counts untested hypotheses");
  // Fisher's method: for p1 = p2 = 0.05 the combined value is 0.0175.
  require(near(fisher_combined(0.05, 0.05), 0.0175, 5e-4), "Fisher combination");
  require(near(fisher_combined(1.0, 1.0), 1.0, 1e-12), "Fisher combination of nulls");
  // Fisher exact test, R fisher.test(matrix(c(3,1,1,3),2)) gives 0.4857143.
  require(near(fisher_exact_two_sided(3U, 1U, 1U, 3U), 0.4857143, 1e-6), "Fisher exact");
}

void test_robust_fit() {
  std::vector<SiteCount> sites(2000U, SiteCount{1U, 1000U});
  const auto clean = fit_beta_binomial(sites, 0.01, 1e-6, 1e-4);
  require(near(clean.mean, 0.001, 1e-6), "fit recovers the rate");
  // Ten true variants among two thousand sites must not move the error model.
  for (std::size_t index = 0U; index < 10U; ++index) {
    sites[index] = SiteCount{200U, 1000U};
  }
  const auto contaminated = fit_beta_binomial(sites, 0.01, 1e-6, 1e-4);
  require(near(contaminated.mean, 0.001, 1e-4), "true variants do not inflate the mean");
  require(contaminated.rho <= 1e-3, "true variants do not inflate the dispersion");
  // In a sparse class, a true variant must not explain itself as error.
  std::vector<SiteCount> sparse(6U, SiteCount{1U, 1000U});
  sparse[0] = SiteCount{120U, 1000U};
  const auto fit = fit_beta_binomial(sparse, 0.01, 1e-6, 1e-4);
  require(fit.mean < 0.005, "sparse class ignores its outlier");
}

void test_empirical_tail() {
  // Pareto rates with alpha = 2: the quantile function by rank, over 1000 of the
  // 50000 sites, so the top one percent used by the fit lies inside them.
  std::vector<double> rates;
  for (std::size_t rank = 1U; rank <= 1000U; ++rank) {
    rates.push_back(0.01 * std::pow(1000.0 / static_cast<double>(rank), 0.5));
  }
  const auto tail = fit_empirical_tail(rates, 49000U, 1e-4);
  require(tail.alpha > 1.5 && tail.alpha < 3.0, "Hill estimate is near the true exponent");
  require(empirical_tail_probability(tail, tail.threshold / 2.0) == 1.0, "rates under the tail are not tested");
  require(empirical_tail_probability(tail, 0.5) < empirical_tail_probability(tail, 0.05),
          "the tail falls with the rate");
  const auto clean = fit_empirical_tail({0.001}, 16000U, 5e-4);
  require(empirical_tail_probability(clean, 0.05) < 1e-6, "a clean context reaches significance");
}

void test_rank_sum() {
  const std::array<std::uint32_t, 4> same_a{5U, 5U, 5U, 5U};
  const std::array<std::uint32_t, 4> same_b{50U, 50U, 50U, 50U};
  require(near(rank_sum_p_value(same_a.data(), same_b.data(), 4U, nullptr), 1.0, 1e-9), "equal distributions");
  const std::array<std::uint32_t, 4> low{20U, 0U, 0U, 0U};
  const std::array<std::uint32_t, 4> high{0U, 0U, 0U, 200U};
  double shift = 0.0;
  require(rank_sum_p_value(low.data(), high.data(), 4U, &shift) < 1e-10, "separated distributions");
  require(near(shift, -3.0, 1e-12), "mean shift in bins");
}

void test_interval() {
  const auto narrow = overdispersed_wilson(100.0, 1000.0, 0.0);
  const auto wide = overdispersed_wilson(100.0, 1000.0, 0.01);
  // Wilson score interval for 100/1000: 0.0829 to 0.1203.
  require(near(narrow.first, 0.0829, 5e-4) && near(narrow.second, 0.1203, 5e-4), "Wilson interval around 10%");
  require(wide.second - wide.first > narrow.second - narrow.first, "overdispersion widens the interval");
  const auto zero = overdispersed_wilson(0.0, 100.0, 0.0);
  require(zero.first == 0.0 && zero.second > 0.0, "the interval at zero is one-sided");
}

void test_repeat_context() {
  const std::string sequence = "GATTTTTGCACACACAG";
  const auto repeats = build_repeat_context(sequence);
  require(repeats.homopolymer_run[3] == 5U, "homopolymer run start");
  require(repeats.homopolymer_run[7] == 5U, "homopolymer run end");
  require(repeats.homopolymer_run[2] == 1U, "no run");
  // The sequence ends and starts with G, which on a circle is a run of two.
  require(repeats.homopolymer_run[1] == 2U && repeats.homopolymer_run[17] == 2U, "run of two across the origin");
  require(repeats.unit[11] == 2U && repeats.span[11] >= 8U, "dinucleotide repeat");
  // A run crossing the origin is one run.
  const auto circular = build_repeat_context("AACGTGCA");
  require(circular.homopolymer_run[1] == 3U, "run across the origin");
  require(circular.homopolymer_run[8] == 3U, "run across the origin, far side");
}

void test_indel_normalization() {
  const std::string sequence = "ACGTTTTTACGCATCATCATGG";
  for (std::size_t position = 4U; position <= 8U; ++position) {
    const auto deletion = normalize_circular_indel(sequence, position, 1U, {});
    require(deletion.left == 4U && deletion.right == 8U && deletion.left_allele == "T",
            "every placement in a run gives one deletion");
  }
  for (std::size_t anchor = 3U; anchor <= 8U; ++anchor) {
    const auto insertion = normalize_circular_indel(sequence, anchor, 0U, "T");
    require(insertion.left == 3U && insertion.right == 8U, "every placement in a run gives one insertion");
  }
  const auto repeat = normalize_circular_indel(sequence, 15U, 3U, {});
  const auto shifted = normalize_circular_indel(sequence, 12U, 3U, {});
  require(repeat.left == shifted.left && repeat.left_allele == shifted.left_allele,
          "tandem repeat placements coincide");
  // Across the origin of a circular sequence.
  const auto around = normalize_circular_indel("AACGTTTA", 1U, 1U, {});
  require(around.left == 8U && around.right == 2U, "normalization crosses the origin");
  // A uniform circle has no boundary; the scan has to stop anyway.
  const auto uniform = normalize_circular_indel("AAAAAA", 3U, 1U, {});
  require(uniform.deletion && uniform.left_allele == "A", "uniform circle terminates");
}

PairHmmParameters flat_hmm() {
  PairHmmParameters parameters;
  parameters.deletion_open.fill(1e-3);
  parameters.insertion_open.fill(1e-3);
  return parameters;
}

void test_pair_hmm() {
  const auto parameters = flat_hmm();
  PairHmmScratch scratch;
  const std::string haplotype = "ACGTTGCATCGGATTACAGGCTTAACGGTCATCAGT";
  const std::string quality(haplotype.size(), 'I');
  const double exact = pair_hmm_log10_likelihood(haplotype, quality, haplotype, parameters, 12U, scratch);
  std::string substituted = haplotype;
  substituted[15] = substituted[15] == 'A' ? 'C' : 'A';
  const double mismatch = pair_hmm_log10_likelihood(substituted, quality, haplotype, parameters, 12U, scratch);
  std::string gapped = haplotype;
  gapped.erase(15U, 3U);
  const double deleted = pair_hmm_log10_likelihood(gapped, quality.substr(0, gapped.size()), haplotype,
                                                   parameters, 12U, scratch);
  require(exact > mismatch, "a mismatch costs likelihood");
  require(exact > deleted, "a gap costs likelihood");
  // A low-quality mismatch costs less than a high-quality one.
  const std::string poor(haplotype.size(), '#');
  const double low_quality = pair_hmm_log10_likelihood(substituted, poor, haplotype, parameters, 12U, scratch);
  const double low_quality_exact = pair_hmm_log10_likelihood(haplotype, poor, haplotype, parameters, 12U, scratch);
  require(low_quality_exact - low_quality < exact - mismatch, "base quality discounts a mismatch");
  // Forward likelihood is deterministic and insensitive to scratch contents.
  PairHmmScratch other;
  require(pair_hmm_log10_likelihood(gapped, quality.substr(0, gapped.size()), haplotype, parameters, 12U, other) ==
              deleted,
          "scratch does not leak state");
  // The right allele explains its own read better than the reference does.
  const double reference_on_gap = deleted;
  const double alternate_on_gap =
      pair_hmm_log10_likelihood(gapped, quality.substr(0, gapped.size()), gapped, parameters, 12U, scratch);
  require(alternate_on_gap > reference_on_gap + 1.0, "a read fits its own haplotype");
  // Long homopolymers tolerate slippage better when the model has learned to.
  auto homopolymer_aware = parameters;
  homopolymer_aware.deletion_open.fill(1e-3);
  homopolymer_aware.deletion_open[8] = 0.1;
  const std::string run_hap = "GATCGTACGCAAAAAAAACGTCGATTGCAG";
  std::string run_read = run_hap;
  run_read.erase(12U, 1U);
  const std::string run_quality(run_read.size(), 'I');
  require(pair_hmm_log10_likelihood(run_read, run_quality, run_hap, homopolymer_aware, 12U, scratch) >
              pair_hmm_log10_likelihood(run_read, run_quality, run_hap, parameters, 12U, scratch),
          "context-dependent gap rates apply");
}

void test_allele_fraction() {
  // 30 reads favour the alternate allele by 10^4, 70 favour the reference.
  std::vector<double> reference;
  std::vector<double> alternate;
  for (int index = 0; index < 100; ++index) {
    const bool carries = index < 30;
    reference.push_back(carries ? -14.0 : -10.0);
    alternate.push_back(carries ? -10.0 : -14.0);
  }
  const std::vector<double> weights(100U, 1.0);
  const auto fit = fit_allele_fraction(reference, alternate, weights);
  require(near(fit.fraction, 0.3, 1e-3), "EM recovers the mixture fraction");
  require(fit.lod > 50.0, "the fitted mixture is much more likely than none");
  // No evidence for the allele: the fraction falls to zero and LOD with it.
  const auto none = fit_allele_fraction(std::vector<double>(100U, -10.0),
                                        std::vector<double>(100U, -14.0), weights);
  require(none.fraction < 1e-3 && none.lod < 0.1, "no evidence gives no allele");
}

void test_rcrs_translation() {
  const auto reference = load_reference("");
  require(reference.size() == 16569U, "bundled rCRS length");
  const std::vector<std::pair<std::string, std::size_t>> proteins = {
      {"MT-ND1", 318U}, {"MT-ND2", 347U}, {"MT-CO1", 513U}, {"MT-CO2", 227U}, {"MT-ATP8", 68U},
      {"MT-ATP6", 226U}, {"MT-CO3", 261U}, {"MT-ND3", 115U}, {"MT-ND4L", 98U}, {"MT-ND4", 459U},
      {"MT-ND5", 603U}, {"MT-ND6", 174U}, {"MT-CYB", 380U}};
  for (const auto &[gene, length] : proteins) {
    const auto protein = translate_rcrs_gene(reference, gene);
    require(protein.size() == length, gene + " has " + std::to_string(length) + " residues");
    require(!protein.empty() && protein.front() == 'M', gene + " starts with Met");
    require(protein.find('X') == std::string::npos, gene + " has no unknown residue");
  }
  // UniProt P00846 (ATP6) and P00156 (CYB) begin with these residues.
  require(translate_rcrs_gene(reference, "MT-ATP6").rfind("MNENLFASFIAPTILGLPAAV", 0) == 0, "ATP6 N-terminus");
  require(translate_rcrs_gene(reference, "MT-CYB").rfind("MTPMRKTNPLMKLINHSFIDLP", 0) == 0, "CYB N-terminus");
}

std::optional<VariantEffect> effect_in(const std::vector<VariantEffect> &effects, const std::string &gene) {
  for (const auto &effect : effects) {
    if (effect.gene == gene) {
      return effect;
    }
  }
  return std::nullopt;
}

void test_known_pathogenic_annotation() {
  const auto reference = load_reference("");
  // m.8993T>G: MT-ATP6 p.Leu156Arg (NARP/Leigh), and the chain in 8H9S.
  const auto narp = annotate_substitution(reference, 8993U, 'G');
  const auto atp6 = effect_in(narp, "MT-ATP6");
  require(atp6 && atp6->hgvs_protein == "p.Leu156Arg", "m.8993T>G protein change");
  require(atp6 && atp6->consequence == "missense_variant", "m.8993T>G consequence");
  require(atp6 && atp6->grantham == 102, "Leu to Arg Grantham distance is 102");
  require(atp6 && atp6->structure && atp6->structure->pdb_id == "8H9S" &&
              atp6->structure->chain == "N" && atp6->structure->residue == 156U,
          "m.8993T>G maps to 8H9S chain N residue 156");
  // m.11778G>A: MT-ND4 p.Arg340His (LHON).
  const auto lhon = annotate_substitution(reference, 11778U, 'A');
  const auto nd4 = effect_in(lhon, "MT-ND4");
  require(nd4 && nd4->hgvs_protein == "p.Arg340His", "m.11778G>A protein change");
  require(nd4 && nd4->grantham == 29, "Arg to His Grantham distance is 29");
  require(nd4 && nd4->structure && nd4->structure->chain == "r", "m.11778G>A maps to chain r");
  // m.3460G>A: MT-ND1 p.Ala52Thr (LHON).
  const auto nd1 = effect_in(annotate_substitution(reference, 3460U, 'A'), "MT-ND1");
  require(nd1 && nd1->hgvs_protein == "p.Ala52Thr", "m.3460G>A protein change");
  // m.14484T>C: MT-ND6 p.Met64Val, on the light strand (LHON).
  const auto nd6 = effect_in(annotate_substitution(reference, 14484U, 'C'), "MT-ND6");
  require(nd6 && nd6->hgvs_protein == "p.Met64Val", "m.14484T>C on the reverse strand");
  // m.3243A>G lies in tRNA-Leu(UUR).
  const auto trna = effect_in(annotate_substitution(reference, 3243U, 'G'), "MT-TL1");
  require(trna && trna->consequence == "tRNA_variant", "m.3243A>G is a tRNA variant");
  // ATP8 and ATP6 overlap: 8572 lies in both reading frames.
  const auto both = annotate_substitution(reference, 8572U, 'A');
  require(effect_in(both, "MT-ATP8").has_value() && effect_in(both, "MT-ATP6").has_value(),
          "overlapping genes are both annotated");
  // Control region and intergenic bases are not protein changes.
  require(annotate_substitution(reference, 16189U, 'C').front().consequence == "control_region_variant",
          "control region");
  require(annotate_substitution(reference, 5580U, 'C').front().consequence == "intergenic_variant", "intergenic");
  // ND1 starts with ATA. ATG is also an initiator, GTA is not.
  const auto start = effect_in(annotate_substitution(reference, 3309U, 'G'), "MT-ND1");
  require(start && start->consequence == "start_retained_variant", "ATA to ATG keeps the start");
  const auto lost = effect_in(annotate_substitution(reference, 3307U, 'G'), "MT-ND1");
  require(lost && lost->consequence == "start_lost", "ATA to GTA loses the start");
  // ND1 ends with TA, completed to TAA by polyadenylation; CAA is a sense codon.
  const auto stop = effect_in(annotate_substitution(reference, 4261U, 'C'), "MT-ND1");
  require(stop && stop->consequence == "stop_lost", "TAA to CAA loses the stop");
}

void test_small_indel_annotation() {
  const auto reference = load_reference("");
  const auto frameshift = annotate_small_indel(reference, 6100U, 5U, 0U, false);
  const auto co1 = effect_in(frameshift, "MT-CO1");
  require(co1 && co1->consequence == "frameshift_variant", "five-base deletion shifts the frame");
  const auto inframe = annotate_small_indel(reference, 6100U, 9U, 0U, false);
  require(effect_in(inframe, "MT-CO1")->consequence == "inframe_deletion", "nine-base deletion keeps the frame");
  const auto insertion = annotate_small_indel(reference, 9500U, 0U, 3U, true);
  require(effect_in(insertion, "MT-CO3")->consequence == "inframe_insertion", "three-base insertion keeps the frame");
  require(annotate_small_indel(reference, 16100U, 1U, 0U, false).front().consequence == "control_region_variant",
          "indel in the control region");
}

void test_large_deletion_impact() {
  // The common deletion, 8483-13459, truncates ATP8 and ND5.
  const auto common = annotate_large_deletion(8483U, 13459U, 16569U);
  const auto has = [](const std::vector<std::string> &list, const std::string &name) {
    return std::find(list.begin(), list.end(), name) != list.end();
  };
  require(has(common.genes_lost, "MT-ATP6") && has(common.genes_lost, "MT-ND4"), "common deletion loses ATP6 and ND4");
  require(has(common.genes_partial, "MT-ATP8") && has(common.genes_partial, "MT-ND5"), "common deletion truncates two genes");
  require(common.transfer_rnas_affected == 5U, "common deletion removes five tRNAs");
  require(common.major_arc && common.contact_zone, "common deletion is in the major arc and contact zone");
  require(!common.heavy_strand_origin_lost && !common.light_strand_origin_lost, "common deletion keeps both origins");
  require(common.fusion.find("MT-ATP8::MT-ND5") == 0U, "common deletion joins ATP8 to ND5");
  // A minor-arc deletion that removes the light-strand origin.
  const auto minor = annotate_large_deletion(5000U, 6000U, 16569U);
  require(minor.light_strand_origin_lost && !minor.major_arc, "light-strand origin loss");
  // A deletion spanning the origin of the linearized reference.
  const auto around = annotate_large_deletion(16500U, 200U, 16569U);
  require(around.heavy_strand_origin_lost, "origin-spanning deletion reaches the heavy-strand origin");
}

// A deletion inside a direct repeat is one event wherever the aligner puts the
// gap. minimap2 puts the common deletion at the 5' end of its 13-bp repeat
// (8470-13446); other aligners may use the 3' end (8483-13459). Both must give
// the call at the 3' placement with the full length.
void test_deletion_placement_in_repeat() {
  std::ifstream in(bundled_rcrs_path());
  const auto reference = parse_reference_fasta(in, bundled_rcrs_path());
  for (const std::size_t first : {std::size_t{8470U}, std::size_t{8483U}}) {
    std::mt19937 random(3U);
    std::string sam = "@SQ\tSN:NC_012920.1\tLN:16569\n";
    const auto add = [&](const std::string &name, const std::size_t index, const std::size_t start,
                         const std::string &cigar, const std::string &sequence) {
      sam += name + std::to_string(index) + "\t" + (index % 2U == 0U ? "0" : "16") +
             "\tNC_012920.1\t" + std::to_string(start) + "\t60\t" + cigar + "\t*\t0\t0\t" +
             sequence + "\t" + std::string(sequence.size(), 'I') + "\n";
    };
    for (std::size_t index = 0U; index < 120U; ++index) {
      const std::size_t start = 7483U - random() % 400U;
      add("d", index, start, std::to_string(first - start) + "M4977D1000M",
          reference.substr(start - 1U, first - start) +
              reference.substr(first - 1U + 4977U, 1000U));
    }
    for (std::size_t index = 0U; index < 300U; ++index) {
      const std::size_t start = 1U + random() % 14000U;
      add("w", index, start, "2500M", reference.substr(start - 1U, 2500U));
    }
    const auto path = std::filesystem::temp_directory_path() /
                      ("mito-precision-repeat-" + std::to_string(first) + ".sam");
    std::ofstream(path, std::ios::binary) << sam;
    std::string json;
    try {
      json = AnalysisEngine().analyze(path.string(), "", {});
    } catch (const std::exception &error) {
      require(false, std::string("deletion placement analysis failed: ") + error.what());
    }
    std::filesystem::remove(path);
    require(json.find("\"pdeletion:8483-13459\"") != std::string::npos,
            "common deletion is called at 8483-13459 with the gap at " + std::to_string(first));
  }
}

} // namespace

int main() {
  test_log_gamma();
  test_beta_binomial();
  test_multiple_testing();
  test_robust_fit();
  test_empirical_tail();
  test_rank_sum();
  test_interval();
  test_repeat_context();
  test_indel_normalization();
  test_pair_hmm();
  test_allele_fraction();
  test_rcrs_translation();
  test_known_pathogenic_annotation();
  test_small_indel_annotation();
  test_large_deletion_impact();
  test_deletion_placement_in_repeat();
  if (failures != 0) {
    std::cerr << failures << " precision test assertion(s) failed\n";
    return 1;
  }
  std::cout << "precision tests passed\n";
  return 0;
}
