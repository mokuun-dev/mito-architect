#include "detail/precision.hpp"

namespace mito::detail {

[[nodiscard]] double log_gamma(double value) noexcept {
  // Lanczos approximation (g = 7, n = 9), accurate to about 15 digits.
  static constexpr std::array<double, 9> coefficients = {
      0.99999999999980993,  676.5203681218851,     -1259.1392167224028,
      771.32342877765313,   -176.61502916214059,   12.507343278686905,
      -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
  constexpr double half_log_two_pi = 0.91893853320467274178;
  if (value < 0.5) {
    constexpr double pi = 3.14159265358979323846;
    return std::log(pi / std::abs(std::sin(pi * value))) -
           log_gamma(1.0 - value);
  }
  value -= 1.0;
  double sum = coefficients[0];
  for (std::size_t index = 1U; index < coefficients.size(); ++index) {
    sum += coefficients[index] / (value + static_cast<double>(index));
  }
  const double t = value + 7.5;
  return half_log_two_pi + (value + 0.5) * std::log(t) - t + std::log(sum);
}

namespace {

[[nodiscard]] double log_beta(const double a, const double b) noexcept {
  return log_gamma(a) + log_gamma(b) - log_gamma(a + b);
}

[[nodiscard]] double log_choose(const double n, const double k) noexcept {
  return log_gamma(n + 1.0) - log_gamma(k + 1.0) - log_gamma(n - k + 1.0);
}

[[nodiscard]] double clamp_probability(const double value) noexcept {
  return value < 0.0 ? 0.0 : (value > 1.0 ? 1.0 : value);
}

} // namespace

[[nodiscard]] double beta_binomial_survival(const std::size_t k,
                                            const std::size_t n, double mean,
                                            const double rho) noexcept {
  if (k == 0U) {
    return 1.0;
  }
  if (k > n) {
    return 0.0;
  }
  mean = std::clamp(mean, 1e-12, 1.0 - 1e-12);
  const auto trials = static_cast<double>(n);
  const auto first = static_cast<double>(k);
  double log_term = 0.0;
  // Successive pmf ratios avoid re-evaluating log-gamma for every term.
  double alpha = 0.0;
  double beta = 0.0;
  const bool binomial = rho <= 1e-12;
  if (binomial) {
    log_term = log_choose(trials, first) + first * std::log(mean) +
               (trials - first) * std::log1p(-mean);
  } else {
    const double concentration = (1.0 - rho) / rho;
    alpha = mean * concentration;
    beta = (1.0 - mean) * concentration;
    log_term = log_choose(trials, first) +
               log_beta(first + alpha, trials - first + beta) -
               log_beta(alpha, beta);
  }
  double term = std::exp(log_term);
  double sum = term;
  for (std::size_t index = k; index < n; ++index) {
    const auto i = static_cast<double>(index);
    const double ratio =
        binomial ? ((trials - i) / (i + 1.0)) * (mean / (1.0 - mean))
                 : ((trials - i) / (i + 1.0)) *
                       ((i + alpha) / (trials - i - 1.0 + beta));
    term *= ratio;
    sum += term;
    if (ratio < 1.0 && term < sum * 1e-14) {
      break;
    }
  }
  return clamp_probability(sum);
}

[[nodiscard]] double fisher_combined(const double first,
                                     const double second) noexcept {
  // For two p-values X = -2 ln(p1 p2) ~ chi2(4), whose survival function is
  // exp(-x / 2) (1 + x / 2), which simplifies to p (1 - ln p) for p = p1 p2.
  const double product = std::max(first * second, 1e-300);
  return clamp_probability(product * (1.0 - std::log(product)));
}

[[nodiscard]] double fisher_exact_two_sided(const std::size_t a,
                                            const std::size_t b,
                                            const std::size_t c,
                                            const std::size_t d) noexcept {
  const auto row_one = static_cast<double>(a + b);
  const auto row_two = static_cast<double>(c + d);
  const auto column_one = static_cast<double>(a + c);
  const auto total = row_one + row_two;
  if (total == 0.0) {
    return 1.0;
  }
  const double log_denominator = log_choose(total, column_one);
  const auto log_probability = [&](const double x) {
    return log_choose(row_one, x) + log_choose(row_two, column_one - x) -
           log_denominator;
  };
  const auto low =
      static_cast<std::size_t>(std::max(0.0, column_one - row_two));
  const auto high = static_cast<std::size_t>(std::min(row_one, column_one));
  const double observed = log_probability(static_cast<double>(a));
  // The relative tolerance treats tables with equal probability as equal.
  const double threshold = observed + 1e-9;
  double probability = 0.0;
  for (std::size_t x = low; x <= high; ++x) {
    const double log_value = log_probability(static_cast<double>(x));
    if (log_value <= threshold) {
      probability += std::exp(log_value);
    }
  }
  return clamp_probability(probability);
}

[[nodiscard]] std::vector<double>
benjamini_hochberg(const std::vector<double> &p_values,
                   const std::size_t total_tests) {
  std::vector<std::size_t> order(p_values.size());
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::stable_sort(order.begin(), order.end(),
                   [&](const std::size_t lhs, const std::size_t rhs) {
                     return p_values[lhs] < p_values[rhs];
                   });
  const auto tests =
      static_cast<double>(std::max(total_tests, p_values.size()));
  std::vector<double> adjusted(p_values.size(), 1.0);
  double running = 1.0;
  for (std::size_t rank = order.size(); rank > 0U; --rank) {
    const auto index = order[rank - 1U];
    running =
        std::min(running, p_values[index] * tests / static_cast<double>(rank));
    adjusted[index] = clamp_probability(running);
  }
  return adjusted;
}

[[nodiscard]] BetaBinomialFit fit_beta_binomial(std::vector<SiteCount> sites,
                                                const double trim_fraction,
                                                const double mean_floor,
                                                const double rho_floor) {
  std::erase_if(sites, [](const SiteCount &site) { return site.n == 0U; });
  BetaBinomialFit fit;
  fit.sites = sites.size();
  fit.mean = mean_floor;
  fit.rho = rho_floor;
  if (sites.empty()) {
    return fit;
  }
  std::stable_sort(sites.begin(), sites.end(),
                   [](const SiteCount &lhs, const SiteCount &rhs) {
                     // Cross-multiplied comparison of k / n, exact in integers.
                     return static_cast<std::uint64_t>(lhs.k) * rhs.n <
                            static_cast<std::uint64_t>(rhs.k) * lhs.n;
                   });
  const auto keep = static_cast<std::size_t>(
      std::ceil(static_cast<double>(sites.size()) * (1.0 - trim_fraction)));
  sites.resize(std::max<std::size_t>(1U, keep));

  // Variance is dominated by its largest terms, and a sparse class can hold a
  // true variant whose own rate would then excuse it. Sites far above a location
  // estimated from the lower 90% are set aside before the moment fit.
  if (sites.size() >= 5U) {
    const std::size_t lower = std::max<std::size_t>(1U, sites.size() * 9U / 10U);
    double lower_n = 0.0;
    double lower_k = 0.0;
    for (std::size_t index = 0U; index < lower; ++index) {
      lower_n += sites[index].n;
      lower_k += sites[index].k;
    }
    const double location = lower_k / lower_n;
    const auto beyond_fence = [&](const SiteCount &site) {
      const double n = site.n;
      const double fence = 6.0 * location +
                           4.0 * std::sqrt(location * (1.0 - location) / n) + 1.0 / n;
      return static_cast<double>(site.k) / n > fence;
    };
    std::erase_if(sites, beyond_fence);
  }

  double total_n = 0.0;
  double total_k = 0.0;
  for (const auto &site : sites) {
    total_n += site.n;
    total_k += site.k;
  }
  const double mean = total_k / total_n;
  fit.mean = std::max(mean, mean_floor);
  if (sites.size() < 2U || mean <= 0.0 || mean >= 1.0) {
    return fit;
  }
  double between = 0.0;
  double binomial_part = 0.0;
  double dispersion_part = 0.0;
  for (const auto &site : sites) {
    const double n = site.n;
    const double weight = n / total_n;
    const double rate = static_cast<double>(site.k) / n;
    between += weight * (rate - mean) * (rate - mean);
    binomial_part += weight * (1.0 - weight) / n;
    dispersion_part += weight * (1.0 - weight) * (1.0 - 1.0 / n);
  }
  const double variance = mean * (1.0 - mean);
  if (dispersion_part > 0.0 && variance > 0.0) {
    const double rho =
        (between - variance * binomial_part) / (variance * dispersion_part);
    fit.rho = std::clamp(rho, rho_floor, 0.5);
  }
  return fit;
}

[[nodiscard]] EmpiricalTail fit_empirical_tail(std::vector<double> nonzero_rates,
                                               const std::size_t zero_sites,
                                               const double minimum_threshold) {
  EmpiricalTail tail;
  tail.threshold = minimum_threshold;
  const std::size_t total = nonzero_rates.size() + zero_sites;
  if (total == 0U || nonzero_rates.empty()) {
    return tail;
  }
  std::sort(nonzero_rates.begin(), nonzero_rates.end(), std::greater<>());
  constexpr double kTailFraction = 0.01;
  constexpr double kTrueVariantFraction = 0.002;
  const auto tail_sites = std::max<std::size_t>(
      10U, static_cast<std::size_t>(std::ceil(static_cast<double>(total) * kTailFraction)));
  if (tail_sites >= nonzero_rates.size()) {
    // Under one percent of sites carry any observation: this is a clean context,
    // and the default steep tail starts at the smallest detectable rate.
    tail.tail_fraction =
        static_cast<double>(nonzero_rates.size()) / static_cast<double>(total);
    return tail;
  }
  tail.tail_fraction = static_cast<double>(tail_sites) / static_cast<double>(total);
  tail.threshold = std::max(minimum_threshold, nonzero_rates[tail_sites]);
  // The largest rates may be true variants, so they only count as exceeding the
  // largest kept rate: the maximum likelihood estimate for a right-censored
  // Pareto sample. Dropping them outright would bias the exponent upward and
  // make the noise tail look lighter than it is.
  const auto censored = std::min(
      tail_sites - 1U,
      static_cast<std::size_t>(std::ceil(static_cast<double>(total) * kTrueVariantFraction)));
  double log_excess = 0.0;
  for (std::size_t index = censored; index < tail_sites; ++index) {
    log_excess += std::log(std::max(nonzero_rates[index], tail.threshold) / tail.threshold);
  }
  log_excess += static_cast<double>(censored) *
                std::log(std::max(nonzero_rates[censored], tail.threshold) / tail.threshold);
  if (log_excess > 0.0) {
    tail.alpha = std::clamp(static_cast<double>(tail_sites - censored) / log_excess, 1.5, 12.0);
  }
  return tail;
}

[[nodiscard]] double empirical_tail_probability(const EmpiricalTail &tail,
                                                const double rate) noexcept {
  if (rate <= tail.threshold || tail.threshold <= 0.0) {
    return rate <= tail.threshold ? 1.0 : 0.0;
  }
  return clamp_probability(tail.tail_fraction * std::pow(rate / tail.threshold, -tail.alpha));
}

[[nodiscard]] std::pair<double, double>
overdispersed_wilson(const double k, const double n,
                     const double rho) noexcept {
  if (n <= 0.0) {
    return {0.0, 1.0};
  }
  constexpr double z = 1.959963984540054;
  // The design effect turns correlated observations into independent ones.
  const double effective = n / (1.0 + (n - 1.0) * std::max(rho, 0.0));
  const double proportion = std::clamp(k / n, 0.0, 1.0);
  const double z2 = z * z;
  const double denominator = 1.0 + z2 / effective;
  const double center = (proportion + z2 / (2.0 * effective)) / denominator;
  const double margin =
      z *
      std::sqrt(proportion * (1.0 - proportion) / effective +
                z2 / (4.0 * effective * effective)) /
      denominator;
  // Without this the bound excludes the observed value at the extremes, where
  // rounding in center - margin can exceed the proportion.
  return {std::max(0.0, std::min(center - margin, proportion)),
          std::min(1.0, std::max(center + margin, proportion))};
}

[[nodiscard]] double rank_sum_p_value(const std::uint32_t *first,
                                      const std::uint32_t *second,
                                      const std::size_t bins,
                                      double *mean_shift) noexcept {
  double count_first = 0.0;
  double count_second = 0.0;
  double sum_first = 0.0;
  double sum_second = 0.0;
  for (std::size_t bin = 0U; bin < bins; ++bin) {
    count_first += first[bin];
    count_second += second[bin];
    sum_first += static_cast<double>(first[bin]) * static_cast<double>(bin);
    sum_second += static_cast<double>(second[bin]) * static_cast<double>(bin);
  }
  if (mean_shift != nullptr) {
    *mean_shift = count_first > 0.0 && count_second > 0.0
                      ? sum_first / count_first - sum_second / count_second
                      : 0.0;
  }
  if (count_first == 0.0 || count_second == 0.0) {
    return 1.0;
  }
  const double total = count_first + count_second;
  double u = 0.0;
  double tie_term = 0.0;
  double below_second = 0.0;
  for (std::size_t bin = 0U; bin < bins; ++bin) {
    const double in_first = first[bin];
    const double in_second = second[bin];
    u += in_first * (below_second + in_second / 2.0);
    below_second += in_second;
    const double tied = in_first + in_second;
    tie_term += tied * tied * tied - tied;
  }
  const double mean = count_first * count_second / 2.0;
  const double variance = count_first * count_second / 12.0 *
                          ((total + 1.0) - tie_term / (total * (total - 1.0)));
  if (variance <= 0.0) {
    return 1.0;
  }
  const double z = std::abs(u - mean) / std::sqrt(variance);
  return clamp_probability(std::erfc(z / std::sqrt(2.0)));
}

[[nodiscard]] RepeatContext build_repeat_context(const std::string &sequence) {
  const std::size_t length = sequence.size();
  RepeatContext context;
  context.homopolymer_run.assign(length + 1U, 1U);
  context.unit.assign(length + 1U, 1U);
  context.span.assign(length + 1U, 1U);
  if (length == 0U) {
    return context;
  }
  const auto base = [&](const std::size_t index) {
    return sequence[index % length];
  };
  constexpr std::size_t maximum_unit = 6U;
  // Scanning from a start that is not inside a repeat keeps a repeat crossing
  // the origin whole; a sequence of one repeated base has no such start.
  for (std::size_t unit = 1U; unit <= maximum_unit && unit < length; ++unit) {
    const std::size_t minimum_span =
        unit == 1U ? std::size_t{2} : std::max<std::size_t>(6U, unit * 2U);
    std::size_t scan_start = 0U;
    while (scan_start < length && base(scan_start) == base(scan_start + unit)) {
      ++scan_start;
    }
    if (scan_start == length) {
      for (std::size_t position = 1U; position <= length; ++position) {
        if (length > context.span[position]) {
          context.unit[position] = static_cast<std::uint8_t>(unit);
          context.span[position] = static_cast<std::uint16_t>(
              std::min<std::size_t>(length, 65535U));
        }
        if (unit == 1U) {
          context.homopolymer_run[position] = context.span[position];
        }
      }
      continue;
    }
    std::size_t offset = 1U;
    while (offset <= length) {
      const std::size_t start = scan_start + offset;
      std::size_t matched = 0U;
      while (matched < length && base(start + matched) == base(start + matched + unit)) {
        ++matched;
      }
      if (matched == 0U) {
        ++offset;
        continue;
      }
      const std::size_t span = matched + unit;
      if (span >= minimum_span) {
        for (std::size_t index = 0U; index < span; ++index) {
          const std::size_t position = (start + index) % length + 1U;
          if (span > context.span[position]) {
            context.unit[position] = static_cast<std::uint8_t>(unit);
            context.span[position] = static_cast<std::uint16_t>(
                std::min<std::size_t>(span, 65535U));
          }
          if (unit == 1U) {
            context.homopolymer_run[position] = static_cast<std::uint16_t>(
                std::min<std::size_t>(span, 65535U));
          }
        }
      }
      offset += matched;
    }
  }
  return context;
}

[[nodiscard]] NormalizedCircularIndel
normalize_circular_indel(const std::string &reference,
                         const std::size_t position,
                         const std::size_t deleted_length,
                         const std::string_view inserted) {
  const std::size_t length = reference.size();
  const auto at = [&](const std::size_t one_based) -> char {
    return reference[(one_based - 1U) % length];
  };
  const auto previous = [&](const std::size_t one_based) {
    return one_based > 1U ? one_based - 1U : length;
  };
  const auto next = [&](const std::size_t one_based) {
    return one_based < length ? one_based + 1U : 1U;
  };
  NormalizedCircularIndel result;
  if (deleted_length > 0U) {
    result.deletion = true;
    std::size_t left = position;
    std::size_t leftward_steps = 0U;
    // A rotation by the whole circle would never stop on a uniform reference.
    while (leftward_steps < length &&
           at(previous(left)) == at(left + deleted_length - 1U)) {
      left = previous(left);
      ++leftward_steps;
    }
    std::size_t right = position;
    std::size_t rightward_steps = 0U;
    while (rightward_steps < length && at(right) == at(right + deleted_length)) {
      right = next(right);
      ++rightward_steps;
    }
    result.left = left;
    result.right = right;
    result.left_allele.reserve(deleted_length);
    result.right_allele.reserve(deleted_length);
    for (std::size_t index = 0U; index < deleted_length; ++index) {
      result.left_allele.push_back(at(left + index));
      result.right_allele.push_back(at(right + index));
    }
    return result;
  }
  std::string left_allele(inserted);
  std::size_t left = position;
  std::size_t leftward_steps = 0U;
  while (leftward_steps < length && !left_allele.empty() &&
         at(left) == left_allele.back()) {
    left_allele.insert(left_allele.begin(), at(left));
    left_allele.pop_back();
    left = previous(left);
    ++leftward_steps;
  }
  std::string right_allele(inserted);
  std::size_t right = position;
  std::size_t rightward_steps = 0U;
  while (rightward_steps < length && !right_allele.empty() &&
         at(next(right)) == right_allele.front()) {
    right = next(right);
    right_allele.erase(right_allele.begin());
    right_allele.push_back(at(right));
    ++rightward_steps;
  }
  result.left = left;
  result.right = right;
  result.left_allele = std::move(left_allele);
  result.right_allele = std::move(right_allele);
  return result;
}

namespace {

// Error probability by Phred score; a power per read base would dominate the HMM.
const std::array<double, 94> kPhredError = [] {
  std::array<double, 94> table{};
  for (std::size_t quality = 0U; quality < table.size(); ++quality) {
    table[quality] = std::pow(10.0, -static_cast<double>(quality) / 10.0);
  }
  return table;
}();

} // namespace

[[nodiscard]] double pair_hmm_log10_likelihood(
    const std::string_view read, const std::string_view qualities,
    const std::string_view haplotype, const PairHmmParameters &parameters,
    const std::size_t band_half_width, PairHmmScratch &scratch) {
  const std::size_t read_length = read.size();
  const std::size_t hap_length = haplotype.size();
  if (read_length == 0U || hap_length == 0U) {
    return 0.0;
  }
  auto &run = scratch.hap_run;
  run.resize(hap_length);
  for (std::size_t index = 0U; index < hap_length;) {
    std::size_t stop = index + 1U;
    while (stop < hap_length && haplotype[stop] == haplotype[index]) {
      ++stop;
    }
    const auto bounded = static_cast<std::uint8_t>(
        std::min(stop - index, PairHmmParameters::kMaximumRun));
    std::fill(run.begin() + static_cast<std::ptrdiff_t>(index),
              run.begin() + static_cast<std::ptrdiff_t>(stop), bounded);
    index = stop;
  }

  const std::size_t columns = hap_length + 1U;
  for (auto *vector :
       {&scratch.match, &scratch.insertion, &scratch.deletion,
        &scratch.next_match, &scratch.next_insertion, &scratch.next_deletion}) {
    vector->assign(columns, 0.0);
  }
  auto &previous_match = scratch.match;
  auto &previous_insertion = scratch.insertion;
  auto &previous_deletion = scratch.deletion;
  auto &current_match = scratch.next_match;
  auto &current_insertion = scratch.next_insertion;
  auto &current_deletion = scratch.next_deletion;

  // Row zero is a uniform free start along the haplotype, as in GATK, because
  // the read window boundary is only approximately placed on the haplotype.
  const double start = 1.0 / static_cast<double>(hap_length);
  for (std::size_t column = 0U; column <= std::min(hap_length, band_half_width);
       ++column) {
    previous_deletion[column] = start;
  }

  const double extend = parameters.gap_extend;
  const double close = 1.0 - extend;
  double log_scale = 0.0;
  std::size_t previous_low = 0U;
  std::size_t previous_high = std::min(hap_length, band_half_width);
  for (std::size_t row = 1U; row <= read_length; ++row) {
    const std::size_t low = row > band_half_width ? row - band_half_width : 1U;
    const std::size_t high = std::min(hap_length, row + band_half_width);
    if (low > high) {
      // The read is longer than the band can reach on this haplotype.
      break;
    }
    const char base = read[row - 1U];
    const auto quality_value =
        row - 1U < qualities.size()
            ? std::clamp(static_cast<int>(qualities[row - 1U]) - 33, 0, 93)
            : 20;
    const double error =
        std::max(parameters.minimum_error, kPhredError[static_cast<std::size_t>(quality_value)]);
    const double equal = 1.0 - error;
    const double unequal = error / 3.0;
    current_deletion[low - 1U] = 0.0;
    current_match[low - 1U] = 0.0;
    current_insertion[low - 1U] = 0.0;
    double row_maximum = 0.0;
    for (std::size_t column = low; column <= high; ++column) {
      const std::size_t hap_index = column - 1U;
      const std::size_t bounded_run = run[hap_index];
      const double deletion_open = parameters.deletion_open[bounded_run];
      const double insertion_open = parameters.insertion_open[bounded_run];
      const double stay = 1.0 - deletion_open - insertion_open;
      const double emission = haplotype[hap_index] == base ? equal : unequal;
      const bool diagonal_in_band =
          column - 1U >= previous_low && column - 1U <= previous_high;
      const bool vertical_in_band =
          column >= previous_low && column <= previous_high;
      const double from_diagonal =
          diagonal_in_band
              ? stay * previous_match[column - 1U] +
                    close * (previous_insertion[column - 1U] +
                             previous_deletion[column - 1U])
              : 0.0;
      const double match = emission * from_diagonal;
      const double insertion =
          vertical_in_band ? insertion_open * previous_match[column] +
                                 extend * previous_insertion[column]
                           : 0.0;
      const double deletion = deletion_open * current_match[column - 1U] +
                              extend * current_deletion[column - 1U];
      current_match[column] = match;
      current_insertion[column] = insertion;
      current_deletion[column] = deletion;
      row_maximum = std::max({row_maximum, match, insertion, deletion});
    }
    if (row_maximum <= 0.0) {
      return log_scale + std::log10(1e-300);
    }
    // Rescaling each row keeps the forward values inside double range.
    if (row_maximum < 1e-100) {
      const double inverse = 1.0 / row_maximum;
      for (std::size_t column = low - 1U; column <= high; ++column) {
        current_match[column] *= inverse;
        current_insertion[column] *= inverse;
        current_deletion[column] *= inverse;
      }
      log_scale += std::log10(row_maximum);
    }
    std::swap(previous_match, current_match);
    std::swap(previous_insertion, current_insertion);
    std::swap(previous_deletion, current_deletion);
    previous_low = low - 1U;
    previous_high = high;
  }
  double total = 0.0;
  for (std::size_t column = previous_low; column <= previous_high; ++column) {
    total += previous_match[column] + previous_insertion[column];
  }
  return log_scale + std::log10(std::max(total, 1e-300));
}

[[nodiscard]] AlleleFractionFit
fit_allele_fraction(const std::vector<double> &log10_reference,
                    const std::vector<double> &log10_alternate,
                    const std::vector<double> &weights) {
  AlleleFractionFit fit;
  const std::size_t reads = log10_reference.size();
  if (reads == 0U) {
    return fit;
  }
  std::vector<double> ratio(reads);
  double total_weight = 0.0;
  for (std::size_t index = 0U; index < reads; ++index) {
    // The ratio is capped so one misplaced read cannot set the likelihood.
    const double log_ratio = std::clamp(
        log10_alternate[index] - log10_reference[index], -30.0, 30.0);
    ratio[index] = std::pow(10.0, log_ratio);
    total_weight += weights[index];
  }
  if (total_weight <= 0.0) {
    return fit;
  }
  double fraction = 0.1;
  for (int iteration = 0; iteration < 200; ++iteration) {
    double expected = 0.0;
    for (std::size_t index = 0U; index < reads; ++index) {
      const double numerator = fraction * ratio[index];
      expected += weights[index] * numerator / (numerator + 1.0 - fraction);
    }
    const double updated = expected / total_weight;
    const bool converged = std::abs(updated - fraction) < 1e-7;
    fraction = updated;
    if (converged) {
      break;
    }
  }
  fit.fraction = fraction;
  for (std::size_t index = 0U; index < reads; ++index) {
    const double numerator = fraction * ratio[index];
    const double mixture = numerator + 1.0 - fraction;
    fit.effective_alternate += weights[index] * numerator / mixture;
    fit.lod += weights[index] * std::log10(std::max(mixture, 1e-300));
  }
  return fit;
}

} // namespace mito::detail
