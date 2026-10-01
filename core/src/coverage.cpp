#include "detail/pipeline.hpp"

namespace mito::detail {

void add_circular_coverage_span(std::vector<std::int64_t> &difference,
                                std::uint64_t &uniform_depth, std::size_t start,
                                std::size_t length) {
  const std::size_t reference_length = difference.size() - 1U;
  if (reference_length == 0 || length == 0) {
    return;
  }

  const auto complete_cycles =
      static_cast<std::uint64_t>(length / reference_length);
  if (complete_cycles >
      std::numeric_limits<std::uint64_t>::max() - uniform_depth) {
    throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                        "coverage depth overflow");
  }
  uniform_depth += complete_cycles;
  const std::size_t remainder = length % reference_length;
  if (remainder == 0) {
    return;
  }

  const std::size_t zero_based_start = (start - 1U) % reference_length;
  const std::size_t linear_end = zero_based_start + remainder;
  ++difference[zero_based_start];
  if (linear_end <= reference_length) {
    --difference[linear_end];
    return;
  }

  --difference[reference_length];
  ++difference[0];
  --difference[linear_end - reference_length];
}

[[nodiscard]] CoverageResult
compute_coverage(const std::vector<ReadRecord> &reads,
                 const std::vector<ReadFeature> &features,
                 const AnalysisConfig &config, std::size_t reference_length,
                 std::size_t bin_count) {
  CoverageResult result;
  if (reference_length == 0 || bin_count == 0) {
    return result;
  }

  std::vector<std::int64_t> difference(reference_length + 1U, 0);
  std::uint64_t uniform_depth = 0;
  for (std::size_t read_index = 0; read_index < reads.size(); ++read_index) {
    throw_if_cancelled(config);
    const auto &read = reads[read_index];
    if (read_index >= features.size() || features[read_index].filtered_numt ||
        (read.flags & 0x4U) != 0U || read.reference_start == 0 ||
        read.reference_name.empty() || read.reference_name == "*" ||
        read.cigar_operations.empty() ||
        looks_like_nuclear_contig(read.reference_name)) {
      continue;
    }

    std::size_t reference_cursor =
        ((read.reference_start - 1U) % reference_length) + 1U;
    for (const auto &operation : read.cigar_operations) {
      if (operation.code == 'M' || operation.code == '=' ||
          operation.code == 'X') {
        add_circular_coverage_span(difference, uniform_depth, reference_cursor,
                                   operation.length);
      }
      if (operation.code == 'M' || operation.code == '=' ||
          operation.code == 'X' || operation.code == 'D' ||
          operation.code == 'N') {
        reference_cursor =
            ((reference_cursor - 1U + (operation.length % reference_length)) %
             reference_length) +
            1U;
      }
    }
  }

  std::vector<std::uint64_t> site_depth(reference_length, 0);
  std::int64_t running_delta = 0;
  long double total_depth = 0.0L;
  std::size_t sites_gt20 = 0;
  for (std::size_t i = 0; i < reference_length; ++i) {
    running_delta += difference[i];
    if (running_delta < 0) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "internal coverage accumulator underflow");
    }
    const auto variable_depth = static_cast<std::uint64_t>(running_delta);
    if (variable_depth >
        std::numeric_limits<std::uint64_t>::max() - uniform_depth) {
      throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                          "coverage depth overflow");
    }
    const auto depth = uniform_depth + variable_depth;
    if (depth > std::numeric_limits<std::size_t>::max()) {
      throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                          "coverage depth exceeds platform size limit");
    }
    site_depth[i] = depth;
    total_depth += static_cast<long double>(depth);
    result.max_depth =
        std::max(result.max_depth, static_cast<std::size_t>(depth));
    if (depth > 20U) {
      ++sites_gt20;
    }
  }

  result.mean_depth = static_cast<double>(
      total_depth / static_cast<long double>(reference_length));
  result.pct_sites_gt20x = (static_cast<double>(sites_gt20) * 100.0) /
                           static_cast<double>(reference_length);

  auto &bins = result.bins;
  bins.reserve(bin_count);
  const std::size_t bin_width =
      std::max<std::size_t>(1, reference_length / bin_count);
  for (std::size_t i = 0; i < bin_count; ++i) {
    const std::size_t start = i * bin_width + 1;
    if (start > reference_length) {
      break;
    }
    const std::size_t end =
        i == bin_count - 1 ? reference_length
                           : std::min(reference_length, start + bin_width - 1U);
    const auto first =
        site_depth.begin() + static_cast<std::ptrdiff_t>(start - 1U);
    const auto last = site_depth.begin() + static_cast<std::ptrdiff_t>(end);
    const long double bin_total = std::accumulate(
        first, last, 0.0L, [](long double sum, std::uint64_t depth) {
          return sum + static_cast<long double>(depth);
        });
    const auto site_count = end - start + 1U;
    const auto mean_depth = static_cast<std::size_t>(
        (bin_total / static_cast<long double>(site_count)) + 0.5L);
    bins.push_back({start, end, mean_depth});
  }

  return result;
}

} // namespace mito::detail
