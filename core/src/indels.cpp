#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] std::vector<NormalizedSmallIndel>
call_small_indels_from_alignment(const ReadRecord &read,
                                 const AnalysisConfig &config,
                                 const std::string &reference) {
  std::vector<NormalizedSmallIndel> calls;
  if (reference.empty() || read.reference_start == 0U ||
      read.cigar_operations.empty() ||
      !passes_snp_alignment_filters(read, config) ||
      config.sv_min_length <= 1U) {
    return calls;
  }

  constexpr std::size_t kMaximumSmallIndelLength = 50U;
  const auto maximum_length =
      std::min(kMaximumSmallIndelLength, config.sv_min_length - 1U);
  std::size_t query_cursor = 0U;
  std::size_t reference_cursor =
      ((read.reference_start - 1U) % reference.size()) + 1U;
  for (const auto &operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;
    if (op == 'M' || op == '=' || op == 'X') {
      query_cursor += len;
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
      continue;
    }
    if (op == 'I') {
      if (len <= maximum_length && query_cursor <= read.sequence.size() &&
          len <= read.sequence.size() - query_cursor) {
        std::string inserted;
        inserted.reserve(len);
        std::optional<std::uint8_t> minimum_quality;
        bool callable = true;
        for (std::size_t offset = 0U; offset < len; ++offset) {
          const char base =
              static_cast<char>(std::toupper(static_cast<unsigned char>(
                  read.sequence[query_cursor + offset])));
          const auto quality = phred_quality_at(read, query_cursor + offset);
          if (!is_base(base) || !quality ||
              *quality < config.min_base_quality) {
            callable = false;
            break;
          }
          inserted.push_back(base);
          minimum_quality =
              minimum_quality ? std::min(*minimum_quality, *quality) : *quality;
        }
        if (callable) {
          std::size_t anchor =
              reference_cursor > 1U ? reference_cursor - 1U : reference.size();
          std::size_t shift = 0U;
          while (anchor < reference.size() && !inserted.empty() &&
                 inserted[shift % inserted.size()] == reference[anchor]) {
            ++shift;
            ++anchor;
          }
          if (!inserted.empty() && shift % inserted.size() != 0U) {
            std::rotate(inserted.begin(),
                        inserted.begin() + static_cast<std::ptrdiff_t>(
                                               shift % inserted.size()),
                        inserted.end());
          }
          NormalizedSmallIndel event;
          event.type = "SMALL_INSERTION";
          event.start = anchor;
          event.end = anchor;
          event.length = len;
          event.reference.assign(1U, reference[anchor - 1U]);
          event.alternate = event.reference + inserted;
          event.normalization = "rcrs_3prime_small_indel_v1";
          event.base_quality = minimum_quality;
          event.id =
              "indel:insertion:" + std::to_string(anchor) + ":" + inserted;
          calls.push_back(std::move(event));
        }
      }
      query_cursor += len;
      continue;
    }
    if (op == 'D') {
      if (len <= maximum_length && len < reference.size()) {
        std::size_t normalized_start = reference_cursor;
        const bool crosses_origin =
            reference_cursor == 1U ||
            len > reference.size() - reference_cursor + 1U;
        while (!crosses_origin && normalized_start <= reference.size() &&
               len <= reference.size() - normalized_start &&
               reference[normalized_start - 1U] ==
                   reference[normalized_start + len - 1U]) {
          ++normalized_start;
        }
        const auto normalized_end =
            ((normalized_start - 1U + (len - 1U)) % reference.size()) + 1U;
        const auto anchor =
            normalized_start > 1U ? normalized_start - 1U : reference.size();
        std::string deleted;
        deleted.reserve(len);
        for (std::size_t offset = 0U; offset < len; ++offset) {
          deleted.push_back(
              reference[(normalized_start - 1U + offset) % reference.size()]);
        }
        NormalizedSmallIndel event;
        event.type = "SMALL_DELETION";
        event.start = normalized_start;
        event.end = normalized_end;
        event.length = len;
        event.reference.assign(1U, reference[anchor - 1U]);
        event.reference += deleted;
        event.alternate.assign(1U, reference[anchor - 1U]);
        event.normalization = crosses_origin ? "rcrs_circular_small_indel_v1"
                                             : "rcrs_3prime_small_indel_v1";
        event.id = "indel:deletion:" + std::to_string(normalized_start) + "-" +
                   std::to_string(normalized_end) + ":" + deleted;
        calls.push_back(std::move(event));
      }
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
      continue;
    }
    if (op == 'N') {
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'S') {
      query_cursor += len;
    }
  }
  std::sort(calls.begin(), calls.end(),
            [](const auto &lhs, const auto &rhs) { return lhs.id < rhs.id; });
  calls.erase(std::unique(calls.begin(), calls.end(),
                          [](const auto &lhs, const auto &rhs) {
                            return lhs.id == rhs.id;
                          }),
              calls.end());
  return calls;
}

} // namespace mito::detail
