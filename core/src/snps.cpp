#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] bool looks_like_nuclear_contig(std::string_view reference_name) {
  std::string name(reference_name);
  std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  if (name.rfind("chr", 0) == 0) {
    name.erase(0, 3);
  }
  if (name == "x" || name == "y") {
    return true;
  }
  const auto chromosome = parse_size(name);
  return chromosome && *chromosome >= 1 && *chromosome <= 22;
}

[[nodiscard]] bool has_nuclear_supplementary_alignment(const ReadRecord &read) {
  return std::any_of(
      read.supplementary_alignments.begin(),
      read.supplementary_alignments.end(), [](const auto &alignment) {
        return looks_like_nuclear_contig(alignment.reference_name);
      });
}

[[nodiscard]] std::vector<std::string>
numt_evidence(const ReadRecord &read, std::size_t reference_length) {
  std::vector<std::string> evidence;
  const std::string id_lower = [&] {
    std::string out = read.id;
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
      return static_cast<char>(std::tolower(c));
    });
    return out;
  }();

  if (id_lower.find("numt") != std::string::npos ||
      id_lower.find("nuclear") != std::string::npos) {
    evidence.emplace_back("read_name_heuristic");
  }
  if (looks_like_nuclear_contig(read.reference_name)) {
    evidence.emplace_back("primary_nuclear_alignment");
  }
  if (has_nuclear_supplementary_alignment(read)) {
    evidence.emplace_back("supplementary_nuclear_alignment");
  }

  const double length_ratio = reference_length == 0
                                  ? 0.0
                                  : static_cast<double>(read.sequence.size()) /
                                        static_cast<double>(reference_length);
  const double gc = gc_fraction(read.sequence);
  if (length_ratio > 1.20) {
    evidence.emplace_back("length_exceeds_mtdna");
  }
  if (gc > 0.62 || gc < 0.25) {
    evidence.emplace_back("atypical_gc_fraction");
  }
  return evidence;
}

[[nodiscard]] double numt_score(const std::vector<std::string> &evidence) {
  double score = 0.05;
  for (const auto &item : evidence) {
    if (item == "primary_nuclear_alignment") {
      score = std::max(score, 0.99);
    } else if (item == "supplementary_nuclear_alignment") {
      score = std::max(score, 0.90);
    } else if (item == "read_name_heuristic") {
      score = std::max(score, 0.82);
    } else if (item == "length_exceeds_mtdna") {
      score += 0.30;
    } else if (item == "atypical_gc_fraction") {
      score += 0.18;
    }
  }
  return std::min(0.99, score);
}

[[nodiscard]] bool is_base(char base) {
  base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
  return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

[[nodiscard]] bool passes_snp_alignment_filters(const ReadRecord &read,
                                                const AnalysisConfig &config) {
  return (read.flags & 0x4U) == 0U &&
         (read.flags & config.excluded_snp_flags) == 0U &&
         read.mapping_quality >= config.min_mapping_quality &&
         read.reference_start != 0 && !read.reference_name.empty() &&
         read.reference_name != "*" &&
         !looks_like_nuclear_contig(read.reference_name);
}

[[nodiscard]] std::optional<std::uint8_t>
phred_quality_at(const ReadRecord &read, std::size_t query_index) {
  if (query_index >= read.sequence.size() ||
      query_index >= read.qualities.size()) {
    return std::nullopt;
  }
  const auto encoded = static_cast<unsigned char>(read.qualities[query_index]);
  if (encoded < 33U || encoded > 126U) {
    return std::nullopt;
  }
  return static_cast<std::uint8_t>(encoded - 33U);
}

[[nodiscard]] std::optional<std::size_t> base_index(char base) {
  switch (static_cast<char>(std::toupper(static_cast<unsigned char>(base)))) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return std::nullopt;
  }
}

[[nodiscard]] std::vector<SnpCall>
call_snps_from_alignment(const ReadRecord &read, const std::string &reference,
                         const AnalysisConfig &config) {
  std::vector<SnpCall> snps;
  if (reference.empty() || read.sequence.empty() ||
      read.cigar_operations.empty() ||
      !passes_snp_alignment_filters(read, config)) {
    return snps;
  }

  std::map<std::size_t, char> molecule_alleles;
  std::size_t read_cursor = 0;
  std::size_t reference_cursor =
      read.reference_start == 0 ? 1 : read.reference_start;

  for (const auto &operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;

    if (op == 'M' || op == '=' || op == 'X') {
      for (std::size_t offset = 0;
           offset < len && read_cursor + offset < read.sequence.size();
           ++offset) {
        const std::size_t position =
            ((reference_cursor - 1 + offset) % reference.size()) + 1;
        const char reference_base = reference[position - 1];
        const char alternate_base = static_cast<char>(std::toupper(
            static_cast<unsigned char>(read.sequence[read_cursor + offset])));
        const auto quality = phred_quality_at(read, read_cursor + offset);
        if (quality && *quality >= config.min_base_quality &&
            is_base(alternate_base) && is_base(reference_base)) {
          auto [allele, inserted] =
              molecule_alleles.try_emplace(position, alternate_base);
          if (!inserted && allele->second != alternate_base) {
            allele->second = 'N';
          }
        }
      }
      read_cursor += len;
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'I' || op == 'S') {
      read_cursor += len;
    } else if (op == 'D' || op == 'N') {
      reference_cursor = ((reference_cursor - 1U + (len % reference.size())) %
                          reference.size()) +
                         1U;
    } else if (op == 'H' || op == 'P') {
      // Does not consume read or reference sequence.
    }
  }

  snps.reserve(molecule_alleles.size());
  for (const auto &[position, alternate] : molecule_alleles) {
    const char reference_base = reference[position - 1U];
    if (is_base(alternate) && alternate != reference_base) {
      SnpCall call;
      call.position = position;
      call.reference = reference_base;
      call.alternate = alternate;
      snps.push_back(std::move(call));
    }
  }

  return snps;
}

[[nodiscard]] std::vector<SnpCall>
call_snps_from_reference_span(const ReadRecord &read,
                              const std::string &reference,
                              const AnalysisConfig &config) {
  std::vector<SnpCall> snps;
  if (reference.empty() || read.sequence.empty() || read.reference_start == 0 ||
      read.reference_name.empty() ||
      !passes_snp_alignment_filters(read, config)) {
    return snps;
  }

  for (std::size_t i = 0; i < read.sequence.size(); ++i) {
    const std::size_t position =
        ((read.reference_start - 1 + i) % reference.size()) + 1;
    const char reference_base = reference[position - 1];
    const char alternate_base = static_cast<char>(
        std::toupper(static_cast<unsigned char>(read.sequence[i])));
    const auto quality = phred_quality_at(read, i);
    if (quality && *quality >= config.min_base_quality &&
        is_base(alternate_base) && is_base(reference_base) &&
        alternate_base != reference_base) {
      SnpCall call;
      call.position = position;
      call.reference = reference_base;
      call.alternate = alternate_base;
      snps.push_back(std::move(call));
    }
  }

  return snps;
}

[[nodiscard]] std::optional<SnpCall> parse_snp_tag(std::string_view token) {
  const auto colon = token.find(':');
  const auto arrow =
      token.find('>', colon == std::string_view::npos ? 0 : colon + 1);
  if (colon == std::string_view::npos || arrow == std::string_view::npos ||
      arrow + 1 >= token.size()) {
    return std::nullopt;
  }

  const auto position = parse_size(token.substr(0, colon));
  if (!position || *position == 0) {
    return std::nullopt;
  }

  const char reference = static_cast<char>(
      std::toupper(static_cast<unsigned char>(token[colon + 1])));
  const char alternate = static_cast<char>(
      std::toupper(static_cast<unsigned char>(token[arrow + 1])));
  const auto valid_base = [](char base) {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
  };
  if (!valid_base(reference) || !valid_base(alternate) ||
      reference == alternate) {
    return std::nullopt;
  }

  SnpCall call;
  call.position = *position;
  call.reference = reference;
  call.alternate = alternate;
  return call;
}

[[nodiscard]] std::vector<SnpCall>
snp_tags_from_header(const ReadRecord &read) {
  std::vector<SnpCall> calls;
  std::istringstream iss(read.id);
  std::string token;
  while (iss >> token) {
    constexpr std::string_view prefix = "snp=";
    if (token.rfind(prefix, 0) == 0) {
      auto call = parse_snp_tag(std::string_view(token).substr(prefix.size()));
      if (call) {
        calls.push_back(std::move(*call));
      }
    }
  }
  return calls;
}

void merge_snps(std::vector<SnpCall> &snps, std::vector<SnpCall> tagged_snps) {
  snps.insert(snps.end(), std::make_move_iterator(tagged_snps.begin()),
              std::make_move_iterator(tagged_snps.end()));
  std::sort(snps.begin(), snps.end(), [](const auto &lhs, const auto &rhs) {
    return snp_key(lhs.position, lhs.reference, lhs.alternate) <
           snp_key(rhs.position, rhs.reference, rhs.alternate);
  });
  snps.erase(std::unique(snps.begin(), snps.end(),
                         [](const auto &lhs, const auto &rhs) {
                           return lhs.position == rhs.position &&
                                  lhs.reference == rhs.reference &&
                                  lhs.alternate == rhs.alternate;
                         }),
             snps.end());
}

[[nodiscard]] bool
observed_phylo_mutation_less(const ObservedPhyloMutation &lhs,
                             const ObservedPhyloMutation &rhs) {
  return std::tie(lhs.position, lhs.kind, lhs.insertion_index, lhs.encoded) <
         std::tie(rhs.position, rhs.kind, rhs.insertion_index, rhs.encoded);
}

[[nodiscard]] ObservedPhyloMutation
observed_phylo_mutation(const PhyloMutation &mutation) {
  ObservedPhyloMutation observed;
  observed.position = mutation.locus.position;
  observed.kind = mutation.kind;
  observed.alternate = mutation.alternate;
  observed.insertion_index = mutation.locus.insertion_index;
  observed.inserted_bases = mutation.inserted_bases;
  observed.encoded = mutation.encoded;
  return observed;
}

} // namespace mito::detail
