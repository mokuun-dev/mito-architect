#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] std::optional<std::size_t> parse_size(std::string_view value) {
  std::size_t parsed = 0;
  const auto *begin = value.data();
  const auto *end = value.data() + value.size();
  const auto result = std::from_chars(begin, end, parsed);
  if (result.ec != std::errc{} || result.ptr != end) {
    return std::nullopt;
  }
  return parsed;
}

[[nodiscard]] std::optional<PhyloMutation>
parse_phylo_mutation(std::string_view raw_token) {
  std::string token = trim_copy(raw_token);
  if (token.empty()) {
    return std::nullopt;
  }

  bool backmutation = false;
  if (token.back() == '!') {
    backmutation = true;
    token.pop_back();
  }
  if (token.empty() || token.find('!') != std::string::npos) {
    return std::nullopt;
  }

  std::size_t digit_count = 0;
  while (digit_count < token.size() &&
         std::isdigit(static_cast<unsigned char>(token[digit_count])) != 0) {
    ++digit_count;
  }
  if (digit_count == 0 || digit_count >= token.size()) {
    return std::nullopt;
  }
  const auto position =
      parse_size(std::string_view(token).substr(0, digit_count));
  if (!position || *position == 0 ||
      *position > static_cast<std::size_t>(kDefaultReferenceLength)) {
    return std::nullopt;
  }

  const auto canonical_base = [](char value) -> std::optional<char> {
    const char upper =
        static_cast<char>(std::toupper(static_cast<unsigned char>(value)));
    if (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T') {
      return upper;
    }
    return std::nullopt;
  };

  PhyloMutation mutation;
  mutation.locus.position = *position;
  mutation.backmutation = backmutation;
  const auto suffix = std::string_view(token).substr(digit_count);
  if (suffix.size() == 1U && (suffix.front() == 'd' || suffix.front() == 'D')) {
    mutation.kind = PhyloMutationKind::deletion;
    mutation.locus.kind = PhyloMutationKind::substitution;
    mutation.encoded = std::to_string(*position) + "d";
    return mutation;
  }

  if (suffix.front() != '.') {
    if (suffix.size() != 1U) {
      return std::nullopt;
    }
    const auto alternate = canonical_base(suffix.front());
    if (!alternate) {
      return std::nullopt;
    }
    mutation.kind = PhyloMutationKind::substitution;
    mutation.locus.kind = PhyloMutationKind::substitution;
    mutation.alternate = *alternate;
    mutation.encoded = std::to_string(*position) + std::string(1, *alternate);
    return mutation;
  }

  std::size_t index_end = 1U;
  if (index_end < suffix.size() && suffix[index_end] == 'X') {
    ++index_end;
    mutation.locus.insertion_index = "X";
  } else {
    const std::size_t index_start = index_end;
    while (index_end < suffix.size() &&
           std::isdigit(static_cast<unsigned char>(suffix[index_end])) != 0) {
      ++index_end;
    }
    if (index_start == index_end) {
      return std::nullopt;
    }
    const auto insertion_index =
        parse_size(suffix.substr(index_start, index_end - index_start));
    if (!insertion_index || *insertion_index == 0) {
      return std::nullopt;
    }
    mutation.locus.insertion_index = std::to_string(*insertion_index);
  }
  if (index_end >= suffix.size()) {
    return std::nullopt;
  }

  mutation.inserted_bases.reserve(suffix.size() - index_end);
  for (const char value : suffix.substr(index_end)) {
    const auto base = canonical_base(value);
    if (!base) {
      return std::nullopt;
    }
    mutation.inserted_bases.push_back(*base);
  }
  mutation.kind = PhyloMutationKind::insertion;
  mutation.locus.kind = PhyloMutationKind::insertion;
  mutation.encoded = std::to_string(*position) + "." +
                     mutation.locus.insertion_index + mutation.inserted_bases;
  return mutation;
}

[[nodiscard]] std::vector<CigarOperation>
parse_cigar(std::string_view cigar, std::string_view read_id) {
  std::vector<CigarOperation> operations;
  if (cigar.empty() || cigar == "*") {
    return operations;
  }

  std::size_t number_start = 0;
  for (std::size_t i = 0; i < cigar.size(); ++i) {
    if (std::isdigit(static_cast<unsigned char>(cigar[i])) != 0) {
      continue;
    }

    const auto length =
        parse_size(cigar.substr(number_start, i - number_start));
    constexpr std::string_view valid_operations = "MIDNSHP=X";
    if (!length || *length == 0 ||
        valid_operations.find(cigar[i]) == std::string_view::npos) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid CIGAR for read '" + std::string(read_id) +
                              "': " + std::string(cigar));
    }
    operations.push_back({*length, cigar[i]});
    number_start = i + 1;
  }

  if (number_start != cigar.size() || operations.empty()) {
    throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                        "invalid CIGAR for read '" + std::string(read_id) +
                            "': " + std::string(cigar));
  }
  return operations;
}

[[nodiscard]] std::vector<std::string>
split_alignment_rule_tokens(std::string_view value) {
  std::vector<std::string> tokens;
  std::size_t cursor = 0U;
  while (cursor < value.size()) {
    while (cursor < value.size() &&
           (std::isspace(static_cast<unsigned char>(value[cursor])) != 0 ||
            value[cursor] == ',')) {
      ++cursor;
    }
    const std::size_t start = cursor;
    while (cursor < value.size() &&
           std::isspace(static_cast<unsigned char>(value[cursor])) == 0 &&
           value[cursor] != ',') {
      ++cursor;
    }
    if (start != cursor) {
      tokens.emplace_back(value.substr(start, cursor - start));
    }
  }
  return tokens;
}

[[nodiscard]] bool is_reference_restore_token(std::string_view token) {
  return !token.empty() &&
         std::all_of(token.begin(), token.end(), [](const char value) {
           return std::isdigit(static_cast<unsigned char>(value)) != 0;
         });
}

[[nodiscard]] std::vector<PhyloAlignmentRule> load_phylo_alignment_rules() {
  const auto path = phylotree_alignment_rules_path();
  std::ifstream input(path);
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "could not open PhyloTree alignment rules: " + path);
  }

  std::vector<PhyloAlignmentRule> rules;
  std::string line;
  std::size_t line_number = 0U;
  while (std::getline(input, line)) {
    ++line_number;
    if (!line.empty() && line.back() == '\r') {
      line.pop_back();
    }
    if (line_number == 1U) {
      if (line != "error,expected") {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "PhyloTree alignment-rule header is invalid");
      }
      continue;
    }
    if (trim_copy(line).empty()) {
      continue;
    }
    const auto delimiter = line.find(',');
    if (delimiter == std::string::npos) {
      throw AnalysisError(AnalysisErrorCode::resource_invalid,
                          "PhyloTree alignment rule has no delimiter at line " +
                              std::to_string(line_number));
    }

    PhyloAlignmentRule rule;
    bool callable_error_side = true;
    for (const auto &token : split_alignment_rule_tokens(
             std::string_view(line).substr(0U, delimiter))) {
      const auto mutation = parse_phylo_mutation(token);
      if (!mutation || mutation->backmutation) {
        // Ambiguity-code rules such as 3106N cannot be emitted by the
        // fail-closed callable-base path and are intentionally inapplicable.
        callable_error_side = false;
        break;
      }
      rule.erroneous_markers.push_back(mutation->encoded);
    }
    if (!callable_error_side) {
      continue;
    }
    if (rule.erroneous_markers.empty()) {
      throw AnalysisError(
          AnalysisErrorCode::resource_invalid,
          "PhyloTree alignment rule has an empty error side at line " +
              std::to_string(line_number));
    }
    std::sort(rule.erroneous_markers.begin(), rule.erroneous_markers.end());
    rule.erroneous_markers.erase(std::unique(rule.erroneous_markers.begin(),
                                             rule.erroneous_markers.end()),
                                 rule.erroneous_markers.end());

    for (const auto &token : split_alignment_rule_tokens(
             std::string_view(line).substr(delimiter + 1U))) {
      if (is_reference_restore_token(token)) {
        continue;
      }
      const auto mutation = parse_phylo_mutation(token);
      if (!mutation || mutation->backmutation) {
        throw AnalysisError(
            AnalysisErrorCode::resource_invalid,
            "unsupported PhyloTree alignment-rule replacement '" + token +
                "' at line " + std::to_string(line_number));
      }
      rule.replacement_markers.push_back(observed_phylo_mutation(*mutation));
    }
    std::sort(rule.replacement_markers.begin(), rule.replacement_markers.end(),
              observed_phylo_mutation_less);
    rule.replacement_markers.erase(
        std::unique(rule.replacement_markers.begin(),
                    rule.replacement_markers.end(),
                    [](const auto &lhs, const auto &rhs) {
                      return lhs.encoded == rhs.encoded;
                    }),
        rule.replacement_markers.end());
    rules.push_back(std::move(rule));
  }
  if (rules.size() < 100U) {
    throw AnalysisError(AnalysisErrorCode::resource_invalid,
                        "PhyloTree alignment-rule resource is incomplete");
  }
  return rules;
}

void apply_phylo_alignment_rules(std::vector<ObservedPhyloMutation> &observed,
                                 const std::vector<PhyloAlignmentRule> &rules) {
  std::map<std::string, ObservedPhyloMutation> markers;
  for (auto &mutation : observed) {
    markers.insert_or_assign(mutation.encoded, std::move(mutation));
  }

  for (std::size_t pass = 0U; pass <= rules.size(); ++pass) {
    bool changed = false;
    for (const auto &rule : rules) {
      if (!std::all_of(
              rule.erroneous_markers.begin(), rule.erroneous_markers.end(),
              [&](const auto &marker) { return markers.contains(marker); })) {
        continue;
      }
      std::map<std::string, ObservedPhyloMutation> candidate = markers;
      for (const auto &marker : rule.erroneous_markers) {
        candidate.erase(marker);
      }
      for (const auto &replacement : rule.replacement_markers) {
        candidate.insert_or_assign(replacement.encoded, replacement);
      }
      const bool same_keys =
          candidate.size() == markers.size() &&
          std::equal(candidate.begin(), candidate.end(), markers.begin(),
                     [](const auto &lhs, const auto &rhs) {
                       return lhs.first == rhs.first;
                     });
      if (!same_keys) {
        markers = std::move(candidate);
        changed = true;
      }
    }
    if (!changed) {
      observed.clear();
      observed.reserve(markers.size());
      for (auto &[_, marker] : markers) {
        observed.push_back(std::move(marker));
      }
      std::sort(observed.begin(), observed.end(), observed_phylo_mutation_less);
      return;
    }
  }
  throw AnalysisError(AnalysisErrorCode::resource_invalid,
                      "PhyloTree alignment rules did not converge");
}

[[nodiscard]] std::vector<ObservedPhyloMutation>
phylo_tags_from_header(const ReadRecord &read) {
  std::vector<ObservedPhyloMutation> calls;
  std::istringstream iss(read.id);
  std::string token;
  while (iss >> token) {
    constexpr std::string_view prefix = "phylo=";
    if (token.rfind(prefix, 0) != 0) {
      continue;
    }
    const auto mutation =
        parse_phylo_mutation(std::string_view(token).substr(prefix.size()));
    if (!mutation || mutation->backmutation) {
      continue;
    }
    calls.push_back(observed_phylo_mutation(*mutation));
  }
  return calls;
}

[[nodiscard]] std::vector<ObservedPhyloMutation>
call_phylo_indels_from_alignment(const ReadRecord &read,
                                 const AnalysisConfig &config,
                                 const std::string &reference) {
  std::vector<ObservedPhyloMutation> calls;
  const std::size_t reference_length = reference.size();
  if (reference_length == 0 || read.sequence.empty() ||
      read.reference_start == 0U || read.cigar_operations.empty() ||
      !passes_snp_alignment_filters(read, config)) {
    return calls;
  }

  constexpr std::size_t kMaxPhylogeneticIndelLength = 50U;
  std::size_t read_cursor = 0;
  std::size_t reference_cursor =
      ((read.reference_start - 1U) % reference_length) + 1U;
  for (const auto &operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;
    if (op == 'M' || op == '=' || op == 'X') {
      read_cursor += len;
      reference_cursor = ((reference_cursor - 1U + (len % reference_length)) %
                          reference_length) +
                         1U;
    } else if (op == 'I') {
      if (len <= kMaxPhylogeneticIndelLength &&
          read_cursor <= read.sequence.size() &&
          len <= read.sequence.size() - read_cursor) {
        std::string inserted;
        inserted.reserve(len);
        bool callable = true;
        for (std::size_t offset = 0; offset < len; ++offset) {
          const char base = static_cast<char>(std::toupper(
              static_cast<unsigned char>(read.sequence[read_cursor + offset])));
          const auto quality = phred_quality_at(read, read_cursor + offset);
          if (!is_base(base) || !quality ||
              *quality < config.min_base_quality) {
            callable = false;
            break;
          }
          inserted.push_back(base);
        }
        if (callable) {
          std::size_t anchor =
              reference_cursor > 1U ? reference_cursor - 1U : reference_length;
          // Mitochondrial variant nomenclature places repeat-equivalent indels
          // at the most 3-prime coordinate. Rotate the inserted sequence while
          // shifting so that non-homopolymer repeats remain
          // haplotype-equivalent. Do not cross the rCRS coordinate boundary:
          // circular-origin events have their own canonical representation and
          // no unique 3-prime endpoint.
          std::size_t shift = 0U;
          while (anchor < reference_length && !inserted.empty() &&
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
          ObservedPhyloMutation mutation;
          mutation.position = anchor;
          mutation.kind = PhyloMutationKind::insertion;
          mutation.insertion_index = "1";
          mutation.inserted_bases = std::move(inserted);
          mutation.encoded = std::to_string(mutation.position) + ".1" +
                             mutation.inserted_bases;
          calls.push_back(std::move(mutation));
        }
      }
      read_cursor += len;
    } else if (op == 'D') {
      if (len <= kMaxPhylogeneticIndelLength) {
        std::size_t normalized_start = reference_cursor;
        const bool crosses_origin =
            len > reference_length - reference_cursor + 1U;
        while (!crosses_origin && normalized_start <= reference_length &&
               len <= reference_length - normalized_start &&
               reference[normalized_start - 1U] ==
                   reference[normalized_start + len - 1U]) {
          ++normalized_start;
        }
        for (std::size_t offset = 0; offset < len; ++offset) {
          ObservedPhyloMutation mutation;
          mutation.position =
              ((normalized_start - 1U + offset) % reference_length) + 1U;
          mutation.kind = PhyloMutationKind::deletion;
          mutation.encoded = std::to_string(mutation.position) + "d";
          calls.push_back(std::move(mutation));
        }
      }
      reference_cursor = ((reference_cursor - 1U + (len % reference_length)) %
                          reference_length) +
                         1U;
    } else if (op == 'N') {
      reference_cursor = ((reference_cursor - 1U + (len % reference_length)) %
                          reference_length) +
                         1U;
    } else if (op == 'S') {
      read_cursor += len;
    } else if (op == 'H' || op == 'P') {
      // Does not consume query or reference sequence.
    }
  }
  return calls;
}

} // namespace mito::detail
