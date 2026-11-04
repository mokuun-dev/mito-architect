#include "detail/pipeline.hpp"
#include "detail/precision.hpp"

// Functional annotation against rCRS. Gene coordinates come from the NCBI
// feature table of NC_012920.1, and codons are read with translation table 2
// (vertebrate mitochondrial). The translation of all 13 proteins is checked by
// scripts/build_structure_catalog.py and core/tests/precision_test.cpp.
namespace mito::detail {
namespace {

enum class Biotype : std::uint8_t { protein_coding, rRNA, tRNA };

struct Feature {
  std::string_view name;
  std::size_t start;
  std::size_t end;
  char strand;
  Biotype biotype;
};

constexpr std::array<Feature, 37> kFeatures = {{
    {"MT-TF", 577U, 647U, '+', Biotype::tRNA},
    {"MT-RNR1", 648U, 1601U, '+', Biotype::rRNA},
    {"MT-TV", 1602U, 1670U, '+', Biotype::tRNA},
    {"MT-RNR2", 1671U, 3229U, '+', Biotype::rRNA},
    {"MT-TL1", 3230U, 3304U, '+', Biotype::tRNA},
    {"MT-ND1", 3307U, 4262U, '+', Biotype::protein_coding},
    {"MT-TI", 4263U, 4331U, '+', Biotype::tRNA},
    {"MT-TQ", 4329U, 4400U, '-', Biotype::tRNA},
    {"MT-TM", 4402U, 4469U, '+', Biotype::tRNA},
    {"MT-ND2", 4470U, 5511U, '+', Biotype::protein_coding},
    {"MT-TW", 5512U, 5579U, '+', Biotype::tRNA},
    {"MT-TA", 5587U, 5655U, '-', Biotype::tRNA},
    {"MT-TN", 5657U, 5729U, '-', Biotype::tRNA},
    {"MT-TC", 5761U, 5826U, '-', Biotype::tRNA},
    {"MT-TY", 5826U, 5891U, '-', Biotype::tRNA},
    {"MT-CO1", 5904U, 7445U, '+', Biotype::protein_coding},
    {"MT-TS1", 7446U, 7514U, '-', Biotype::tRNA},
    {"MT-TD", 7518U, 7585U, '+', Biotype::tRNA},
    {"MT-CO2", 7586U, 8269U, '+', Biotype::protein_coding},
    {"MT-TK", 8295U, 8364U, '+', Biotype::tRNA},
    {"MT-ATP8", 8366U, 8572U, '+', Biotype::protein_coding},
    {"MT-ATP6", 8527U, 9207U, '+', Biotype::protein_coding},
    {"MT-CO3", 9207U, 9990U, '+', Biotype::protein_coding},
    {"MT-TG", 9991U, 10058U, '+', Biotype::tRNA},
    {"MT-ND3", 10059U, 10404U, '+', Biotype::protein_coding},
    {"MT-TR", 10405U, 10469U, '+', Biotype::tRNA},
    {"MT-ND4L", 10470U, 10766U, '+', Biotype::protein_coding},
    {"MT-ND4", 10760U, 12137U, '+', Biotype::protein_coding},
    {"MT-TH", 12138U, 12206U, '+', Biotype::tRNA},
    {"MT-TS2", 12207U, 12265U, '+', Biotype::tRNA},
    {"MT-TL2", 12266U, 12336U, '+', Biotype::tRNA},
    {"MT-ND5", 12337U, 14148U, '+', Biotype::protein_coding},
    {"MT-ND6", 14149U, 14673U, '-', Biotype::protein_coding},
    {"MT-TE", 14674U, 14742U, '-', Biotype::tRNA},
    {"MT-CYB", 14747U, 15887U, '+', Biotype::protein_coding},
    {"MT-TT", 15888U, 15953U, '+', Biotype::tRNA},
    {"MT-TP", 15956U, 16023U, '-', Biotype::tRNA},
}};

// Replication-origin regions and the control region as listed by MITOMAP.
constexpr std::size_t kHeavyOriginStart = 110U;
constexpr std::size_t kHeavyOriginEnd = 441U;
constexpr std::size_t kLightOriginStart = 5721U;
constexpr std::size_t kLightOriginEnd = 5798U;
constexpr std::size_t kControlRegionStart = 16024U;
constexpr std::size_t kControlRegionEnd = 576U;
constexpr std::size_t kMajorArcStart = 5799U;
constexpr std::size_t kMajorArcEnd = 16569U;

[[nodiscard]] constexpr std::string_view biotype_name(const Biotype biotype) noexcept {
  switch (biotype) {
  case Biotype::protein_coding:
    return "protein_coding";
  case Biotype::rRNA:
    return "rRNA";
  case Biotype::tRNA:
    return "tRNA";
  }
  return "unknown";
}

[[nodiscard]] constexpr char complement(const char base) noexcept {
  switch (base) {
  case 'A':
    return 'T';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  case 'T':
    return 'A';
  default:
    return 'N';
  }
}

[[nodiscard]] constexpr int code_index(const char base) noexcept {
  switch (base) {
  case 'T':
    return 0;
  case 'C':
    return 1;
  case 'A':
    return 2;
  case 'G':
    return 3;
  default:
    return -1;
  }
}

/** NCBI translation table 2; returns 'X' for a codon with an unknown base. */
[[nodiscard]] char translate_codon(const std::string_view codon) noexcept {
  constexpr std::string_view table =
      "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
  const int first = code_index(codon[0]);
  const int second = code_index(codon[1]);
  const int third = code_index(codon[2]);
  if (first < 0 || second < 0 || third < 0) {
    return 'X';
  }
  return table[static_cast<std::size_t>(16 * first + 4 * second + third)];
}

[[nodiscard]] bool is_start_codon(const std::string_view codon) noexcept {
  return codon == "ATG" || codon == "ATA" || codon == "ATT" || codon == "ATC" ||
         codon == "GTG";
}

struct CodonLocation {
  std::size_t codon_index = 0;
  std::size_t phase = 0;
};

[[nodiscard]] CodonLocation codon_location(const Feature &gene,
                                           const std::size_t position) noexcept {
  const std::size_t offset =
      gene.strand == '+' ? position - gene.start : gene.end - position;
  return {offset / 3U, offset % 3U};
}

/** Codon on the coding strand. A truncated final codon is completed with A,
 * since polyadenylation supplies the rest of the stop codon in these genes. */
[[nodiscard]] std::string gene_codon(const std::string &reference,
                                     const Feature &gene,
                                     const std::size_t codon_index) {
  std::string codon;
  codon.reserve(3U);
  for (std::size_t offset = codon_index * 3U; offset < codon_index * 3U + 3U;
       ++offset) {
    if (gene.strand == '+') {
      const auto position = gene.start + offset;
      codon.push_back(position <= gene.end ? reference[position - 1U] : 'A');
    } else {
      codon.push_back(gene.end >= gene.start + offset
                          ? complement(reference[gene.end - offset - 1U])
                          : 'A');
    }
  }
  return codon;
}

[[nodiscard]] std::string three_letter(const char code) {
  return std::string(amino_acid_three_letter(code));
}

struct CatalogEntry {
  std::string pdb_id;
  std::string chain;
  std::string complex_name;
  std::vector<std::pair<std::size_t, std::size_t>> modelled;
  std::set<std::size_t> differs;
};

[[nodiscard]] std::map<std::string, CatalogEntry> load_structure_catalog() {
  std::ifstream input(structure_catalog_path());
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "structure catalog could not be opened: " +
                            structure_catalog_path());
  }
  std::map<std::string, CatalogEntry> catalog;
  std::string line;
  std::size_t line_number = 0U;
  while (std::getline(input, line)) {
    ++line_number;
    strip_trailing_carriage_return(line);
    if (line_number == 1U || line.empty()) {
      continue;
    }
    const auto fields = split_tab(line);
    if (fields.size() != 8U) {
      throw AnalysisError(AnalysisErrorCode::resource_invalid,
                          "structure catalog line " +
                              std::to_string(line_number) +
                              " must have 8 fields");
    }
    CatalogEntry entry;
    entry.pdb_id = fields[2];
    entry.chain = fields[3];
    entry.complex_name = fields[4];
    std::stringstream ranges(fields[6]);
    std::string range;
    while (std::getline(ranges, range, ',')) {
      const auto dash = range.find('-');
      const auto low = parse_size(range.substr(0, dash));
      const auto high =
          dash == std::string::npos ? low : parse_size(range.substr(dash + 1));
      if (!low || !high || *low > *high) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "structure catalog line " +
                                std::to_string(line_number) +
                                " has an invalid residue range");
      }
      entry.modelled.emplace_back(*low, *high);
    }
    std::stringstream differs(fields[7]);
    std::string change;
    while (std::getline(differs, change, ',')) {
      // Format: rCRS residue, number, model residue (for example A52T).
      const auto number = change.size() > 2U
                              ? parse_size(std::string_view(change).substr(
                                    1U, change.size() - 2U))
                              : std::nullopt;
      if (!number) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "structure catalog line " +
                                std::to_string(line_number) +
                                " has an invalid residue difference");
      }
      entry.differs.insert(*number);
    }
    catalog.emplace(fields[0], std::move(entry));
  }
  return catalog;
}

/** A residue maps only where a verified model has exactly the rCRS residue. */
[[nodiscard]] std::optional<StructureMapping>
map_structure(const std::string_view gene, const std::size_t residue) {
  static const auto catalog = load_structure_catalog();
  const auto entry = catalog.find(std::string(gene));
  if (entry == catalog.end() || entry->second.differs.contains(residue)) {
    return std::nullopt;
  }
  const bool modelled = std::any_of(
      entry->second.modelled.begin(), entry->second.modelled.end(),
      [&](const auto &range) {
        return range.first <= residue && residue <= range.second;
      });
  if (!modelled) {
    return std::nullopt;
  }
  return StructureMapping{entry->second.pdb_id, entry->second.chain, residue,
                          entry->second.complex_name};
}

[[nodiscard]] bool in_control_region(const std::size_t position) noexcept {
  return position >= kControlRegionStart || position <= kControlRegionEnd;
}

[[nodiscard]] VariantEffect noncoding_effect(const Feature &gene) {
  VariantEffect effect;
  effect.gene = std::string(gene.name);
  effect.biotype = std::string(biotype_name(gene.biotype));
  effect.consequence =
      gene.biotype == Biotype::tRNA ? "tRNA_variant" : "rRNA_variant";
  return effect;
}

[[nodiscard]] VariantEffect regional_effect(const std::size_t position) {
  VariantEffect effect;
  if (in_control_region(position)) {
    effect.gene = "MT-CR";
    effect.biotype = "control_region";
    effect.consequence = "control_region_variant";
  } else {
    effect.biotype = "intergenic";
    effect.consequence = "intergenic_variant";
  }
  return effect;
}

/** Fraction of circular interval [start, end] overlapping one feature. */
[[nodiscard]] std::size_t circular_overlap(const std::size_t start,
                                           const std::size_t end,
                                           const std::size_t feature_start,
                                           const std::size_t feature_end,
                                           const std::size_t length) noexcept {
  const auto linear = [](const std::size_t a, const std::size_t b,
                         const std::size_t c, const std::size_t d) {
    const auto low = std::max(a, c);
    const auto high = std::min(b, d);
    return low <= high ? high - low + 1U : std::size_t{0};
  };
  if (start <= end) {
    return linear(start, end, feature_start, feature_end);
  }
  return linear(start, length, feature_start, feature_end) +
         linear(1U, end, feature_start, feature_end);
}

} // namespace

[[nodiscard]] std::string_view amino_acid_three_letter(const char code) noexcept {
  switch (code) {
  case 'A':
    return "Ala";
  case 'R':
    return "Arg";
  case 'N':
    return "Asn";
  case 'D':
    return "Asp";
  case 'C':
    return "Cys";
  case 'Q':
    return "Gln";
  case 'E':
    return "Glu";
  case 'G':
    return "Gly";
  case 'H':
    return "His";
  case 'I':
    return "Ile";
  case 'L':
    return "Leu";
  case 'K':
    return "Lys";
  case 'M':
    return "Met";
  case 'F':
    return "Phe";
  case 'P':
    return "Pro";
  case 'S':
    return "Ser";
  case 'T':
    return "Thr";
  case 'W':
    return "Trp";
  case 'Y':
    return "Tyr";
  case 'V':
    return "Val";
  case '*':
    return "Ter";
  default:
    return "Xaa";
  }
}

[[nodiscard]] std::optional<int> grantham_distance(const char first,
                                                   const char second) noexcept {
  struct Properties {
    char code;
    double composition;
    double polarity;
    double volume;
  };
  // Composition, polarity and molecular volume from Grantham (1974), Table 1.
  static constexpr std::array<Properties, 20> properties = {{
      {'S', 1.42, 9.2, 32.0},  {'R', 0.65, 10.5, 124.0},
      {'L', 0.0, 4.9, 111.0},  {'P', 0.39, 8.0, 32.5},
      {'T', 0.71, 8.6, 61.0},  {'A', 0.0, 8.1, 31.0},
      {'V', 0.0, 5.9, 84.0},   {'G', 0.74, 9.0, 3.0},
      {'I', 0.0, 5.2, 111.0},  {'F', 0.0, 5.2, 132.0},
      {'Y', 0.20, 6.2, 136.0}, {'C', 2.75, 5.5, 55.0},
      {'H', 0.58, 10.4, 96.0}, {'Q', 0.89, 10.5, 85.0},
      {'N', 1.33, 11.6, 56.0}, {'K', 0.33, 11.3, 119.0},
      {'D', 1.38, 13.0, 54.0}, {'E', 0.92, 12.3, 83.0},
      {'M', 0.0, 5.7, 105.0},  {'W', 0.13, 5.4, 170.0},
  }};
  const auto find = [](const char code) -> const Properties * {
    for (const auto &item : properties) {
      if (item.code == code) {
        return &item;
      }
    }
    return nullptr;
  };
  const auto *lhs = find(first);
  const auto *rhs = find(second);
  if (lhs == nullptr || rhs == nullptr) {
    return std::nullopt;
  }
  const double composition = lhs->composition - rhs->composition;
  const double polarity = lhs->polarity - rhs->polarity;
  const double volume = lhs->volume - rhs->volume;
  const double distance =
      50.723 * std::sqrt(1.833 * composition * composition +
                         0.1018 * polarity * polarity +
                         0.000399 * volume * volume);
  return static_cast<int>(std::lround(distance));
}

[[nodiscard]] std::string translate_rcrs_gene(const std::string &reference,
                                              const std::string_view gene) {
  for (const auto &feature : kFeatures) {
    if (feature.name != gene || feature.biotype != Biotype::protein_coding) {
      continue;
    }
    std::string protein;
    const std::size_t codons = (feature.end - feature.start + 3U) / 3U;
    for (std::size_t index = 0U; index < codons; ++index) {
      const char amino = index == 0U
                             ? 'M'
                             : translate_codon(gene_codon(reference, feature, index));
      if (amino == '*') {
        break;
      }
      protein.push_back(amino);
    }
    return protein;
  }
  return {};
}

[[nodiscard]] std::vector<VariantEffect>
annotate_substitution(const std::string &reference, const std::size_t position,
                      const char alternate) {
  std::vector<VariantEffect> effects;
  for (const auto &gene : kFeatures) {
    if (position < gene.start || position > gene.end) {
      continue;
    }
    if (gene.biotype != Biotype::protein_coding) {
      effects.push_back(noncoding_effect(gene));
      continue;
    }
    const auto location = codon_location(gene, position);
    const auto reference_codon =
        gene_codon(reference, gene, location.codon_index);
    auto alternate_codon = reference_codon;
    alternate_codon[location.phase] =
        gene.strand == '+' ? alternate : complement(alternate);
    const bool first_codon = location.codon_index == 0U;
    char reference_amino = translate_codon(reference_codon);
    char alternate_amino = translate_codon(alternate_codon);
    VariantEffect effect;
    effect.gene = std::string(gene.name);
    effect.biotype = "protein_coding";
    effect.amino_acid_position = location.codon_index + 1U;
    effect.codon_change = reference_codon + ">" + alternate_codon;
    const auto residue_number = std::to_string(effect.amino_acid_position);
    if (first_codon) {
      // Every mitochondrial initiator codon encodes methionine.
      reference_amino = 'M';
      alternate_amino = is_start_codon(alternate_codon) ? 'M' : alternate_amino;
      if (alternate_amino == 'M') {
        effect.consequence = "start_retained_variant";
        effect.hgvs_protein = "p.Met1=";
      } else {
        effect.consequence = "start_lost";
        effect.hgvs_protein = "p.Met1?";
      }
    } else if (reference_amino == alternate_amino) {
      effect.consequence =
          reference_amino == '*' ? "stop_retained_variant" : "synonymous_variant";
      effect.hgvs_protein =
          "p." + three_letter(reference_amino) + residue_number + "=";
    } else if (alternate_amino == '*') {
      effect.consequence = "stop_gained";
      effect.hgvs_protein =
          "p." + three_letter(reference_amino) + residue_number + "Ter";
    } else if (reference_amino == '*') {
      effect.consequence = "stop_lost";
      effect.hgvs_protein = "p.Ter" + residue_number +
                            three_letter(alternate_amino) + "extTer?";
    } else {
      effect.consequence = "missense_variant";
      effect.hgvs_protein = "p." + three_letter(reference_amino) +
                            residue_number + three_letter(alternate_amino);
      effect.grantham = grantham_distance(reference_amino, alternate_amino);
    }
    effect.reference_amino_acid = reference_amino;
    effect.alternate_amino_acid = alternate_amino;
    if (reference_amino != '*') {
      effect.structure = map_structure(gene.name, effect.amino_acid_position);
    }
    effects.push_back(std::move(effect));
  }
  if (effects.empty()) {
    effects.push_back(regional_effect(position));
  }
  return effects;
}

[[nodiscard]] std::vector<VariantEffect>
annotate_small_indel(const std::string &reference, const std::size_t position,
                     const std::size_t reference_length_removed,
                     const std::size_t inserted_length, const bool insertion) {
  std::vector<VariantEffect> effects;
  const std::size_t event_end =
      insertion ? position : position + reference_length_removed - 1U;
  for (const auto &gene : kFeatures) {
    // An insertion between the two last bases of a gene still falls inside it;
    // one after the last base does not, hence the strict comparison at the end.
    const bool overlaps = insertion ? position >= gene.start && position < gene.end
                                    : position <= gene.end && event_end >= gene.start;
    if (!overlaps) {
      continue;
    }
    if (gene.biotype != Biotype::protein_coding) {
      effects.push_back(noncoding_effect(gene));
      continue;
    }
    const std::size_t affected =
        insertion ? position
                  : std::clamp(gene.strand == '+' ? position : event_end,
                               gene.start, gene.end);
    const auto location = codon_location(gene, affected);
    const auto reference_codon =
        gene_codon(reference, gene, location.codon_index);
    const char amino = location.codon_index == 0U
                           ? 'M'
                           : translate_codon(reference_codon);
    const std::size_t change_length =
        insertion ? inserted_length : reference_length_removed;
    VariantEffect effect;
    effect.gene = std::string(gene.name);
    effect.biotype = "protein_coding";
    effect.amino_acid_position = location.codon_index + 1U;
    effect.reference_amino_acid = amino;
    const auto residue = three_letter(amino) +
                         std::to_string(effect.amino_acid_position);
    if (change_length % 3U != 0U) {
      effect.consequence = "frameshift_variant";
      effect.hgvs_protein = "p." + residue + "fs";
    } else if (insertion) {
      effect.consequence = "inframe_insertion";
      effect.hgvs_protein = "p." + residue + "ins";
    } else {
      effect.consequence = "inframe_deletion";
      effect.hgvs_protein = "p." + residue + "del";
    }
    if (amino != '*') {
      effect.structure = map_structure(gene.name, effect.amino_acid_position);
    }
    effects.push_back(std::move(effect));
  }
  if (effects.empty()) {
    effects.push_back(regional_effect(position));
  }
  return effects;
}

[[nodiscard]] LargeDeletionImpact
annotate_large_deletion(const std::size_t start, const std::size_t end,
                        const std::size_t reference_length) {
  LargeDeletionImpact impact;
  const auto overlaps = [&](const std::size_t low, const std::size_t high) {
    return circular_overlap(start, end, low, high, reference_length);
  };
  const Feature *left_partial = nullptr;
  const Feature *right_partial = nullptr;
  for (const auto &gene : kFeatures) {
    const auto overlap = overlaps(gene.start, gene.end);
    if (overlap == 0U) {
      continue;
    }
    if (gene.biotype == Biotype::tRNA) {
      ++impact.transfer_rnas_affected;
    }
    if (overlap == gene.end - gene.start + 1U) {
      impact.genes_lost.emplace_back(gene.name);
      continue;
    }
    impact.genes_partial.emplace_back(gene.name);
    if (gene.biotype == Biotype::protein_coding && gene.strand == '+') {
      const std::size_t before = start > 1U ? start - 1U : reference_length;
      const std::size_t after = end < reference_length ? end + 1U : 1U;
      if (before >= gene.start && before <= gene.end) {
        left_partial = &gene;
      }
      if (after >= gene.start && after <= gene.end && right_partial == nullptr) {
        right_partial = &gene;
      }
    }
  }
  if (left_partial != nullptr && right_partial != nullptr) {
    const std::size_t left_phase = (start - left_partial->start) % 3U;
    const std::size_t right_phase = (end + 1U - right_partial->start) % 3U;
    impact.fusion = std::string(left_partial->name) +
                    (left_partial == right_partial
                         ? " internal deletion"
                         : "::" + std::string(right_partial->name)) +
                    (left_phase == right_phase ? " (in frame)" : " (out of frame)");
  }
  impact.heavy_strand_origin_lost =
      overlaps(kHeavyOriginStart, kHeavyOriginEnd) > 0U;
  impact.light_strand_origin_lost =
      overlaps(kLightOriginStart, kLightOriginEnd) > 0U;
  impact.major_arc = start <= end && start >= kMajorArcStart && end <= kMajorArcEnd;
  impact.contact_zone =
      start >= 6000U && start <= 9000U && end >= 13000U && end <= 16000U;
  return impact;
}

} // namespace mito::detail
