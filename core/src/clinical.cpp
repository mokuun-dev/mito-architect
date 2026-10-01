#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] std::string lowercase_token(std::string_view value) {
  std::string normalized;
  normalized.reserve(value.size());
  for (const char c : value) {
    const auto byte = static_cast<unsigned char>(c);
    if (std::isspace(byte) != 0 || c == '-') {
      if (normalized.empty() || normalized.back() != '_') {
        normalized.push_back('_');
      }
    } else {
      normalized.push_back(static_cast<char>(std::tolower(byte)));
    }
  }
  while (!normalized.empty() && normalized.back() == '_') {
    normalized.pop_back();
  }
  return normalized;
}

[[nodiscard]] std::string
normalize_clinical_significance(std::string_view significance) {
  const auto value = lowercase_token(trim_copy(significance));
  if (value.empty() || value == "not_provided" || value == "not_specified") {
    return "not_provided";
  }
  if (value.find("conflict") != std::string::npos ||
      (value.find("pathogenic") != std::string::npos &&
       value.find("benign") != std::string::npos)) {
    return "conflicting";
  }
  if (value.find("pathogenic/likely_pathogenic") != std::string::npos ||
      value.find("pathogenic,_likely_pathogenic") != std::string::npos) {
    return "pathogenic_or_likely_pathogenic";
  }
  if (value.find("likely_pathogenic") != std::string::npos) {
    return "likely_pathogenic";
  }
  if (value.find("pathogenic") != std::string::npos) {
    return "pathogenic";
  }
  if (value.find("benign/likely_benign") != std::string::npos ||
      value.find("benign,_likely_benign") != std::string::npos) {
    return "benign_or_likely_benign";
  }
  if (value.find("likely_benign") != std::string::npos) {
    return "likely_benign";
  }
  if (value.find("benign") != std::string::npos) {
    return "benign";
  }
  if (value.find("uncertain") != std::string::npos || value == "vus") {
    return "uncertain_significance";
  }
  for (const auto category : {"risk_factor", "association", "drug_response",
                              "protective", "affects"}) {
    if (value.find(category) != std::string::npos) {
      return category;
    }
  }
  return "other";
}

[[nodiscard]] std::string
clinical_significance_group(std::string_view normalized) {
  if (normalized == "pathogenic" || normalized == "likely_pathogenic" ||
      normalized == "pathogenic_or_likely_pathogenic") {
    return "pathogenic";
  }
  if (normalized == "benign" || normalized == "likely_benign" ||
      normalized == "benign_or_likely_benign") {
    return "benign";
  }
  if (normalized == "uncertain_significance") {
    return "uncertain";
  }
  if (normalized == "conflicting") {
    return "conflicting";
  }
  if (normalized == "not_provided" || normalized == "other") {
    return "unclassified";
  }
  return "other_assertion";
}

[[nodiscard]] std::string join_values(const std::vector<std::string> &values,
                                      std::string_view delimiter) {
  std::string joined;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i != 0U) {
      joined += delimiter;
    }
    joined += values[i];
  }
  return joined;
}

void finalize_clinical_annotation(SnpCall &call) {
  auto &assertions = call.clinical_assertions;
  std::sort(assertions.begin(), assertions.end());
  assertions.erase(std::unique(assertions.begin(), assertions.end()),
                   assertions.end());

  std::vector<std::string> diseases;
  std::vector<std::string> normalized_significances;
  std::vector<std::string> significance_groups;
  call.references.clear();
  call.sources.clear();
  call.clinvar_allele_id.clear();
  call.mitomap_url.clear();
  for (const auto &assertion : assertions) {
    if (!assertion.disease.empty()) {
      diseases.push_back(assertion.disease);
    }
    if (assertion.normalized_significance != "not_provided") {
      normalized_significances.push_back(assertion.normalized_significance);
    }
    const auto group =
        clinical_significance_group(assertion.normalized_significance);
    if (group != "unclassified") {
      significance_groups.push_back(group);
    }
    call.sources.push_back(assertion.source);
    call.references.insert(call.references.end(), assertion.references.begin(),
                           assertion.references.end());
    if (assertion.source == "ClinVar" && call.clinvar_allele_id.empty()) {
      call.clinvar_allele_id = assertion.allele_id;
    }
    if (assertion.source == "MITOMAP" && call.mitomap_url.empty()) {
      call.mitomap_url = assertion.source_url;
    }
  }
  for (auto *values : {&diseases, &normalized_significances,
                       &significance_groups, &call.references, &call.sources}) {
    std::sort(values->begin(), values->end());
    values->erase(std::unique(values->begin(), values->end()), values->end());
  }

  const bool explicit_conflict =
      std::find(significance_groups.begin(), significance_groups.end(),
                "conflicting") != significance_groups.end();
  const std::size_t classified_group_count =
      static_cast<std::size_t>(std::count_if(
          significance_groups.begin(), significance_groups.end(),
          [](const std::string &group) { return group != "conflicting"; }));
  const bool incompatible_groups = classified_group_count > 1U;
  if (explicit_conflict || incompatible_groups) {
    call.clinical_conflict_status = "conflicting";
    call.clinical_consensus_significance = "conflicting";
  } else if (assertions.size() == 1U) {
    call.clinical_conflict_status = "single_assertion";
  } else {
    call.clinical_conflict_status = "consistent";
  }

  if (call.clinical_consensus_significance.empty()) {
    if (normalized_significances.empty()) {
      call.clinical_consensus_significance = "not_provided";
    } else if (normalized_significances.size() == 1U) {
      call.clinical_consensus_significance = normalized_significances.front();
    } else if (significance_groups.size() == 1U &&
               significance_groups.front() == "pathogenic") {
      call.clinical_consensus_significance = "pathogenic_or_likely_pathogenic";
    } else if (significance_groups.size() == 1U &&
               significance_groups.front() == "benign") {
      call.clinical_consensus_significance = "benign_or_likely_benign";
    } else {
      call.clinical_consensus_significance = normalized_significances.front();
    }
  }
  call.phenotype = join_values(diseases, "; ");
  call.pathogenicity = call.clinical_consensus_significance;
}

[[nodiscard]] std::map<std::string, SnpCall> load_clinical_annotations() {
  std::map<std::string, SnpCall> annotations;
  std::ifstream input(clinical_annotations_path());
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "could not open clinical annotations: " +
                            clinical_annotations_path());
  }

  constexpr std::string_view expected_header =
      "position\tref\talt\tgene\tconsequence\tprotein\tresidue\tstructure_id\t"
      "structure_chain\tstructure_residue\tstructure_"
      "complex\tsource\tassertion_id\t"
      "allele_id\tdisease\tclinical_significance\treview_status\tassertion_"
      "date\t"
      "source_url\treferences\tresource_version\tretrieved_at";
  std::string line;
  std::size_t line_number = 0U;
  bool header_seen = false;
  while (std::getline(input, line)) {
    ++line_number;
    strip_trailing_carriage_return(line);
    if (line.empty()) {
      continue;
    }
    if (!header_seen) {
      header_seen = true;
      if (line != expected_header) {
        throw AnalysisError(
            AnalysisErrorCode::resource_invalid,
            "clinical annotation header does not match schema 1.0");
      }
      continue;
    }

    auto fields = split_tab(line);
    if (fields.size() != 22U) {
      throw AnalysisError(
          AnalysisErrorCode::resource_invalid,
          "clinical annotation row must contain 22 fields at line " +
              std::to_string(line_number));
    }
    for (auto &field : fields) {
      field = trim_copy(field);
    }
    const auto position = parse_size(fields[0]);
    const auto valid_allele = [](const std::string &allele) {
      if (allele.size() != 1U) {
        return false;
      }
      const char base = static_cast<char>(
          std::toupper(static_cast<unsigned char>(allele.front())));
      return base == 'A' || base == 'C' || base == 'G' || base == 'T';
    };
    if (!position || *position == 0U ||
        *position > static_cast<std::size_t>(kDefaultReferenceLength) ||
        !valid_allele(fields[1]) || !valid_allele(fields[2]) ||
        std::toupper(static_cast<unsigned char>(fields[1].front())) ==
            std::toupper(static_cast<unsigned char>(fields[2].front())) ||
        fields[11].empty()) {
      throw AnalysisError(AnalysisErrorCode::resource_invalid,
                          "invalid clinical variant/assertion key at line " +
                              std::to_string(line_number));
    }
    if (!fields[18].empty() && !fields[18].starts_with("https://") &&
        !fields[18].starts_with("http://")) {
      throw AnalysisError(AnalysisErrorCode::resource_invalid,
                          "clinical source_url must be HTTP(S) at line " +
                              std::to_string(line_number));
    }

    SnpCall call;
    call.position = *position;
    call.reference = static_cast<char>(
        std::toupper(static_cast<unsigned char>(fields[1].front())));
    call.alternate = static_cast<char>(
        std::toupper(static_cast<unsigned char>(fields[2].front())));
    call.gene = fields[3];
    call.consequence = fields[4];
    call.protein = fields[5];
    call.residue = fields[6];
    call.structure_id = fields[7];
    call.structure_chain = fields[8];
    if (!fields[9].empty()) {
      const auto residue = parse_size(fields[9]);
      if (!residue || *residue == 0U) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "invalid clinical structure residue at line " +
                                std::to_string(line_number));
      }
      call.structure_residue = *residue;
    }
    call.structure_complex = fields[10];
    ClinicalAssertion assertion;
    assertion.source = fields[11];
    assertion.assertion_id = fields[12];
    assertion.allele_id = fields[13];
    assertion.disease = fields[14];
    assertion.clinical_significance = fields[15];
    assertion.normalized_significance =
        normalize_clinical_significance(fields[15]);
    assertion.review_status = fields[16];
    assertion.assertion_date = fields[17];
    assertion.source_url = fields[18];
    assertion.references = split_semicolon(fields[19]);
    assertion.resource_version = fields[20];
    assertion.retrieved_at = fields[21];
    call.clinical_assertions.push_back(std::move(assertion));

    const auto key = snp_key(call.position, call.reference, call.alternate);
    auto [it, inserted] = annotations.try_emplace(key, call);
    if (!inserted) {
      auto &existing = it->second;
      const auto merge_scalar = [&](std::string &target,
                                    const std::string &incoming,
                                    std::string_view name) {
        if (target.empty()) {
          target = incoming;
        } else if (!incoming.empty() && target != incoming) {
          throw AnalysisError(AnalysisErrorCode::resource_invalid,
                              "conflicting clinical " + std::string(name) +
                                  " at line " + std::to_string(line_number));
        }
      };
      merge_scalar(existing.gene, call.gene, "gene");
      merge_scalar(existing.consequence, call.consequence, "consequence");
      merge_scalar(existing.protein, call.protein, "protein");
      merge_scalar(existing.residue, call.residue, "residue");
      merge_scalar(existing.structure_id, call.structure_id, "structure_id");
      merge_scalar(existing.structure_chain, call.structure_chain,
                   "structure_chain");
      merge_scalar(existing.structure_complex, call.structure_complex,
                   "structure_complex");
      if (existing.structure_residue == 0U) {
        existing.structure_residue = call.structure_residue;
      } else if (call.structure_residue != 0U &&
                 existing.structure_residue != call.structure_residue) {
        throw AnalysisError(AnalysisErrorCode::resource_invalid,
                            "conflicting clinical structure_residue at line " +
                                std::to_string(line_number));
      }
      existing.clinical_assertions.insert(existing.clinical_assertions.end(),
                                          call.clinical_assertions.begin(),
                                          call.clinical_assertions.end());
    }
  }

  if (!header_seen || annotations.empty()) {
    throw AnalysisError(AnalysisErrorCode::resource_invalid,
                        "clinical annotation resource contains no assertions");
  }
  for (auto &[_, call] : annotations) {
    finalize_clinical_annotation(call);
  }
  return annotations;
}

void apply_clinical_annotations(
    std::vector<SnpCall> &snps,
    const std::map<std::string, SnpCall> &annotations) {
  for (auto &snp : snps) {
    const auto it =
        annotations.find(snp_key(snp.position, snp.reference, snp.alternate));
    if (it == annotations.end()) {
      continue;
    }
    const auto position = snp.position;
    const auto reference = snp.reference;
    const auto alternate = snp.alternate;
    snp = it->second;
    snp.position = position;
    snp.reference = reference;
    snp.alternate = alternate;
  }
}

} // namespace mito::detail
