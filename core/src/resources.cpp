#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] std::vector<std::string> split_tab(std::string_view line) {
  std::vector<std::string> fields;
  std::size_t start = 0;
  while (start <= line.size()) {
    const std::size_t next = line.find('\t', start);
    if (next == std::string_view::npos) {
      fields.emplace_back(line.substr(start));
      break;
    }
    fields.emplace_back(line.substr(start, next - start));
    start = next + 1;
  }
  return fields;
}

[[nodiscard]] std::string trim_copy(std::string_view value) {
  std::size_t first = 0;
  while (first < value.size() &&
         std::isspace(static_cast<unsigned char>(value[first])) != 0) {
    ++first;
  }
  std::size_t last = value.size();
  while (last > first &&
         std::isspace(static_cast<unsigned char>(value[last - 1])) != 0) {
    --last;
  }
  return std::string(value.substr(first, last - first));
}

void strip_trailing_carriage_return(std::string &value) {
  if (!value.empty() && value.back() == '\r') {
    value.pop_back();
  }
}

[[nodiscard]] std::vector<std::string> split_semicolon(std::string_view line) {
  std::vector<std::string> fields;
  std::size_t start = 0;
  while (start <= line.size()) {
    const std::size_t next = line.find(';', start);
    const auto token = trim_copy(
        line.substr(start, next == std::string_view::npos ? line.size() - start
                                                          : next - start));
    if (!token.empty()) {
      fields.push_back(token);
    }
    if (next == std::string_view::npos) {
      break;
    }
    start = next + 1;
  }
  return fields;
}

[[nodiscard]] std::string parse_reference_fasta(std::istream &in,
                                                const std::string &source) {
  std::string reference;
  std::string line;
  std::size_t header_count = 0;
  while (std::getline(in, line)) {
    if (!line.empty() && line.front() == '>') {
      ++header_count;
      if (header_count > 1) {
        throw AnalysisError(
            AnalysisErrorCode::reference_invalid,
            "reference FASTA must contain exactly one sequence: " + source);
      }
      continue;
    }
    for (const char base : line) {
      if (std::isspace(static_cast<unsigned char>(base)) != 0) {
        continue;
      }
      const char upper =
          static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
      if (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T' ||
          upper == 'N') {
        reference.push_back(upper);
      } else {
        throw AnalysisError(AnalysisErrorCode::reference_invalid,
                            "reference FASTA contains an unsupported base: " +
                                source);
      }
    }
  }

  if (header_count == 0) {
    throw AnalysisError(AnalysisErrorCode::reference_invalid,
                        "reference FASTA is missing a header: " + source);
  }
  if (reference.empty()) {
    throw AnalysisError(AnalysisErrorCode::reference_invalid,
                        "reference FASTA is empty: " + source);
  }
  return reference;
}

[[nodiscard]] std::string bundled_rcrs_path() {
#ifdef MITO_DEFAULT_RCRS_FASTA_PATH
  return MITO_DEFAULT_RCRS_FASTA_PATH;
#else
  return "core/data/rcrs.fasta";
#endif
}

[[nodiscard]] std::string bundled_clinical_annotations_path() {
#ifdef MITO_CLINICAL_ANNOTATIONS_PATH
  return MITO_CLINICAL_ANNOTATIONS_PATH;
#else
  return "core/data/clinical_annotations.tsv";
#endif
}

[[nodiscard]] std::string clinical_annotations_path() {
  if (const char *override_path = std::getenv("MITO_CLINICAL_ANNOTATIONS")) {
    if (override_path[0] != '\0') {
      return override_path;
    }
  }
  return bundled_clinical_annotations_path();
}

[[nodiscard]] std::string phylotree_path() {
  if (const char *override_path = std::getenv("MITO_PHYLOTREE")) {
    if (override_path[0] != '\0') {
      return override_path;
    }
  }
#ifdef MITO_PHYLOTREE_PATH
  return MITO_PHYLOTREE_PATH;
#else
  return "core/data/phylotree-rcrs-17.3.xml";
#endif
}

[[nodiscard]] std::string phylotree_weights_path() {
  if (const char *override_path = std::getenv("MITO_PHYLOTREE_WEIGHTS")) {
    if (override_path[0] != '\0') {
      return override_path;
    }
  }
#ifdef MITO_PHYLOTREE_WEIGHTS_PATH
  return MITO_PHYLOTREE_WEIGHTS_PATH;
#else
  return "core/data/phylotree-rcrs-17.3-weights.txt";
#endif
}

[[nodiscard]] std::string phylotree_alignment_rules_path() {
  if (const char *override_path =
          std::getenv("MITO_PHYLOTREE_ALIGNMENT_RULES")) {
    if (override_path[0] != '\0') {
      return override_path;
    }
  }
#ifdef MITO_PHYLOTREE_ALIGNMENT_RULES_PATH
  return MITO_PHYLOTREE_ALIGNMENT_RULES_PATH;
#else
  return "core/data/phylotree-rcrs-17.3-rules.csv";
#endif
}

[[nodiscard]] std::string resource_manifest_path() {
#ifdef MITO_RESOURCE_MANIFEST_PATH
  return MITO_RESOURCE_MANIFEST_PATH;
#else
  return "core/data/resource_manifest.tsv";
#endif
}

[[nodiscard]] std::vector<ResourceRecord> load_resource_manifest() {
  std::vector<ResourceRecord> records;
  std::ifstream input(resource_manifest_path());
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::resource_open_failed,
                        "could not open resource manifest: " +
                            resource_manifest_path());
  }
  std::string line;
  bool header = true;
  while (std::getline(input, line)) {
    strip_trailing_carriage_return(line);
    if (line.empty() || header) {
      header = false;
      continue;
    }
    const auto fields = split_tab(line);
    if (fields.size() != 7) {
      throw AnalysisError(AnalysisErrorCode::resource_invalid,
                          "resource manifest row must contain seven fields");
    }
    records.push_back({fields[0], fields[1], fields[2], fields[3], fields[4],
                       fields[5], fields[6]});
  }
  if (records.empty()) {
    throw AnalysisError(AnalysisErrorCode::resource_invalid,
                        "resource manifest contains no resources");
  }
  return records;
}

void apply_runtime_resource_overrides(std::vector<ResourceRecord> &records) {
  const auto clinical_path = clinical_annotations_path();
  if (clinical_path == bundled_clinical_annotations_path()) {
    return;
  }
  const auto clinical =
      std::find_if(records.begin(), records.end(), [](const auto &resource) {
        return resource.name == "clinical-curated";
      });
  if (clinical == records.end()) {
    throw AnalysisError(AnalysisErrorCode::resource_invalid,
                        "resource manifest is missing clinical-curated");
  }
  clinical->version = "external-unpinned";
  clinical->path = clinical_path;
  clinical->sha256 = "not-verified";
  clinical->source = "MITO_CLINICAL_ANNOTATIONS override";
  clinical->license = "not-recorded";
  clinical->retrieved = "not-recorded";
}

void throw_if_cancelled(const AnalysisConfig &config) {
  if (config.should_cancel && config.should_cancel()) {
    throw AnalysisError(AnalysisErrorCode::analysis_cancelled,
                        "analysis cancelled");
  }
}

[[nodiscard]] std::string load_reference(const std::string &reference_path) {
  const std::string effective_path =
      reference_path.empty() ? bundled_rcrs_path() : reference_path;
  std::ifstream in(effective_path);
  if (!in) {
    throw AnalysisError(AnalysisErrorCode::reference_open_failed,
                        "could not open reference FASTA: " + effective_path);
  }

  auto reference = parse_reference_fasta(in, effective_path);
  if (reference.size() != static_cast<std::size_t>(kDefaultReferenceLength) &&
      reference_path.empty()) {
    throw AnalysisError(AnalysisErrorCode::reference_invalid,
                        "bundled rCRS reference length is not 16569 bp: " +
                            std::to_string(reference.size()));
  }
  return reference;
}

[[nodiscard]] std::vector<GeneAnnotation> default_genes() {
  return {
      {"MT-RNR1", 648, 1601, "+", "rRNA"},
      {"MT-RNR2", 1671, 3229, "+", "rRNA"},
      {"MT-ND1", 3307, 4262, "+", "protein_coding"},
      {"MT-ND2", 4470, 5511, "+", "protein_coding"},
      {"MT-CO1", 5904, 7445, "+", "protein_coding"},
      {"MT-CO2", 7586, 8269, "+", "protein_coding"},
      {"MT-ATP8", 8366, 8572, "+", "protein_coding"},
      {"MT-ATP6", 8527, 9207, "+", "protein_coding"},
      {"MT-CO3", 9207, 9990, "+", "protein_coding"},
      {"MT-ND3", 10059, 10404, "+", "protein_coding"},
      {"MT-ND4L", 10470, 10766, "+", "protein_coding"},
      {"MT-ND4", 10760, 12137, "+", "protein_coding"},
      {"MT-ND5", 12337, 14148, "+", "protein_coding"},
      {"MT-ND6", 14149, 14673, "-", "protein_coding"},
      {"MT-CYB", 14747, 15887, "+", "protein_coding"},
  };
}

} // namespace mito::detail
