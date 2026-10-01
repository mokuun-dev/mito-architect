#include "mito_c_api.h"

#include "mito/analysis_engine.hpp"
#include "mito/version.hpp"

#include <cstdlib>
#include <cstring>
#include <cstddef>
#include <exception>
#include <limits>
#include <new>
#include <string>

#include "native/error_state.hpp"

using mito::ffi::clear_error;
using mito::ffi::set_error;

namespace {
[[nodiscard]] bool valid_utf8(const char *value) noexcept {
  if (value == nullptr) {
    return true;
  }
  const auto *bytes = reinterpret_cast<const unsigned char *>(value);
  for (std::size_t index = 0U; bytes[index] != 0U;) {
    const auto byte = bytes[index];
    if (byte < 0x80U) {
      ++index;
      continue;
    }
    const auto continuation = [&](const std::size_t offset) {
      return bytes[index + offset] != 0U &&
             (bytes[index + offset] & 0xc0U) == 0x80U;
    };
    std::size_t length = 0U;
    if (byte >= 0xc2U && byte <= 0xdfU && continuation(1U)) length = 2U;
    else if (byte == 0xe0U && continuation(1U) && continuation(2U) && bytes[index + 1U] >= 0xa0U) length = 3U;
    else if (byte >= 0xe1U && byte <= 0xecU && continuation(1U) && continuation(2U)) length = 3U;
    else if (byte == 0xedU && continuation(1U) && continuation(2U) && bytes[index + 1U] <= 0x9fU) length = 3U;
    else if (byte >= 0xeeU && byte <= 0xefU && continuation(1U) && continuation(2U)) length = 3U;
    else if (byte == 0xf0U && continuation(1U) && continuation(2U) && continuation(3U) && bytes[index + 1U] >= 0x90U) length = 4U;
    else if (byte >= 0xf1U && byte <= 0xf3U && continuation(1U) && continuation(2U) && continuation(3U)) length = 4U;
    else if (byte == 0xf4U && continuation(1U) && continuation(2U) && continuation(3U) && bytes[index + 1U] <= 0x8fU) length = 4U;
    if (length == 0U) return false;
    index += length;
  }
  return true;
}

[[nodiscard]] const char *copy_to_c_string(const std::string &value) noexcept {
  if (value.size() == std::numeric_limits<std::size_t>::max()) {
    set_error("MITO-E1601", "analysis JSON length overflows C allocation");
    return nullptr;
  }
  auto *result = static_cast<char *>(std::malloc(value.size() + 1U));
  if (result == nullptr) {
    set_error("MITO-E1601", "malloc failed while returning analysis JSON");
    return nullptr;
  }
  std::memcpy(result, value.c_str(), value.size() + 1U);
  return result;
}
} // namespace

extern "C" void *mito_engine_new(void) noexcept {
  try {
    clear_error();
    return new mito::AnalysisEngine();
  } catch (const std::bad_alloc &) {
    set_error("MITO-E1601", "allocation failed while creating mito engine");
  } catch (const std::exception &error) {
    set_error("MITO-E9001", error.what());
  } catch (...) {
    set_error("MITO-E9001", "unknown error while creating mito engine");
  }
  return nullptr;
}

extern "C" void mito_engine_delete(void *engine) noexcept {
  delete static_cast<mito::AnalysisEngine *>(engine);
}

extern "C" bool mito_engine_has_htslib(void) noexcept {
#ifdef MITO_HAS_HTSLIB
  return true;
#else
  return false;
#endif
}

extern "C" const char *mito_engine_version(void) noexcept {
  return mito::kEngineVersion;
}

extern "C" const char *mito_engine_schema_version(void) noexcept {
  return mito::kResultSchemaVersion;
}

extern "C" const char *mito_engine_error_schema_version(void) noexcept {
  return mito::kErrorSchemaVersion;
}

extern "C" const char *mito_engine_analyze_with_config_v6(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, std::size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool emit_evidence_graph, std::size_t max_evidence_observations,
    std::size_t max_phase_links, std::size_t evidence_page_size,
    const char *molecule_id_tag, const char *umi_tag, const char *duplex_tag,
    std::size_t min_architecture_molecules,
    std::size_t max_candidate_architectures, double architecture_max_distance,
    double architecture_ambiguity_margin,
    double architecture_min_overlap_fraction,
    double architecture_consensus_fraction,
    double architecture_optional_fraction,
    std::size_t architecture_stability_replicates,
    std::uint64_t architecture_seed, bool (*should_cancel)(void *),
    void *cancel_user_data) noexcept {
  if (engine == nullptr) {
    set_error("MITO-E1001", "mito_engine_analyze received a null engine");
    return nullptr;
  }
  if (input_path == nullptr) {
    set_error("MITO-E1001", "mito_engine_analyze received a null input path");
    return nullptr;
  }

  try {
    clear_error();
    mito::AnalysisConfig config;
    config.filter_numt = filter_numt;
    config.threads = threads == 0 ? 1 : threads;
    config.min_mapping_quality = min_mapping_quality;
    config.min_base_quality = min_base_quality;
    config.excluded_snp_flags = excluded_snp_flags;
    config.numt_threshold = numt_threshold;
    config.allow_development_tags = allow_development_tags;
    config.result_schema = emit_evidence_graph ? mito::ResultSchema::v0_6
                                               : mito::ResultSchema::v0_5;
    config.max_evidence_observations = max_evidence_observations;
    config.max_phase_links = max_phase_links;
    config.evidence_page_size = evidence_page_size;
    config.min_architecture_molecules = min_architecture_molecules;
    config.max_candidate_architectures = max_candidate_architectures;
    config.architecture_max_distance = architecture_max_distance;
    config.architecture_ambiguity_margin = architecture_ambiguity_margin;
    config.architecture_min_overlap_fraction =
        architecture_min_overlap_fraction;
    config.architecture_consensus_fraction = architecture_consensus_fraction;
    config.architecture_optional_fraction = architecture_optional_fraction;
    config.architecture_stability_replicates =
        architecture_stability_replicates;
    config.architecture_seed = architecture_seed;
    config.molecule_id_tag = molecule_id_tag == nullptr
                                 ? std::string{}
                                 : std::string(molecule_id_tag);
    config.umi_tag = umi_tag == nullptr ? std::string{} : std::string(umi_tag);
    config.duplex_tag =
        duplex_tag == nullptr ? std::string{} : std::string(duplex_tag);
    if (should_cancel != nullptr) {
      config.should_cancel = [should_cancel, cancel_user_data] {
        return should_cancel(cancel_user_data);
      };
    }
    const std::string reference_path =
        ref_path == nullptr ? std::string{} : std::string(ref_path);
    const auto json = static_cast<mito::AnalysisEngine *>(engine)->analyze(
        input_path, reference_path, config);
    return copy_to_c_string(json);
  } catch (const mito::AnalysisError &error) {
    const auto code = mito::analysis_error_code_name(error.code());
    set_error(code, error.what());
  } catch (const std::bad_alloc &) {
    set_error("MITO-E1601", "allocation failed during mtDNA analysis");
  } catch (const std::exception &error) {
    set_error("MITO-E9001", error.what());
  } catch (...) {
    set_error("MITO-E9001", "unknown error during mtDNA analysis");
  }
  return nullptr;
}

extern "C" void mito_engine_analyze_options_v1_init(
    mito_engine_analyze_options_v1 *options) noexcept {
  if (options == nullptr) return;
  *options = {};
  options->struct_size = sizeof(*options);
  options->abi_version = MITO_ENGINE_ANALYZE_OPTIONS_ABI_V1;
  options->filter_numt = true;
  options->threads = 1U;
  options->min_mapping_quality = 20U;
  options->min_base_quality = 10U;
  options->excluded_snp_flags = 0xF00U;
  options->numt_threshold = 0.30;
  options->max_evidence_observations = 5'000'000U;
  options->max_phase_links = 1'000'000U;
  options->max_phase_work = 20'000'000U;
  options->max_phase_molecule_references = 5'000'000U;
  options->max_result_bytes = 128U * 1024U * 1024U;
  options->evidence_page_size = 4096U;
  options->min_architecture_molecules = 2U;
  options->max_candidate_architectures = 256U;
  options->architecture_max_distance = 0.20;
  options->architecture_ambiguity_margin = 0.05;
  options->architecture_min_overlap_fraction = 0.50;
  options->architecture_consensus_fraction = 0.80;
  options->architecture_optional_fraction = 0.20;
  options->architecture_stability_replicates = 32U;
  options->architecture_seed = 0x4d49544f41524348ULL;
}

extern "C" const char *mito_engine_analyze_with_options_v1(
    void *engine, const char *input_path, const char *ref_path,
    const mito_engine_analyze_options_v1 *options) noexcept {
  if (engine == nullptr || input_path == nullptr || options == nullptr) {
    set_error("MITO-E1001", "mito_engine options call received a null required argument");
    return nullptr;
  }
  if (options->abi_version != MITO_ENGINE_ANALYZE_OPTIONS_ABI_V1 ||
      options->struct_size != sizeof(*options)) {
    set_error("MITO-E1001", "unsupported mito_engine_analyze_options_v1 layout");
    return nullptr;
  }
  if (!valid_utf8(input_path) || !valid_utf8(ref_path) ||
      !valid_utf8(options->molecule_id_tag) || !valid_utf8(options->umi_tag) ||
      !valid_utf8(options->duplex_tag)) {
    set_error("MITO-E1001", "C ABI string arguments must be valid UTF-8");
    return nullptr;
  }
  try {
    clear_error();
    mito::AnalysisConfig config;
    config.filter_numt = options->filter_numt;
    config.threads = options->threads == 0U ? 1U : options->threads;
    config.min_mapping_quality = options->min_mapping_quality;
    config.min_base_quality = options->min_base_quality;
    config.excluded_snp_flags = options->excluded_snp_flags;
    config.numt_threshold = options->numt_threshold;
    config.allow_development_tags = options->allow_development_tags;
    config.result_schema = options->emit_evidence_graph ? mito::ResultSchema::v0_6
                                                         : mito::ResultSchema::v0_5;
    config.max_evidence_observations = options->max_evidence_observations;
    config.max_phase_links = options->max_phase_links;
    config.max_phase_work = options->max_phase_work;
    config.max_phase_molecule_references = options->max_phase_molecule_references;
    config.max_result_bytes = options->max_result_bytes;
    config.evidence_page_size = options->evidence_page_size;
    config.molecule_id_tag = options->molecule_id_tag == nullptr ? "" : options->molecule_id_tag;
    config.umi_tag = options->umi_tag == nullptr ? "" : options->umi_tag;
    config.duplex_tag = options->duplex_tag == nullptr ? "" : options->duplex_tag;
    config.min_architecture_molecules = options->min_architecture_molecules;
    config.max_candidate_architectures = options->max_candidate_architectures;
    config.architecture_max_distance = options->architecture_max_distance;
    config.architecture_ambiguity_margin = options->architecture_ambiguity_margin;
    config.architecture_min_overlap_fraction = options->architecture_min_overlap_fraction;
    config.architecture_consensus_fraction = options->architecture_consensus_fraction;
    config.architecture_optional_fraction = options->architecture_optional_fraction;
    config.architecture_stability_replicates = options->architecture_stability_replicates;
    config.architecture_seed = options->architecture_seed;
    if (options->should_cancel != nullptr) {
      config.should_cancel = [callback = options->should_cancel,
                              user_data = options->cancel_user_data] {
        return callback(user_data);
      };
    }
    const auto json = static_cast<mito::AnalysisEngine *>(engine)->analyze(
        input_path, ref_path == nullptr ? std::string{} : std::string(ref_path), config);
    return copy_to_c_string(json);
  } catch (const mito::AnalysisError &error) {
    set_error(mito::analysis_error_code_name(error.code()), error.what());
  } catch (const std::bad_alloc &) {
    set_error("MITO-E1601", "allocation failed during mtDNA analysis");
  } catch (const std::exception &error) {
    set_error("MITO-E9001", error.what());
  } catch (...) {
    set_error("MITO-E9001", "unknown error during mtDNA analysis");
  }
  return nullptr;
}

extern "C" const char *mito_engine_get_last_error(void) noexcept {
  return mito::ffi::last_error();
}

extern "C" const char *mito_engine_get_last_error_code(void) noexcept {
  return mito::ffi::last_error_code();
}

extern "C" void mito_engine_free_string(const char *value) noexcept {
  std::free(const_cast<char *>(value));
}
