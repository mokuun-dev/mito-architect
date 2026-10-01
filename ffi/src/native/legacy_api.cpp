#include "mito_c_api.h"
#include <cstddef>

extern "C" const char *mito_engine_analyze(void *engine, const char *input_path,
                                           const char *ref_path) noexcept {
  return mito_engine_analyze_with_options(engine, input_path, ref_path, true,
                                          1);
}

extern "C" const char *
mito_engine_analyze_with_options(void *engine, const char *input_path,
                                 const char *ref_path, bool filter_numt,
                                 std::size_t threads) noexcept {
  return mito_engine_analyze_with_cancel(
      engine, input_path, ref_path, filter_numt, threads, nullptr, nullptr);
}

extern "C" const char *mito_engine_analyze_with_cancel(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, std::size_t threads, bool (*should_cancel)(void *),
    void *cancel_user_data) noexcept {
  return mito_engine_analyze_with_config(engine, input_path, ref_path,
                                         filter_numt, threads, 20, 10, 0xF00,
                                         should_cancel, cancel_user_data);
}

extern "C" const char *mito_engine_analyze_with_config(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, std::size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    bool (*should_cancel)(void *), void *cancel_user_data) noexcept {
  return mito_engine_analyze_with_config_v2(
      engine, input_path, ref_path, filter_numt, threads, min_mapping_quality,
      min_base_quality, excluded_snp_flags, 0.30, false, should_cancel,
      cancel_user_data);
}

extern "C" const char *mito_engine_analyze_with_config_v2(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, std::size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool (*should_cancel)(void *), void *cancel_user_data) noexcept {
  return mito_engine_analyze_with_config_v3(
      engine, input_path, ref_path, filter_numt, threads, min_mapping_quality,
      min_base_quality, excluded_snp_flags, numt_threshold,
      allow_development_tags, false, 5'000'000U, should_cancel,
      cancel_user_data);
}

extern "C" const char *mito_engine_analyze_with_config_v3(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, std::size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool emit_evidence_graph, std::size_t max_evidence_observations,
    bool (*should_cancel)(void *), void *cancel_user_data) noexcept {
  return mito_engine_analyze_with_config_v4(
      engine, input_path, ref_path, filter_numt, threads, min_mapping_quality,
      min_base_quality, excluded_snp_flags, numt_threshold,
      allow_development_tags, emit_evidence_graph, max_evidence_observations,
      1'000'000U, should_cancel, cancel_user_data);
}

extern "C" const char *mito_engine_analyze_with_config_v4(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, std::size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool emit_evidence_graph, std::size_t max_evidence_observations,
    std::size_t max_phase_links, bool (*should_cancel)(void *),
    void *cancel_user_data) noexcept {
  return mito_engine_analyze_with_config_v5(
      engine, input_path, ref_path, filter_numt, threads, min_mapping_quality,
      min_base_quality, excluded_snp_flags, numt_threshold,
      allow_development_tags, emit_evidence_graph, max_evidence_observations,
      max_phase_links, 4096U, nullptr, nullptr, nullptr, should_cancel,
      cancel_user_data);
}

extern "C" const char *mito_engine_analyze_with_config_v5(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, std::size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool emit_evidence_graph, std::size_t max_evidence_observations,
    std::size_t max_phase_links, std::size_t evidence_page_size,
    const char *molecule_id_tag, const char *umi_tag, const char *duplex_tag,
    bool (*should_cancel)(void *), void *cancel_user_data) noexcept {
  return mito_engine_analyze_with_config_v6(
      engine, input_path, ref_path, filter_numt, threads, min_mapping_quality,
      min_base_quality, excluded_snp_flags, numt_threshold,
      allow_development_tags, emit_evidence_graph, max_evidence_observations,
      max_phase_links, evidence_page_size, molecule_id_tag, umi_tag, duplex_tag,
      2U, 256U, 0.20, 0.05, 0.50, 0.80, 0.20, 32U, 0x4d49544f41524348ULL,
      should_cancel, cancel_user_data);
}
