#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#define MITO_NOEXCEPT noexcept
#else
#define MITO_NOEXCEPT
#endif

/* Ownership and lifetime contract:
 * - engine must be a live handle from new; delete(NULL) is valid. Delete once,
 *   after all calls have returned. Foreign/stale pointers are caller errors.
 * - input_path is required. All strings are NUL-terminated and borrowed for
 *   this synchronous call; NULL reference/tags select defaults.
 * - successful analysis strings are independent owned allocations: free each
 *   exactly once with mito_engine_free_string (NULL is accepted), even after
 *   deleting the engine. Version and error strings are borrowed, never freed.
 * - error state is per calling thread; copy it before the next new/analyze on
 *   that thread. Diagnostics truncate at 2047 bytes; error codes do not change.
 * - callbacks may run on worker threads; user_data must remain alive and be
 *   thread-safe until analysis returns. Foreign callbacks must never unwind.
 * - concurrent calls must not mutate process environment/resource files.
 * All exported functions contain native exceptions. */
void *mito_engine_new(void) MITO_NOEXCEPT;
void mito_engine_delete(void *engine) MITO_NOEXCEPT;
bool mito_engine_has_htslib(void) MITO_NOEXCEPT;
const char *mito_engine_version(void) MITO_NOEXCEPT;
const char *mito_engine_schema_version(void) MITO_NOEXCEPT;
const char *mito_engine_error_schema_version(void) MITO_NOEXCEPT;
const char *mito_engine_analyze(void *engine, const char *input_path,
                                const char *ref_path) MITO_NOEXCEPT;
const char *mito_engine_analyze_with_options(void *engine,
                                             const char *input_path,
                                             const char *ref_path,
                                             bool filter_numt,
                                             size_t threads) MITO_NOEXCEPT;
const char *
mito_engine_analyze_with_cancel(void *engine, const char *input_path,
                                const char *ref_path, bool filter_numt,
                                size_t threads, bool (*should_cancel)(void *),
                                void *cancel_user_data) MITO_NOEXCEPT;
const char *mito_engine_analyze_with_config(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    bool (*should_cancel)(void *), void *cancel_user_data) MITO_NOEXCEPT;
const char *mito_engine_analyze_with_config_v2(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool (*should_cancel)(void *), void *cancel_user_data) MITO_NOEXCEPT;
const char *mito_engine_analyze_with_config_v3(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool emit_evidence_graph, size_t max_evidence_observations,
    bool (*should_cancel)(void *), void *cancel_user_data) MITO_NOEXCEPT;
const char *mito_engine_analyze_with_config_v4(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool emit_evidence_graph, size_t max_evidence_observations,
    size_t max_phase_links, bool (*should_cancel)(void *),
    void *cancel_user_data) MITO_NOEXCEPT;
const char *mito_engine_analyze_with_config_v5(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool emit_evidence_graph, size_t max_evidence_observations,
    size_t max_phase_links, size_t evidence_page_size,
    const char *molecule_id_tag, const char *umi_tag, const char *duplex_tag,
    bool (*should_cancel)(void *), void *cancel_user_data) MITO_NOEXCEPT;
const char *mito_engine_analyze_with_config_v6(
    void *engine, const char *input_path, const char *ref_path,
    bool filter_numt, size_t threads, unsigned char min_mapping_quality,
    unsigned char min_base_quality, unsigned short excluded_snp_flags,
    double numt_threshold, bool allow_development_tags,
    bool emit_evidence_graph, size_t max_evidence_observations,
    size_t max_phase_links, size_t evidence_page_size,
    const char *molecule_id_tag, const char *umi_tag, const char *duplex_tag,
    size_t min_architecture_molecules, size_t max_candidate_architectures,
    double architecture_max_distance, double architecture_ambiguity_margin,
    double architecture_min_overlap_fraction,
    double architecture_consensus_fraction,
    double architecture_optional_fraction,
    size_t architecture_stability_replicates, uint64_t architecture_seed,
    bool (*should_cancel)(void *), void *cancel_user_data) MITO_NOEXCEPT;
/* Extensible replacement for positional config_v* calls. Initialize with the
 * function below; the ABI version and exact structure size are checked before
 * any option is read. */
#define MITO_ENGINE_ANALYZE_OPTIONS_ABI_V1 UINT32_C(1)
typedef struct mito_engine_analyze_options_v1 {
  uint32_t struct_size;
  uint32_t abi_version;
  bool filter_numt;
  size_t threads;
  unsigned char min_mapping_quality;
  unsigned char min_base_quality;
  unsigned short excluded_snp_flags;
  double numt_threshold;
  bool allow_development_tags;
  bool emit_evidence_graph;
  size_t max_evidence_observations;
  size_t max_phase_links;
  size_t max_phase_work;
  size_t max_phase_molecule_references;
  size_t max_result_bytes;
  size_t evidence_page_size;
  const char *molecule_id_tag;
  const char *umi_tag;
  const char *duplex_tag;
  size_t min_architecture_molecules;
  size_t max_candidate_architectures;
  double architecture_max_distance;
  double architecture_ambiguity_margin;
  double architecture_min_overlap_fraction;
  double architecture_consensus_fraction;
  double architecture_optional_fraction;
  size_t architecture_stability_replicates;
  uint64_t architecture_seed;
  bool (*should_cancel)(void *);
  void *cancel_user_data;
} mito_engine_analyze_options_v1;
void mito_engine_analyze_options_v1_init(
    mito_engine_analyze_options_v1 *options) MITO_NOEXCEPT;
const char *mito_engine_analyze_with_options_v1(
    void *engine, const char *input_path, const char *ref_path,
    const mito_engine_analyze_options_v1 *options) MITO_NOEXCEPT;
const char *mito_engine_get_last_error(void) MITO_NOEXCEPT;
const char *mito_engine_get_last_error_code(void) MITO_NOEXCEPT;
void mito_engine_free_string(const char *value) MITO_NOEXCEPT;

#ifdef __cplusplus
}
#endif

#undef MITO_NOEXCEPT
