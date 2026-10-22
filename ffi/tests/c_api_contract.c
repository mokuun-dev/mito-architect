#include "mito_c_api.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define REQUIRE(x)                                                             \
  do {                                                                         \
    if (!(x)) {                                                                \
      fprintf(stderr, "line %d: %s\n", __LINE__, #x);                          \
      return 1;                                                                \
    }                                                                          \
  } while (0)
static bool cancelled(void *data) {
  unsigned *calls = data;
  ++*calls;
  return true;
}
/* The Rust binding (ffi/src/tests.rs) asserts the same LP64 offsets, so a
 * field reordering on either side fails a test even when sizeof is unchanged. */
static int check_options_layout(void) {
#if SIZE_MAX == UINT64_MAX && UINTPTR_MAX == UINT64_MAX
  REQUIRE(offsetof(mito_engine_analyze_options_v1, struct_size) == 0);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, abi_version) == 4);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, filter_numt) == 8);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, threads) == 16);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, min_mapping_quality) == 24);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, min_base_quality) == 25);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, excluded_snp_flags) == 26);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, numt_threshold) == 32);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, allow_development_tags) == 40);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, emit_evidence_graph) == 41);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, max_evidence_observations) == 48);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, max_phase_links) == 56);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, max_phase_work) == 64);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, max_phase_molecule_references) == 72);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, max_result_bytes) == 80);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, evidence_page_size) == 88);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, molecule_id_tag) == 96);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, umi_tag) == 104);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, duplex_tag) == 112);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, min_architecture_molecules) == 120);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, max_candidate_architectures) == 128);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, architecture_max_distance) == 136);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, architecture_ambiguity_margin) == 144);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, architecture_min_overlap_fraction) == 152);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, architecture_consensus_fraction) == 160);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, architecture_optional_fraction) == 168);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, architecture_stability_replicates) == 176);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, architecture_seed) == 184);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, should_cancel) == 192);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, cancel_user_data) == 200);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, max_input_records) == 208);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, max_input_bases) == 216);
  REQUIRE(offsetof(mito_engine_analyze_options_v1, allow_external_cram_references) == 224);
  REQUIRE(MITO_ENGINE_ANALYZE_OPTIONS_V1_BASE_SIZE == 208);
  REQUIRE(sizeof(mito_engine_analyze_options_v1) == 232);
#endif
  return 0;
}

int main(int argc, char **argv) {
  REQUIRE(argc == 2);
  REQUIRE(check_options_layout() == 0);
  mito_engine_delete(NULL);
  mito_engine_free_string(NULL);
  REQUIRE(mito_engine_analyze(NULL, argv[1], NULL) == NULL);
  REQUIRE(strcmp(mito_engine_get_last_error_code(), "MITO-E1001") == 0);
  void *engine = mito_engine_new();
  REQUIRE(engine != NULL);
  REQUIRE(strlen(mito_engine_get_last_error()) == 0);
  REQUIRE(mito_engine_analyze(engine, NULL, NULL) == NULL);
  REQUIRE(strcmp(mito_engine_get_last_error_code(), "MITO-E1001") == 0);
  unsigned calls = 0;
  REQUIRE(mito_engine_analyze_with_cancel(engine, argv[1], NULL, true, 1,
                                          cancelled, &calls) == NULL);
  REQUIRE(calls > 0);
  REQUIRE(strcmp(mito_engine_get_last_error_code(), "MITO-E1501") == 0);
  const char *first = mito_engine_analyze(engine, argv[1], NULL);
  REQUIRE(first != NULL);
  REQUIRE(strlen(mito_engine_get_last_error()) == 0);
  const char *second =
      mito_engine_analyze_with_options(engine, argv[1], NULL, true, 0);
  REQUIRE(second != NULL && second != first);
  REQUIRE(strcmp(first, second) == 0);
  mito_engine_analyze_options_v1 config;
  mito_engine_analyze_options_v1_init(&config);
  const char *latest = mito_engine_analyze_with_options_v1(
      engine, argv[1], NULL, &config);
  REQUIRE(latest != NULL && strcmp(first, latest) == 0);
  config.struct_size -= 1;
  REQUIRE(mito_engine_analyze_with_options_v1(engine, argv[1], NULL, &config) == NULL);
  REQUIRE(strcmp(mito_engine_get_last_error_code(), "MITO-E1001") == 0);
  mito_engine_analyze_options_v1_init(&config);
  const char invalid_utf8[] = { (char)0xff, 0 };
  REQUIRE(mito_engine_analyze_with_options_v1(engine, invalid_utf8, NULL, &config) == NULL);
  REQUIRE(strcmp(mito_engine_get_last_error_code(), "MITO-E1001") == 0);
  /* A caller built against the 1.0 header passes the base size; appended
   * fields are ignored and take their defaults. */
  mito_engine_analyze_options_v1_init(&config);
  config.struct_size = (uint32_t)MITO_ENGINE_ANALYZE_OPTIONS_V1_BASE_SIZE;
  config.max_input_records = 1;
  const char *base_layout =
      mito_engine_analyze_with_options_v1(engine, argv[1], NULL, &config);
  REQUIRE(base_layout != NULL && strcmp(first, base_layout) == 0);
  mito_engine_free_string(base_layout);
  mito_engine_analyze_options_v1_init(&config);
  config.max_input_records = 1;
  REQUIRE(mito_engine_analyze_with_options_v1(engine, argv[1], NULL, &config) == NULL);
  REQUIRE(strcmp(mito_engine_get_last_error_code(), "MITO-E1601") == 0);
  mito_engine_analyze_options_v1_init(&config);
  config.max_result_bytes = 1;
  REQUIRE(mito_engine_analyze_with_options_v1(engine, argv[1], NULL, &config) == NULL);
  REQUIRE(strcmp(mito_engine_get_last_error_code(), "MITO-E1601") == 0);
  mito_engine_delete(engine);
  REQUIRE(strstr(first, "\"schema_version\":\"0.5\"") != NULL);
  mito_engine_free_string(second);
  mito_engine_free_string(first);
  mito_engine_free_string(latest);
  REQUIRE(strlen(mito_engine_version()) > 0);
  return 0;
}
