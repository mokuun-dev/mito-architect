#include "mito_c_api.h"
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
int main(int argc, char **argv) {
  REQUIRE(argc == 2);
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
