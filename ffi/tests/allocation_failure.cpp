#include "mito_c_api.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>

namespace {
thread_local bool fail_allocations = false;
}
void *operator new(std::size_t size) {
  if (fail_allocations)
    throw std::bad_alloc();
  if (void *value = std::malloc(size == 0 ? 1 : size))
    return value;
  throw std::bad_alloc();
}
void operator delete(void *value) noexcept { std::free(value); }
void operator delete(void *value, std::size_t) noexcept { std::free(value); }

int main() {
  bool escaped = false;
  const char *result = nullptr;
  fail_allocations = true;
  try {
    result = mito_engine_analyze(nullptr, nullptr, nullptr);
  } catch (...) {
    escaped = true;
  }
  fail_allocations = false;
  if (escaped || result != nullptr ||
      std::strcmp(mito_engine_get_last_error_code(), "MITO-E1001") != 0) {
    std::fputs("null-argument diagnostic allocated or escaped the C ABI\n",
               stderr);
    return 1;
  }
  fail_allocations = true;
  void *engine = mito_engine_new();
  fail_allocations = false;
  if (engine != nullptr ||
      std::strcmp(mito_engine_get_last_error_code(), "MITO-E1601") != 0 ||
      std::strlen(mito_engine_get_last_error()) == 0)
    return 2;
  // The thread's error state and allocator recover for the next call.
  engine = mito_engine_new();
  if (engine == nullptr || std::strlen(mito_engine_get_last_error()) != 0)
    return 3;
  mito_engine_delete(engine);
  return 0;
}
