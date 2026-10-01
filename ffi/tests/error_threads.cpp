#include "mito_c_api.h"
#include <atomic>
#include <barrier>
#include <cstring>
#include <thread>

int main() {
  std::barrier ready(2);
  std::atomic<bool> good{true};
  auto worker = [&](bool null_engine) {
    void *engine = null_engine ? nullptr : mito_engine_new();
    const char *expected =
        null_engine ? "mito_engine_analyze received a null engine"
                    : "mito_engine_analyze received a null input path";
    for (int i = 0; i < 100; ++i) {
      if (mito_engine_analyze(engine, nullptr, nullptr) != nullptr)
        good = false;
      ready.arrive_and_wait();
      if (std::strcmp(mito_engine_get_last_error(), expected) != 0)
        good = false;
      ready.arrive_and_wait();
    }
    mito_engine_delete(engine);
  };
  std::thread first(worker, true), second(worker, false);
  first.join();
  second.join();
  return good ? 0 : 1;
}
