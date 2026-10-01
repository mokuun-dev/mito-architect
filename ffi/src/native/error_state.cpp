#include "error_state.hpp"
#include <algorithm>
#include <array>
#include <cstring>

namespace mito::ffi {
namespace {
thread_local std::array<char, 2048> error_message{};
thread_local std::array<char, 32> error_code{};
} // namespace
void clear_error() noexcept {
  error_message[0] = '\0';
  error_code[0] = '\0';
}
void set_error(std::string_view code, const char *message) noexcept {
  const auto size = std::min(code.size(), error_code.size() - 1U);
  std::memcpy(error_code.data(), code.data(), size);
  error_code[size] = '\0';
  if (message == nullptr)
    message = "unknown error";
  std::size_t i = 0;
  while (i + 1U < error_message.size() && message[i] != '\0') {
    error_message[i] = message[i];
    ++i;
  }
  error_message[i] = '\0';
}
const char *last_error() noexcept { return error_message.data(); }
const char *last_error_code() noexcept {
  return error_code[0] == '\0' ? "MITO-E9001" : error_code.data();
}
} // namespace mito::ffi
