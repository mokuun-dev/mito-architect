#pragma once
#include <string_view>

namespace mito::ffi {
// No dynamic allocation: usable even while handling std::bad_alloc.
void clear_error() noexcept;
void set_error(std::string_view code, const char *message) noexcept;
const char *last_error() noexcept;
const char *last_error_code() noexcept;
} // namespace mito::ffi
