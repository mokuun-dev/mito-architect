add_library(mito_ffi STATIC
  src/mito_c_api.cpp
  src/native/error_state.cpp
  src/native/legacy_api.cpp
)

target_include_directories(mito_ffi
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(mito_ffi PUBLIC mito_core)
target_compile_features(mito_ffi PUBLIC cxx_std_20)
mito_enable_warnings(mito_ffi)

install(TARGETS mito_ffi ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)


if(MITO_BUILD_TESTS AND BUILD_TESTING)
  add_executable(mito_c_api_contract tests/c_api_contract.c)
  target_link_libraries(mito_c_api_contract PRIVATE mito_ffi)
  set_target_properties(mito_c_api_contract PROPERTIES LINKER_LANGUAGE CXX)
  add_test(NAME mito_c_api_contract COMMAND mito_c_api_contract
    "${PROJECT_SOURCE_DIR}/fixtures/tiny.fastq")
  foreach(test_name allocation_failure error_threads)
    add_executable(mito_ffi_${test_name} tests/${test_name}.cpp)
    target_link_libraries(mito_ffi_${test_name} PRIVATE mito_ffi)
    target_compile_features(mito_ffi_${test_name} PRIVATE cxx_std_20)
    add_test(NAME mito_ffi_${test_name} COMMAND mito_ffi_${test_name})
  endforeach()
endif()
