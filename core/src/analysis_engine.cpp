#include "detail/pipeline.hpp"

namespace mito {

std::string_view
analysis_error_code_name(const AnalysisErrorCode code) noexcept {
  switch (code) {
  case AnalysisErrorCode::invalid_configuration:
    return "MITO-E1001";
  case AnalysisErrorCode::input_open_failed:
    return "MITO-E1101";
  case AnalysisErrorCode::input_format_unsupported:
    return "MITO-E1102";
  case AnalysisErrorCode::input_parse_failed:
    return "MITO-E1103";
  case AnalysisErrorCode::input_empty:
    return "MITO-E1104";
  case AnalysisErrorCode::reference_open_failed:
    return "MITO-E1201";
  case AnalysisErrorCode::reference_invalid:
    return "MITO-E1202";
  case AnalysisErrorCode::resource_open_failed:
    return "MITO-E1301";
  case AnalysisErrorCode::resource_invalid:
    return "MITO-E1302";
  case AnalysisErrorCode::dependency_unavailable:
    return "MITO-E1401";
  case AnalysisErrorCode::analysis_cancelled:
    return "MITO-E1501";
  case AnalysisErrorCode::resource_exhausted:
    return "MITO-E1601";
  case AnalysisErrorCode::internal_error:
    return "MITO-E9001";
  }
  return "MITO-E9001";
}

std::string AnalysisEngine::analyze(const std::string &input_path,
                                    const std::string &reference_path,
                                    const AnalysisConfig &config) const {
  return detail::analyze_impl(input_path, reference_path, config, nullptr);
}

ProfiledAnalysis
AnalysisEngine::analyze_profiled(const std::string &input_path,
                                 const std::string &reference_path,
                                 const AnalysisConfig &config) const {
  ProfiledAnalysis result;
  result.json =
      detail::analyze_impl(input_path, reference_path, config, &result.timings);
  return result;
}

} // namespace mito
