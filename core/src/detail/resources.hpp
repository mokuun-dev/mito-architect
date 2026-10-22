#pragma once
#include "common.hpp"
#include "evidence.hpp"
#include "reads.hpp"
#include "variants.hpp"

namespace mito::detail {

struct ResourceRecord {
  std::string name;
  std::string version;
  std::string path;
  std::string sha256;
  std::string source;
  std::string license;
  std::string retrieved;
};

struct HaplogroupDefinition {
  std::string name;
  std::vector<PhyloMutation> mutations;
  std::size_t source_order = 0U;
};

struct HaplogroupCandidate {
  std::string name;
  double score = 0.0;
  std::vector<std::string> matched;
  std::vector<std::string> missing;
  std::vector<std::string> extra;
};

struct ClusterHaplogroupAssignment {
  std::string best = "unassigned";
  double quality = 0.0;
  bool callable_ranges_known = false;
  bool contamination_warning = false;
  std::vector<std::string> observed_markers;
  std::vector<std::pair<std::size_t, std::size_t>> callable_ranges;
  std::vector<HaplogroupCandidate> candidates;
};

} // namespace mito::detail
