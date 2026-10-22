#pragma once
#include "common.hpp"
#include "evidence.hpp"
#include "reads.hpp"
#include "resources.hpp"
#include "variants.hpp"

namespace mito::detail {

struct ProtocolTagValues {
  std::set<std::string> values;
  bool missing = false;
};

struct CallablePhyloRanges {
  bool known = false;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
};

struct AlignmentSegment {
  std::string reference_name;
  std::size_t reference_start = 0;
  std::size_t reference_end = 0;
  std::size_t query_start = 0;
  std::size_t query_end = 0;
  char strand = '+';
  std::uint8_t mapping_quality = 0;
};

struct CanonicalJunctionPath {
  std::vector<std::string> ids;
  std::vector<std::string> orientations;
};

using FeatureTokens = std::vector<std::string>;

/** Interned feature tokens, ranked by ascending profile frequency. */
using TokenRanks = std::vector<std::uint32_t>;

struct TokenProfile {
  /** Sorted global token ranks; rarer tokens come first. */
  TokenRanks tokens;
  /** Leading tokens indexed for exact Jaccard prefix filtering. */
  std::size_t prefix_length = 0;
  std::vector<std::size_t> feature_indices;
};

struct Neighborhood {
  std::vector<std::size_t> profiles;
  std::size_t read_count = 0;
};

struct TokenPosting {
  std::size_t profile = 0;
  /** Position of the token within the profile's sorted ranks. */
  std::uint32_t position = 0;
};

/** Token rank -> postings, in ascending profile order, for every profile
 * whose indexed prefix contains the rank. */
using InvertedTokenIndex = std::vector<std::vector<TokenPosting>>;

/** Reusable per-query scratch for region_query, sized to the profile count. */
struct RegionQueryScratch {
  std::vector<std::size_t> marks;
  std::vector<std::size_t> overlaps;
  std::vector<unsigned char> pruned;
  std::size_t epoch = 0;
};

struct MoleculeArchitectureEvidence {
  bool eligible = false;
  double callable_fraction = 0.0;
  std::vector<std::pair<std::size_t, ObservationState>> observations;
  std::vector<std::size_t> alternate_signature;
};

struct ArchitectureScore {
  std::size_t architecture_index = 0U;
  double distance = 1.0;
  double score = 1.0;
  double overlap_fraction = 0.0;
};

} // namespace mito::detail
