#pragma once
#include "common.hpp"
#include "reads.hpp"

namespace mito::detail {

struct ClinicalAssertion {
  std::string source;
  std::string assertion_id;
  std::string allele_id;
  std::string disease;
  std::string clinical_significance;
  std::string normalized_significance;
  std::string review_status;
  std::string assertion_date;
  std::string source_url;
  std::vector<std::string> references;
  std::string resource_version;
  std::string retrieved_at;

  auto operator<=>(const ClinicalAssertion &) const = default;
};

struct SnpCall {
  std::size_t position = 0;
  char reference = 'N';
  char alternate = 'N';
  std::string gene;
  std::string consequence;
  std::string protein;
  std::string residue;
  std::string phenotype;
  std::string pathogenicity;
  std::vector<std::string> references;
  std::vector<std::string> sources;
  std::string structure_id;
  std::string structure_chain;
  std::size_t structure_residue = 0;
  std::string structure_complex;
  std::string clinvar_allele_id;
  std::string mitomap_url;
  std::string clinical_conflict_status;
  std::string clinical_consensus_significance;
  std::vector<ClinicalAssertion> clinical_assertions;
};

enum class PhyloMutationKind : std::uint8_t {
  substitution,
  deletion,
  insertion,
};

struct PhyloMutationLocus {
  std::size_t position = 0;
  PhyloMutationKind kind = PhyloMutationKind::substitution;
  std::string insertion_index;

  auto operator<=>(const PhyloMutationLocus &) const = default;
};

struct PhyloMutation {
  PhyloMutationLocus locus;
  PhyloMutationKind kind = PhyloMutationKind::substitution;
  char alternate = 'N';
  std::string inserted_bases;
  std::string encoded;
  double weight = 1.0;
  bool backmutation = false;
};

struct ObservedPhyloMutation {
  std::size_t position = 0;
  PhyloMutationKind kind = PhyloMutationKind::substitution;
  char alternate = 'N';
  std::string insertion_index;
  std::string inserted_bases;
  std::string encoded;
};

struct PhyloAlignmentRule {
  std::vector<std::string> erroneous_markers;
  std::vector<ObservedPhyloMutation> replacement_markers;
};

struct SvCall {
  std::string id;
  std::string type;
  std::size_t start = 0;
  std::size_t end = 0;
  std::size_t length = 0;
  std::vector<std::string> supporting_reads;
  bool known_event = false;
  std::vector<std::string> evidence_sources;
  std::vector<std::string> orientations;
  std::size_t segment_count = 1;
};

struct ComplexSvCall {
  std::string id;
  std::vector<std::string> junction_ids;
  std::vector<std::string> junction_orientations;
  std::vector<std::string> supporting_reads;
  std::size_t segment_count = 0;
};

struct ReadFeature {
  std::string id;
  std::size_t length = 0;
  double mean_quality = 0.0;
  double numt_score = 0.0;
  bool filtered_numt = false;
  std::vector<std::string> numt_evidence;
  int cluster_id = -1;
  bool outlier = false;
  std::uint8_t mapping_quality = 0;
  std::uint16_t flags = 0;
  std::string reference_name;
  std::map<std::string, std::string> aux_tags;
  std::vector<SnpCall> snps;
  std::vector<ObservedPhyloMutation> haplogroup_markers;
  std::vector<std::pair<std::size_t, std::size_t>> haplogroup_ranges;
  bool haplogroup_range_known = false;
  std::vector<std::string> sv_ids;
  std::vector<std::string> complex_event_ids;
};

struct CoverageBin {
  std::size_t start = 0;
  std::size_t end = 0;
  std::size_t depth = 0;
};

struct CoverageResult {
  std::vector<CoverageBin> bins;
  double mean_depth = 0.0;
  double pct_sites_gt20x = 0.0;
  std::size_t max_depth = 0;
};

struct FeatureExtractionResult {
  ReadFeature feature;
  std::vector<SvCall> svs;
  std::vector<ComplexSvCall> complex_events;
};

struct SnpAggregate {
  SnpCall call;
  std::size_t alternate_depth = 0;
  std::size_t reference_depth = 0;
  std::size_t other_depth = 0;
  std::size_t callable_depth = 0;
  double heteroplasmy = 0.0;
  double ci95_low = 0.0;
  double ci95_high = 0.0;
  std::vector<std::string> supporting_reads;
  std::array<std::size_t, 2> alternate_strand_depths{};
  std::array<std::size_t, 2> reference_strand_depths{};
  std::array<std::size_t, 2> other_strand_depths{};
  double alternate_quality_sum = 0.0;
  double reference_quality_sum = 0.0;
  double other_quality_sum = 0.0;
  std::uint8_t alternate_quality_min = std::numeric_limits<std::uint8_t>::max();
  std::uint8_t reference_quality_min = std::numeric_limits<std::uint8_t>::max();
  std::uint8_t other_quality_min = std::numeric_limits<std::uint8_t>::max();
  std::uint8_t alternate_quality_max = 0;
  std::uint8_t reference_quality_max = 0;
  std::uint8_t other_quality_max = 0;
  double alternate_read_position_sum = 0.0;
  double reference_read_position_sum = 0.0;
  double other_read_position_sum = 0.0;
};

struct LocusAlleleAccumulator {
  std::array<std::size_t, 4> depths{};
  std::array<std::array<std::size_t, 2>, 4> strand_depths{};
  std::array<double, 4> quality_sums{};
  std::array<std::uint8_t, 4> quality_mins{};
  std::array<std::uint8_t, 4> quality_maxes{};
  std::array<double, 4> read_position_sums{};

  LocusAlleleAccumulator() {
    quality_mins.fill(std::numeric_limits<std::uint8_t>::max());
  }
};

struct MoleculeAlleleObservation {
  char base = 'N';
  std::uint8_t quality = 0;
  double center_proximity = 0.0;
  std::size_t query_index = 0;
  AlignmentFragmentId alignment_fragment_id;
  std::uint8_t mapping_quality = 0;
  char strand = '+';
  bool has_passing_observation = false;
  bool covered = false;
  bool conflicted = false;
};

} // namespace mito::detail
