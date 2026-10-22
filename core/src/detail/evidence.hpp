#pragma once
#include "common.hpp"
#include "reads.hpp"
#include "variants.hpp"

namespace mito::detail {

enum class ObservationState : std::uint8_t {
  reference,
  alternate,
  event_absent,
  not_callable,
  low_quality,
  conflict,
};

struct EvidenceEvent {
  std::string id;
  std::string type;
  std::size_t start = 0;
  std::size_t end = 0;
  std::size_t length = 0;
  std::string reference;
  std::string alternate;
  std::string normalization;
  std::string source_projection;
  std::string negative_evidence_rule;
  bool absence_assessable = false;
  std::vector<std::string> component_event_ids;
  std::vector<std::string> supporting_molecules;
};

struct EvidenceObservation {
  std::size_t id = 0;
  MoleculeIndex molecule_index;
  std::size_t event_index = 0;
  AlignmentFragmentId alignment_fragment_id;
  ObservationState state = ObservationState::not_callable;
  std::string observed_allele;
  std::optional<std::uint8_t> base_quality;
  std::uint8_t mapping_quality = 0;
  char strand = '+';
  std::optional<double> center_proximity;
  std::string evidence_source;
};

struct AlignmentCallabilityEvidence {
  AlignmentFragmentId alignment_fragment_id;
  bool eligible = false;
  std::string status;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  std::size_t callable_bases = 0;
  std::map<std::string, std::size_t> reference_exclusion_counts;
  std::vector<std::size_t> disrupted_adjacency_anchors;
  std::size_t inserted_query_bases = 0;
  std::size_t soft_clipped_query_bases = 0;
};

struct MoleculeCallabilityEvidence {
  MoleculeIndex molecule_index;
  bool known = false;
  std::string status;
  std::vector<std::pair<std::size_t, std::size_t>> ranges;
  std::size_t callable_bases = 0;
  std::vector<AlignmentCallabilityEvidence> alignments;
};

struct CallabilityResult {
  std::vector<MoleculeCallabilityEvidence> molecules;
};

struct PhaseLink {
  std::size_t event_a_index = 0;
  std::size_t event_b_index = 0;
  std::size_t jointly_callable = 0;
  std::size_t both_alternate = 0;
  std::size_t a_alternate_b_absent = 0;
  std::size_t a_absent_b_alternate = 0;
  std::size_t neither_alternate = 0;
  std::size_t jointly_uncertain = 0;
  double co_alternate_fraction = 0.0;
  double ci95_low = 0.0;
  double ci95_high = 0.0;
  double expected_co_alternate_fraction = 0.0;
  double linkage_delta = 0.0;
  bool complete_callability = false;
  std::vector<std::size_t> supporting_molecule_indices;
  std::vector<std::size_t> uncertain_molecule_indices;
};

enum class ArchitectureAssignmentStatus : std::uint8_t {
  assigned,
  ambiguous,
  unassigned,
  ineligible,
};

struct ArchitectureEventProfile {
  std::size_t event_index = 0;
  std::size_t alternate = 0;
  std::size_t absent = 0;
  std::size_t uncertain = 0;
  std::size_t not_callable = 0;
  double alternate_fraction = 0.0;
  bool expected_alternate = false;
  bool expected_absent = false;
};

struct ArchitectureAssignment {
  ArchitectureAssignmentStatus status =
      ArchitectureAssignmentStatus::ineligible;
  std::optional<std::size_t> architecture_index;
  double distance = 1.0;
  double score = 1.0;
  double overlap_fraction = 0.0;
  double confidence = 0.0;
  std::vector<std::size_t> candidate_indices;
};

struct CandidateArchitecture {
  std::string id;
  std::vector<std::size_t> seed_signature;
  std::vector<std::size_t> defining_event_indices;
  std::vector<std::size_t> optional_event_indices;
  std::vector<ArchitectureEventProfile> event_profile;
  std::vector<std::size_t> member_molecule_indices;
  std::vector<std::size_t> ambiguous_molecule_indices;
  std::vector<std::size_t> representative_molecule_indices;
  double estimated_fraction = 0.0;
  double ci95_low = 0.0;
  double ci95_high = 0.0;
  double median_callable_fraction = 0.0;
  double mean_assignment_confidence = 0.0;
  double min_assignment_confidence = 0.0;
  double signature_stability = 0.0;
  double recovery_rate = 0.0;
  double abundance_standard_deviation = 0.0;
};

struct ArchitectureInferenceResult {
  std::vector<CandidateArchitecture> architectures;
  std::vector<ArchitectureAssignment> assignments;
  std::size_t eligible_molecules = 0;
  std::size_t assigned_molecules = 0;
  std::size_t ambiguous_molecules = 0;
  std::size_t unassigned_molecules = 0;
  std::string status = "NOT_RUN";
  std::vector<std::string> qc_flags;
};

struct NormalizedSmallIndel {
  std::string id;
  std::string type;
  std::size_t start = 0;
  std::size_t end = 0;
  std::size_t length = 0;
  std::string reference;
  std::string alternate;
  std::string normalization;
  std::optional<std::uint8_t> base_quality;
};

struct SparseEvidenceStore {
  std::vector<EvidenceEvent> events;
  std::vector<EvidenceObservation> observations;
};

struct EvidenceAggregationResult {
  std::map<std::string, SnpAggregate> variants;
  SparseEvidenceStore store;
};

struct VariantAlleleEvidenceSummary {
  std::size_t count = 0U;
  std::array<std::size_t, 2> strand_depths{};
  double base_quality_sum = 0.0;
  std::size_t base_quality_count = 0U;
  std::uint8_t base_quality_min = std::numeric_limits<std::uint8_t>::max();
  std::uint8_t base_quality_max = 0U;
  double mapping_quality_sum = 0.0;
  std::uint8_t mapping_quality_min = std::numeric_limits<std::uint8_t>::max();
  std::uint8_t mapping_quality_max = 0U;
  double read_position_sum = 0.0;
  std::size_t read_position_count = 0U;
};

struct UnifiedVariantEvidenceSummary {
  VariantAlleleEvidenceSummary alternate;
  VariantAlleleEvidenceSummary reference;
  VariantAlleleEvidenceSummary other;
  std::size_t event_absent = 0U;
  std::size_t low_quality = 0U;
  std::size_t conflict = 0U;
  bool multi_allelic = false;
};

} // namespace mito::detail
