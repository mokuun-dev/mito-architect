#pragma once
#include "common.hpp"

namespace mito::detail {

struct GeneAnnotation {
  std::string name;
  std::size_t start;
  std::size_t end;
  std::string strand;
  std::string biotype;
};

struct CigarOperation {
  std::size_t length = 0;
  char code = 'M';
};

struct SupplementaryAlignment {
  std::string reference_name;
  std::size_t reference_start = 0;
  char strand = '+';
  std::string cigar;
  std::vector<CigarOperation> cigar_operations;
  std::uint8_t mapping_quality = 0;
  std::size_t edit_distance = 0;
};

struct ReadRecord {
  std::string id;
  std::string sequence;
  std::string qualities;
  std::size_t reference_start = 0;
  std::string reference_name;
  std::string cigar;
  std::vector<CigarOperation> cigar_operations;
  std::uint16_t flags = 0;
  std::uint8_t mapping_quality = 0;
  std::map<std::string, std::string> aux_tags;
  std::vector<SupplementaryAlignment> supplementary_alignments;
};

struct InputRecords {
  std::vector<ReadRecord> reads;
  bool alignment_input = false;
  bool has_nuclear_contigs = false;
};

/**
 * A nuclear contig in an alignment header is useful context, but does not by
 * itself prove that the input was aligned against a versioned competitive
 * reference.  This record keeps those two facts separate in the result
 * contract.
 */
struct NumtAssessment {
  std::string mode = "unaligned_fastq";
  bool nuclear_contigs_present = false;
  bool specificity_assessable = false;
  std::string provenance_status = "NOT_APPLICABLE";
  std::string provenance_reason;
};

struct AlignmentFragmentId {
  std::size_t value = 0;

  auto operator<=>(const AlignmentFragmentId &) const = default;
};

struct MoleculeIndex {
  std::size_t value = 0;

  auto operator<=>(const MoleculeIndex &) const = default;
};

struct AlignmentFragmentEvidence {
  AlignmentFragmentId id;
  MoleculeIndex molecule_index;
  std::size_t source_record_index = 0;
  std::string molecule_id;
  std::string role;
  bool selected_representative = false;
  std::uint16_t flags = 0;
  std::uint8_t mapping_quality = 0;
  std::string reference_name;
  std::size_t reference_start = 0;
  std::string cigar;
};

struct MoleculeAssemblyEvidence {
  MoleculeIndex index;
  std::string id;
  std::string identity_policy;
  std::string assembly_status;
  std::vector<AlignmentFragmentId> fragment_ids;
  AlignmentFragmentId representative_fragment_id;
  std::size_t primary_candidate_count = 0;
  bool ambiguous = false;
  bool analysis_eligible = true;
  std::vector<std::string> source_qnames;
  std::map<std::string, std::string> protocol_metadata;
  std::vector<std::string> protocol_flags;
  std::vector<std::string> exclusion_reasons;
  std::vector<std::string> warnings;
};

struct MoleculeAssemblyResult {
  std::vector<ReadRecord> representatives;
  /** Original alignment records retained for provenance and future
   * multi-fragment evidence. */
  std::vector<ReadRecord> source_alignment_records;
  std::vector<AlignmentFragmentEvidence> fragments;
  std::vector<MoleculeAssemblyEvidence> molecules;
};

} // namespace mito::detail
