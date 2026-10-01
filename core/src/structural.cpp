#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] std::optional<SvCall> parse_sv_tag(std::string_view read_id,
                                                 std::string_view token) {
  const auto colon = token.find(':');
  const auto dash =
      token.find('-', colon == std::string_view::npos ? 0 : colon + 1);
  if (colon == std::string_view::npos || dash == std::string_view::npos) {
    return std::nullopt;
  }

  const auto type = token.substr(0, colon);
  const auto start = parse_size(token.substr(colon + 1, dash - colon - 1));
  const auto end = parse_size(token.substr(dash + 1));
  if (!start || !end || *start == 0 || *end < *start) {
    return std::nullopt;
  }

  SvCall sv;
  sv.type = std::string(type);
  std::transform(
      sv.type.begin(), sv.type.end(), sv.type.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (sv.type == "del") {
    sv.type = "deletion";
  } else if (sv.type == "ins") {
    sv.type = "insertion";
  }
  sv.start = *start;
  sv.end = *end;
  sv.length = sv.end - sv.start + 1;
  sv.supporting_reads.emplace_back(read_id);
  sv.evidence_sources.emplace_back("development_tag");
  sv.known_event = sv.type == "deletion" && sv.start >= 8400 &&
                   sv.start <= 8500 && sv.end >= 13400 && sv.end <= 13550;
  sv.id =
      sv.type + ":" + std::to_string(sv.start) + "-" + std::to_string(sv.end);
  return sv;
}

[[nodiscard]] std::vector<SvCall> sv_tags_from_header(const ReadRecord &read) {
  std::vector<SvCall> calls;
  std::istringstream iss(read.id);
  std::string token;
  while (iss >> token) {
    constexpr std::string_view prefix = "sv=";
    if (token.rfind(prefix, 0) == 0) {
      auto call =
          parse_sv_tag(read.id, std::string_view(token).substr(prefix.size()));
      if (call) {
        calls.push_back(std::move(*call));
      }
    }
  }
  return calls;
}

[[nodiscard]] std::vector<SvCall>
parse_cigar_svs(const ReadRecord &read, std::size_t sv_min_length,
                std::size_t reference_length) {
  std::vector<SvCall> calls;
  if (read.cigar_operations.empty() || (read.flags & 0x4U) != 0U ||
      read.reference_start == 0 || read.reference_name.empty() ||
      read.reference_name == "*" ||
      looks_like_nuclear_contig(read.reference_name)) {
    return calls;
  }

  std::size_t reference_cursor =
      read.reference_start == 0 ? 1 : read.reference_start;
  for (const auto &operation : read.cigar_operations) {
    const auto len = operation.length;
    const char op = operation.code;
    if ((op == 'D' || op == 'N') && len >= sv_min_length) {
      if (len - 1U >
          std::numeric_limits<std::size_t>::max() - reference_cursor) {
        throw AnalysisError(
            AnalysisErrorCode::input_parse_failed,
            "CIGAR structural-variant coordinate overflow for read '" +
                read.id + "'");
      }
      SvCall sv;
      sv.type = "deletion";
      sv.start = reference_cursor;
      sv.end = reference_cursor + len - 1;
      sv.length = len;
      sv.supporting_reads.emplace_back(read.id);
      sv.evidence_sources.emplace_back("cigar");
      sv.known_event = sv.start >= 8400 && sv.start <= 8500 &&
                       sv.end >= 13400 && sv.end <= 13550;
      sv.id =
          "deletion:" + std::to_string(sv.start) + "-" + std::to_string(sv.end);
      calls.push_back(std::move(sv));
    } else if (op == 'I' && len >= sv_min_length) {
      SvCall sv;
      sv.type = "insertion";
      sv.start =
          reference_cursor > 1U ? reference_cursor - 1U : reference_length;
      sv.end = sv.start;
      sv.length = len;
      sv.supporting_reads.emplace_back(read.id);
      sv.evidence_sources.emplace_back("cigar");
      sv.id =
          "insertion:" + std::to_string(sv.start) + "+" + std::to_string(len);
      calls.push_back(std::move(sv));
    } else if (op == 'S' && len >= sv_min_length) {
      SvCall sv;
      const bool leading =
          reference_cursor ==
          (read.reference_start == 0 ? 1 : read.reference_start);
      sv.type = leading ? "soft_clip_left" : "soft_clip_right";
      sv.start = leading ? reference_cursor
                         : (reference_cursor > 1U ? reference_cursor - 1U
                                                  : reference_length);
      sv.end = reference_cursor;
      sv.length = len;
      sv.supporting_reads.emplace_back(read.id);
      sv.evidence_sources.emplace_back("cigar");
      sv.id =
          sv.type + ":" + std::to_string(sv.start) + "+" + std::to_string(len);
      calls.push_back(std::move(sv));
    }

    if (op == 'M' || op == '=' || op == 'X' || op == 'D' || op == 'N') {
      if (len > std::numeric_limits<std::size_t>::max() - reference_cursor) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "CIGAR reference coordinate overflow for read '" +
                                read.id + "'");
      }
      reference_cursor += len;
    }
  }

  return calls;
}

[[nodiscard]] bool is_mitochondrial_contig(std::string_view name) {
  std::string normalized(name);
  std::transform(
      normalized.begin(), normalized.end(), normalized.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return normalized == "mt" || normalized == "chrm" || normalized == "m" ||
         normalized == "nc_012920.1" || normalized == "nc_012920";
}

[[nodiscard]] std::optional<AlignmentSegment>
make_alignment_segment(std::string reference_name, std::size_t reference_start,
                       char strand, std::uint8_t mapping_quality,
                       const std::vector<CigarOperation> &operations) {
  if (!is_mitochondrial_contig(reference_name) || reference_start == 0 ||
      operations.empty()) {
    return std::nullopt;
  }
  std::size_t leading_clip = 0;
  std::size_t query_extent = 0;
  std::size_t query_aligned = 0;
  std::size_t reference_span = 0;
  bool seen_alignment = false;
  for (const auto &operation : operations) {
    if (operation.code == 'M' || operation.code == 'I' ||
        operation.code == 'S' || operation.code == 'H' ||
        operation.code == '=' || operation.code == 'X') {
      if (operation.length >
          std::numeric_limits<std::size_t>::max() - query_extent) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment query extent overflow");
      }
      query_extent += operation.length;
    }
    if (!seen_alignment && (operation.code == 'S' || operation.code == 'H')) {
      if (operation.length >
          std::numeric_limits<std::size_t>::max() - leading_clip) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment query coordinate overflow");
      }
      leading_clip += operation.length;
      continue;
    }
    if (operation.code == 'M' || operation.code == '=' ||
        operation.code == 'X') {
      seen_alignment = true;
      if (operation.length >
              std::numeric_limits<std::size_t>::max() - query_aligned ||
          operation.length >
              std::numeric_limits<std::size_t>::max() - reference_span) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment segment length overflow");
      }
      query_aligned += operation.length;
      reference_span += operation.length;
    } else if (operation.code == 'I') {
      seen_alignment = true;
      if (operation.length >
          std::numeric_limits<std::size_t>::max() - query_aligned) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment query length overflow");
      }
      query_aligned += operation.length;
    } else if (operation.code == 'D' || operation.code == 'N') {
      seen_alignment = true;
      if (operation.length >
          std::numeric_limits<std::size_t>::max() - reference_span) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "alignment reference length overflow");
      }
      reference_span += operation.length;
    }
  }
  if (query_aligned == 0 || reference_span == 0) {
    return std::nullopt;
  }
  if (leading_clip > query_extent ||
      query_aligned > query_extent - leading_clip) {
    throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                        "alignment clipping exceeds query extent");
  }
  const auto query_start = strand == '+'
                               ? leading_clip
                               : query_extent - leading_clip - query_aligned;
  if (reference_span - 1U >
          std::numeric_limits<std::size_t>::max() - reference_start ||
      query_aligned - 1U >
          std::numeric_limits<std::size_t>::max() - query_start) {
    throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                        "alignment segment coordinate overflow");
  }
  return AlignmentSegment{std::move(reference_name),
                          reference_start,
                          reference_start + reference_span - 1U,
                          query_start,
                          query_start + query_aligned - 1U,
                          strand,
                          mapping_quality};
}

[[nodiscard]] std::vector<AlignmentSegment>
alignment_segments(const ReadRecord &read) {
  std::vector<AlignmentSegment> segments;
  if (const auto primary =
          make_alignment_segment(read.reference_name, read.reference_start,
                                 (read.flags & 0x10U) != 0U ? '-' : '+',
                                 read.mapping_quality, read.cigar_operations)) {
    segments.push_back(*primary);
  }
  for (const auto &alignment : read.supplementary_alignments) {
    if (const auto segment = make_alignment_segment(
            alignment.reference_name, alignment.reference_start,
            alignment.strand, alignment.mapping_quality,
            alignment.cigar_operations)) {
      segments.push_back(*segment);
    }
  }
  std::sort(
      segments.begin(), segments.end(), [](const auto &lhs, const auto &rhs) {
        return std::tie(lhs.query_start, lhs.query_end, lhs.reference_start,
                        lhs.strand) < std::tie(rhs.query_start, rhs.query_end,
                                               rhs.reference_start, rhs.strand);
      });
  segments.erase(std::unique(segments.begin(), segments.end(),
                             [](const auto &lhs, const auto &rhs) {
                               return lhs.query_start == rhs.query_start &&
                                      lhs.query_end == rhs.query_end &&
                                      lhs.reference_start ==
                                          rhs.reference_start &&
                                      lhs.reference_end == rhs.reference_end &&
                                      lhs.strand == rhs.strand;
                             }),
                 segments.end());
  return segments;
}

[[nodiscard]] std::vector<SvCall>
split_alignment_svs(const ReadRecord &read, std::size_t reference_length,
                    std::size_t sv_min_length) {
  std::vector<SvCall> calls;
  const auto segments = alignment_segments(read);
  if (segments.size() < 2 || reference_length == 0) {
    return calls;
  }
  for (std::size_t i = 1; i < segments.size(); ++i) {
    const auto &previous = segments[i - 1];
    const auto &current = segments[i];
    const auto query_gap = static_cast<std::int64_t>(current.query_start) -
                           static_cast<std::int64_t>(previous.query_end) - 1;
    SvCall event;
    event.supporting_reads.push_back(read.id);
    event.evidence_sources.emplace_back("split_alignment");
    event.segment_count = segments.size();
    event.orientations.push_back(std::string(1, previous.strand) + "/" +
                                 current.strand);

    const auto leaving_position = previous.strand == '+'
                                      ? previous.reference_end
                                      : previous.reference_start;
    const auto entering_position =
        current.strand == '+' ? current.reference_start : current.reference_end;
    if (leaving_position > reference_length ||
        entering_position > reference_length) {
      throw AnalysisError(
          AnalysisErrorCode::input_parse_failed,
          "split alignment coordinate exceeds the mitochondrial reference for "
          "read '" +
              read.id + "'");
    }

    if (previous.strand != current.strand) {
      event.type = "inversion";
      event.start = std::min(leaving_position, entering_position);
      event.end = std::max(leaving_position, entering_position);
      event.length = event.end - event.start + 1U;
    } else {
      const bool forward = previous.strand == '+';
      const bool follows_linear_reference =
          forward ? entering_position > leaving_position
                  : entering_position < leaving_position;
      const auto linear_gap =
          follows_linear_reference
              ? (forward ? entering_position - leaving_position - 1U
                         : leaving_position - entering_position - 1U)
              : 0U;
      const auto overlap =
          follows_linear_reference
              ? 0U
              : (forward ? leaving_position - entering_position + 1U
                         : entering_position - leaving_position + 1U);
      const auto circular_gap =
          follows_linear_reference
              ? reference_length
              : (forward ? reference_length - leaving_position +
                               entering_position - 1U
                         : leaving_position - 1U + reference_length -
                               entering_position);
      const bool crosses_origin =
          !follows_linear_reference && circular_gap < overlap;
      const bool ambiguous_topology =
          !follows_linear_reference && circular_gap == overlap;

      if (crosses_origin) {
        event.type = "circular_origin";
        event.start = std::max(leaving_position, entering_position);
        event.end = std::min(leaving_position, entering_position);
        event.length = circular_gap;
      } else if (ambiguous_topology && overlap >= sv_min_length) {
        event.type = "ambiguous_adjacency";
        event.start = std::min(leaving_position, entering_position);
        event.end = std::max(leaving_position, entering_position);
        event.length = overlap;
      } else if (follows_linear_reference && linear_gap >= sv_min_length) {
        event.type = "deletion";
        event.start = std::min(leaving_position, entering_position) + 1U;
        event.end = std::max(leaving_position, entering_position) - 1U;
        event.length = linear_gap;
      } else if (!follows_linear_reference && overlap >= sv_min_length) {
        event.type = "duplication";
        event.start = std::min(leaving_position, entering_position);
        event.end = std::max(leaving_position, entering_position);
        event.length = overlap;
      } else if (query_gap >= static_cast<std::int64_t>(sv_min_length)) {
        event.type = "insertion";
        event.start = std::min(leaving_position, entering_position);
        event.end = event.start;
        event.length = static_cast<std::size_t>(query_gap);
      } else {
        continue;
      }
    }
    event.known_event = event.type == "deletion" && event.start >= 8400 &&
                        event.start <= 8500 && event.end >= 13400 &&
                        event.end <= 13550;
    if (event.type == "insertion") {
      event.id = event.type + ":" + std::to_string(event.start) + "+" +
                 std::to_string(event.length);
    } else {
      event.id = event.type + ":" + std::to_string(event.start) + "-" +
                 std::to_string(event.end);
    }
    calls.push_back(std::move(event));
  }
  return calls;
}

[[nodiscard]] char flip_strand(char strand) {
  if (strand == '+') {
    return '-';
  }
  if (strand == '-') {
    return '+';
  }
  throw AnalysisError(AnalysisErrorCode::internal_error,
                      "invalid strand in split-junction orientation");
}

[[nodiscard]] std::string
reverse_complement_orientation(std::string_view orientation) {
  if (orientation.size() != 3U || orientation[1] != '/') {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "invalid split-junction orientation '" +
                            std::string(orientation) + "'");
  }
  return std::string{flip_strand(orientation[2]), '/',
                     flip_strand(orientation[0])};
}

[[nodiscard]] CanonicalJunctionPath
canonical_junction_path(const std::vector<SvCall> &junctions) {
  CanonicalJunctionPath forward;
  std::vector<std::string> forward_tokens;
  forward.ids.reserve(junctions.size());
  forward.orientations.reserve(junctions.size());
  forward_tokens.reserve(junctions.size());
  for (const auto &junction : junctions) {
    if (junction.orientations.size() != 1U) {
      throw AnalysisError(
          AnalysisErrorCode::internal_error,
          "split junction does not have exactly one orientation");
    }
    forward.ids.push_back(junction.id);
    forward.orientations.push_back(junction.orientations.front());
    forward_tokens.push_back(junction.id + "@" + junction.orientations.front());
  }

  CanonicalJunctionPath reverse;
  std::vector<std::string> reverse_tokens;
  reverse.ids.reserve(junctions.size());
  reverse.orientations.reserve(junctions.size());
  reverse_tokens.reserve(junctions.size());
  for (auto it = junctions.rbegin(); it != junctions.rend(); ++it) {
    const auto orientation =
        reverse_complement_orientation(it->orientations.front());
    reverse.ids.push_back(it->id);
    reverse.orientations.push_back(orientation);
    reverse_tokens.push_back(it->id + "@" + orientation);
  }
  return reverse_tokens < forward_tokens ? std::move(reverse)
                                         : std::move(forward);
}

[[nodiscard]] std::optional<ComplexSvCall>
coalesce_complex_sv(std::string_view read_id,
                    const std::vector<SvCall> &junctions) {
  if (junctions.size() < 2U) {
    return std::nullopt;
  }

  ComplexSvCall event;
  auto canonical_path = canonical_junction_path(junctions);
  event.junction_ids = std::move(canonical_path.ids);
  event.junction_orientations = std::move(canonical_path.orientations);
  event.supporting_reads.emplace_back(read_id);
  for (const auto &junction : junctions) {
    event.segment_count = std::max(event.segment_count, junction.segment_count);
  }

  // The exact canonical path is the identifier. This is deliberately not a
  // truncated hash: an event ID remains collision-free and independently
  // reviewable, while reversal canonicalization merges opposite-strand reads.
  event.id = "complex:";
  for (std::size_t i = 0; i < event.junction_ids.size(); ++i) {
    if (i != 0U) {
      event.id.push_back('|');
    }
    event.id += event.junction_ids[i];
    event.id.push_back('@');
    event.id += event.junction_orientations[i];
  }
  return event;
}

void merge_sv(std::map<std::string, SvCall> &svs, SvCall sv) {
  const auto existing = svs.find(sv.id);
  if (existing == svs.end()) {
    for (auto *values :
         {&sv.supporting_reads, &sv.evidence_sources, &sv.orientations}) {
      std::sort(values->begin(), values->end());
      values->erase(std::unique(values->begin(), values->end()), values->end());
    }
    svs.emplace(sv.id, std::move(sv));
    return;
  }

  auto &merged = existing->second;
  if (merged.type != sv.type || merged.start != sv.start ||
      merged.end != sv.end || merged.length != sv.length) {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "canonical SV ID collision for '" + sv.id + "'");
  }
  const auto merge_unique = [](std::vector<std::string> &target,
                               const std::vector<std::string> &incoming) {
    target.insert(target.end(), incoming.begin(), incoming.end());
    std::sort(target.begin(), target.end());
    target.erase(std::unique(target.begin(), target.end()), target.end());
  };
  merge_unique(merged.supporting_reads, sv.supporting_reads);
  merge_unique(merged.evidence_sources, sv.evidence_sources);
  merge_unique(merged.orientations, sv.orientations);
  merged.known_event = merged.known_event || sv.known_event;
  merged.segment_count = std::max(merged.segment_count, sv.segment_count);
}

void merge_complex_sv(std::map<std::string, ComplexSvCall> &events,
                      ComplexSvCall event) {
  const auto existing = events.find(event.id);
  if (existing == events.end()) {
    std::sort(event.supporting_reads.begin(), event.supporting_reads.end());
    event.supporting_reads.erase(std::unique(event.supporting_reads.begin(),
                                             event.supporting_reads.end()),
                                 event.supporting_reads.end());
    events.emplace(event.id, std::move(event));
    return;
  }

  auto &merged = existing->second;
  if (merged.junction_ids != event.junction_ids ||
      merged.junction_orientations != event.junction_orientations) {
    throw AnalysisError(AnalysisErrorCode::internal_error,
                        "canonical complex-SV ID collision for '" + event.id +
                            "'");
  }
  merged.supporting_reads.insert(merged.supporting_reads.end(),
                                 event.supporting_reads.begin(),
                                 event.supporting_reads.end());
  std::sort(merged.supporting_reads.begin(), merged.supporting_reads.end());
  merged.supporting_reads.erase(std::unique(merged.supporting_reads.begin(),
                                            merged.supporting_reads.end()),
                                merged.supporting_reads.end());
  merged.segment_count = std::max(merged.segment_count, event.segment_count);
}

} // namespace mito::detail
