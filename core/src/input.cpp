#include "detail/pipeline.hpp"
#ifdef MITO_HAS_HTSLIB
#include <htslib/sam.h>
#endif

namespace mito::detail {

[[nodiscard]] std::vector<SupplementaryAlignment>
parse_supplementary_alignments(std::string_view encoded_tag,
                               std::string_view read_id) {
  if (encoded_tag.empty()) {
    throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                        "empty SA tag for read '" + std::string(read_id) + "'");
  }

  std::vector<SupplementaryAlignment> alignments;
  std::size_t entry_start = 0;
  std::size_t entry_index = 1;
  while (entry_start < encoded_tag.size()) {
    const auto entry_end = encoded_tag.find(';', entry_start);
    if (entry_end == std::string_view::npos) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "unterminated SA entry " +
                              std::to_string(entry_index) + " for read '" +
                              std::string(read_id) + "'");
    }
    const auto entry = encoded_tag.substr(entry_start, entry_end - entry_start);
    if (entry.empty()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "empty SA entry " + std::to_string(entry_index) +
                              " for read '" + std::string(read_id) + "'");
    }

    std::array<std::string_view, 6> fields{};
    std::size_t field_start = 0;
    std::size_t field_count = 0;
    bool has_extra_field = false;
    while (field_count < fields.size()) {
      const auto comma = entry.find(',', field_start);
      fields[field_count++] =
          entry.substr(field_start, comma == std::string_view::npos
                                        ? entry.size() - field_start
                                        : comma - field_start);
      if (comma == std::string_view::npos) {
        field_start = entry.size();
        break;
      }
      field_start = comma + 1U;
      if (field_count == fields.size()) {
        has_extra_field = true;
      }
    }
    if (field_count != fields.size() || has_extra_field ||
        field_start != entry.size()) {
      throw AnalysisError(
          AnalysisErrorCode::input_parse_failed,
          "SA entry " + std::to_string(entry_index) +
              " must contain exactly 6 comma-separated fields for read '" +
              std::string(read_id) + "'");
    }

    const auto position = parse_size(fields[1]);
    const auto mapq = parse_size(fields[4]);
    const auto edit_distance = parse_size(fields[5]);
    if (fields[0].empty() || fields[0] == "*" || !position || *position == 0 ||
        (fields[2] != "+" && fields[2] != "-") || !mapq ||
        *mapq > std::numeric_limits<std::uint8_t>::max() || !edit_distance) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid SA entry " + std::to_string(entry_index) +
                              " for read '" + std::string(read_id) + "'");
    }

    SupplementaryAlignment alignment;
    alignment.reference_name = std::string(fields[0]);
    alignment.reference_start = *position;
    alignment.strand = fields[2].front();
    alignment.cigar = std::string(fields[3]);
    alignment.cigar_operations =
        parse_cigar(fields[3], std::string(read_id) + " SA");
    if (alignment.cigar_operations.empty()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "SA entry " + std::to_string(entry_index) +
                              " has no alignment CIGAR for read '" +
                              std::string(read_id) + "'");
    }
    alignment.mapping_quality = static_cast<std::uint8_t>(*mapq);
    alignment.edit_distance = *edit_distance;
    alignments.push_back(std::move(alignment));

    entry_start = entry_end + 1U;
    ++entry_index;
  }
  return alignments;
}

[[nodiscard]] std::map<std::string, std::string>
parse_sam_aux_tags(const std::vector<std::string> &fields,
                   std::size_t line_number) {
  std::map<std::string, std::string> tags;
  for (std::size_t i = 11; i < fields.size(); ++i) {
    if (fields[i].size() < 5 || fields[i][2] != ':' || fields[i][4] != ':') {
      if (fields[i].rfind("SA", 0) == 0) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "malformed SA auxiliary field at SAM line " +
                                std::to_string(line_number));
      }
      continue;
    }
    if (fields[i].compare(0, 2, "SA") == 0 && fields[i][3] != 'Z') {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "SA auxiliary field must use SAM type Z at line " +
                              std::to_string(line_number));
    }
    tags.emplace(fields[i].substr(0, 2), fields[i].substr(5));
  }
  return tags;
}

[[nodiscard]] std::string snp_key(std::size_t position, char reference,
                                  char alternate) {
  std::string key = std::to_string(position);
  key.push_back(':');
  key.push_back(
      static_cast<char>(std::toupper(static_cast<unsigned char>(reference))));
  key.push_back(':');
  key.push_back(
      static_cast<char>(std::toupper(static_cast<unsigned char>(alternate))));
  return key;
}

[[nodiscard]] InputRecords parse_fastq(std::istream &input,
                                       const AnalysisConfig &config) {
  InputRecords records;
  std::string header;
  std::string sequence;
  std::string plus;
  std::string quality;

  while (std::getline(input, header)) {
    throw_if_cancelled(config);
    strip_trailing_carriage_return(header);
    if (header.empty()) {
      continue;
    }
    if (header.front() != '@') {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "FASTQ parser expected a header starting with '@'");
    }
    if (!std::getline(input, sequence) || !std::getline(input, plus) ||
        !std::getline(input, quality)) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "FASTQ record is truncated");
    }
    strip_trailing_carriage_return(sequence);
    strip_trailing_carriage_return(plus);
    strip_trailing_carriage_return(quality);
    if (plus.empty() || plus.front() != '+') {
      throw AnalysisError(
          AnalysisErrorCode::input_parse_failed,
          "FASTQ parser expected a separator starting with '+' for read '" +
              header.substr(1) + "'");
    }
    if (sequence.empty()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "FASTQ sequence is empty for read '" +
                              header.substr(1) + "'");
    }
    if (quality.size() != sequence.size()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "FASTQ sequence/quality length mismatch for read '" +
                              header.substr(1) + "'");
    }
    ReadRecord record;
    record.id = header.substr(1);
    record.sequence = sequence;
    record.qualities = quality;
    records.reads.push_back(std::move(record));
  }

  return records;
}

[[nodiscard]] InputRecords parse_sam(std::istream &input,
                                     const AnalysisConfig &config) {
  InputRecords records;
  records.alignment_input = true;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    throw_if_cancelled(config);
    ++line_number;
    strip_trailing_carriage_return(line);
    if (line.rfind("@SQ\t", 0) == 0) {
      const auto header_fields = split_tab(line);
      for (const auto &field : header_fields) {
        if (field.rfind("SN:", 0) == 0 &&
            looks_like_nuclear_contig(field.substr(3))) {
          records.has_nuclear_contigs = true;
        }
      }
      continue;
    }
    if (line.empty() || line.front() == '@') {
      continue;
    }

    const auto fields = split_tab(line);
    if (fields.size() < 11) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "SAM record has fewer than 11 fields at line " +
                              std::to_string(line_number));
    }

    const auto flags = parse_size(fields[1]);
    const auto reference_start = parse_size(fields[3]);
    const auto mapping_quality = parse_size(fields[4]);
    if (!flags || *flags > std::numeric_limits<std::uint16_t>::max()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid SAM FLAG at line " +
                              std::to_string(line_number));
    }
    if (!reference_start) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid SAM POS at line " +
                              std::to_string(line_number));
    }
    if (!mapping_quality ||
        *mapping_quality > std::numeric_limits<std::uint8_t>::max()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid SAM MAPQ at line " +
                              std::to_string(line_number));
    }
    if (fields[9] != "*" && fields[10] != "*" &&
        fields[9].size() != fields[10].size()) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "SAM sequence/quality length mismatch at line " +
                              std::to_string(line_number));
    }

    ReadRecord record;
    record.id = fields[0];
    record.flags = static_cast<std::uint16_t>(*flags);
    record.reference_name = fields[2];
    record.reference_start = *reference_start;
    record.mapping_quality = static_cast<std::uint8_t>(*mapping_quality);
    record.cigar = fields[5];
    record.cigar_operations = parse_cigar(record.cigar, record.id);
    record.sequence = fields[9] == "*" ? std::string{} : fields[9];
    record.qualities = fields[10] == "*" ? std::string{} : fields[10];
    if (!record.sequence.empty() && !record.cigar_operations.empty()) {
      std::size_t query_length = 0;
      for (const auto &operation : record.cigar_operations) {
        if (operation.code == 'M' || operation.code == 'I' ||
            operation.code == 'S' || operation.code == '=' ||
            operation.code == 'X') {
          if (operation.length >
              std::numeric_limits<std::size_t>::max() - query_length) {
            throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                                "SAM CIGAR query length overflow at line " +
                                    std::to_string(line_number));
          }
          query_length += operation.length;
        }
      }
      if (query_length != record.sequence.size()) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "SAM CIGAR/query length mismatch at line " +
                                std::to_string(line_number));
      }
    }
    record.aux_tags = parse_sam_aux_tags(fields, line_number);
    records.reads.push_back(std::move(record));
  }
  return records;
}

#ifdef MITO_HAS_HTSLIB
class MitoFileReader {
public:
  explicit MitoFileReader(std::string path) : path_(std::move(path)) {}

  [[nodiscard]] InputRecords read_all(const AnalysisConfig &config) const {
    htsFile *raw_file = sam_open(path_.c_str(), "r");
    if (raw_file == nullptr) {
      throw AnalysisError(AnalysisErrorCode::input_open_failed,
                          "htslib could not open alignment file: " + path_);
    }
    HtsFileHandle file(raw_file);

    bam_hdr_t *raw_header = sam_hdr_read(file.get());
    if (raw_header == nullptr) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "htslib could not read alignment header: " + path_);
    }
    BamHeaderHandle header(raw_header);

    bam1_t *raw_record = bam_init1();
    if (raw_record == nullptr) {
      throw AnalysisError(AnalysisErrorCode::resource_exhausted,
                          "htslib could not allocate BAM record");
    }
    BamRecordHandle record(raw_record);

    InputRecords records;
    records.alignment_input = true;
    for (int i = 0; i < header.get()->n_targets; ++i) {
      if (looks_like_nuclear_contig(header.get()->target_name[i])) {
        records.has_nuclear_contigs = true;
        break;
      }
    }
    int read_status = 0;
    while ((read_status = sam_read1(file.get(), header.get(), record.get())) >=
           0) {
      throw_if_cancelled(config);
      records.reads.push_back(convert_record(*header.get(), *record.get()));
    }
    if (read_status < -1) {
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "htslib failed while reading alignment records: " +
                              path_);
    }
    return records;
  }

private:
  struct HtsFileHandle {
    explicit HtsFileHandle(htsFile *value) : value_(value) {}
    ~HtsFileHandle() { hts_close(value_); }
    HtsFileHandle(const HtsFileHandle &) = delete;
    HtsFileHandle &operator=(const HtsFileHandle &) = delete;
    [[nodiscard]] htsFile *get() const { return value_; }
    htsFile *value_;
  };

  struct BamHeaderHandle {
    explicit BamHeaderHandle(bam_hdr_t *value) : value_(value) {}
    ~BamHeaderHandle() { bam_hdr_destroy(value_); }
    BamHeaderHandle(const BamHeaderHandle &) = delete;
    BamHeaderHandle &operator=(const BamHeaderHandle &) = delete;
    [[nodiscard]] bam_hdr_t *get() const { return value_; }
    bam_hdr_t *value_;
  };

  struct BamRecordHandle {
    explicit BamRecordHandle(bam1_t *value) : value_(value) {}
    ~BamRecordHandle() { bam_destroy1(value_); }
    BamRecordHandle(const BamRecordHandle &) = delete;
    BamRecordHandle &operator=(const BamRecordHandle &) = delete;
    [[nodiscard]] bam1_t *get() const { return value_; }
    bam1_t *value_;
  };

  [[nodiscard]] static ReadRecord convert_record(const bam_hdr_t &header,
                                                 const bam1_t &record) {
    ReadRecord out;
    out.id = bam_get_qname(&record);
    out.flags = record.core.flag;
    out.mapping_quality = record.core.qual;
    out.reference_start =
        record.core.pos < 0 ? 1 : static_cast<std::size_t>(record.core.pos) + 1;
    if (record.core.tid >= 0 && record.core.tid < header.n_targets) {
      out.reference_name = header.target_name[record.core.tid];
    }
    out.cigar = cigar_string(record);
    out.cigar_operations = parse_cigar(out.cigar, out.id);
    out.sequence = sequence_string(record);
    out.qualities = quality_string(record);
    out.aux_tags = aux_tags(record);
    return out;
  }

  [[nodiscard]] static std::string cigar_string(const bam1_t &record) {
    std::string cigar;
    const auto *raw_cigar = bam_get_cigar(&record);
    for (std::uint32_t i = 0; i < record.core.n_cigar; ++i) {
      cigar += std::to_string(bam_cigar_oplen(raw_cigar[i]));
      cigar.push_back(bam_cigar_opchr(raw_cigar[i]));
    }
    return cigar;
  }

  [[nodiscard]] static std::string sequence_string(const bam1_t &record) {
    std::string sequence;
    sequence.reserve(static_cast<std::size_t>(record.core.l_qseq));
    const auto *seq = bam_get_seq(&record);
    for (int i = 0; i < record.core.l_qseq; ++i) {
      sequence.push_back(seq_nt16_str[bam_seqi(seq, i)]);
    }
    return sequence;
  }

  [[nodiscard]] static std::string quality_string(const bam1_t &record) {
    std::string qualities;
    qualities.reserve(static_cast<std::size_t>(record.core.l_qseq));
    const auto *qual = bam_get_qual(&record);
    for (int i = 0; i < record.core.l_qseq; ++i) {
      qualities.push_back(qual[i] == 0xFFU ? '!'
                                           : static_cast<char>(qual[i] + 33U));
    }
    return qualities;
  }

  [[nodiscard]] static std::map<std::string, std::string>
  aux_tags(const bam1_t &record) {
    std::map<std::string, std::string> tags;
    errno = 0;
    const std::uint8_t *value = bam_aux_first(&record);
    if (value == nullptr) {
      if (errno == ENOENT) {
        return tags;
      }
      throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                          "invalid auxiliary fields for read '" +
                              std::string(bam_get_qname(&record)) + "'");
    }
    while (value != nullptr) {
      const auto *raw_tag = bam_aux_tag(value);
      const std::string tag(raw_tag, raw_tag + 2);
      const char type = bam_aux_type(value);
      if (tag == "SA" && type != 'Z') {
        throw AnalysisError(
            AnalysisErrorCode::input_parse_failed,
            "SA auxiliary field must use SAM type Z for read '" +
                std::string(bam_get_qname(&record)) + "'");
      }
      std::string encoded;
      if (type == 'Z' || type == 'H') {
        if (const char *text = bam_aux2Z(value)) {
          encoded = text;
        } else {
          throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                              "invalid string auxiliary field " + tag +
                                  " for read '" +
                                  std::string(bam_get_qname(&record)) + "'");
        }
      } else if (type == 'A') {
        encoded.assign(1U, bam_aux2A(value));
      } else if (type == 'c' || type == 'C' || type == 's' || type == 'S' ||
                 type == 'i' || type == 'I') {
        encoded = std::to_string(bam_aux2i(value));
      } else if (type == 'f' || type == 'd') {
        std::ostringstream out;
        out << std::setprecision(std::numeric_limits<double>::max_digits10)
            << bam_aux2f(value);
        encoded = std::move(out).str();
      } else if (type == 'B') {
        const char element_type = static_cast<char>(value[1]);
        const auto length = bam_auxB_len(value);
        std::ostringstream out;
        out << element_type;
        for (std::uint32_t index = 0U; index < length; ++index) {
          out << ',';
          if (element_type == 'f') {
            out << std::setprecision(std::numeric_limits<double>::max_digits10)
                << bam_auxB2f(value, index);
          } else {
            out << bam_auxB2i(value, index);
          }
        }
        encoded = std::move(out).str();
      } else {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "unsupported auxiliary field type for " + tag +
                                " on read '" +
                                std::string(bam_get_qname(&record)) + "'");
      }
      if (!tags.emplace(tag, std::move(encoded)).second) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "duplicate auxiliary field " + tag + " for read '" +
                                std::string(bam_get_qname(&record)) + "'");
      }
      errno = 0;
      value = bam_aux_next(&record, value);
      if (value == nullptr && errno != ENOENT) {
        throw AnalysisError(AnalysisErrorCode::input_parse_failed,
                            "invalid auxiliary fields for read '" +
                                std::string(bam_get_qname(&record)) + "'");
      }
    }
    return tags;
  }

  std::string path_;
};
#endif

[[nodiscard]] bool
has_extension(const std::filesystem::path &path,
              std::initializer_list<std::string_view> extensions) {
  auto ext = path.extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char value) {
    return static_cast<char>(std::tolower(value));
  });
  return std::any_of(extensions.begin(), extensions.end(),
                     [&](std::string_view value) { return ext == value; });
}

[[nodiscard]] InputRecords read_input(const std::string &input_path,
                                      const AnalysisConfig &config) {
  const std::filesystem::path path(input_path);
  if (has_extension(path, {".fa", ".fas", ".fasta", ".fna"})) {
    throw AnalysisError(
        AnalysisErrorCode::input_format_unsupported,
        "FASTA input is not supported; use FASTQ, SAM, BAM, or CRAM");
  }

#ifdef MITO_HAS_HTSLIB
  if (has_extension(path, {".sam", ".bam", ".cram"})) {
    return MitoFileReader(input_path).read_all(config);
  }
#endif

  std::ifstream input(input_path, std::ios::binary);
  if (!input) {
    throw AnalysisError(AnalysisErrorCode::input_open_failed,
                        "could not open input file: " + input_path);
  }

  if (has_extension(path, {".sam"})) {
    return parse_sam(input, config);
  }
  if (has_extension(path, {".bam", ".cram"})) {
    throw AnalysisError(AnalysisErrorCode::dependency_unavailable,
                        "BAM/CRAM support requires htslib; install htslib "
                        "development headers and rebuild");
  }
  return parse_fastq(input, config);
}

} // namespace mito::detail
