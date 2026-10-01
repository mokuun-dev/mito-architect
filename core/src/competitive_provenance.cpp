#include "detail/pipeline.hpp"

namespace mito::detail {
namespace {

[[nodiscard]] bool is_sha256(std::string_view value) {
  return value.size() == 64U && std::all_of(value.begin(), value.end(),
                                             [](unsigned char byte) {
    return std::isxdigit(byte) != 0;
  });
}

[[nodiscard]] std::string lower_ascii(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char byte) {
                   return static_cast<char>(std::tolower(byte));
                 });
  return value;
}

[[nodiscard]] std::uint32_t sha256_word(const std::array<std::uint8_t, 64U> &b,
                                        const std::size_t offset) {
  return (static_cast<std::uint32_t>(b[offset]) << 24U) |
         (static_cast<std::uint32_t>(b[offset + 1U]) << 16U) |
         (static_cast<std::uint32_t>(b[offset + 2U]) << 8U) |
         static_cast<std::uint32_t>(b[offset + 3U]);
}

[[nodiscard]] std::string sha256_file(const std::filesystem::path &path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    return {};
  }
  constexpr std::array<std::uint32_t, 64U> constants = {
      0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU,
      0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U, 0xd807aa98U, 0x12835b01U,
      0x243185beU, 0x550c7dc3U, 0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U,
      0xc19bf174U, 0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
      0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU, 0x983e5152U,
      0xa831c66dU, 0xb00327c8U, 0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U,
      0x06ca6351U, 0x14292967U, 0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU,
      0x53380d13U, 0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
      0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U, 0xd192e819U,
      0xd6990624U, 0xf40e3585U, 0x106aa070U, 0x19a4c116U, 0x1e376c08U,
      0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU,
      0x682e6ff3U, 0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
      0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U};
  std::array<std::uint32_t, 8U> state = {0x6a09e667U, 0xbb67ae85U,
                                          0x3c6ef372U, 0xa54ff53aU,
                                          0x510e527fU, 0x9b05688cU,
                                          0x1f83d9abU, 0x5be0cd19U};
  std::array<std::uint8_t, 64U> block{};
  std::size_t used = 0U;
  std::uint64_t bytes = 0U;
  const auto compress = [&] {
    std::array<std::uint32_t, 64U> words{};
    for (std::size_t index = 0U; index < 16U; ++index) {
      words[index] = sha256_word(block, index * 4U);
    }
    for (std::size_t index = 16U; index < words.size(); ++index) {
      const auto s0 = std::rotr(words[index - 15U], 7U) ^
                      std::rotr(words[index - 15U], 18U) ^
                      (words[index - 15U] >> 3U);
      const auto s1 = std::rotr(words[index - 2U], 17U) ^
                      std::rotr(words[index - 2U], 19U) ^
                      (words[index - 2U] >> 10U);
      words[index] = words[index - 16U] + s0 + words[index - 7U] + s1;
    }
    auto [a, b, c, d, e, f, g, h] = state;
    for (std::size_t index = 0U; index < words.size(); ++index) {
      const auto s1 = std::rotr(e, 6U) ^ std::rotr(e, 11U) ^ std::rotr(e, 25U);
      const auto choose = (e & f) ^ (~e & g);
      const auto temp1 = h + s1 + choose + constants[index] + words[index];
      const auto s0 = std::rotr(a, 2U) ^ std::rotr(a, 13U) ^ std::rotr(a, 22U);
      const auto majority = (a & b) ^ (a & c) ^ (b & c);
      const auto temp2 = s0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + temp1;
      d = c;
      c = b;
      b = a;
      a = temp1 + temp2;
    }
    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
    state[4] += e;
    state[5] += f;
    state[6] += g;
    state[7] += h;
  };
  char byte = 0;
  while (input.get(byte)) {
    if (bytes == std::numeric_limits<std::uint64_t>::max()) {
      return {};
    }
    block[used++] = static_cast<std::uint8_t>(byte);
    ++bytes;
    if (used == block.size()) {
      compress();
      used = 0U;
    }
  }
  if (!input.eof()) {
    return {};
  }
  const auto bit_length = bytes * 8U;
  block[used++] = 0x80U;
  if (used > 56U) {
    std::fill(block.begin() + static_cast<std::ptrdiff_t>(used), block.end(), 0U);
    compress();
    used = 0U;
  }
  std::fill(block.begin() + static_cast<std::ptrdiff_t>(used),
            block.begin() + 56, 0U);
  for (std::size_t index = 0U; index < 8U; ++index) {
    block[63U - index] = static_cast<std::uint8_t>(bit_length >> (index * 8U));
  }
  compress();

  std::ostringstream out;
  out << std::hex << std::setfill('0');
  for (const auto word : state) {
    out << std::setw(8) << word;
  }
  return out.str();
}

[[nodiscard]] std::optional<std::map<std::string, std::string>>
read_manifest(const std::filesystem::path &path) {
  std::ifstream input(path);
  std::string line;
  if (!std::getline(input, line) || line != "field\tvalue") {
    return std::nullopt;
  }
  std::map<std::string, std::string> fields;
  while (std::getline(input, line)) {
    const auto tab = line.find('\t');
    if (tab == std::string::npos || tab == 0U ||
        line.find('\t', tab + 1U) != std::string::npos) {
      return std::nullopt;
    }
    auto [_, inserted] = fields.emplace(line.substr(0U, tab),
                                        line.substr(tab + 1U));
    if (!inserted) {
      return std::nullopt;
    }
  }
  return input.eof() ? std::optional{std::move(fields)} : std::nullopt;
}

[[nodiscard]] const std::string *field(
    const std::map<std::string, std::string> &fields, std::string_view name) {
  const auto found = fields.find(std::string(name));
  return found == fields.end() ? nullptr : &found->second;
}

[[nodiscard]] bool same_path(const std::filesystem::path &lhs,
                             const std::filesystem::path &rhs) {
  std::error_code error;
  return std::filesystem::weakly_canonical(lhs, error) ==
             std::filesystem::weakly_canonical(rhs, error) &&
         !error;
}

[[nodiscard]] bool valid_competitive_reference(
    const std::filesystem::path &manifest_path,
    const std::string &expected_hash, const std::string &expected_index_path,
    const std::string &expected_index_hash) {
  if (!is_sha256(expected_hash) || !is_sha256(expected_index_hash) ||
      lower_ascii(sha256_file(manifest_path)) != lower_ascii(expected_hash)) {
    return false;
  }
  const auto manifest = read_manifest(manifest_path);
  if (!manifest || field(*manifest, "schema_version") == nullptr ||
      *field(*manifest, "schema_version") != "competitive-reference-1.0") {
    return false;
  }
  const auto *nuclear_hash = field(*manifest, "nuclear_fasta_sha256");
  const auto *mito_hash = field(*manifest, "mitochondrial_fasta_sha256");
  const auto *index_path = field(*manifest, "reference_index_path");
  const auto *index_hash = field(*manifest, "reference_index_sha256");
  if (nuclear_hash == nullptr || mito_hash == nullptr || index_path == nullptr ||
      index_hash == nullptr || !is_sha256(*nuclear_hash) ||
      !is_sha256(*mito_hash) || !is_sha256(*index_hash) ||
      !same_path(*index_path, expected_index_path) ||
      lower_ascii(*index_hash) != lower_ascii(expected_index_hash)) {
    return false;
  }
  return lower_ascii(sha256_file(*index_path)) == lower_ascii(*index_hash);
}

} // namespace

[[nodiscard]] NumtAssessment assess_numt_provenance(
    const std::string &input_path, const InputRecords &input) {
  NumtAssessment result;
  result.nuclear_contigs_present = input.has_nuclear_contigs;
  if (!input.alignment_input) {
    result.provenance_reason = "input is not an alignment";
    return result;
  }
  result.mode = "mt_only_or_unknown";
  result.provenance_status = "NOT_VERIFIED";
  if (!input.has_nuclear_contigs) {
    result.provenance_reason = "alignment header has no nuclear contigs";
    return result;
  }
  const std::filesystem::path bam_path(input_path);
  const auto manifest_path = bam_path.parent_path() / "competitive-mapping-manifest.tsv";
  const auto manifest = read_manifest(manifest_path);
  if (!manifest) {
    result.provenance_reason = "competitive mapping manifest is absent or malformed";
    return result;
  }
  const auto *schema = field(*manifest, "schema_version");
  const auto *assessability = field(*manifest, "numt_assessability");
  const auto *output_bam = field(*manifest, "output_bam_sha256");
  const auto *output_bai = field(*manifest, "output_bai_sha256");
  const auto *reference_manifest = field(*manifest, "reference_manifest_path");
  const auto *reference_manifest_hash = field(*manifest, "reference_manifest_sha256");
  const auto *reference_index = field(*manifest, "reference_index_path");
  const auto *reference_index_hash = field(*manifest, "reference_index_sha256");
  if (schema == nullptr || *schema != "competitive-mapping-1.1" ||
      assessability == nullptr ||
      *assessability != "ASSESSABLE_BY_VERSIONED_COMPETITIVE_ALIGNMENT_INPUT" ||
      output_bam == nullptr || output_bai == nullptr || reference_manifest == nullptr ||
      reference_manifest_hash == nullptr || reference_index == nullptr ||
      reference_index_hash == nullptr || !is_sha256(*output_bam) ||
      !is_sha256(*output_bai) || !same_path(bam_path.parent_path() / "competitive.bam", bam_path)) {
    result.provenance_reason = "competitive mapping manifest does not match the input contract";
    return result;
  }
  const auto bai_path = std::filesystem::path(input_path + ".bai");
  if (lower_ascii(sha256_file(bam_path)) != lower_ascii(*output_bam) ||
      lower_ascii(sha256_file(bai_path)) != lower_ascii(*output_bai)) {
    result.provenance_reason = "competitive mapping artifact checksum does not match its manifest";
    return result;
  }
  if (!valid_competitive_reference(*reference_manifest, *reference_manifest_hash,
                                   *reference_index, *reference_index_hash)) {
    result.provenance_reason = "competitive reference provenance cannot be verified";
    return result;
  }
  result.mode = "competitive_alignment";
  result.specificity_assessable = true;
  result.provenance_status = "VERIFIED";
  result.provenance_reason = "versioned competitive mapping and reference checksums verified";
  return result;
}

} // namespace mito::detail
