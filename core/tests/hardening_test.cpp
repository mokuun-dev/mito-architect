#include "detail/pipeline.hpp"

#ifdef MITO_HAS_HTSLIB
#include <htslib/faidx.h>
#include <htslib/sam.h>

#include <stdlib.h> // setenv: htslib builds are POSIX
#endif

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <random>

namespace {

using namespace mito;
using namespace mito::detail;

int failures = 0;

// Unique per test run so parallel ctest invocations do not collide.
const std::string run_suffix = std::to_string(std::random_device{}());

void require(const bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAILED: " << message << "\n";
    ++failures;
  }
}

// Reference DBSCAN: the pre-optimization algorithm, which compared every
// profile with every other using string tokens.
std::vector<int> brute_force_dbscan(const std::vector<ReadFeature> &features,
                                    const double epsilon,
                                    const std::size_t min_cluster_size) {
  std::vector<std::size_t> active;
  for (std::size_t i = 0; i < features.size(); ++i) {
    if (!features[i].filtered_numt) {
      active.push_back(i);
    }
  }
  std::vector<FeatureTokens> tokens;
  std::vector<std::vector<std::size_t>> members;
  std::map<FeatureTokens, std::size_t> lookup;
  for (const auto index : active) {
    auto profile = feature_tokens(features[index]);
    const auto found = lookup.find(profile);
    if (found != lookup.end()) {
      members[found->second].push_back(index);
      continue;
    }
    lookup.emplace(profile, tokens.size());
    tokens.push_back(std::move(profile));
    members.push_back({index});
  }
  const auto distance = [&](const std::size_t lhs, const std::size_t rhs) {
    std::vector<std::string> common;
    std::set_intersection(tokens[lhs].begin(), tokens[lhs].end(),
                          tokens[rhs].begin(), tokens[rhs].end(),
                          std::back_inserter(common));
    const auto union_size =
        tokens[lhs].size() + tokens[rhs].size() - common.size();
    return union_size == 0U ? 0.0
                            : 1.0 - (static_cast<double>(common.size()) /
                                     static_cast<double>(union_size));
  };
  const auto region = [&](const std::size_t point) {
    std::vector<std::size_t> neighbors;
    std::size_t reads = 0;
    for (std::size_t candidate = 0; candidate < tokens.size(); ++candidate) {
      const auto lhs = tokens[point].size();
      const auto rhs = tokens[candidate].size();
      if (static_cast<double>(std::min(lhs, rhs)) /
              static_cast<double>(std::max(lhs, rhs)) <
          1.0 - epsilon) {
        continue;
      }
      if (distance(point, candidate) <= epsilon) {
        neighbors.push_back(candidate);
        reads += members[candidate].size();
      }
    }
    return std::pair{neighbors, reads};
  };
  constexpr int kUnvisited = -2;
  std::vector<int> labels(tokens.size(), kUnvisited);
  int next_cluster = 0;
  for (std::size_t point = 0; point < tokens.size(); ++point) {
    if (labels[point] != kUnvisited) {
      continue;
    }
    auto [neighbors, reads] = region(point);
    if (reads < min_cluster_size) {
      labels[point] = -1;
      continue;
    }
    const int cluster = next_cluster++;
    labels[point] = cluster;
    std::vector<unsigned char> queued(tokens.size(), 0U);
    for (const auto neighbor : neighbors) {
      queued[neighbor] = 1U;
    }
    for (std::size_t cursor = 0; cursor < neighbors.size(); ++cursor) {
      const auto neighbor = neighbors[cursor];
      if (labels[neighbor] == -1) {
        labels[neighbor] = cluster;
      }
      if (labels[neighbor] != kUnvisited) {
        continue;
      }
      labels[neighbor] = cluster;
      auto [expanded, expanded_reads] = region(neighbor);
      if (expanded_reads >= min_cluster_size) {
        for (const auto candidate : expanded) {
          if (queued[candidate] == 0U) {
            queued[candidate] = 1U;
            neighbors.push_back(candidate);
          }
        }
      }
    }
  }
  std::vector<int> result(features.size(), -1);
  for (std::size_t profile = 0; profile < tokens.size(); ++profile) {
    for (const auto feature : members[profile]) {
      result[feature] = labels[profile];
    }
  }
  return result;
}

void clustering_matches_brute_force() {
  std::mt19937_64 random(20261022U);
  for (int round = 0; round < 12; ++round) {
    // A small marker pool makes overlaps frequent; a few shared "true"
    // markers mimic homoplasmic variants on top of per-read noise.
    const std::size_t pool = 8U + static_cast<std::size_t>(round) * 6U;
    std::vector<ReadFeature> features(160);
    for (std::size_t index = 0; index < features.size(); ++index) {
      auto &feature = features[index];
      feature.id = "read" + std::to_string(index);
      feature.length = 1000U + random() % 5000U;
      feature.filtered_numt = random() % 17U == 0U;
      const auto marker_count = random() % 12U;
      for (std::size_t marker = 0; marker < marker_count; ++marker) {
        SnpCall snp;
        snp.position = 1U + random() % pool;
        snp.reference = 'A';
        snp.alternate = "CGT"[random() % 3U];
        feature.snps.push_back(snp);
      }
      if (random() % 4U == 0U) {
        feature.sv_ids.push_back("deletion:" + std::to_string(random() % 4U));
      }
      if (index % 3U == 0U) {
        SnpCall shared;
        shared.position = 263U;
        shared.reference = 'A';
        shared.alternate = 'G';
        feature.snps.push_back(shared);
      }
    }
    for (const double epsilon :
         {0.0, 0.1, 0.25, 0.5, 0.72, 0.9, 0.99, 0.9999999999, 1.0}) {
      for (const std::size_t min_size : {std::size_t{1}, std::size_t{2},
                                         std::size_t{3}, std::size_t{6}}) {
        auto clustered = features;
        AnalysisConfig config;
        assign_dbscan_clusters(clustered, epsilon, min_size, config);
        const auto expected = brute_force_dbscan(features, epsilon, min_size);
        for (std::size_t index = 0; index < features.size(); ++index) {
          if (features[index].filtered_numt) {
            continue;
          }
          if (clustered[index].cluster_id != expected[index] ||
              clustered[index].outlier != (expected[index] < 0)) {
            require(false, "DBSCAN differs from brute force (round " +
                               std::to_string(round) + ", epsilon " +
                               std::to_string(epsilon) + ", min " +
                               std::to_string(min_size) + ", read " +
                               std::to_string(index) + ")");
            break;
          }
        }
      }
    }
  }
}

std::string write_temp(const std::string &name, const std::string &content) {
  const auto path = std::filesystem::temp_directory_path() /
                    ("mito-hardening-" + run_suffix + "-" + name);
  std::ofstream(path, std::ios::binary) << content;
  return path.string();
}

AnalysisErrorCode error_code_of(const std::string &input,
                                const AnalysisConfig &config) {
  try {
    (void)AnalysisEngine().analyze(input, "", config);
  } catch (const AnalysisError &error) {
    return error.code();
  }
  return AnalysisErrorCode::internal_error; // sentinel: no error
}

void input_budget_fails_closed() {
  const auto fastq = write_temp(
      "budget.fastq", "@a\nACGTACGT\n+\nIIIIIIII\n@b\nACGTACGT\n+\nIIIIIIII\n"
                      "@c\nACGTACGT\n+\nIIIIIIII\n");
  AnalysisConfig config;
  config.max_input_records = 2U;
  require(error_code_of(fastq, config) == AnalysisErrorCode::resource_exhausted,
          "FASTQ record budget");
  config.max_input_records = 3U;
  config.max_input_bases = 23U;
  require(error_code_of(fastq, config) == AnalysisErrorCode::resource_exhausted,
          "FASTQ base budget");
  config.max_input_bases = 24U;
  require(error_code_of(fastq, config) == AnalysisErrorCode::internal_error,
          "FASTQ within budget succeeds");
  config = {};
  config.max_input_records = 0U;
  require(error_code_of(fastq, config) ==
              AnalysisErrorCode::invalid_configuration,
          "zero record budget is rejected");
  std::filesystem::remove(fastq);
}

std::string rcrs() {
  std::ifstream in(bundled_rcrs_path());
  return parse_reference_fasta(in, bundled_rcrs_path());
}

void origin_spanning_deletion_is_reported() {
  const auto reference = rcrs();
  const auto sequence = reference.substr(16549U, 10U) + reference.substr(0U, 10U);
  const auto sam = write_temp(
      "origin.sam", "@SQ\tSN:NC_012920.1\tLN:16569\nwrap\t0\tNC_012920.1\t16550"
                    "\t60\t10M15D10M\t*\t0\t0\t" +
                        sequence + "\t" + std::string(20U, 'I') + "\n");
  try {
    const auto json = AnalysisEngine().analyze(sam, "", {});
    require(json.find("\"deletion:16560-5\"") != std::string::npos,
            "circular deletion ID is reported");
  } catch (const std::exception &error) {
    require(false, std::string("origin-spanning deletion aborted: ") +
                       error.what());
  }
  std::filesystem::remove(sam);
}

#ifdef MITO_HAS_HTSLIB
std::string md5_of(const std::string &sequence) {
  auto *context = hts_md5_init();
  hts_md5_update(context, sequence.data(),
                 static_cast<unsigned long>(sequence.size()));
  unsigned char digest[16];
  hts_md5_final(digest, context);
  hts_md5_destroy(context);
  char hex[33];
  hts_md5_hex(hex, digest);
  return hex;
}

// Writes a CRAM against `fasta` whose header names the given contigs with
// explicit M5/UR tags, holding one read per contig.
void write_cram(const std::string &path, const std::string &fasta,
                const std::vector<std::pair<std::string, std::string>> &contigs) {
  std::string header_text;
  std::string reads;
  for (const auto &[name, sequence] : contigs) {
    header_text += "@SQ\tSN:" + name + "\tLN:" +
                   std::to_string(sequence.size()) + "\tM5:" + md5_of(sequence) +
                   "\tUR:" + fasta + "\n";
    reads += "r_" + name + "\t0\t" + name + "\t101\t60\t50M\t*\t0\t0\t" +
             sequence.substr(100U, 50U) + "\t" + std::string(50U, 'I') + "\n";
  }
  htsFile *out = sam_open(path.c_str(), "wc");
  hts_set_fai_filename(out, fasta.c_str());
  sam_hdr_t *header = sam_hdr_parse(header_text.size(), header_text.c_str());
  require(sam_hdr_write(out, header) == 0, "CRAM header written");
  bam1_t *record = bam_init1();
  std::size_t start = 0U;
  while (start < reads.size()) {
    const auto end = reads.find('\n', start);
    kstring_t line = KS_INITIALIZE;
    kputsn(reads.data() + start, end - start, &line);
    require(sam_parse1(&line, header, record) >= 0, "CRAM record parsed");
    require(sam_write1(out, header, record) >= 0, "CRAM record written");
    ks_free(&line);
    start = end + 1U;
  }
  bam_destroy1(record);
  sam_hdr_destroy(header);
  require(sam_close(out) == 0, "CRAM closed");
}

void cram_reference_policy() {
  // No remote MD5 service or shared cache during the test.
  const auto empty = std::filesystem::temp_directory_path() /
                     ("mito-hardening-refpath-" + run_suffix);
  std::filesystem::create_directories(empty);
  ::setenv("REF_PATH", (empty.string() + "/%s").c_str(), 1);
  ::setenv("REF_CACHE", (empty.string() + "/%s").c_str(), 1);

  const auto mito = rcrs();
  std::string nuclear(4000U, 'A');
  std::mt19937 random(7U);
  for (auto &base : nuclear) {
    base = "ACGT"[random() % 4U];
  }
  const auto fasta = write_temp(
      "refs.fa", ">chrM\n" + mito + "\n>chr1\n" + nuclear + "\n");
  require(fai_build(fasta.c_str()) == 0, "FASTA index built");

  // chrM matches the analysis reference by MD5 despite the different name.
  const auto mito_cram = write_temp("mito.cram", "");
  write_cram(mito_cram, fasta, {{"chrM", mito}});
  // Reads on chr1 need a reference the analysis does not have.
  const auto nuclear_cram = write_temp("nuclear.cram", "");
  write_cram(nuclear_cram, fasta, {{"chrM", mito}, {"chr1", nuclear}});

  AnalysisConfig strict;
  require(error_code_of(mito_cram, strict) == AnalysisErrorCode::internal_error,
          "MD5-matched CRAM decodes with external lookups disabled");
  require(error_code_of(nuclear_cram, strict) ==
              AnalysisErrorCode::input_parse_failed,
          "header UR is not followed with external lookups disabled");
  AnalysisConfig permissive;
  permissive.allow_external_cram_references = true;
  require(error_code_of(nuclear_cram, permissive) ==
              AnalysisErrorCode::internal_error,
          "header UR is followed when external lookups are allowed");

  for (const auto &path :
       {mito_cram, nuclear_cram, fasta, fasta + ".fai"}) {
    std::filesystem::remove(path);
  }
  std::filesystem::remove_all(empty);
}
#endif

} // namespace

int main() {
  clustering_matches_brute_force();
  input_budget_fails_closed();
  origin_spanning_deletion_is_reported();
#ifdef MITO_HAS_HTSLIB
  cram_reference_policy();
#endif
  if (failures != 0) {
    std::cerr << failures << " hardening check(s) failed\n";
    return 1;
  }
  return 0;
}
