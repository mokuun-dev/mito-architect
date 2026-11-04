#include "detail/pipeline.hpp"
#include "detail/precision.hpp"

// Sample-calibrated calling of SNVs, small indels and large deletions.
//
//  1. Pileup per molecule (primary and supplementary alignments together, so
//     reads crossing the origin count once), rebuilt without molecules whose
//     divergence from the sample consensus marks them as NUMT-like outliers.
//  2. SNVs: a beta-binomial error model per read-oriented trinucleotide context
//     and homopolymer class, tested on each strand and combined by Fisher's
//     method, with Benjamini-Hochberg control over all hypotheses
//     (deepSNV/Shearwater; Gerstung 2012, 2014).
//  3. Small indels: placements normalized on the circular reference (Tan 2015),
//     per-repeat-length error rates measured from the sample itself (the idea
//     behind DRAGstr), then a per-read pair HMM and an EM allele-fraction fit
//     (Mutect2), which also recovers reads whose gap the aligner fragmented.
//  4. Large deletions: junctions from CIGAR gaps and split alignments, corrected
//     for query overlap, clustered in (placement, length) because microhomology
//     shifts placement but not length, refined by maximum likelihood, and
//     quantified against wild-type molecules spanning both breakpoints
//     (LostArc, mitoSAlt, NanoDel).
namespace mito::detail {
namespace {

constexpr double kFalseDiscoveryRate = 0.01;
constexpr std::size_t kQualityBins = 8U;
constexpr std::size_t kProximityBins = 5U;
constexpr std::size_t kMinimumStrandDepth = 10U;
constexpr std::size_t kMaximumSmallIndel = 49U;
constexpr std::size_t kMinimumLargeDeletion = 50U;
constexpr std::size_t kIndelEdgeGuard = 8U;
constexpr std::size_t kJunctionAnchor = 100U;
constexpr std::size_t kSpanningFlank = 25U;
constexpr std::size_t kMaximumReportedFiltered = 200U;
constexpr std::size_t kMaximumHmmOtherReads = 150U;
constexpr std::size_t kJunctionPlacementTolerance = 30U;
constexpr std::size_t kJunctionLengthTolerance = 8U;
constexpr double kBackgroundEffectMultiple = 2.0;
constexpr double kWeakerStrandLimit = 0.01;
constexpr std::size_t kMinimumLowLevelSupport = 5U;
constexpr std::size_t kGapWindow = 5U;
constexpr double kLowLevelFraction = 0.25;
constexpr std::size_t kUnreliableRun = 6U;
constexpr double kSlipArtifactMultiple = 2.0;
constexpr std::size_t kFragmentedEventLength = 7U;
// Realignment flank cap. Human mtDNA repeats are tens of bases, but the span
// comes from the sample consensus, which a crafted input controls.
constexpr std::size_t kMaximumRescoreFlank = 100U;
// Work budgets, in banded DP cells (plus intervals scanned for deletions).
// Measured peaks on real and simulated runs (2026-11-04): realignment 6.6e5
// cells per molecule (GM12878 nanopore, 1.5e9 in all), breakpoint search 1.8e4
// per molecule and 3.2e7 in all. The budgets sit 40x or more above those and
// cap a crafted input at about a minute of CPU.
constexpr std::uint64_t kRescoreWorkBase = 20'000'000'000ULL;
constexpr std::uint64_t kRescoreWorkPerMolecule = 20'000'000ULL;
constexpr std::uint64_t kRefineWorkBase = 1'000'000'000ULL;
constexpr std::uint64_t kRefineWorkPerMolecule = 100'000ULL;

struct PileupCell {
  std::array<std::array<std::uint32_t, 2>, 4> counts{};
  std::array<std::array<std::uint32_t, kQualityBins>, 4> quality{};
  std::array<std::array<std::uint32_t, kProximityBins>, 4> proximity{};
  /** Bases within kGapWindow of a small gap in the same alignment. */
  std::array<std::uint32_t, 4> near_gap{};
};

struct RawIndel {
  bool deletion = false;
  std::size_t left = 0;
  std::size_t right = 0;
  std::string left_allele;
  std::string right_allele;
  std::array<std::uint32_t, 2> strands{};
  std::vector<std::uint32_t> molecules;
};

struct AlignedInterval {
  std::uint32_t start = 0;
  std::uint32_t end = 0;
  std::uint32_t molecule = 0;
};

struct Junction {
  std::uint32_t left = 0;  // Last retained base before the gap.
  std::uint32_t right = 0; // First retained base after the gap.
  std::uint32_t length = 0;
  std::uint32_t molecule = 0;
  std::uint8_t strand = 0;
  bool from_split = false;
  /** Index in the representative's SEQ just after the left part. */
  std::uint32_t query_index = 0;
  /** SEQ runs against the reference for this junction. */
  bool reverse_sequence = false;
};

struct MoleculeView {
  std::uint32_t index = 0;
  const ReadRecord *representative = nullptr;
  std::vector<const ReadRecord *> records;
};

struct WorkerPileup {
  std::vector<PileupCell> cells;
  std::array<std::vector<std::int64_t>, 2> span_delta;
  std::unordered_map<std::string, RawIndel> indels;
  /** Evidence entry by raw deletion placement. Reads repeat the same placements,
   * and normalizing and keying each one is most of the cost of a noisy read. */
  std::unordered_map<std::uint64_t, RawIndel *> deletion_entries;
  std::vector<AlignedInterval> intervals;
  std::vector<std::uint32_t> stamp;
};

struct Pileup {
  std::vector<PileupCell> cells;
  std::array<std::vector<std::uint32_t>, 2> span;
  std::vector<RawIndel> indels;
  std::vector<AlignedInterval> intervals;
};

[[nodiscard]] constexpr std::size_t fold(const std::size_t position,
                                         const std::size_t length) noexcept {
  return (position - 1U) % length + 1U;
}

[[nodiscard]] constexpr std::size_t quality_bin(const std::uint8_t quality) noexcept {
  const std::size_t bin = quality / 5U;
  return bin < kQualityBins ? bin : kQualityBins - 1U;
}

[[nodiscard]] constexpr char complement_base(const char base) noexcept {
  switch (base) {
  case 'A':
    return 'T';
  case 'C':
    return 'G';
  case 'G':
    return 'C';
  case 'T':
    return 'A';
  default:
    return 'N';
  }
}

[[nodiscard]] std::string reverse_complement_string(const std::string_view sequence) {
  std::string result(sequence.size(), 'N');
  for (std::size_t index = 0U; index < sequence.size(); ++index) {
    result[sequence.size() - 1U - index] =
        complement_base(ascii_upper(sequence[index]));
  }
  return result;
}

[[nodiscard]] std::size_t aligned_query_bases(const ReadRecord &read) noexcept {
  std::size_t total = 0U;
  for (const auto &operation : read.cigar_operations) {
    if (operation.code == 'M' || operation.code == '=' || operation.code == 'X') {
      total += operation.length;
    }
  }
  return total;
}

[[nodiscard]] bool usable_record(const ReadRecord &read,
                                 const AnalysisConfig &fragment_config,
                                 const std::size_t length) {
  return !read.cigar_operations.empty() && !read.sequence.empty() &&
         read.reference_start != 0U && read.reference_start <= length &&
         passes_snp_alignment_filters(read, fragment_config);
}

[[nodiscard]] std::vector<MoleculeView>
molecule_views(const MoleculeAssemblyResult &assembly,
               const std::vector<ReadFeature> &features,
               const AnalysisConfig &fragment_config,
               const std::size_t length) {
  std::vector<MoleculeView> views;
  views.reserve(assembly.molecules.size());
  for (std::size_t index = 0U; index < assembly.molecules.size(); ++index) {
    if (index >= features.size() || features[index].filtered_numt ||
        !assembly.molecules[index].analysis_eligible ||
        index >= assembly.representatives.size()) {
      continue;
    }
    MoleculeView view;
    view.index = static_cast<std::uint32_t>(index);
    view.representative = &assembly.representatives[index];
    for (const auto fragment_id : assembly.molecules[index].fragment_ids) {
      const auto &record = source_record_for_fragment(
          assembly, assembly.fragments[fragment_id.value]);
      if (usable_record(record, fragment_config, length)) {
        view.records.push_back(&record);
      }
    }
    if (!view.records.empty()) {
      views.push_back(std::move(view));
    }
  }
  return views;
}

/** Calls body(position, query_index) for each aligned base of one record. */
template <typename Body>
void for_each_aligned_base(const ReadRecord &read, const std::size_t length,
                           Body &&body) {
  std::size_t query = 0U;
  std::size_t cursor = read.reference_start;
  for (const auto &operation : read.cigar_operations) {
    const char op = operation.code;
    if (op == 'M' || op == '=' || op == 'X') {
      for (std::size_t offset = 0U;
           offset < operation.length && query + offset < read.sequence.size();
           ++offset) {
        body(fold(cursor + offset, length), query + offset);
      }
      query += operation.length;
      cursor += operation.length;
    } else if (op == 'I' || op == 'S') {
      query += operation.length;
    } else if (op == 'D' || op == 'N') {
      cursor += operation.length;
    }
  }
}

void add_record_to_pileup(const ReadRecord &read, const std::uint32_t molecule,
                          const std::string &reference,
                          const AnalysisConfig &config, WorkerPileup &worker) {
  const std::size_t length = reference.size();
  const std::size_t strand = (read.flags & 0x10U) != 0U ? 1U : 0U;
  const std::size_t aligned_total = aligned_query_bases(read);
  std::size_t query = 0U;
  std::size_t cursor = read.reference_start;
  std::size_t aligned_so_far = 0U;
  std::size_t interval_start = read.reference_start;
  const auto read_length = read.sequence.size();

  const auto close_interval = [&](const std::size_t end) {
    if (end < interval_start) {
      return;
    }
    if (end - interval_start + 1U >= length) {
      worker.intervals.push_back({1U, static_cast<std::uint32_t>(length), molecule});
      worker.span_delta[strand][1U] += 1;
      worker.span_delta[strand][length + 1U] -= 1;
      return;
    }
    const auto first = fold(interval_start, length);
    const auto last = fold(end, length);
    if (first <= last) {
      worker.intervals.push_back({static_cast<std::uint32_t>(first),
                                  static_cast<std::uint32_t>(last), molecule});
      worker.span_delta[strand][first] += 1;
      worker.span_delta[strand][last + 1U] -= 1;
    } else {
      worker.intervals.push_back({static_cast<std::uint32_t>(first),
                                  static_cast<std::uint32_t>(length), molecule});
      worker.intervals.push_back(
          {1U, static_cast<std::uint32_t>(last), molecule});
      worker.span_delta[strand][first] += 1;
      worker.span_delta[strand][length + 1U] -= 1;
      worker.span_delta[strand][1U] += 1;
      worker.span_delta[strand][last + 1U] -= 1;
    }
  };

  const auto add_indel = [&](const bool deletion, const std::size_t position,
                             const std::size_t size,
                             const std::string_view inserted) {
    const std::uint64_t placement =
        (static_cast<std::uint64_t>(position) << 16U) | static_cast<std::uint64_t>(size);
    if (deletion) {
      const auto known = worker.deletion_entries.find(placement);
      if (known != worker.deletion_entries.end()) {
        if (known->second != nullptr) {
          known->second->strands[strand] += 1U;
          known->second->molecules.push_back(molecule);
        }
        return;
      }
    }
    auto normalized = normalize_circular_indel(reference, position,
                                               deletion ? size : 0U, inserted);
    if (normalized.left_allele.find('N') != std::string::npos) {
      if (deletion) {
        worker.deletion_entries.emplace(placement, nullptr);
      }
      return;
    }
    std::string key;
    key.reserve(normalized.left_allele.size() + 16U);
    key.push_back(deletion ? 'D' : 'I');
    key += std::to_string(normalized.left);
    key.push_back(':');
    key += normalized.left_allele;
    auto [entry, inserted_entry] = worker.indels.try_emplace(std::move(key));
    auto &evidence = entry->second;
    if (inserted_entry) {
      evidence.deletion = deletion;
      evidence.left = normalized.left;
      evidence.right = normalized.right;
      evidence.left_allele = std::move(normalized.left_allele);
      evidence.right_allele = std::move(normalized.right_allele);
    }
    evidence.strands[strand] += 1U;
    evidence.molecules.push_back(molecule);
    if (deletion) {
      worker.deletion_entries.emplace(placement, &evidence);
    }
  };

  const auto is_small_gap = [](const CigarOperation &item) {
    return item.code == 'I' || (item.code == 'D' && item.length <= kMaximumSmallIndel);
  };
  for (std::size_t operation_index = 0U; operation_index < read.cigar_operations.size();
       ++operation_index) {
    const auto &operation = read.cigar_operations[operation_index];
    const char op = operation.code;
    const auto size = operation.length;
    if (op == 'M' || op == '=' || op == 'X') {
      const bool gap_before =
          operation_index > 0U && is_small_gap(read.cigar_operations[operation_index - 1U]);
      const bool gap_after = operation_index + 1U < read.cigar_operations.size() &&
                             is_small_gap(read.cigar_operations[operation_index + 1U]);
      for (std::size_t offset = 0U;
           offset < size && query + offset < read_length; ++offset) {
        const auto position = fold(cursor + offset, length);
        if (worker.stamp[position] == molecule + 1U) {
          continue;
        }
        const auto index = query + offset;
        const auto quality = phred_quality_at(read, index);
        const auto base = base_index(read.sequence[index]);
        if (!quality || *quality < config.min_base_quality || !base) {
          continue;
        }
        worker.stamp[position] = molecule + 1U;
        auto &cell = worker.cells[position];
        cell.counts[*base][strand] += 1U;
        cell.quality[*base][quality_bin(*quality)] += 1U;
        const std::size_t distance = std::min(index, read_length - 1U - index);
        const std::size_t proximity_bin = std::min<std::size_t>(
            kProximityBins - 1U,
            (2U * distance * kProximityBins) / read_length);
        cell.proximity[*base][proximity_bin] += 1U;
        if ((gap_before && offset < kGapWindow) || (gap_after && size - 1U - offset < kGapWindow)) {
          cell.near_gap[*base] += 1U;
        }
      }
      query += size;
      cursor += size;
      aligned_so_far += size;
    } else if (op == 'I') {
      if (size <= kMaximumSmallIndel && aligned_so_far >= kIndelEdgeGuard &&
          aligned_total - aligned_so_far >= kIndelEdgeGuard &&
          query + size <= read_length) {
        // Inserted bases carry low quality because they are often errors. They are
        // counted anyway: this is the noise the background has to describe, and
        // read rescoring sees them too. Only unknown bases are left out.
        std::string inserted;
        inserted.reserve(size);
        bool callable = true;
        for (std::size_t offset = 0U; offset < size; ++offset) {
          const char base = ascii_upper(read.sequence[query + offset]);
          if (!is_base(base)) {
            callable = false;
            break;
          }
          inserted.push_back(base);
        }
        if (callable) {
          add_indel(false, cursor > 1U ? fold(cursor - 1U, length) : length, size,
                    inserted);
        }
      }
      query += size;
    } else if (op == 'S') {
      query += size;
    } else if (op == 'D' || op == 'N') {
      if (op == 'D' && size <= kMaximumSmallIndel &&
          aligned_so_far >= kIndelEdgeGuard &&
          aligned_total - aligned_so_far >= kIndelEdgeGuard) {
        add_indel(true, fold(cursor, length), size, {});
      }
      if (size >= kMinimumLargeDeletion) {
        // A long gap splits the alignment: nothing here spans the deleted arc.
        close_interval(cursor - 1U);
        interval_start = cursor + size;
      }
      cursor += size;
    }
  }
  close_interval(cursor - 1U);
}

[[nodiscard]] Pileup build_pileup(const std::vector<MoleculeView> &molecules,
                                  const std::vector<char> &excluded,
                                  const std::string &reference,
                                  const AnalysisConfig &config) {
  const std::size_t length = reference.size();
  std::vector<WorkerPileup> workers(parallel_worker_limit(config));
  parallel_for_index(
      molecules.size(), config, 16U,
      [&](const std::size_t index, const std::size_t worker_index) {
        if (excluded[index] != 0) {
          return;
        }
        auto &worker = workers[worker_index];
        if (worker.cells.empty()) {
          worker.cells.resize(length + 1U);
          worker.span_delta[0].assign(length + 2U, 0);
          worker.span_delta[1].assign(length + 2U, 0);
          worker.stamp.assign(length + 1U, 0U);
        }
        for (const auto *record : molecules[index].records) {
          add_record_to_pileup(*record, molecules[index].index, reference, config,
                               worker);
        }
      });

  Pileup pileup;
  pileup.cells.resize(length + 1U);
  std::array<std::vector<std::int64_t>, 2> delta;
  delta[0].assign(length + 2U, 0);
  delta[1].assign(length + 2U, 0);
  std::map<std::string, RawIndel> indels;
  for (auto &worker : workers) {
    if (worker.cells.empty()) {
      continue;
    }
    for (std::size_t position = 1U; position <= length; ++position) {
      auto &target = pileup.cells[position];
      const auto &source = worker.cells[position];
      for (std::size_t base = 0U; base < 4U; ++base) {
        target.counts[base][0] += source.counts[base][0];
        target.counts[base][1] += source.counts[base][1];
        for (std::size_t bin = 0U; bin < kQualityBins; ++bin) {
          target.quality[base][bin] += source.quality[base][bin];
        }
        for (std::size_t bin = 0U; bin < kProximityBins; ++bin) {
          target.proximity[base][bin] += source.proximity[base][bin];
        }
        target.near_gap[base] += source.near_gap[base];
      }
    }
    for (std::size_t strand = 0U; strand < 2U; ++strand) {
      for (std::size_t position = 0U; position < length + 2U; ++position) {
        delta[strand][position] += worker.span_delta[strand][position];
      }
    }
    for (auto &[key, evidence] : worker.indels) {
      auto [entry, inserted] = indels.try_emplace(key, std::move(evidence));
      if (!inserted) {
        entry->second.strands[0] += evidence.strands[0];
        entry->second.strands[1] += evidence.strands[1];
        entry->second.molecules.insert(entry->second.molecules.end(),
                                       evidence.molecules.begin(),
                                       evidence.molecules.end());
      }
    }
    pileup.intervals.insert(pileup.intervals.end(), worker.intervals.begin(),
                            worker.intervals.end());
  }
  for (std::size_t strand = 0U; strand < 2U; ++strand) {
    pileup.span[strand].assign(length + 1U, 0U);
    std::int64_t running = 0;
    for (std::size_t position = 1U; position <= length; ++position) {
      running += delta[strand][position];
      pileup.span[strand][position] =
          static_cast<std::uint32_t>(std::max<std::int64_t>(0, running));
    }
  }
  pileup.indels.reserve(indels.size());
  for (auto &[_, evidence] : indels) {
    std::sort(evidence.molecules.begin(), evidence.molecules.end());
    evidence.molecules.erase(
        std::unique(evidence.molecules.begin(), evidence.molecules.end()),
        evidence.molecules.end());
    pileup.indels.push_back(std::move(evidence));
  }
  std::sort(pileup.intervals.begin(), pileup.intervals.end(),
            [](const AlignedInterval &lhs, const AlignedInterval &rhs) {
              return std::tie(lhs.start, lhs.end, lhs.molecule) <
                     std::tie(rhs.start, rhs.end, rhs.molecule);
            });
  return pileup;
}

/** Removes what build_pileup counted for a set of molecules from the total. The
 * contribution of a molecule does not depend on the others, so this equals a
 * rebuild without them, at the cost of those few molecules only. */
void subtract_pileup(Pileup &total, const Pileup &removed, const std::vector<char> &removed_molecule) {
  for (std::size_t position = 0U; position < total.cells.size(); ++position) {
    auto &target = total.cells[position];
    const auto &source = removed.cells[position];
    for (std::size_t base = 0U; base < 4U; ++base) {
      target.counts[base][0] -= source.counts[base][0];
      target.counts[base][1] -= source.counts[base][1];
      for (std::size_t bin = 0U; bin < kQualityBins; ++bin) {
        target.quality[base][bin] -= source.quality[base][bin];
      }
      for (std::size_t bin = 0U; bin < kProximityBins; ++bin) {
        target.proximity[base][bin] -= source.proximity[base][bin];
      }
      target.near_gap[base] -= source.near_gap[base];
    }
  }
  for (std::size_t strand = 0U; strand < 2U; ++strand) {
    for (std::size_t position = 0U; position < total.span[strand].size(); ++position) {
      total.span[strand][position] -= removed.span[strand][position];
    }
  }
  // Both lists are sorted by the same order of (kind, placement, allele) keys.
  const auto key = [](const RawIndel &indel) {
    return std::tie(indel.deletion, indel.left, indel.left_allele);
  };
  std::map<std::tuple<bool, std::size_t, std::string>, const RawIndel *> removed_indels;
  for (const auto &indel : removed.indels) {
    removed_indels.emplace(std::tuple<bool, std::size_t, std::string>(key(indel)), &indel);
  }
  std::erase_if(total.indels, [&](RawIndel &indel) {
    const auto found =
        removed_indels.find(std::tuple<bool, std::size_t, std::string>(key(indel)));
    if (found == removed_indels.end()) {
      return false;
    }
    indel.strands[0] -= found->second->strands[0];
    indel.strands[1] -= found->second->strands[1];
    std::erase_if(indel.molecules, [&](const std::uint32_t molecule) {
      return removed_molecule[molecule] != 0;
    });
    return indel.molecules.empty();
  });
  std::erase_if(total.intervals, [&](const AlignedInterval &interval) {
    return removed_molecule[interval.molecule] != 0;
  });
}

[[nodiscard]] std::string consensus_from_pileup(const Pileup &pileup,
                                                const std::string &reference) {
  std::string consensus = reference;
  for (std::size_t position = 1U; position <= reference.size(); ++position) {
    const auto &cell = pileup.cells[position];
    const auto reference_index = base_index(reference[position - 1U]);
    std::size_t best = reference_index.value_or(0U);
    std::uint64_t best_count = cell.counts[best][0] + cell.counts[best][1];
    for (std::size_t base = 0U; base < 4U; ++base) {
      const std::uint64_t count = cell.counts[base][0] + cell.counts[base][1];
      if (count > best_count) {
        best = base;
        best_count = count;
      }
    }
    if (best_count > 0U) {
      consensus[position - 1U] = "ACGT"[best];
    }
  }
  return consensus;
}

/** Substitution rate of each molecule against the sample consensus. */
[[nodiscard]] std::vector<double>
molecule_divergence(const std::vector<MoleculeView> &molecules,
                    const std::string &consensus, const AnalysisConfig &config) {
  const std::size_t length = consensus.size();
  std::vector<double> divergence(molecules.size(), -1.0);
  parallel_for_index(molecules.size(), config, 32U,
                     [&](const std::size_t index, std::size_t) {
                       std::uint64_t mismatches = 0U;
                       std::uint64_t compared = 0U;
                       for (const auto *record : molecules[index].records) {
                         for_each_aligned_base(
                             *record, length,
                             [&](const std::size_t position, const std::size_t query) {
                               const auto quality = phred_quality_at(*record, query);
                               const char base = ascii_upper(record->sequence[query]);
                               if (!quality || *quality < config.min_base_quality ||
                                   !is_base(base)) {
                                 return;
                               }
                               ++compared;
                               if (base != consensus[position - 1U]) {
                                 ++mismatches;
                               }
                             });
                       }
                       if (compared >= 200U) {
                         divergence[index] = static_cast<double>(mismatches) /
                                             static_cast<double>(compared);
                       }
                     });
  return divergence;
}

[[nodiscard]] double median_of(std::vector<double> values) {
  if (values.empty()) {
    return 0.0;
  }
  const auto middle = values.begin() + static_cast<std::ptrdiff_t>(values.size() / 2U);
  std::nth_element(values.begin(), middle, values.end());
  return *middle;
}

// ------------------------------------------------------------------- SNVs

struct ContextIndex {
  static constexpr std::size_t kClasses = 64U * 3U * 3U;

  [[nodiscard]] static std::optional<std::size_t>
  of(const char previous, const char center, const char next, const char alternate,
     const std::size_t homopolymer) noexcept {
    const auto p = base_index(previous);
    const auto c = base_index(center);
    const auto n = base_index(next);
    const auto a = base_index(alternate);
    if (!p || !c || !n || !a || *a == *c) {
      return std::nullopt;
    }
    const std::size_t slot = *a < *c ? *a : *a - 1U;
    const std::size_t bin = homopolymer <= 3U ? 0U : (homopolymer <= 5U ? 1U : 2U);
    return ((*p * 16U + *c * 4U + *n) * 3U + slot) * 3U + bin;
  }
};

struct SnvHypothesis {
  std::size_t position = 0;
  std::size_t allele = 0;
  std::size_t alternate_forward = 0;
  std::size_t alternate_reverse = 0;
  std::size_t depth_forward = 0;
  std::size_t depth_reverse = 0;
  double p_value = 1.0;
  /** The weaker strand's own significance; 1 when a strand lacks depth. */
  double weaker_strand_p = 1.0;
  double background = 0.0;
  double rho = 0.0;
};

[[nodiscard]] std::size_t context_class_for(const std::string &consensus,
                                            const RepeatContext &repeats,
                                            const std::size_t position,
                                            const std::size_t allele,
                                            const std::size_t strand) {
  const std::size_t length = consensus.size();
  const char previous = consensus[(position + length - 2U) % length];
  const char center = consensus[position - 1U];
  const char next = consensus[position % length];
  const std::size_t run = std::max({repeats.homopolymer_run[fold(position + length - 1U, length)],
                                    repeats.homopolymer_run[position],
                                    repeats.homopolymer_run[fold(position + 1U, length)]});
  const char alternate = "ACGT"[allele];
  const auto index =
      strand == 0U
          ? ContextIndex::of(previous, center, next, alternate, run)
          : ContextIndex::of(complement_base(next), complement_base(center),
                             complement_base(previous),
                             complement_base(alternate), run);
  return index.value_or(ContextIndex::kClasses);
}

[[nodiscard]] PrecisionRepeat repeat_at(const RepeatContext &repeats,
                                        const std::size_t position) {
  return {repeats.homopolymer_run[position], repeats.unit[position],
          repeats.span[position]};
}

struct StrandBackground {
  std::array<BetaBinomialFit, 2> strands;
};

[[nodiscard]] std::vector<PrecisionSnv>
call_snvs(const Pileup &pileup, const std::string &reference,
          const std::string &consensus, const RepeatContext &repeats,
          PrecisionDiagnostics &diagnostics) {
  const std::size_t length = reference.size();
  std::vector<std::vector<SiteCount>> class_sites(ContextIndex::kClasses + 1U);
  std::vector<SiteCount> all_sites;
  all_sites.reserve(length * 6U);
  std::uint64_t total_minor = 0U;
  std::uint64_t total_depth = 0U;
  for (std::size_t position = 1U; position <= length; ++position) {
    const auto major = base_index(consensus[position - 1U]);
    if (!major) {
      continue;
    }
    const auto &cell = pileup.cells[position];
    for (std::size_t strand = 0U; strand < 2U; ++strand) {
      std::uint32_t depth = 0U;
      for (std::size_t base = 0U; base < 4U; ++base) {
        depth += cell.counts[base][strand];
      }
      if (depth == 0U) {
        continue;
      }
      total_depth += depth;
      for (std::size_t base = 0U; base < 4U; ++base) {
        if (base == *major) {
          continue;
        }
        total_minor += cell.counts[base][strand];
        const SiteCount site{cell.counts[base][strand], depth};
        class_sites[context_class_for(consensus, repeats, position, base, strand)]
            .push_back(site);
        all_sites.push_back(site);
      }
    }
  }
  diagnostics.mean_substitution_rate =
      total_depth == 0U ? 0.0
                        : static_cast<double>(total_minor) /
                              (3.0 * static_cast<double>(total_depth));

  constexpr double kTrim = 0.01;
  constexpr double kRhoFloor = 2e-4;
  constexpr double kMeanFloor = 2e-5;
  const auto global_fit = fit_beta_binomial(all_sites, kTrim, kMeanFloor, kRhoFloor);
  std::vector<BetaBinomialFit> fits(ContextIndex::kClasses + 1U, global_fit);
  // Classes with few sites shrink to the genome-wide rate instead of trusting
  // a rate estimated from a handful of positions.
  constexpr double kShrinkageSites = 30.0;
  for (std::size_t index = 0U; index < class_sites.size(); ++index) {
    if (class_sites[index].empty()) {
      continue;
    }
    auto fit = fit_beta_binomial(std::move(class_sites[index]), kTrim, kMeanFloor,
                                 kRhoFloor);
    const double weight = static_cast<double>(fit.sites) /
                          (static_cast<double>(fit.sites) + kShrinkageSites);
    fit.mean = weight * fit.mean + (1.0 - weight) * global_fit.mean;
    fit.rho = std::max(weight * fit.rho + (1.0 - weight) * global_fit.rho, kRhoFloor);
    fits[index] = fit;
  }

  std::vector<SnvHypothesis> hypotheses;
  std::size_t covered_positions = 0U;
  for (std::size_t position = 1U; position <= length; ++position) {
    const auto &cell = pileup.cells[position];
    const auto major = base_index(consensus[position - 1U]);
    const auto reference_base = base_index(reference[position - 1U]);
    const std::size_t forward = cell.counts[0][0] + cell.counts[1][0] +
                                cell.counts[2][0] + cell.counts[3][0];
    const std::size_t reverse = cell.counts[0][1] + cell.counts[1][1] +
                                cell.counts[2][1] + cell.counts[3][1];
    if (forward + reverse == 0U || !major) {
      continue;
    }
    ++covered_positions;
    // A placeholder reference base (rCRS 3107) cannot have an alternate allele.
    if (!reference_base) {
      continue;
    }
    for (std::size_t base = 0U; base < 4U; ++base) {
      if (base == *reference_base) {
        continue;
      }
      const std::size_t alternate_forward = cell.counts[base][0];
      const std::size_t alternate_reverse = cell.counts[base][1];
      if (alternate_forward + alternate_reverse == 0U) {
        continue;
      }
      SnvHypothesis hypothesis;
      hypothesis.position = position;
      hypothesis.allele = base;
      hypothesis.alternate_forward = alternate_forward;
      hypothesis.alternate_reverse = alternate_reverse;
      hypothesis.depth_forward = forward;
      hypothesis.depth_reverse = reverse;
      // The major allele has no error class of its own; the genome-wide rate
      // is enough, since a true major allele dwarfs any error model.
      const auto &forward_fit =
          base == *major
              ? global_fit
              : fits[context_class_for(consensus, repeats, position, base, 0U)];
      const auto &reverse_fit =
          base == *major
              ? global_fit
              : fits[context_class_for(consensus, repeats, position, base, 1U)];
      const double p_forward = beta_binomial_survival(
          alternate_forward, forward, forward_fit.mean, forward_fit.rho);
      const double p_reverse = beta_binomial_survival(
          alternate_reverse, reverse, reverse_fit.mean, reverse_fit.rho);
      const bool both_strands = forward >= kMinimumStrandDepth &&
                                reverse >= kMinimumStrandDepth;
      if (both_strands) {
        hypothesis.p_value = fisher_combined(p_forward, p_reverse);
        hypothesis.weaker_strand_p = std::max(p_forward, p_reverse);
      } else {
        hypothesis.weaker_strand_p = 0.0;
        hypothesis.p_value = forward >= reverse ? p_forward : p_reverse;
      }
      const double total = static_cast<double>(forward + reverse);
      hypothesis.background = (forward_fit.mean * static_cast<double>(forward) +
                               reverse_fit.mean * static_cast<double>(reverse)) /
                              total;
      hypothesis.rho = (forward_fit.rho * static_cast<double>(forward) +
                        reverse_fit.rho * static_cast<double>(reverse)) /
                       total;
      hypotheses.push_back(hypothesis);
    }
  }
  diagnostics.reference_positions_covered = covered_positions;
  diagnostics.snv_hypotheses = covered_positions * 3U;

  std::vector<double> p_values;
  p_values.reserve(hypotheses.size());
  for (const auto &hypothesis : hypotheses) {
    p_values.push_back(hypothesis.p_value);
  }
  const auto q_values = benjamini_hochberg(p_values, diagnostics.snv_hypotheses);

  std::vector<PrecisionSnv> calls;
  for (std::size_t index = 0U; index < hypotheses.size(); ++index) {
    const auto &hypothesis = hypotheses[index];
    if (hypothesis.alternate_forward + hypothesis.alternate_reverse < 2U) {
      continue;
    }
    const std::size_t depth = hypothesis.depth_forward + hypothesis.depth_reverse;
    const std::size_t alternate_depth =
        hypothesis.alternate_forward + hypothesis.alternate_reverse;
    const double raw = static_cast<double>(alternate_depth) / static_cast<double>(depth);
    PrecisionSnv call;
    call.position = hypothesis.position;
    call.reference = reference[hypothesis.position - 1U];
    call.alternate = "ACGT"[hypothesis.allele];
    call.depth = depth;
    call.alternate_depth = alternate_depth;
    call.strands = {hypothesis.alternate_forward, hypothesis.alternate_reverse,
                    hypothesis.depth_forward, hypothesis.depth_reverse};
    call.p_value = hypothesis.p_value;
    call.q_value = q_values[index];
    call.background_rate = hypothesis.background;
    call.repeat = repeat_at(repeats, hypothesis.position);
    // The observed rate is f (1 - 3 e) + (1 - f) e, with e the background.
    const double correction = std::max(0.25, 1.0 - 4.0 * hypothesis.background);
    const auto corrected = [&](const double value) {
      return std::clamp((value - hypothesis.background) / correction, 0.0, 1.0);
    };
    call.fraction = corrected(raw);
    const auto interval = overdispersed_wilson(static_cast<double>(alternate_depth),
                                               static_cast<double>(depth), hypothesis.rho);
    call.ci95_low = std::min(corrected(interval.first), call.fraction);
    call.ci95_high = std::max(corrected(interval.second), call.fraction);

    const auto &cell = pileup.cells[hypothesis.position];
    double weighted_quality = 0.0;
    double quality_count = 0.0;
    for (std::size_t bin = 0U; bin < kQualityBins; ++bin) {
      weighted_quality += (static_cast<double>(bin) * 5.0 + 2.5) *
                          cell.quality[hypothesis.allele][bin];
      quality_count += cell.quality[hypothesis.allele][bin];
    }
    call.alternate_mean_quality =
        quality_count > 0.0 ? weighted_quality / quality_count : 0.0;

    bool pass = call.q_value <= kFalseDiscoveryRate;
    if (raw < kBackgroundEffectMultiple * hypothesis.background) {
      call.filters.emplace_back("below_error_floor");
    }
    if (raw < 0.5 && alternate_depth < kMinimumLowLevelSupport) {
      call.filters.emplace_back("low_support");
    }
    if (hypothesis.depth_forward >= kMinimumStrandDepth &&
        hypothesis.depth_reverse >= kMinimumStrandDepth) {
      const double rate_forward = static_cast<double>(hypothesis.alternate_forward) /
                                  static_cast<double>(hypothesis.depth_forward);
      const double rate_reverse = static_cast<double>(hypothesis.alternate_reverse) /
                                  static_cast<double>(hypothesis.depth_reverse);
      const double strand_p = fisher_exact_two_sided(
          hypothesis.alternate_forward,
          hypothesis.depth_forward - hypothesis.alternate_forward,
          hypothesis.alternate_reverse,
          hypothesis.depth_reverse - hypothesis.alternate_reverse);
      if (strand_p < 1e-3 &&
          std::min(rate_forward, rate_reverse) <
              0.25 * std::max(rate_forward, rate_reverse)) {
        call.filters.emplace_back("strand_bias");
      }
      // Sequencer errors depend on the strand read, so a real allele has to
      // stand out on both (the AND model of Shearwater). Combining the strands
      // alone would let one strong artifact carry a silent opposite strand.
      if (hypothesis.weaker_strand_p > kWeakerStrandLimit && raw < 0.2) {
        call.filters.emplace_back("one_strand_support");
      }
    } else {
      call.flags.emplace_back("single_strand_depth");
    }
    const auto major = base_index(consensus[hypothesis.position - 1U]);
    if (major && *major != hypothesis.allele && alternate_depth >= 5U) {
      double quality_shift = 0.0;
      const double quality_p =
          rank_sum_p_value(cell.quality[hypothesis.allele].data(),
                           cell.quality[*major].data(), kQualityBins, &quality_shift);
      // True alleles read like the rest of the site, and errors read worse. At low
      // levels there are few alternate bases, so the test has less power and needs
      // a looser threshold, at the cost of a small loss of real low-level alleles.
      const double quality_limit = raw < 0.05 ? 0.02 : 1e-3;
      if (quality_p < quality_limit && quality_shift <= -0.6) {
        call.filters.emplace_back("low_alternate_base_quality");
      }
      double proximity_shift = 0.0;
      const double proximity_p = rank_sum_p_value(
          cell.proximity[hypothesis.allele].data(), cell.proximity[*major].data(),
          kProximityBins, &proximity_shift);
      if (proximity_p < 1e-4 && proximity_shift <= -0.75) {
        call.filters.emplace_back("read_end_bias");
      }
      // A length error is often written as a gap plus a mismatch beside it, so a
      // "variant" seen mostly next to gaps is the aligner's spelling of an indel.
      const double alternate_near_gap =
          static_cast<double>(cell.near_gap[hypothesis.allele]) / static_cast<double>(alternate_depth);
      const double major_total = cell.counts[*major][0] + cell.counts[*major][1];
      const double major_near_gap =
          major_total > 0.0 ? static_cast<double>(cell.near_gap[*major]) / major_total : 0.0;
      if (raw < kLowLevelFraction && alternate_near_gap >= 0.4 &&
          alternate_near_gap >= 3.0 * major_near_gap) {
        call.filters.emplace_back("gap_adjacent");
      }
    }
    if (call.repeat.span >= 8U) {
      call.flags.emplace_back("repeat_context");
    }
    pass = pass && call.filters.empty();
    if (!pass && call.q_value > kFalseDiscoveryRate) {
      call.filters.emplace_back("not_significant");
    }
    call.pass = pass;
    calls.push_back(std::move(call));
  }
  return calls;
}

// ------------------------------------------------------------- small indels

struct IndelClass {
  bool deletion = false;
  std::size_t length_bin = 1;
  std::size_t repeat_bin = 1;
  std::size_t strand = 0;

  auto operator<=>(const IndelClass &) const = default;
};

[[nodiscard]] bool allele_is_unit_of(const std::string_view allele,
                                     const char base) noexcept {
  return std::all_of(allele.begin(), allele.end(),
                     [&](const char value) { return value == base; });
}

/** In-repeat events are binned by run length; everything else is bin 1. */
[[nodiscard]] std::size_t indel_repeat_bin(const RawIndel &indel,
                                           const std::string &consensus,
                                           const RepeatContext &repeats) {
  const std::size_t length = consensus.size();
  const std::size_t site =
      indel.deletion ? indel.left : fold(indel.left + 1U, length);
  const std::size_t run = repeats.homopolymer_run[site];
  if (run >= 2U && allele_is_unit_of(indel.left_allele, consensus[site - 1U])) {
    return std::min<std::size_t>(run, 12U);
  }
  // Placements are normalized on rCRS and context comes from this sample, so a
  // fixed variant can leave the allele's run one base to the side of its placement.
  if (allele_is_unit_of(indel.left_allele, indel.left_allele.front())) {
    for (const std::size_t neighbor : {fold(site + length - 1U, length), fold(site + 1U, length)}) {
      if (consensus[neighbor - 1U] == indel.left_allele.front() &&
          repeats.homopolymer_run[neighbor] >= 2U) {
        return std::min<std::size_t>(repeats.homopolymer_run[neighbor], 12U);
      }
    }
  }
  // Beside a long run, a slip often aligns as a gap on the neighboring base, so
  // such an event has that run's noise whatever its allele looks like.
  const std::size_t event_end = fold(site + indel.left_allele.size() - 1U, length);
  const std::size_t neighbor_run =
      std::max({static_cast<std::size_t>(repeats.homopolymer_run[fold(site + length - 1U, length)]),
                static_cast<std::size_t>(repeats.homopolymer_run[fold(event_end + 1U, length)]),
                run});
  constexpr std::size_t kInfluentialRun = 4U;
  if (neighbor_run >= kInfluentialRun) {
    return std::min<std::size_t>(neighbor_run, 12U);
  }
  const std::size_t unit = repeats.unit[site];
  const std::size_t span = repeats.span[site];
  if (unit >= 2U && span >= 2U * unit && indel.left_allele.size() % unit == 0U) {
    // Dinucleotide and longer repeats behave like a run of their copy number.
    return std::min<std::size_t>(span / unit + 1U, 12U);
  }
  return 1U;
}

/** Groups indel rates for the empirical tail: type, length and run length. */
struct TailGroup {
  bool deletion = false;
  std::size_t length_bin = 1;
  std::size_t run_group = 0;

  auto operator<=>(const TailGroup &) const = default;
};

struct IndelModel {
  std::map<IndelClass, BetaBinomialFit> fits;
  std::map<TailGroup, EmpiricalTail> tails;
  PairHmmParameters hmm;

  [[nodiscard]] const EmpiricalTail &tail_for(const bool deletion, const std::size_t length_bin,
                                              const std::size_t repeat_bin) const {
    const auto own = tails.find({deletion, length_bin, std::min<std::size_t>(repeat_bin, 7U)});
    if (own != tails.end()) {
      return own->second;
    }
    return tails.at({deletion, length_bin, 0U});
  }
};

/** Per-repeat-length error model, measured from this sample's own reads. */
[[nodiscard]] IndelModel fit_indel_model(const std::vector<RawIndel> &indels,
                                         const Pileup &pileup,
                                         const std::string &consensus,
                                         const RepeatContext &repeats,
                                         PrecisionDiagnostics &diagnostics) {
  const std::size_t length = consensus.size();
  std::map<IndelClass, std::vector<SiteCount>> sites;
  std::map<IndelClass, std::size_t> loci;

  // Every homopolymer run is one locus for a one-base length change.
  double depth_sum = 0.0;
  std::size_t positions = 0U;
  for (std::size_t position = 1U; position <= length; ++position) {
    depth_sum += pileup.span[0][position] + pileup.span[1][position];
    ++positions;
    const char base = consensus[position - 1U];
    const char before = consensus[(position + length - 2U) % length];
    if (base == before) {
      continue;
    }
    const std::size_t run = std::min<std::size_t>(repeats.homopolymer_run[position], 12U);
    for (std::size_t strand = 0U; strand < 2U; ++strand) {
      ++loci[{true, 1U, run, strand}];
      if (run >= 2U) {
        ++loci[{false, 1U, run, strand}];
      } else {
        loci[{false, 1U, 1U, strand}] += 3U;
      }
      for (std::size_t bin = 2U; bin <= 3U; ++bin) {
        ++loci[{true, bin, 1U, strand}];
        ++loci[{false, bin, 1U, strand}];
      }
    }
  }
  const double mean_strand_depth = positions == 0U
                                       ? 0.0
                                       : depth_sum / (2.0 * static_cast<double>(positions));

  std::vector<std::array<IndelClass, 2>> indel_classes;
  indel_classes.reserve(indels.size());
  for (const auto &indel : indels) {
    const std::size_t length_bin = std::min<std::size_t>(indel.left_allele.size(), 3U);
    const std::size_t repeat_bin = indel_repeat_bin(indel, consensus, repeats);
    for (std::size_t strand = 0U; strand < 2U; ++strand) {
      const IndelClass key{indel.deletion, length_bin, repeat_bin, strand};
      const auto depth = pileup.span[strand][indel.left];
      const auto count = std::min(indel.strands[strand], depth);
      // A majority allele is a fixed variant of this sample, not an error.
      if (depth > 0U && count * 2U < depth) {
        sites[key].push_back({count, depth});
      }
    }
  }
  IndelModel model;
  const auto zero_depth = static_cast<std::uint32_t>(std::max(1.0, mean_strand_depth));
  constexpr double kTrim = 0.01;
  constexpr double kMeanFloor = 1e-5;
  constexpr double kRhoFloor = 5e-4;
  constexpr std::size_t kMinimumClassSites = 20U;

  // Every class also has a parent over all contexts of its type and length. A
  // class with few sites takes the parent's rate: with only a handful of sites,
  // a true variant would be most of the data and would explain itself as error.
  std::map<IndelClass, std::vector<SiteCount>> parent_sites;
  std::map<IndelClass, BetaBinomialFit> parents;
  for (const auto &[key, observed] : sites) {
    auto &parent = parent_sites[{key.deletion, key.length_bin, 0U, key.strand}];
    parent.insert(parent.end(), observed.begin(), observed.end());
  }
  for (auto &[key, observed] : parent_sites) {
    const std::size_t alleles = key.deletion ? 1U : (key.length_bin == 1U ? 4U : 2U);
    const std::size_t all_sites = positions * alleles;
    if (all_sites > observed.size()) {
      observed.insert(observed.end(), all_sites - observed.size(), SiteCount{0U, zero_depth});
    }
    parents[key] = fit_beta_binomial(std::move(observed), kTrim, kMeanFloor, kRhoFloor);
  }
  for (auto &[key, observed] : sites) {
    const auto known_loci = loci.find(key);
    const std::size_t zero_sites =
        known_loci != loci.end() && known_loci->second > observed.size()
            ? known_loci->second - observed.size()
            : 0U;
    observed.insert(observed.end(), zero_sites, SiteCount{0U, zero_depth});
    const auto parent = parents.find({key.deletion, key.length_bin, 0U, key.strand});
    // One-base slips stay on their own class: their rate depends strongly on run
    // length, and the trend fitted below covers the runs with few sites.
    if (key.length_bin >= 2U && observed.size() < kMinimumClassSites &&
        parent != parents.end()) {
      auto fit = parent->second;
      fit.sites = observed.size();
      model.fits[key] = fit;
    } else {
      model.fits[key] = fit_beta_binomial(std::move(observed), kTrim, kMeanFloor, kRhoFloor);
    }
  }

  // The empirical noise tail, per run length and, for groups too small to fit,
  // per event type and length across all contexts.
  {
    constexpr std::size_t kMinimumTailSites = 200U;
    const double detectable =
        1.0 / std::max(2.0, 2.0 * mean_strand_depth);
    std::map<TailGroup, std::vector<double>> group_rates;
    std::map<TailGroup, std::size_t> group_loci;
    for (const auto &[key, count] : loci) {
      // loci counts each strand separately; rates below combine both.
      if (key.strand == 0U) {
        group_loci[{key.deletion, key.length_bin, std::min<std::size_t>(key.repeat_bin, 7U)}] += count;
      }
    }
    for (const auto &indel : indels) {
      const auto depth = pileup.span[0][indel.left] + pileup.span[1][indel.left];
      if (depth == 0U) {
        continue;
      }
      const double rate = static_cast<double>(indel.molecules.size()) / static_cast<double>(depth);
      const std::size_t length_bin = std::min<std::size_t>(indel.left_allele.size(), 3U);
      const std::size_t repeat_bin = indel_repeat_bin(indel, consensus, repeats);
      group_rates[{indel.deletion, length_bin, std::min<std::size_t>(repeat_bin, 7U)}].push_back(rate);
      group_rates[{indel.deletion, length_bin, 0U}].push_back(rate);
    }
    for (const bool deletion : {true, false}) {
      for (std::size_t length_bin = 1U; length_bin <= 3U; ++length_bin) {
        const std::size_t alleles = deletion ? 1U : (length_bin == 1U ? 4U : 2U);
        group_loci[{deletion, length_bin, 0U}] = positions * alleles;
        group_rates.try_emplace({deletion, length_bin, 0U});
      }
    }
    for (auto &[group, rates] : group_rates) {
      const auto known = group_loci.find(group);
      const std::size_t total = known == group_loci.end() ? rates.size() : known->second;
      if (group.run_group != 0U && std::max(total, rates.size()) < kMinimumTailSites) {
        continue;
      }
      const std::size_t zero_sites = total > rates.size() ? total - rates.size() : 0U;
      model.tails[group] = fit_empirical_tail(std::move(rates), zero_sites, detectable);
    }
  }

  std::map<IndelClass, std::size_t> sites_before_trend;
  for (const auto &[key, count] : loci) {
    sites_before_trend[key] = count;
  }

  // Long runs are rare: a log-linear trend over run length keeps their error
  // estimate from falling below what shorter runs already establish.
  for (const bool deletion : {true, false}) {
    for (std::size_t strand = 0U; strand < 2U; ++strand) {
      double weight_sum = 0.0;
      double x_sum = 0.0;
      double y_sum = 0.0;
      double xx_sum = 0.0;
      double xy_sum = 0.0;
      for (std::size_t run = 2U; run <= 12U; ++run) {
        const auto fit = model.fits.find({deletion, 1U, run, strand});
        if (fit == model.fits.end() || fit->second.sites < 2U ||
            fit->second.mean <= kMeanFloor) {
          continue;
        }
        const double x = static_cast<double>(run);
        const double y = std::log(fit->second.mean);
        const double weight = static_cast<double>(fit->second.sites);
        weight_sum += weight;
        x_sum += weight * x;
        y_sum += weight * y;
        xx_sum += weight * x * x;
        xy_sum += weight * x * y;
      }
      const double denominator = weight_sum * xx_sum - x_sum * x_sum;
      if (weight_sum <= 0.0 || denominator <= 1e-9) {
        continue;
      }
      const double slope =
          std::max(0.0, (weight_sum * xy_sum - x_sum * y_sum) / denominator);
      const double intercept = (y_sum - slope * x_sum) / weight_sum;
      double shorter_run_mean = 0.0;
      for (std::size_t run = 2U; run <= 12U; ++run) {
        const double trend =
            std::min(0.8, std::exp(intercept + slope * static_cast<double>(run)));
        auto &fit = model.fits[{deletion, 1U, run, strand}];
        if (fit.sites < 8U) {
          fit.mean = std::max(fit.mean, trend);
          fit.rho = std::max(fit.rho, kRhoFloor);
        }
        // A longer run does not slip less. Dips come from few sites, and from
        // setting aside the run's own outliers, so the rate is kept monotone.
        fit.mean = std::max(fit.mean, shorter_run_mean);
        shorter_run_mean = fit.mean;
        // Slips of two or more bases are about ten times rarer per extra base,
        // but they still grow with run length, unlike the parent class's rate.
        for (std::size_t length_bin = 2U; length_bin <= 3U; ++length_bin) {
          const auto multi = model.fits.find({deletion, length_bin, run, strand});
          if (multi != model.fits.end()) {
            multi->second.mean = std::max(
                multi->second.mean,
                trend * std::pow(0.1, static_cast<double>(length_bin - 1U)));
          }
        }
      }
    }
  }

  // Gap-open rates per position: a run offers run-many places to open a gap.
  const auto strand_mean = [&](const bool deletion, const std::size_t run) {
    const auto lookup = [&](const std::size_t strand) {
      const auto fit = model.fits.find(
          {deletion, 1U, std::min<std::size_t>(run, 12U), strand});
      return fit == model.fits.end() ? -1.0 : fit->second.mean;
    };
    const double forward = lookup(0U);
    const double reverse = lookup(1U);
    if (forward < 0.0 && reverse < 0.0) {
      return -1.0;
    }
    return forward < 0.0 ? reverse : (reverse < 0.0 ? forward : (forward + reverse) / 2.0);
  };
  // Runs of a length rCRS barely has are extrapolated; reports show measurements.
  constexpr std::size_t kMinimumMeasuredSites = 3U;
  const auto measured_sites = [&](const std::size_t run) {
    const auto fit = sites_before_trend.find({true, 1U, std::min<std::size_t>(run, 12U), 0U});
    return fit == sites_before_trend.end() ? std::size_t{0} : fit->second;
  };
  double last_deletion = 1e-3;
  double last_insertion = 1e-3;
  for (std::size_t run = 1U; run <= PairHmmParameters::kMaximumRun; ++run) {
    const double deletion_rate = strand_mean(true, run);
    const double insertion_rate = strand_mean(false, run);
    last_deletion = deletion_rate > 0.0 ? deletion_rate : last_deletion;
    last_insertion = insertion_rate > 0.0 ? insertion_rate : last_insertion;
    const double per_position = 1.0 / static_cast<double>(run);
    model.hmm.deletion_open[run] =
        std::clamp(last_deletion * per_position, 2e-5, 0.25);
    model.hmm.insertion_open[run] =
        std::clamp(last_insertion * per_position, 2e-5, 0.25);
    const bool measured_run = run <= 12U && measured_sites(run) >= kMinimumMeasuredSites;
    if (measured_run) {
      if (deletion_rate > 0.0) {
        diagnostics.homopolymer_deletion_rate.emplace_back(run, deletion_rate);
      }
      if (insertion_rate > 0.0) {
        diagnostics.homopolymer_insertion_rate.emplace_back(run, insertion_rate);
      }
    }
    if (run < diagnostics.slip_rate_by_run.size()) {
      diagnostics.slip_rate_by_run[run] =
          std::max(deletion_rate, 0.0) + std::max(insertion_rate, 0.0);
    }
  }
  model.hmm.deletion_open[0] = model.hmm.deletion_open[1];
  model.hmm.insertion_open[0] = model.hmm.insertion_open[1];
  return model;
}

struct WindowBounds {
  std::size_t query_start = 0;
  std::size_t query_end = 0; // Exclusive.
};

/** Where a read's alignment stands at every few CIGAR operations. A noisy long
 * read has thousands of operations, and each candidate looks up each read. */
struct CigarCheckpoint {
  std::uint32_t reference_cursor = 0;
  std::uint32_t query = 0;
  std::uint32_t operation = 0;
};

constexpr std::size_t kCigarCheckpointStride = 32U;

[[nodiscard]] std::vector<CigarCheckpoint> cigar_checkpoints(const ReadRecord &read) {
  std::vector<CigarCheckpoint> checkpoints;
  checkpoints.reserve(read.cigar_operations.size() / kCigarCheckpointStride + 1U);
  std::size_t cursor = read.reference_start;
  std::size_t query = 0U;
  for (std::size_t index = 0U; index < read.cigar_operations.size(); ++index) {
    if (index % kCigarCheckpointStride == 0U) {
      checkpoints.push_back({static_cast<std::uint32_t>(cursor), static_cast<std::uint32_t>(query),
                             static_cast<std::uint32_t>(index)});
    }
    const auto &operation = read.cigar_operations[index];
    const char op = operation.code;
    if (op == 'M' || op == '=' || op == 'X') {
      cursor += operation.length;
      query += operation.length;
    } else if (op == 'I' || op == 'S') {
      query += operation.length;
    } else if (op == 'D' || op == 'N') {
      cursor += operation.length;
    }
  }
  return checkpoints;
}

/** Query interval for reference window [first, last]; nullopt if not spanned. */
[[nodiscard]] std::optional<WindowBounds>
read_window(const ReadRecord &read, const std::vector<CigarCheckpoint> &checkpoints,
            const std::size_t first, const std::size_t last) {
  // The walk starts at the last checkpoint at or before the window; operations
  // before it cannot contain the window's first base.
  const auto after = std::upper_bound(
      checkpoints.begin(), checkpoints.end(), first,
      [](const std::size_t value, const CigarCheckpoint &item) { return value < item.reference_cursor; });
  const CigarCheckpoint origin =
      after == checkpoints.begin() ? CigarCheckpoint{} : *(after - 1);
  std::size_t cursor = after == checkpoints.begin() ? read.reference_start : origin.reference_cursor;
  std::size_t query = origin.query;
  std::optional<std::size_t> start;
  for (std::size_t index = origin.operation; index < read.cigar_operations.size(); ++index) {
    const auto &operation = read.cigar_operations[index];
    const auto size = operation.length;
    const char op = operation.code;
    if (op == 'M' || op == '=' || op == 'X') {
      if (!start && first >= cursor && first < cursor + size) {
        start = query + (first - cursor);
      }
      if (start && last >= cursor && last < cursor + size) {
        return WindowBounds{*start, query + (last - cursor) + 1U};
      }
      query += size;
      cursor += size;
    } else if (op == 'I' || op == 'S') {
      query += size;
    } else if (op == 'D' || op == 'N') {
      if (size >= kMinimumLargeDeletion && start) {
        return std::nullopt;
      }
      if (!start && first >= cursor && first < cursor + size) {
        start = query;
      }
      if (start && last >= cursor && last < cursor + size) {
        return WindowBounds{*start, query};
      }
      cursor += size;
    }
  }
  return std::nullopt;
}

[[nodiscard]] std::size_t reference_end_of(const ReadRecord &read) noexcept {
  return read.reference_start + reference_consuming_length(read) - 1U;
}

struct IndelRead {
  const ReadRecord *record = nullptr;
  std::uint32_t molecule = 0;
  std::size_t reference_end = 0;
  std::vector<CigarCheckpoint> checkpoints;
};

struct RescoreResult {
  bool rescored = false;
  double fraction = 0.0;
  double lod = 0.0;
  double effective_alternate = 0.0;
  double effective_depth = 0.0;
  std::array<double, 2> strand_alternate{};
  std::array<double, 2> strand_depth{};
};

[[nodiscard]] std::string apply_indel_to_window(std::string window,
                                                const std::size_t window_first,
                                                const RawIndel &indel) {
  if (indel.deletion) {
    window.erase(indel.left - window_first, indel.left_allele.size());
  } else {
    window.insert(indel.left - window_first + 1U, indel.left_allele);
  }
  return window;
}

[[nodiscard]] RescoreResult
rescore_indel(const RawIndel &indel, const std::vector<const RawIndel *> &competitors,
              const std::vector<IndelRead> &reads, const std::string &consensus,
              const PairHmmParameters &hmm, const std::size_t flank,
              PairHmmScratch &scratch, WorkBudget &budget) {
  RescoreResult result;
  const std::size_t length = consensus.size();
  const std::size_t allele_length = indel.left_allele.size();
  const std::size_t event_end =
      indel.deletion ? indel.right + allele_length - 1U : indel.right + 1U;
  if (indel.left <= flank + 1U || event_end + flank >= length || indel.right < indel.left) {
    return result;
  }
  const std::size_t window_first = indel.left - flank;
  const std::size_t window_last = event_end + flank;

  // Haplotypes use the sample consensus, so fixed variants do not penalize
  // either allele. The rCRS placeholder base never occurs in a molecule.
  std::string reference_haplotype =
      consensus.substr(window_first - 1U, window_last - window_first + 1U);
  std::string alternate_haplotype =
      apply_indel_to_window(reference_haplotype, window_first, indel);
  // Other common alleles at this locus compete with the reference allele. Without
  // them a read with one extra base fits "two extra" better than "none", and a
  // true +1 allele would be credited to a +2 hypothesis.
  std::vector<std::string> competing_haplotypes;
  for (const auto *competitor : competitors) {
    if (competitor->left < window_first ||
        competitor->left + competitor->left_allele.size() > window_last) {
      continue;
    }
    competing_haplotypes.push_back(
        apply_indel_to_window(reference_haplotype, window_first, *competitor));
    std::erase(competing_haplotypes.back(), 'N');
  }
  std::erase(reference_haplotype, 'N');
  std::erase(alternate_haplotype, 'N');

  std::vector<const IndelRead *> carriers;
  std::vector<const IndelRead *> others;
  for (const auto &read : reads) {
    if (read.record->reference_start > window_first || read.reference_end < window_last) {
      continue;
    }
    if (std::binary_search(indel.molecules.begin(), indel.molecules.end(),
                           read.molecule)) {
      carriers.push_back(&read);
    } else {
      others.push_back(&read);
    }
  }
  if (carriers.empty()) {
    return result;
  }
  // Every carrier is scored; the rest are a systematic sample, reweighted.
  const std::size_t other_stride =
      others.size() > kMaximumHmmOtherReads ? others.size() / kMaximumHmmOtherReads : 1U;
  std::size_t others_taken = 0U;
  for (std::size_t index = 0U; index < others.size(); index += other_stride) {
    ++others_taken;
  }
  const double other_weight =
      others_taken == 0U ? 1.0
                         : static_cast<double>(others.size()) /
                               static_cast<double>(others_taken);

  std::vector<double> log_reference;
  std::vector<double> log_alternate;
  std::vector<double> weights;
  std::vector<std::uint8_t> strands;
  const std::size_t band = allele_length + 10U;
  const auto score = [&](const IndelRead &read, const double weight) {
    const auto window = read_window(*read.record, read.checkpoints, window_first, window_last);
    if (!window || window->query_end <= window->query_start ||
        window->query_end > read.record->sequence.size()) {
      return;
    }
    const auto span = window->query_end - window->query_start;
    // Banded DP cells for the reference, competing and alternate haplotypes.
    budget.charge(static_cast<std::uint64_t>(span) * (2U * band + 1U) *
                  (2U + competing_haplotypes.size()));
    const auto bases =
        std::string_view(read.record->sequence).substr(window->query_start, span);
    const auto qualities =
        window->query_end <= read.record->qualities.size()
            ? std::string_view(read.record->qualities).substr(window->query_start, span)
            : std::string_view{};
    double best_other = pair_hmm_log10_likelihood(bases, qualities, reference_haplotype,
                                                  hmm, band, scratch);
    for (const auto &haplotype : competing_haplotypes) {
      best_other = std::max(best_other, pair_hmm_log10_likelihood(bases, qualities, haplotype,
                                                                  hmm, band, scratch));
    }
    log_reference.push_back(best_other);
    log_alternate.push_back(pair_hmm_log10_likelihood(
        bases, qualities, alternate_haplotype, hmm, band, scratch));
    weights.push_back(weight);
    strands.push_back((read.record->flags & 0x10U) != 0U ? 1U : 0U);
  };
  for (const auto *carrier : carriers) {
    score(*carrier, 1.0);
  }
  for (std::size_t index = 0U; index < others.size(); index += other_stride) {
    score(*others[index], other_weight);
  }
  if (log_reference.empty()) {
    return result;
  }
  const auto fit = fit_allele_fraction(log_reference, log_alternate, weights);
  result.rescored = true;
  result.fraction = fit.fraction;
  result.lod = fit.lod;
  constexpr double kInformativeLog10Ratio = 1.0;
  for (std::size_t index = 0U; index < log_reference.size(); ++index) {
    // Each read counts for the allele it clearly fits, with no prior, so these
    // counts can be compared with the error rate; a posterior at the fitted
    // fraction would already have removed the gap noise. Reads that fit neither
    // allele by 10:1 say nothing and leave the denominator as well.
    const double log_ratio = log_alternate[index] - log_reference[index];
    if (std::abs(log_ratio) < kInformativeLog10Ratio) {
      continue;
    }
    if (log_ratio > 0.0) {
      result.strand_alternate[strands[index]] += weights[index];
      result.effective_alternate += weights[index];
    }
    result.strand_depth[strands[index]] += weights[index];
    result.effective_depth += weights[index];
  }
  return result;
}

[[nodiscard]] std::string indel_hgvs(const RawIndel &indel,
                                     const std::string &reference) {
  const std::size_t length = reference.size();
  const std::size_t size = indel.right_allele.size();
  if (indel.deletion) {
    const auto end = fold(indel.right + size - 1U, length);
    return size == 1U ? "m." + std::to_string(indel.right) + "del"
                      : "m." + std::to_string(indel.right) + "_" +
                            std::to_string(end) + "del";
  }
  // An insertion copying the bases before it is a duplication in HGVS.
  bool duplication = true;
  for (std::size_t offset = 0U; offset < size; ++offset) {
    const auto position = fold(indel.right + length - size + 1U + offset, length);
    if (reference[position - 1U] != indel.right_allele[offset]) {
      duplication = false;
      break;
    }
  }
  if (duplication) {
    const auto first = fold(indel.right + length - size + 1U, length);
    return size == 1U ? "m." + std::to_string(indel.right) + "dup"
                      : "m." + std::to_string(first) + "_" +
                            std::to_string(indel.right) + "dup";
  }
  return "m." + std::to_string(indel.right) + "_" +
         std::to_string(fold(indel.right + 1U, length)) + "ins" + indel.right_allele;
}

/** rCRS 3107 is a placeholder (N): reads lack that base by construction, and an
 * aligner may put the resulting gap on either side of it, sometimes with a
 * compensating change a few bases away. */
[[nodiscard]] bool is_placeholder_artifact(const RawIndel &indel,
                                           const std::string &reference) {
  if (indel.left_allele.size() > 3U) {
    return false;
  }
  constexpr std::size_t kPlaceholderReach = 4U;
  for (const auto position : {indel.left, indel.right}) {
    for (std::size_t neighbor = position > kPlaceholderReach ? position - kPlaceholderReach : 1U;
         neighbor <= std::min(reference.size(), position + kPlaceholderReach); ++neighbor) {
      if (reference[neighbor - 1U] == 'N') {
        return true;
      }
    }
  }
  return false;
}

[[nodiscard]] std::vector<PrecisionIndel>
call_indels(const Pileup &pileup, const std::vector<MoleculeView> &molecules,
            const std::vector<char> &excluded, const std::string &reference,
            const std::string &consensus, const RepeatContext &repeats,
            const AnalysisConfig &config, const std::vector<std::string> &molecule_ids,
            PrecisionDiagnostics &diagnostics) {
  const std::size_t length = reference.size();
  std::vector<RawIndel> indels;
  indels.reserve(pileup.indels.size());
  for (const auto &indel : pileup.indels) {
    if (!is_placeholder_artifact(indel, reference)) {
      indels.push_back(indel);
    }
  }
  diagnostics.indel_candidates = indels.size();
  const auto model = fit_indel_model(indels, pileup, consensus, repeats, diagnostics);

  std::vector<IndelRead> reads;
  for (std::size_t index = 0U; index < molecules.size(); ++index) {
    if (excluded[index] != 0) {
      continue;
    }
    for (const auto *record : molecules[index].records) {
      reads.push_back({record, molecules[index].index, reference_end_of(*record), {}});
    }
  }

  const auto fit_for = [&](const IndelClass &key) {
    const auto fit = model.fits.find(key);
    return fit == model.fits.end() ? BetaBinomialFit{1e-5, 5e-4, 0U} : fit->second;
  };

  struct Scored {
    std::size_t index = 0;
    std::size_t bin = 1;
    std::size_t repeat_bin = 1;
    std::array<std::size_t, 2> depth{};
    double raw_rate = 0.0;
    double background = 0.0;
  };
  std::vector<Scored> scored;
  for (std::size_t index = 0U; index < indels.size(); ++index) {
    const auto &indel = indels[index];
    const std::size_t total = indel.molecules.size();
    Scored item;
    item.index = index;
    item.bin = std::min<std::size_t>(indel.left_allele.size(), 3U);
    item.repeat_bin = indel_repeat_bin(indel, consensus, repeats);
    item.depth = {pileup.span[0][indel.left], pileup.span[1][indel.left]};
    const double depth = static_cast<double>(item.depth[0] + item.depth[1]);
    if (depth == 0.0) {
      continue;
    }
    item.raw_rate = static_cast<double>(total) / depth;
    const auto forward = fit_for({indel.deletion, item.bin, item.repeat_bin, 0U});
    const auto reverse = fit_for({indel.deletion, item.bin, item.repeat_bin, 1U});
    item.background = (forward.mean * static_cast<double>(item.depth[0]) +
                       reverse.mean * static_cast<double>(item.depth[1])) /
                      depth;
    scored.push_back(item);
  }

  // Rescore what could plausibly be real: events above the error model on raw
  // counts, and longer events, whose gaps an aligner often splits into pieces
  // so that raw counts understate them.
  std::vector<std::size_t> to_rescore;
  std::vector<const RawIndel *> plausible;
  for (std::size_t index = 0U; index < scored.size(); ++index) {
    const auto &item = scored[index];
    const auto &indel = indels[item.index];
    const std::size_t total = indel.molecules.size();
    if (total < 2U) {
      continue;
    }
    const auto forward = fit_for({indel.deletion, item.bin, item.repeat_bin, 0U});
    const auto reverse = fit_for({indel.deletion, item.bin, item.repeat_bin, 1U});
    const double raw_p = fisher_combined(
        beta_binomial_survival(indel.strands[0], item.depth[0], forward.mean, forward.rho),
        beta_binomial_survival(indel.strands[1], item.depth[1], reverse.mean, reverse.rho));
    // Realignment rarely gains much for a short event. If even three times its
    // raw rate would stay inside this sample's noise tail, it cannot pass, and
    // the HMM would be spent on noise. Long events keep the HMM, since aligners
    // fragment them and their raw counts say little.
    const double best_case_rate = std::min(1.0, 3.0 * item.raw_rate);
    constexpr double kOptimisticDiscoveries = 100.0;
    const double reachable_p =
        kFalseDiscoveryRate * kOptimisticDiscoveries / static_cast<double>(length * 2U);
    const bool within_noise =
        empirical_tail_probability(model.tail_for(indel.deletion, item.bin, item.repeat_bin),
                                   best_case_rate) > reachable_p;
    const bool long_event = indel.left_allele.size() >= 4U && item.raw_rate >= 0.002;
    if ((raw_p < 1e-3 && !within_noise) || long_event) {
      to_rescore.push_back(index);
    }
    if (item.raw_rate >= std::max(0.02, 2.0 * item.background)) {
      plausible.push_back(&indel);
    }
  }
  std::sort(to_rescore.begin(), to_rescore.end(),
            [&](const std::size_t lhs, const std::size_t rhs) {
              if (scored[lhs].raw_rate != scored[rhs].raw_rate) {
                return scored[lhs].raw_rate > scored[rhs].raw_rate;
              }
              return scored[lhs].index < scored[rhs].index;
            });
  constexpr std::size_t kMaximumRescored = 2000U;
  if (to_rescore.size() > kMaximumRescored) {
    to_rescore.resize(kMaximumRescored);
  }
  if (!to_rescore.empty()) {
    parallel_for_index(reads.size(), config, 64U, [&](const std::size_t index, std::size_t) {
      reads[index].checkpoints = cigar_checkpoints(*reads[index].record);
    });
  }
  std::vector<RescoreResult> rescored(scored.size());
  std::vector<PairHmmScratch> scratches(parallel_worker_limit(config));
  WorkBudget budget(config, kRescoreWorkBase + kRescoreWorkPerMolecule * molecules.size(),
                    "indel realignment");
  parallel_for_index(
      to_rescore.size(), config, 1U,
      [&](const std::size_t index, const std::size_t worker) {
        const auto scored_index = to_rescore[index];
        const auto &indel = indels[scored[scored_index].index];
        // The repeat span comes from the sample consensus, which a crafted input
        // controls; human mtDNA repeats are far shorter than this cap.
        const std::size_t flank = std::min<std::size_t>(
            kMaximumRescoreFlank,
            std::max<std::size_t>(30U, static_cast<std::size_t>(repeats.span[indel.left]) + 15U));
        std::vector<const RawIndel *> competitors;
        for (const auto *other : plausible) {
          // Alleles of one locus overlap once placed at their outermost positions.
          if (other != &indel && other->left <= indel.right + indel.left_allele.size() &&
              indel.left <= other->right + other->left_allele.size() &&
              competitors.size() < 4U) {
            competitors.push_back(other);
          }
        }
        rescored[scored_index] = rescore_indel(indel, competitors, reads, consensus, model.hmm,
                                               flank, scratches[worker], budget);
      });
  diagnostics.indels_rescored = to_rescore.size();

  std::vector<double> p_values(scored.size(), 1.0);
  struct Counts {
    std::array<std::size_t, 2> alternate{};
    std::array<std::size_t, 2> depth{};
  };
  std::vector<Counts> counts(scored.size());
  for (std::size_t index = 0U; index < scored.size(); ++index) {
    const auto &item = scored[index];
    const auto &indel = indels[item.index];
    auto &count = counts[index];
    // The error model is fitted on molecules that show the gap in their alignment,
    // so a short event is tested on that same count. Rescored counts run higher on
    // noisy reads, and testing them against it would mistake that for signal.
    // Long events do use them: an aligner fragments those, so raw counts mislead.
    if (rescored[index].rescored && indel.left_allele.size() >= kFragmentedEventLength) {
      for (std::size_t strand = 0U; strand < 2U; ++strand) {
        count.alternate[strand] =
            static_cast<std::size_t>(std::lround(rescored[index].strand_alternate[strand]));
        count.depth[strand] =
            static_cast<std::size_t>(std::lround(rescored[index].strand_depth[strand]));
      }
    } else {
      for (std::size_t strand = 0U; strand < 2U; ++strand) {
        count.alternate[strand] = indel.strands[strand];
        count.depth[strand] = std::max<std::size_t>(item.depth[strand], indel.strands[strand]);
      }
    }
    const auto forward = fit_for({indel.deletion, item.bin, item.repeat_bin, 0U});
    const auto reverse = fit_for({indel.deletion, item.bin, item.repeat_bin, 1U});
    const double p_forward = beta_binomial_survival(count.alternate[0], count.depth[0],
                                                    forward.mean, forward.rho);
    const double p_reverse = beta_binomial_survival(count.alternate[1], count.depth[1],
                                                    reverse.mean, reverse.rho);
    const double sampling_p = count.depth[0] >= kMinimumStrandDepth &&
                                      count.depth[1] >= kMinimumStrandDepth
                                  ? fisher_combined(p_forward, p_reverse)
                                  : (count.depth[0] >= count.depth[1] ? p_forward : p_reverse);
    // A count model explains sampling. Motif-specific errors recur at the same
    // site on every read, so a call also has to exceed what the noise itself
    // reaches across this sample's genome.
    const auto total_depth = count.depth[0] + count.depth[1];
    const double tested_rate =
        total_depth == 0U ? 0.0
                          : static_cast<double>(count.alternate[0] + count.alternate[1]) /
                                static_cast<double>(total_depth);
    const double noise_p = empirical_tail_probability(
        model.tail_for(indel.deletion, item.bin, item.repeat_bin), tested_rate);
    p_values[index] = std::max(sampling_p, noise_p);
  }
  const std::size_t indel_hypotheses = std::max(scored.size(), length * 2U);
  const auto q_values = benjamini_hochberg(p_values, indel_hypotheses);
  // How high one-base indel noise reaches outside repeats in this sample. It is a
  // measured rate, reported so that a reader can judge any low-level indel call.
  diagnostics.indel_noise_p99 = std::max(model.tail_for(true, 1U, 1U).threshold,
                                         model.tail_for(false, 1U, 1U).threshold);

  std::vector<PrecisionIndel> calls;
  for (std::size_t index = 0U; index < scored.size(); ++index) {
    const auto &item = scored[index];
    const auto &indel = indels[item.index];
    const auto &count = counts[index];
    const auto &score = rescored[index];
    const std::size_t alternate_total = count.alternate[0] + count.alternate[1];
    const std::size_t depth_total = count.depth[0] + count.depth[1];
    if (depth_total == 0U || indel.molecules.size() < 2U) {
      continue;
    }
    PrecisionIndel call;
    call.deletion = indel.deletion;
    call.left = indel.left;
    call.right = indel.right;
    call.length = indel.left_allele.size();
    call.allele = indel.left_allele;
    call.id = std::string(indel.deletion ? "pdel:" : "pins:") +
              std::to_string(indel.left) + ":" + indel.left_allele;
    const std::size_t anchor = indel.deletion
                                   ? (indel.left > 1U ? indel.left - 1U : length)
                                   : indel.left;
    call.vcf_position = anchor;
    const char anchor_base = reference[anchor - 1U];
    if (indel.deletion) {
      call.vcf_reference = std::string(1U, anchor_base) + indel.left_allele;
      call.vcf_alternate = std::string(1U, anchor_base);
    } else {
      call.vcf_reference = std::string(1U, anchor_base);
      call.vcf_alternate = std::string(1U, anchor_base) + indel.left_allele;
    }
    call.hgvs = indel_hgvs(indel, reference);
    call.depth = item.depth[0] + item.depth[1];
    call.raw_alternate_depth = indel.molecules.size();
    call.rescored = score.rescored;
    call.lod = score.lod;
    call.effective_alternate_depth = score.rescored
                                         ? score.effective_alternate
                                         : static_cast<double>(alternate_total);
    call.effective_depth =
        score.rescored ? score.effective_depth : static_cast<double>(depth_total);
    call.strands = {count.alternate[0], count.alternate[1], count.depth[0], count.depth[1]};
    call.p_value = p_values[index];
    call.q_value = q_values[index];
    call.background_rate = item.background;
    call.repeat = {repeats.homopolymer_run[indel.left], repeats.unit[indel.left],
                   repeats.span[indel.left]};
    const auto forward = fit_for({indel.deletion, item.bin, item.repeat_bin, 0U});
    const auto reverse = fit_for({indel.deletion, item.bin, item.repeat_bin, 1U});
    const double rho = (forward.rho + reverse.rho) / 2.0;
    const double observed =
        static_cast<double>(alternate_total) / static_cast<double>(depth_total);
    // The EM fit already accounts for gap noise through its gap-open rates.
    // Counted events have no such model, so their background is removed here.
    const double correction = std::max(0.25, 1.0 - 2.0 * item.background);
    const auto corrected = [&](const double value) {
      return score.rescored
                 ? std::clamp(value, 0.0, 1.0)
                 : std::clamp((value - item.background) / correction, 0.0, 1.0);
    };
    call.fraction = corrected(score.rescored ? score.fraction : observed);
    const double interval_count = score.rescored
                                      ? score.fraction * call.effective_depth
                                      : call.effective_alternate_depth;
    const auto interval = overdispersed_wilson(interval_count, call.effective_depth, rho);
    call.ci95_low = std::min(corrected(interval.first), call.fraction);
    call.ci95_high = std::max(corrected(interval.second), call.fraction);
    if (alternate_total < 3U) {
      call.filters.emplace_back("low_support");
    }
    if (observed < kBackgroundEffectMultiple * item.background) {
      call.filters.emplace_back("below_error_floor");
    }
    if (score.rescored && score.lod < 3.0) {
      call.filters.emplace_back("low_lod");
    }
    // Most reads carrying a real short indel also show it in their alignment. A
    // call resting mainly on realignment fits noisy bases better than it fits
    // this allele. Long events are exempt, because aligners do fragment them.
    constexpr double kMinimumAlignmentSupport = 0.6;
    if (score.rescored && call.length < kFragmentedEventLength &&
        static_cast<double>(call.raw_alternate_depth) <
            kMinimumAlignmentSupport * call.effective_alternate_depth) {
      call.filters.emplace_back("realignment_only_support");
    }
    if (count.depth[0] >= kMinimumStrandDepth && count.depth[1] >= kMinimumStrandDepth) {
      const double rate_forward =
          static_cast<double>(count.alternate[0]) / static_cast<double>(count.depth[0]);
      const double rate_reverse =
          static_cast<double>(count.alternate[1]) / static_cast<double>(count.depth[1]);
      const double strand_p = fisher_exact_two_sided(
          count.alternate[0], count.depth[0] - std::min(count.alternate[0], count.depth[0]),
          count.alternate[1], count.depth[1] - std::min(count.alternate[1], count.depth[1]));
      if (strand_p < 1e-3 &&
          std::min(rate_forward, rate_reverse) < 0.25 * std::max(rate_forward, rate_reverse)) {
        call.filters.emplace_back("strand_bias");
      }
      const double p_forward = beta_binomial_survival(count.alternate[0], count.depth[0],
                                                      forward.mean, forward.rho);
      const double p_reverse = beta_binomial_survival(count.alternate[1], count.depth[1],
                                                      reverse.mean, reverse.rho);
      if (std::max(p_forward, p_reverse) > kWeakerStrandLimit && observed < 0.2) {
        call.filters.emplace_back("one_strand_support");
      }
    } else {
      call.flags.emplace_back("single_strand_depth");
    }
    if (item.repeat_bin >= 7U) {
      call.flags.emplace_back("long_repeat_uncertain");
    }
    if (!score.rescored) {
      call.flags.emplace_back("counted_without_rescoring");
    }
    if (call.q_value > kFalseDiscoveryRate) {
      call.filters.emplace_back("not_significant");
    }
    call.pass = call.filters.empty();
    if (call.pass) {
      for (const auto molecule : indel.molecules) {
        call.supporting_molecules.push_back(molecule_ids[molecule]);
      }
    }
    calls.push_back(std::move(call));
  }
  return calls;
}

// ---------------------------------------------------------- large deletions

void collect_cigar_junctions(const ReadRecord &read, const std::uint32_t molecule,
                             const bool sequence_is_representative,
                             const std::size_t length,
                             std::vector<Junction> &junctions) {
  std::size_t query = 0U;
  std::size_t cursor = read.reference_start;
  std::size_t aligned_before = 0U;
  const std::size_t aligned_total = aligned_query_bases(read);
  for (const auto &operation : read.cigar_operations) {
    const auto size = operation.length;
    const char op = operation.code;
    if (op == 'M' || op == '=' || op == 'X') {
      query += size;
      cursor += size;
      aligned_before += size;
    } else if (op == 'I' || op == 'S') {
      query += size;
    } else if (op == 'D' || op == 'N') {
      if (size >= kMinimumLargeDeletion && size < length &&
          aligned_before >= kJunctionAnchor &&
          aligned_total - aligned_before >= kJunctionAnchor && cursor > 1U) {
        Junction junction;
        junction.left = static_cast<std::uint32_t>(fold(cursor - 1U, length));
        junction.right = static_cast<std::uint32_t>(fold(cursor + size, length));
        junction.length = static_cast<std::uint32_t>(size);
        junction.molecule = molecule;
        junction.strand = (read.flags & 0x10U) != 0U ? 1U : 0U;
        // Refinement reads the junction from this record's SEQ only when it
        // is the representative, whose coordinates the refinement step uses.
        junction.query_index = sequence_is_representative
                                   ? static_cast<std::uint32_t>(query)
                                   : std::numeric_limits<std::uint32_t>::max();
        junction.reverse_sequence = false;
        junctions.push_back(junction);
      }
      cursor += size;
    }
  }
}

void collect_split_junctions(const ReadRecord &representative,
                             const std::uint32_t molecule, const std::size_t length,
                             const std::uint8_t minimum_mapq,
                             std::vector<Junction> &junctions) {
  const auto segments = alignment_segments(representative);
  if (segments.size() < 2U) {
    return;
  }
  const char representative_strand =
      (representative.flags & 0x10U) != 0U ? '-' : '+';
  for (std::size_t index = 1U; index < segments.size(); ++index) {
    const auto &previous = segments[index - 1U];
    const auto &current = segments[index];
    if (previous.strand != current.strand ||
        previous.mapping_quality < minimum_mapq ||
        current.mapping_quality < minimum_mapq ||
        previous.reference_end - previous.reference_start + 1U < kJunctionAnchor ||
        current.reference_end - current.reference_start + 1U < kJunctionAnchor) {
      continue;
    }
    const bool forward = previous.strand == '+';
    // In reference orientation, a reverse-strand path runs current then previous.
    std::int64_t left = forward ? static_cast<std::int64_t>(previous.reference_end)
                                : static_cast<std::int64_t>(current.reference_end);
    std::int64_t right = forward ? static_cast<std::int64_t>(current.reference_start)
                                 : static_cast<std::int64_t>(previous.reference_start);
    const auto query_gap = static_cast<std::int64_t>(current.query_start) -
                           static_cast<std::int64_t>(previous.query_end) - 1;
    // Two alignments that both claim the same read bases are microhomology:
    // keep those bases once, or the deletion looks shorter than it is.
    if (query_gap < 0) {
      left += query_gap;
    }
    if (query_gap > 50) {
      continue;
    }
    std::int64_t gap = right - left - 1;
    const auto reference_length = static_cast<std::int64_t>(length);
    if (gap < 0) {
      // The path leaves the end of the rCRS linearization and re-enters at its
      // start; on the circle the two segments are separated by this much.
      gap += reference_length;
    }
    if (gap < static_cast<std::int64_t>(kMinimumLargeDeletion) ||
        gap > reference_length - static_cast<std::int64_t>(kMinimumLargeDeletion)) {
      continue;
    }
    Junction junction;
    const auto normalized_left =
        left < 1 ? static_cast<std::size_t>(left + reference_length)
                 : static_cast<std::size_t>(left);
    junction.left = static_cast<std::uint32_t>(fold(normalized_left, length));
    junction.right = static_cast<std::uint32_t>(fold(static_cast<std::size_t>(right), length));
    junction.length = static_cast<std::uint32_t>(gap);
    junction.molecule = molecule;
    junction.strand = forward ? 0U : 1U;
    junction.from_split = true;
    // Index in SEQ just after the left part, in read order.
    const std::size_t query_after_left =
        query_gap < 0
            ? static_cast<std::size_t>(static_cast<std::int64_t>(previous.query_end) + 1 + query_gap)
            : previous.query_end + 1U;
    const bool sequence_reversed = representative_strand == '-';
    const auto read_length = representative.sequence.size();
    if (query_after_left <= read_length) {
      junction.query_index = static_cast<std::uint32_t>(
          sequence_reversed ? read_length - query_after_left : query_after_left);
    } else {
      junction.query_index = std::numeric_limits<std::uint32_t>::max();
    }
    junction.reverse_sequence = (representative_strand == '-') != !forward;
    junctions.push_back(junction);
  }
}

struct NormalizedJunction {
  std::size_t left = 0;  // Last retained base, leftmost placement.
  std::size_t right = 0; // First retained base, leftmost placement.
  std::size_t shift = 0; // Equivalent placements to the right.
};

[[nodiscard]] NormalizedJunction normalize_junction(const std::string &reference,
                                                    const std::size_t left,
                                                    const std::size_t right) {
  const std::size_t length = reference.size();
  const auto at = [&](const std::size_t position) {
    return reference[(position + length - 1U) % length];
  };
  std::size_t new_left = left;
  std::size_t new_right = right;
  std::size_t steps = 0U;
  while (steps < 200U && at(new_left) == at(new_right + length - 1U) &&
         at(new_left) != 'N') {
    new_left = new_left > 1U ? new_left - 1U : length;
    new_right = new_right > 1U ? new_right - 1U : length;
    ++steps;
  }
  std::size_t shift = 0U;
  std::size_t probe_left = new_left;
  std::size_t probe_right = new_right;
  while (shift < 400U && at(probe_left + 1U) == at(probe_right) &&
         at(probe_right) != 'N') {
    probe_left = fold(probe_left + 1U, length);
    probe_right = fold(probe_right + 1U, length);
    ++shift;
  }
  return {new_left, new_right, shift};
}

[[nodiscard]] std::string junction_haplotype(const std::string &reference,
                                             const std::size_t left,
                                             const std::size_t right,
                                             const std::size_t flank) {
  const std::size_t length = reference.size();
  std::string haplotype;
  haplotype.reserve(2U * flank);
  for (std::size_t offset = 0U; offset < flank; ++offset) {
    haplotype.push_back(reference[fold(left + length - flank + 1U + offset, length) - 1U]);
  }
  for (std::size_t offset = 0U; offset < flank; ++offset) {
    haplotype.push_back(reference[fold(right + offset, length) - 1U]);
  }
  std::erase(haplotype, 'N');
  return haplotype;
}

/**
 * Links junctions whose placements differ by at most the placement tolerance
 * and whose lengths differ by at most the length tolerance; the labels are
 * union-find roots.
 *
 * Unions run in the order of a pairwise scan over the sorted junctions, so the
 * roots match that scan exactly, but only matching pairs are visited. Sorting
 * by (left, length) makes the matching seconds for one left value a contiguous
 * length range, found by binary search. A junction with the same (left, length)
 * as the one before it is skipped: every pair it would link was linked by its
 * predecessor, so its unions are no-ops. Work is linear in the junction count
 * times the tolerance window, not quadratic in junctions sharing a placement.
 */
[[nodiscard]] std::vector<std::size_t> cluster_junctions(std::vector<Junction> &junctions,
                                                         const AnalysisConfig &config) {
  std::sort(junctions.begin(), junctions.end(), [](const Junction &lhs, const Junction &rhs) {
    return std::tie(lhs.left, lhs.length, lhs.molecule, lhs.right) <
           std::tie(rhs.left, rhs.length, rhs.molecule, rhs.right);
  });
  const std::size_t count = junctions.size();
  std::vector<std::size_t> parent(count);
  std::iota(parent.begin(), parent.end(), std::size_t{0});
  const auto find = [&](std::size_t item) {
    while (parent[item] != item) {
      parent[item] = parent[parent[item]];
      item = parent[item];
    }
    return item;
  };
  // End of the run of equal left values that starts at or contains each index.
  std::vector<std::size_t> block_end(count);
  for (std::size_t index = count; index-- > 0U;) {
    block_end[index] = index + 1U < count && junctions[index + 1U].left == junctions[index].left
                           ? block_end[index + 1U]
                           : index + 1U;
  }
  const auto by_length = [](const Junction &junction, const std::uint64_t value) {
    return junction.length < value;
  };
  const auto length_below = [](const std::uint64_t value, const Junction &junction) {
    return value < junction.length;
  };
  for (std::size_t first = 0U; first < count; ++first) {
    if ((first & 4095U) == 0U) {
      throw_if_cancelled(config);
    }
    const auto &anchor = junctions[first];
    if (first > 0U && junctions[first - 1U].left == anchor.left &&
        junctions[first - 1U].length == anchor.length) {
      continue;
    }
    const std::uint64_t shortest =
        anchor.length > kJunctionLengthTolerance ? anchor.length - kJunctionLengthTolerance : 0U;
    const std::uint64_t longest =
        static_cast<std::uint64_t>(anchor.length) + kJunctionLengthTolerance;
    std::size_t block = first + 1U;
    while (block < count && junctions[block].left - anchor.left <= kJunctionPlacementTolerance) {
      const auto block_last = junctions.begin() + static_cast<std::ptrdiff_t>(block_end[block]);
      const auto low = std::lower_bound(junctions.begin() + static_cast<std::ptrdiff_t>(block),
                                        block_last, shortest, by_length);
      const auto high = std::upper_bound(low, block_last, longest, length_below);
      for (auto second = low; second != high; ++second) {
        const auto index = static_cast<std::size_t>(second - junctions.begin());
        parent[find(index)] = find(first);
      }
      block = block_end[block];
    }
  }
  std::vector<std::size_t> labels(junctions.size());
  for (std::size_t index = 0U; index < junctions.size(); ++index) {
    labels[index] = find(index);
  }
  return labels;
}

/** Interval containing [first, last] on the circle, split at the origin. */
[[nodiscard]] std::size_t count_spanning(const std::vector<AlignedInterval> &intervals,
                                         const std::size_t position,
                                         const std::size_t length,
                                         const std::vector<char> &ignored,
                                         std::vector<std::uint32_t> &stamp,
                                         const std::uint32_t generation) {
  if (position <= kSpanningFlank || position + kSpanningFlank > length) {
    return 0U;
  }
  const std::size_t first = position - kSpanningFlank;
  const std::size_t last = position + kSpanningFlank;
  std::size_t count = 0U;
  for (const auto &interval : intervals) {
    if (interval.start > first) {
      break;
    }
    if (interval.end >= last && ignored[interval.molecule] == 0 &&
        stamp[interval.molecule] != generation) {
      stamp[interval.molecule] = generation;
      ++count;
    }
  }
  return count;
}

[[nodiscard]] std::vector<PrecisionDeletion>
call_large_deletions(const std::vector<MoleculeView> &molecules,
                     const std::vector<char> &excluded, const Pileup &pileup,
                     const std::string &reference, const std::string &consensus,
                     const PairHmmParameters &hmm, const AnalysisConfig &config,
                     const std::vector<std::string> &molecule_ids,
                     PrecisionDiagnostics &diagnostics) {
  const std::size_t length = reference.size();
  std::vector<std::vector<Junction>> per_molecule(molecules.size());
  parallel_for_index(
      molecules.size(), config, 32U, [&](const std::size_t index, std::size_t) {
        if (excluded[index] != 0) {
          return;
        }
        const auto &view = molecules[index];
        auto &found = per_molecule[index];
        for (const auto *record : view.records) {
          collect_cigar_junctions(*record, view.index, record == view.representative ||
                                      record->sequence == view.representative->sequence,
                                  length, found);
        }
        collect_split_junctions(*view.representative, view.index, length,
                                config.min_mapping_quality, found);
      });
  std::vector<Junction> junctions;
  for (auto &found : per_molecule) {
    junctions.insert(junctions.end(), found.begin(), found.end());
  }
  // Each junction goes to its leftmost placement before clustering, because
  // microhomology moves a placement without changing the event.
  // Placement is judged on this sample's own sequence: a fixed variant inside a
  // direct repeat, such as 8473C of haplogroup D4 in the common-deletion repeat,
  // limits how far a junction can shift.
  for (auto &junction : junctions) {
    const auto normalized = normalize_junction(consensus, junction.left, junction.right);
    // The gap moves in the read with the placement, so the refinement window
    // stays centred on it. Without this, a gap that an aligner put at the 3' end
    // of a direct repeat leaves the window off-centre, and the search picks a
    // nearby placement of the wrong length.
    const std::size_t moved = (junction.left + length - normalized.left) % length;
    if (junction.query_index != std::numeric_limits<std::uint32_t>::max() && moved != 0U) {
      if (junction.reverse_sequence) {
        junction.query_index += static_cast<std::uint32_t>(moved);
      } else if (junction.query_index >= moved) {
        junction.query_index -= static_cast<std::uint32_t>(moved);
      } else {
        junction.query_index = std::numeric_limits<std::uint32_t>::max();
      }
    }
    junction.left = static_cast<std::uint32_t>(normalized.left);
    junction.right = static_cast<std::uint32_t>(normalized.right);
  }
  std::vector<PrecisionDeletion> calls;
  if (junctions.empty()) {
    return calls;
  }
  const auto labels = cluster_junctions(junctions, config);
  std::map<std::size_t, std::vector<std::size_t>> clusters;
  for (std::size_t index = 0U; index < junctions.size(); ++index) {
    clusters[labels[index]].push_back(index);
  }

  std::unordered_map<std::uint32_t, std::size_t> view_of_molecule;
  view_of_molecule.reserve(molecules.size());
  for (std::size_t index = 0U; index < molecules.size(); ++index) {
    view_of_molecule.emplace(molecules[index].index, index);
  }

  std::size_t singletons = 0U;
  std::vector<std::vector<std::size_t>> kept_clusters;
  std::set<std::uint32_t> junction_molecules;
  for (auto &[_, members] : clusters) {
    std::set<std::uint32_t> distinct;
    for (const auto member : members) {
      distinct.insert(junctions[member].molecule);
      junction_molecules.insert(junctions[member].molecule);
    }
    if (distinct.size() < 2U) {
      ++singletons;
    } else {
      kept_clusters.push_back(members);
    }
  }
  diagnostics.junction_molecules = junction_molecules.size();
  diagnostics.singleton_junctions = singletons;
  const double used_molecules = std::max(1.0, static_cast<double>(diagnostics.molecules_used));
  const double chimera_rate = static_cast<double>(singletons) / used_molecules;
  // A random chimera lands in one cluster's window with this probability.
  const double window_area = static_cast<double>(4U * kJunctionPlacementTolerance *
                                                 kJunctionLengthTolerance);
  const double arrangement_area =
      static_cast<double>(length) * static_cast<double>(length) / 2.0;
  const double background = std::max(chimera_rate * window_area / arrangement_area, 1e-7);

  std::vector<std::uint32_t> stamp(molecule_ids.size(), 0U);
  std::vector<char> ignored(stamp.size(), 0);
  std::uint32_t generation = 0U;
  PairHmmScratch scratch;
  std::vector<double> p_values;
  // Each cluster scans the aligned intervals and may run the breakpoint search,
  // so the cluster count, which crafted junctions can inflate, is budgeted.
  WorkBudget budget(config, kRefineWorkBase + kRefineWorkPerMolecule * molecules.size(),
                    "deletion breakpoint search");
  for (auto &members : kept_clusters) {
    budget.charge(2U * static_cast<std::uint64_t>(pileup.intervals.size()) + members.size());
    // Candidate placements: the most frequent, by exact (left, right).
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::size_t> placements;
    for (const auto member : members) {
      ++placements[{junctions[member].left, junctions[member].right}];
    }
    std::vector<std::pair<std::pair<std::uint32_t, std::uint32_t>, std::size_t>> ranked(
        placements.begin(), placements.end());
    std::stable_sort(ranked.begin(), ranked.end(),
                     [](const auto &lhs, const auto &rhs) { return lhs.second > rhs.second; });
    constexpr std::size_t kMaximumPlacements = 6U;
    if (ranked.size() > kMaximumPlacements) {
      ranked.resize(kMaximumPlacements);
    }
    // Maximum likelihood over reads that carry the junction in their sequence.
    constexpr std::size_t kRefineFlank = 40U;
    constexpr std::size_t kMaximumRefineReads = 40U;
    constexpr std::size_t kRefineBand = 16U;
    constexpr std::size_t kPlacementSearchWidth = 9U;
    std::vector<std::pair<std::string, std::string>> windows;
    for (const auto member : members) {
      if (windows.size() >= kMaximumRefineReads) {
        break;
      }
      const auto &junction = junctions[member];
      if (junction.query_index == std::numeric_limits<std::uint32_t>::max()) {
        continue;
      }
      const auto view = view_of_molecule.find(junction.molecule);
      if (view == view_of_molecule.end()) {
        continue;
      }
      const auto &sequence = molecules[view->second].representative->sequence;
      const auto &qualities = molecules[view->second].representative->qualities;
      if (junction.query_index < kRefineFlank ||
          junction.query_index + kRefineFlank > sequence.size()) {
        continue;
      }
      const auto window_start = junction.query_index - kRefineFlank;
      auto bases = sequence.substr(window_start, 2U * kRefineFlank);
      auto bases_quality = qualities.size() >= window_start + 2U * kRefineFlank
                               ? qualities.substr(window_start, 2U * kRefineFlank)
                               : std::string{};
      if (junction.reverse_sequence) {
        bases = reverse_complement_string(bases);
        std::reverse(bases_quality.begin(), bases_quality.end());
      }
      windows.emplace_back(std::move(bases), std::move(bases_quality));
    }
    std::pair<std::uint32_t, std::uint32_t> best = ranked.front().first;
    if (!windows.empty()) {
      // Reads rarely place a break exactly: a variant in a direct repeat can leave
      // the true placement almost unobserved. So the search also covers the area
      // around the most supported placement. Placements that spell the same junction
      // are one hypothesis, scored once, and the likelihood decides among the rest.
      std::vector<std::pair<std::uint32_t, std::uint32_t>> candidates;
      candidates.reserve(ranked.size() + kPlacementSearchWidth * kPlacementSearchWidth);
      for (const auto &[placement, support] : ranked) {
        candidates.push_back(placement);
      }
      const auto center = ranked.front().first;
      const auto half_width = static_cast<std::int64_t>(kPlacementSearchWidth / 2U);
      for (std::int64_t left_shift = -half_width; left_shift <= half_width; ++left_shift) {
        for (std::int64_t right_shift = -half_width; right_shift <= half_width; ++right_shift) {
          const auto left = static_cast<std::int64_t>(center.first) + left_shift;
          const auto right = static_cast<std::int64_t>(center.second) + right_shift;
          if (left >= 1 && right >= 1 && left <= static_cast<std::int64_t>(length) &&
              right <= static_cast<std::int64_t>(length)) {
            candidates.emplace_back(static_cast<std::uint32_t>(left), static_cast<std::uint32_t>(right));
          }
        }
      }
      std::unordered_set<std::string> scored_haplotypes;
      double best_score = -std::numeric_limits<double>::infinity();
      for (const auto &placement : candidates) {
        auto haplotype = junction_haplotype(consensus, placement.first, placement.second,
                                             kRefineFlank + 10U);
        if (!scored_haplotypes.insert(haplotype).second) {
          continue;
        }
        budget.charge(static_cast<std::uint64_t>(windows.size()) * 2U * kRefineFlank *
                      (2U * kRefineBand + 1U));
        double score = 0.0;
        for (const auto &[bases, quality] : windows) {
          score += pair_hmm_log10_likelihood(bases, quality, haplotype, hmm, kRefineBand, scratch);
        }
        if (score > best_score + 1e-9) {
          best_score = score;
          best = placement;
        }
      }
    }
    const auto normalized = normalize_junction(consensus, best.first, best.second);
    const std::size_t left_retained = normalized.left;
    const std::size_t right_retained = normalized.right;
    const std::size_t deletion_length =
        right_retained > left_retained ? right_retained - left_retained - 1U
                                       : length - left_retained + right_retained - 1U;
    if (deletion_length < kMinimumLargeDeletion) {
      continue;
    }

    PrecisionDeletion call;
    call.length = deletion_length;
    call.left_start = fold(left_retained + 1U, length);
    call.left_end = fold(right_retained + length - 1U, length);
    call.start = fold(left_retained + 1U + normalized.shift, length);
    call.end = fold(right_retained + normalized.shift + length - 1U, length);
    call.id = "pdeletion:" + std::to_string(call.start) + "-" + std::to_string(call.end);
    for (std::size_t offset = 0U; offset < normalized.shift; ++offset) {
      call.microhomology.push_back(consensus[fold(left_retained + 1U + offset, length) - 1U]);
    }
    call.breakpoint_variants = placements.size();
    std::set<std::uint32_t> supporting;
    std::array<std::size_t, 2> strands{};
    std::fill(ignored.begin(), ignored.end(), 0);
    for (const auto member : members) {
      const auto molecule = junctions[member].molecule;
      if (supporting.insert(molecule).second) {
        strands[junctions[member].strand] += 1U;
        ignored[molecule] = 1;
      }
    }
    call.junction_molecules = supporting.size();
    call.junction_forward = strands[0];
    call.junction_reverse = strands[1];
    // Breakpoints sit at the 3' placement here: any placement within a repeat
    // gives the same wild-type spanning count to within the repeat length.
    ++generation;
    call.wild_type_left = count_spanning(pileup.intervals, call.left_start, length,
                                         ignored, stamp, generation);
    ++generation;
    call.wild_type_right = count_spanning(pileup.intervals, call.end, length, ignored,
                                          stamp, generation);
    const double wild_type =
        (static_cast<double>(call.wild_type_left) + static_cast<double>(call.wild_type_right)) / 2.0;
    const double denominator = static_cast<double>(call.junction_molecules) + wild_type;
    if (call.wild_type_left + call.wild_type_right > 0U) {
      call.fraction = static_cast<double>(call.junction_molecules) / denominator;
      const auto interval = overdispersed_wilson(
          static_cast<double>(call.junction_molecules), denominator, 0.0);
      call.ci95_low = interval.first;
      call.ci95_high = interval.second;
    } else {
      call.flags.emplace_back("breakpoints_not_spanned");
    }
    // Median depth inside the arc against outside is independent evidence: a
    // real high-level deletion lowers it, a chimera does not.
    std::vector<double> inside;
    std::vector<double> outside;
    for (std::size_t position = 1U; position <= length; position += 20U) {
      const bool within = call.left_start <= call.left_end
                              ? position >= call.left_start && position <= call.left_end
                              : position >= call.left_start || position <= call.left_end;
      (within ? inside : outside)
          .push_back(static_cast<double>(pileup.span[0][position] + pileup.span[1][position]));
    }
    const double median_inside = median_of(inside);
    const double median_outside = median_of(outside);
    if (!inside.empty() && !outside.empty() && median_outside > 0.0) {
      call.depth_ratio = median_inside / median_outside;
      if (call.fraction && *call.fraction >= 0.2 &&
          *call.depth_ratio > 1.0 - *call.fraction / 2.0) {
        call.flags.emplace_back("depth_discordant");
      }
    }
    if (call.junction_molecules < 3U) {
      call.filters.emplace_back("low_support");
    }
    if (call.junction_molecules >= 8U &&
        (call.junction_forward == 0U || call.junction_reverse == 0U)) {
      call.filters.emplace_back("single_strand_support");
    }
    if (call.microhomology.size() >= 5U) {
      call.flags.emplace_back("direct_repeat_flanked");
    }
    // Reads are linearized at the rCRS origin, so a read crossing it is split and
    // short anchors there produce spurious gaps. Wild-type support cannot be
    // measured at these placements either.
    constexpr std::size_t kOriginMargin = 100U;
    if (call.left_start <= kOriginMargin || call.end <= kOriginMargin ||
        call.left_start + kOriginMargin > length || call.end + kOriginMargin > length) {
      call.filters.emplace_back("linearization_origin");
    }
    // At its 3' placement, so a variant inside the repeat does not hide the event.
    call.known_common_deletion = call.start == 8483U && call.length == 4977U;
    call.impact = annotate_large_deletion(call.start, call.end, length);
    const auto trials = static_cast<std::size_t>(std::llround(std::max(
        denominator, static_cast<double>(call.junction_molecules))));
    call.p_value =
        beta_binomial_survival(call.junction_molecules, trials, background, 0.0);
    for (const auto molecule : supporting) {
      call.supporting_molecules.push_back(molecule_ids[molecule]);
    }
    p_values.push_back(call.p_value);
    calls.push_back(std::move(call));
  }
  const auto q_values = benjamini_hochberg(p_values, calls.size());
  for (std::size_t index = 0U; index < calls.size(); ++index) {
    calls[index].q_value = q_values[index];
    if (calls[index].q_value > kFalseDiscoveryRate) {
      calls[index].filters.emplace_back("not_significant");
    }
    calls[index].pass = calls[index].filters.empty();
  }
  std::sort(calls.begin(), calls.end(), [](const auto &lhs, const auto &rhs) {
    return std::tie(lhs.start, lhs.end) < std::tie(rhs.start, rhs.end);
  });
  return calls;
}

// ----------------------------------------------------- haplotype deconvolution

struct PhasedObservation {
  std::uint16_t site = 0;
  bool alternate = false;
};

struct MixtureFit {
  std::vector<double> weights;
  std::vector<std::vector<double>> allele;
  double log_likelihood = 0.0;
};

[[nodiscard]] MixtureFit
fit_bernoulli_mixture(const std::vector<std::vector<PhasedObservation>> &observations,
                      const std::size_t sites, MixtureFit fit, const double error) {
  const std::size_t components = fit.weights.size();
  std::vector<double> log_joint(components);
  for (int iteration = 0; iteration < 40; ++iteration) {
    std::vector<double> weight_sum(components, 0.0);
    std::vector<std::vector<double>> alternate_sum(components, std::vector<double>(sites, 0.0));
    std::vector<std::vector<double>> total_sum(components, std::vector<double>(sites, 0.0));
    double log_likelihood = 0.0;
    for (const auto &molecule : observations) {
      for (std::size_t k = 0U; k < components; ++k) {
        double value = std::log(std::max(fit.weights[k], 1e-12));
        for (const auto &observation : molecule) {
          const double theta = fit.allele[k][observation.site];
          value += std::log(observation.alternate ? theta : 1.0 - theta);
        }
        log_joint[k] = value;
      }
      const double best = *std::max_element(log_joint.begin(), log_joint.end());
      double normalizer = 0.0;
      for (auto &value : log_joint) {
        value = std::exp(value - best);
        normalizer += value;
      }
      log_likelihood += best + std::log(normalizer);
      for (std::size_t k = 0U; k < components; ++k) {
        const double posterior = log_joint[k] / normalizer;
        weight_sum[k] += posterior;
        for (const auto &observation : molecule) {
          total_sum[k][observation.site] += posterior;
          if (observation.alternate) {
            alternate_sum[k][observation.site] += posterior;
          }
        }
      }
    }
    const double total = static_cast<double>(observations.size());
    for (std::size_t k = 0U; k < components; ++k) {
      fit.weights[k] = weight_sum[k] / total;
      for (std::size_t site = 0U; site < sites; ++site) {
        // A site a component never covers keeps its value: no observation there
        // is not evidence about that component's allele.
        if (total_sum[k][site] > 0.0) {
          fit.allele[k][site] = std::clamp(alternate_sum[k][site] / total_sum[k][site],
                                           error, 1.0 - error);
        }
      }
    }
    const bool converged = std::abs(log_likelihood - fit.log_likelihood) < 1e-4 * total;
    fit.log_likelihood = log_likelihood;
    if (converged) {
      break;
    }
  }
  return fit;
}

/** Number of sites where two components carry different alleles. */
[[nodiscard]] std::size_t distinguishing_sites(const std::vector<double> &first,
                                               const std::vector<double> &second) {
  std::size_t count = 0U;
  for (std::size_t site = 0U; site < first.size(); ++site) {
    if (std::abs(first[site] - second[site]) >= 0.6) {
      ++count;
    }
  }
  return count;
}

void deconvolve_haplotypes(const std::vector<MoleculeView> &molecules,
                           const std::vector<char> &excluded, const std::size_t length,
                           const std::size_t molecule_count, const AnalysisConfig &config,
                           PrecisionCallResult &result) {
  result.molecule_haplotype.assign(molecule_count, -1);
  std::vector<const PrecisionSnv *> sites;
  std::vector<int> site_of_position(length + 1U, -1);
  constexpr std::size_t kMaximumSites = 400U;
  for (const auto &snv : result.snvs) {
    // Repeat-context alleles cannot tell molecules apart reliably.
    if (snv.pass && snv.fraction >= 0.02 && snv.fraction <= 0.98 && snv.repeat.span < 8U &&
        site_of_position[snv.position] < 0 && sites.size() < kMaximumSites) {
      site_of_position[snv.position] = static_cast<int>(sites.size());
      sites.push_back(&snv);
    }
  }

  std::vector<std::vector<PhasedObservation>> phased(molecules.size());
  if (!sites.empty()) {
    parallel_for_index(molecules.size(), config, 32U, [&](const std::size_t index, std::size_t) {
      if (excluded[index] != 0) {
        return;
      }
      auto &observations = phased[index];
      for (const auto *record : molecules[index].records) {
        for_each_aligned_base(*record, length, [&](const std::size_t position, const std::size_t query) {
          const int site = site_of_position[position];
          if (site < 0) {
            return;
          }
          const auto quality = phred_quality_at(*record, query);
          const char base = ascii_upper(record->sequence[query]);
          if (quality && *quality >= config.min_base_quality && is_base(base)) {
            observations.push_back({static_cast<std::uint16_t>(site),
                                    base == sites[static_cast<std::size_t>(site)]->alternate});
          }
        });
      }
      std::sort(observations.begin(), observations.end(),
                [](const PhasedObservation &lhs, const PhasedObservation &rhs) {
                  return lhs.site < rhs.site;
                });
      observations.erase(std::unique(observations.begin(), observations.end(),
                                     [](const PhasedObservation &lhs, const PhasedObservation &rhs) {
                                       return lhs.site == rhs.site;
                                     }),
                         observations.end());
    });
  }

  std::vector<std::vector<PhasedObservation>> informative;
  std::vector<std::size_t> informative_index;
  for (std::size_t index = 0U; index < phased.size(); ++index) {
    if (!phased[index].empty()) {
      informative.push_back(phased[index]);
      informative_index.push_back(index);
    }
  }

  constexpr double kAlleleError = 0.02;
  constexpr std::size_t kMaximumComponents = 4U;
  constexpr std::size_t kMinimumDefiningSites = 3U;
  const std::size_t site_count = sites.size();
  MixtureFit best;
  best.weights = {1.0};
  best.allele.assign(1U, std::vector<double>(site_count));
  for (std::size_t site = 0U; site < site_count; ++site) {
    best.allele[0][site] = std::clamp(sites[site]->fraction, kAlleleError, 1.0 - kAlleleError);
  }
  if (!informative.empty() && site_count > 0U) {
    best = fit_bernoulli_mixture(informative, site_count, std::move(best), kAlleleError);
    const double log_samples = std::log(static_cast<double>(informative.size()));
    const auto bic = [&](const MixtureFit &fit) {
      const double parameters = static_cast<double>(fit.weights.size() * site_count +
                                                    fit.weights.size() - 1U);
      return -2.0 * fit.log_likelihood + parameters * log_samples;
    };
    while (best.weights.size() < kMaximumComponents) {
      // Split the most mixed site of the most mixed component.
      std::size_t split_component = 0U;
      std::size_t split_site = 0U;
      double split_score = -1.0;
      for (std::size_t k = 0U; k < best.weights.size(); ++k) {
        for (std::size_t site = 0U; site < site_count; ++site) {
          const double theta = best.allele[k][site];
          const double score = best.weights[k] * theta * (1.0 - theta);
          if (score > split_score) {
            split_score = score;
            split_component = k;
            split_site = site;
          }
        }
      }
      MixtureFit candidate = best;
      candidate.weights[split_component] /= 2.0;
      candidate.weights.push_back(candidate.weights[split_component]);
      candidate.allele.push_back(candidate.allele[split_component]);
      candidate.allele[split_component][split_site] = kAlleleError;
      candidate.allele.back()[split_site] = 1.0 - kAlleleError;
      candidate = fit_bernoulli_mixture(informative, site_count, std::move(candidate), kAlleleError);

      // A new component must beat the penalty, hold real molecules, and differ
      // from every other component at several sites. One recurrent error on its
      // own must not look like a second lineage.
      const std::size_t added = candidate.weights.size() - 1U;
      bool distinct = candidate.weights[added] * static_cast<double>(informative.size()) >= 5.0;
      for (std::size_t k = 0U; k < added && distinct; ++k) {
        distinct = distinguishing_sites(candidate.allele[k], candidate.allele[added]) >=
                   kMinimumDefiningSites;
      }
      if (!distinct || bic(candidate) >= bic(best) - 10.0) {
        break;
      }
      best = std::move(candidate);
    }
  }

  const std::size_t components = best.weights.size();
  // The largest component is haplotype zero, in a stable order.
  std::vector<std::size_t> order(components);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::stable_sort(order.begin(), order.end(), [&](const std::size_t lhs, const std::size_t rhs) {
    return best.weights[lhs] > best.weights[rhs];
  });
  std::vector<std::size_t> rank(components);
  for (std::size_t position = 0U; position < components; ++position) {
    rank[order[position]] = position;
  }
  std::vector<std::size_t> assigned(components, 0U);
  std::vector<double> log_joint(components);
  for (std::size_t index = 0U; index < informative.size(); ++index) {
    for (std::size_t k = 0U; k < components; ++k) {
      double value = std::log(std::max(best.weights[k], 1e-12));
      for (const auto &observation : informative[index]) {
        const double theta = best.allele[k][observation.site];
        value += std::log(observation.alternate ? theta : 1.0 - theta);
      }
      log_joint[k] = value;
    }
    const double top = *std::max_element(log_joint.begin(), log_joint.end());
    double normalizer = 0.0;
    std::size_t winner = 0U;
    for (std::size_t k = 0U; k < components; ++k) {
      const double value = std::exp(log_joint[k] - top);
      normalizer += value;
      if (log_joint[k] == top) {
        winner = k;
      }
    }
    if (1.0 / normalizer >= 0.9) {
      ++assigned[winner];
      result.molecule_haplotype[molecules[informative_index[index]].index] =
          static_cast<int>(rank[winner]);
    }
  }
  if (components == 1U) {
    // With one lineage every usable molecule belongs to it, informative or not.
    for (std::size_t index = 0U; index < molecules.size(); ++index) {
      if (excluded[index] == 0) {
        result.molecule_haplotype[molecules[index].index] = 0;
      }
    }
    assigned[0] = result.diagnostics.molecules_used;
  }

  const auto add_marker = [](ReadFeature &feature, const PrecisionSnv &snv) {
    ObservedPhyloMutation marker;
    marker.position = snv.position;
    marker.kind = PhyloMutationKind::substitution;
    marker.alternate = snv.alternate;
    marker.encoded = std::to_string(snv.position) + std::string(1U, snv.alternate);
    feature.haplogroup_markers.push_back(std::move(marker));
  };
  for (std::size_t position = 0U; position < components; ++position) {
    const std::size_t k = order[position];
    PrecisionHaplotype haplotype;
    haplotype.id = static_cast<int>(position);
    haplotype.fraction = best.weights[k];
    haplotype.molecules = assigned[k];
    haplotype.feature.id = "haplotype-" + std::to_string(position);
    haplotype.feature.cluster_id = static_cast<int>(position);
    haplotype.feature.haplogroup_range_known = true;
    haplotype.feature.haplogroup_ranges = result.callable_ranges;
    for (const auto &snv : result.snvs) {
      if (!snv.pass) {
        continue;
      }
      const int site = site_of_position[snv.position];
      const bool phased_site =
          site >= 0 && sites[static_cast<std::size_t>(site)] == &snv;
      const bool carries = phased_site ? best.allele[k][static_cast<std::size_t>(site)] >= 0.5
                                       : snv.fraction >= 0.5;
      if (carries) {
        add_marker(haplotype.feature, snv);
      }
      if (phased_site && position > 0U) {
        const double in_major = best.allele[order[0]][static_cast<std::size_t>(site)];
        const double here = best.allele[k][static_cast<std::size_t>(site)];
        if (std::abs(here - in_major) >= 0.6) {
          haplotype.defining_variants.push_back(
              "snv:" + std::to_string(snv.position) + ":" + std::string(1U, snv.reference) +
              ":" + std::string(1U, snv.alternate));
        }
      }
    }
    std::sort(haplotype.feature.haplogroup_markers.begin(),
              haplotype.feature.haplogroup_markers.end(), observed_phylo_mutation_less);
    result.haplotypes.push_back(std::move(haplotype));
  }
}

template <typename Call>
void retain_reportable(std::vector<Call> &calls) {
  std::vector<Call> kept;
  std::vector<Call> filtered;
  for (auto &call : calls) {
    (call.pass ? kept : filtered).push_back(std::move(call));
  }
  std::stable_sort(filtered.begin(), filtered.end(),
                   [](const auto &lhs, const auto &rhs) { return lhs.p_value < rhs.p_value; });
  if (filtered.size() > kMaximumReportedFiltered) {
    filtered.resize(kMaximumReportedFiltered);
  }
  kept.insert(kept.end(), std::make_move_iterator(filtered.begin()),
              std::make_move_iterator(filtered.end()));
  calls = std::move(kept);
}

} // namespace

[[nodiscard]] PrecisionCallResult
call_precision_variants(const MoleculeAssemblyResult &assembly,
                        const std::vector<ReadFeature> &features,
                        const std::string &reference,
                        const AnalysisConfig &config) {
  PrecisionCallResult result;
  if (reference.empty()) {
    return result;
  }
  AnalysisConfig fragment_config = config;
  fragment_config.excluded_snp_flags =
      static_cast<std::uint16_t>(fragment_config.excluded_snp_flags & ~0x800U);
  const auto molecules = molecule_views(assembly, features, fragment_config, reference.size());
  if (molecules.empty()) {
    return result;
  }
  throw_if_cancelled(config);

  std::vector<char> excluded(molecules.size(), 0);
  auto pileup = build_pileup(molecules, excluded, reference, config);
  auto consensus = consensus_from_pileup(pileup, reference);

  const auto divergence = molecule_divergence(molecules, consensus, config);
  std::vector<double> known;
  known.reserve(divergence.size());
  for (const auto value : divergence) {
    if (value >= 0.0) {
      known.push_back(value);
    }
  }
  const double median_divergence = median_of(known);
  std::vector<double> absolute_deviation;
  absolute_deviation.reserve(known.size());
  for (const auto value : known) {
    absolute_deviation.push_back(std::abs(value - median_divergence));
  }
  const double scale = 1.4826 * median_of(std::move(absolute_deviation));
  // True haplotype mixtures differ by well under one percent, while a NUMT copy
  // or a chimeric contaminant differs by several; the floor keeps a clean, tight
  // read population from excluding ordinary reads.
  const double divergence_limit = median_divergence + std::max(6.0 * scale, 0.02);
  std::size_t excluded_count = 0U;
  for (std::size_t index = 0U; index < molecules.size(); ++index) {
    if (divergence[index] > divergence_limit) {
      excluded[index] = 1;
      ++excluded_count;
    }
  }
  if (excluded_count > 0U) {
    std::vector<char> included(molecules.size(), 1);
    std::vector<char> removed_molecule(assembly.molecules.size(), 0);
    for (std::size_t index = 0U; index < molecules.size(); ++index) {
      if (excluded[index] != 0) {
        included[index] = 0;
        removed_molecule[molecules[index].index] = 1;
      }
    }
    subtract_pileup(pileup, build_pileup(molecules, included, reference, config), removed_molecule);
    consensus = consensus_from_pileup(pileup, reference);
  }
  throw_if_cancelled(config);

  auto &diagnostics = result.diagnostics;
  diagnostics.molecules_used = molecules.size() - excluded_count;
  diagnostics.molecules_excluded_divergent = excluded_count;
  diagnostics.median_molecule_divergence = median_divergence;
  diagnostics.divergence_limit = divergence_limit;
  std::vector<double> depths;
  depths.reserve(reference.size());
  for (std::size_t position = 1U; position <= reference.size(); ++position) {
    const auto depth = pileup.span[0][position] + pileup.span[1][position];
    // Uncovered positions are not low depth: an amplicon covers part of the genome.
    if (depth > 0U) {
      depths.push_back(static_cast<double>(depth));
    }
  }
  diagnostics.median_depth = median_of(std::move(depths));

  std::vector<std::string> molecule_ids(assembly.molecules.size());
  for (std::size_t index = 0U; index < assembly.molecules.size(); ++index) {
    molecule_ids[index] = assembly.molecules[index].id;
  }

  const auto repeats = build_repeat_context(consensus);
  result.snvs = call_snvs(pileup, reference, consensus, repeats, diagnostics);
  throw_if_cancelled(config);
  result.indels = call_indels(pileup, molecules, excluded, reference, consensus, repeats,
                              config, molecule_ids, diagnostics);
  throw_if_cancelled(config);

  // Beside a long run, a slipped read can also align as a pair of mismatches with
  // no gap at all, as at rCRS 309-310. Up to about the run's measured slip rate,
  // a low-level allele here cannot be told from that.
  for (auto &call : result.snvs) {
    const std::size_t length = reference.size();
    const std::size_t touching_run = std::max(
        {static_cast<std::size_t>(repeats.homopolymer_run[fold(call.position + length - 1U, length)]),
         static_cast<std::size_t>(repeats.homopolymer_run[call.position]),
         static_cast<std::size_t>(repeats.homopolymer_run[fold(call.position + 1U, length)])});
    if (touching_run < kUnreliableRun) {
      continue;
    }
    const double slip = diagnostics.slip_rate_by_run[std::min(
        touching_run, diagnostics.slip_rate_by_run.size() - 1U)];
    const double raw = static_cast<double>(call.alternate_depth) / static_cast<double>(call.depth);
    if (raw < std::min(0.5, kSlipArtifactMultiple * slip)) {
      call.filters.emplace_back("repeat_alignment_artifact");
      call.pass = false;
    }
  }

  // The deletion refinement uses the same gap model as the small-indel scoring.
  PairHmmParameters refinement_hmm;
  refinement_hmm.deletion_open.fill(2e-3);
  refinement_hmm.insertion_open.fill(2e-3);
  result.deletions =
      call_large_deletions(molecules, excluded, pileup, reference, consensus, refinement_hmm,
                           config, molecule_ids, diagnostics);

  retain_reportable(result.snvs);
  retain_reportable(result.indels);
  retain_reportable(result.deletions);
  const auto by_snv_position = [](const PrecisionSnv &lhs, const PrecisionSnv &rhs) {
    return std::tie(lhs.position, lhs.alternate) < std::tie(rhs.position, rhs.alternate);
  };
  std::sort(result.snvs.begin(), result.snvs.end(), by_snv_position);
  std::sort(result.indels.begin(), result.indels.end(),
            [](const PrecisionIndel &lhs, const PrecisionIndel &rhs) {
              return std::tie(lhs.left, lhs.deletion, lhs.allele) <
                     std::tie(rhs.left, rhs.deletion, rhs.allele);
            });
  for (auto &snv : result.snvs) {
    if (snv.pass) {
      snv.effects = annotate_substitution(reference, snv.position, snv.alternate);
    }
  }
  for (auto &indel : result.indels) {
    if (indel.pass) {
      indel.effects = annotate_small_indel(reference, indel.right, indel.deletion ? indel.length : 0U,
                                           indel.deletion ? 0U : indel.length, !indel.deletion);
    }
  }
  result.depth_by_position.assign(reference.size() + 1U, 0U);
  for (std::size_t position = 1U; position <= reference.size(); ++position) {
    result.depth_by_position[position] = pileup.span[0][position] + pileup.span[1][position];
  }
  // A position needs a tenth of the typical depth to count as assessable, and
  // always a few molecules; below that an absent variant means nothing.
  const double callable_depth = std::max(3.0, 0.1 * diagnostics.median_depth);
  std::optional<std::size_t> range_start;
  for (std::size_t position = 1U; position <= reference.size(); ++position) {
    const double depth =
        static_cast<double>(pileup.span[0][position] + pileup.span[1][position]);
    if (depth >= callable_depth) {
      if (!range_start) {
        range_start = position;
      }
    } else if (range_start) {
      result.callable_ranges.emplace_back(*range_start, position - 1U);
      range_start.reset();
    }
  }
  if (range_start) {
    result.callable_ranges.emplace_back(*range_start, reference.size());
  }
  deconvolve_haplotypes(molecules, excluded, reference.size(), assembly.molecules.size(), config,
                        result);
  result.consensus = std::move(consensus);
  result.assessed = true;
  return result;
}

[[nodiscard]] ReadFeature
consensus_haplogroup_feature(const PrecisionCallResult &calls,
                             const std::string &reference,
                             const MoleculeAssemblyResult &assembly) {
  static_cast<void>(assembly);
  ReadFeature feature;
  feature.id = "consensus";
  feature.cluster_id = 0;
  feature.haplogroup_range_known = true;
  feature.haplogroup_ranges = calls.callable_ranges;
  const std::size_t length = reference.size();
  for (const auto &snv : calls.snvs) {
    if (!snv.pass || snv.fraction < 0.5) {
      continue;
    }
    ObservedPhyloMutation marker;
    marker.position = snv.position;
    marker.kind = PhyloMutationKind::substitution;
    marker.alternate = snv.alternate;
    marker.encoded = std::to_string(snv.position) + std::string(1U, snv.alternate);
    feature.haplogroup_markers.push_back(std::move(marker));
  }
  for (const auto &indel : calls.indels) {
    if (!indel.pass || indel.fraction < 0.5) {
      continue;
    }
    if (indel.deletion) {
      for (std::size_t offset = 0U; offset < indel.length; ++offset) {
        ObservedPhyloMutation marker;
        marker.position = fold(indel.right + offset, length);
        marker.kind = PhyloMutationKind::deletion;
        marker.encoded = std::to_string(marker.position) + "d";
        feature.haplogroup_markers.push_back(std::move(marker));
      }
    } else {
      const auto normalized =
          normalize_circular_indel(reference, indel.right, 0U, indel.allele);
      ObservedPhyloMutation marker;
      marker.position = normalized.right;
      marker.kind = PhyloMutationKind::insertion;
      marker.insertion_index = "1";
      marker.inserted_bases = normalized.right_allele;
      marker.encoded = std::to_string(normalized.right) + ".1" + normalized.right_allele;
      feature.haplogroup_markers.push_back(std::move(marker));
    }
  }
  std::sort(feature.haplogroup_markers.begin(), feature.haplogroup_markers.end(),
            observed_phylo_mutation_less);
  return feature;
}

[[nodiscard]] CoverageResult
coverage_from_molecule_depth(const std::vector<std::uint32_t> &depth, const std::size_t bin_count) {
  CoverageResult coverage;
  const std::size_t length = depth.empty() ? 0U : depth.size() - 1U;
  if (length == 0U || bin_count == 0U) {
    return coverage;
  }
  std::uint64_t total = 0U;
  std::size_t above_twenty = 0U;
  for (std::size_t position = 1U; position <= length; ++position) {
    total += depth[position];
    coverage.max_depth = std::max<std::size_t>(coverage.max_depth, depth[position]);
    above_twenty += depth[position] > 20U ? 1U : 0U;
  }
  coverage.mean_depth = static_cast<double>(total) / static_cast<double>(length);
  coverage.pct_sites_gt20x = 100.0 * static_cast<double>(above_twenty) / static_cast<double>(length);
  coverage.bins.reserve(bin_count);
  for (std::size_t bin = 0U; bin < bin_count; ++bin) {
    const std::size_t start = bin * length / bin_count + 1U;
    const std::size_t end = (bin + 1U) * length / bin_count;
    std::uint64_t sum = 0U;
    for (std::size_t position = start; position <= end; ++position) {
      sum += depth[position];
    }
    const std::size_t span = end >= start ? end - start + 1U : 1U;
    coverage.bins.push_back({start, end, static_cast<std::size_t>((sum + span / 2U) / span)});
  }
  return coverage;
}

void assign_haplotype_clusters(const PrecisionCallResult &calls,
                               std::vector<ReadFeature> &features) {
  for (std::size_t index = 0U; index < features.size(); ++index) {
    const int haplotype =
        index < calls.molecule_haplotype.size() ? calls.molecule_haplotype[index] : -1;
    features[index].cluster_id = haplotype;
    features[index].outlier = haplotype < 0;
  }
}

void apply_compact_projection(const PrecisionCallResult &calls,
                              const std::string &reference,
                              std::vector<ReadFeature> &features,
                              std::map<std::string, SvCall> &svs,
                              std::map<std::string, ComplexSvCall> &complex_events) {
  static_cast<void>(reference);
  std::unordered_set<std::string> allowed_snps;
  std::unordered_set<std::string> allowed_markers;
  for (const auto &snv : calls.snvs) {
    if (snv.pass) {
      allowed_snps.insert(snp_key(snv.position, snv.reference, snv.alternate));
      allowed_markers.insert(std::to_string(snv.position) + std::string(1U, snv.alternate));
    }
  }
  for (const auto &indel : calls.indels) {
    if (!indel.pass) {
      continue;
    }
    if (indel.deletion) {
      for (std::size_t offset = 0U; offset < indel.length; ++offset) {
        allowed_markers.insert(
            std::to_string(fold(indel.right + offset, reference.size())) + "d");
      }
    } else {
      const auto normalized =
          normalize_circular_indel(reference, indel.right, 0U, indel.allele);
      allowed_markers.insert(std::to_string(normalized.right) + ".1" + normalized.right_allele);
    }
  }

  std::map<std::string, SvCall> kept;
  std::unordered_map<std::string, std::vector<std::string>> molecule_events;
  std::unordered_set<std::string> kept_ids;
  for (const auto &deletion : calls.deletions) {
    if (!deletion.pass) {
      continue;
    }
    SvCall sv;
    sv.id = "deletion:" + std::to_string(deletion.start) + "-" + std::to_string(deletion.end);
    sv.type = "deletion";
    sv.start = deletion.start;
    sv.end = deletion.end;
    sv.length = deletion.length;
    sv.supporting_reads = deletion.supporting_molecules;
    sv.known_event = deletion.known_common_deletion;
    sv.evidence_sources.emplace_back("consensus_junction");
    sv.orientations.emplace_back("+/+");
    sv.segment_count = 2U;
    for (const auto &molecule : deletion.supporting_molecules) {
      molecule_events[molecule].push_back(sv.id);
    }
    kept_ids.insert(sv.id);
    kept.emplace(sv.id, std::move(sv));
  }
  for (auto &[id, sv] : svs) {
    const bool rearrangement = sv.type == "inversion" || sv.type == "duplication" ||
                               sv.type == "ambiguous_adjacency" ||
                               (sv.type == "insertion" && sv.length >= kMinimumLargeDeletion);
    if (rearrangement && sv.supporting_reads.size() >= 3U && !kept.contains(id)) {
      kept_ids.insert(id);
      kept.emplace(id, std::move(sv));
    }
  }
  svs = std::move(kept);
  std::erase_if(complex_events, [](const auto &entry) {
    return entry.second.supporting_reads.size() < 3U;
  });

  for (auto &feature : features) {
    // MD and cs tags restate the read's errors and would dominate the result.
    std::erase_if(feature.aux_tags, [](const auto &entry) {
      constexpr std::size_t kLongestKeptTag = 32U;
      return entry.second.size() > kLongestKeptTag;
    });
    // Noisy reads have hundreds of tiny callable gaps. A marker falling in one is
    // simply not observed, so bridging them loses nothing at the cluster level.
    constexpr std::size_t kBridgedCallableGap = 25U;
    auto &ranges = feature.haplogroup_ranges;
    std::size_t written = 0U;
    for (std::size_t index = 0U; index < ranges.size(); ++index) {
      if (written != 0U && ranges[index].first <= ranges[written - 1U].second + kBridgedCallableGap) {
        ranges[written - 1U].second = std::max(ranges[written - 1U].second, ranges[index].second);
      } else {
        ranges[written++] = ranges[index];
      }
    }
    ranges.resize(written);
    std::erase_if(feature.snps, [&](const SnpCall &snp) {
      return !allowed_snps.contains(snp_key(snp.position, snp.reference, snp.alternate));
    });
    std::erase_if(feature.haplogroup_markers, [&](const ObservedPhyloMutation &marker) {
      return !allowed_markers.contains(marker.encoded);
    });
    std::erase_if(feature.sv_ids, [&](const std::string &id) { return !kept_ids.contains(id); });
    const auto events = molecule_events.find(feature.id);
    if (events != molecule_events.end()) {
      feature.sv_ids.insert(feature.sv_ids.end(), events->second.begin(), events->second.end());
      sort_unique_strings(feature.sv_ids);
    }
    std::erase_if(feature.complex_event_ids, [&](const std::string &id) {
      return !complex_events.contains(id);
    });
  }
}

} // namespace mito::detail
