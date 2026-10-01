#include "detail/pipeline.hpp"

namespace mito::detail {

[[nodiscard]] std::size_t effective_thread_count(std::size_t requested_threads,
                                                 std::size_t work_items) {
  if (work_items <= 1) {
    return 1;
  }
  const auto hardware_threads =
      std::max(1U, std::thread::hardware_concurrency());
  const std::size_t requested = requested_threads == 0 ? 1 : requested_threads;
  return std::max<std::size_t>(
      1, std::min<std::size_t>({requested, work_items,
                                static_cast<std::size_t>(hardware_threads)}));
}

[[nodiscard]] FeatureExtractionResult extract_read_feature(
    const ReadRecord &read, const std::string &reference,
    const AnalysisConfig &config,
    const std::map<std::string, SnpCall> &clinical_annotations) {
  throw_if_cancelled(config);
  FeatureExtractionResult result;
  auto &feature = result.feature;
  feature.id = read.id;
  feature.length = read.sequence.size();
  feature.mean_quality = mean_quality(read.qualities);
  feature.numt_evidence = numt_evidence(read, reference.size());
  feature.numt_score = numt_score(feature.numt_evidence);
  feature.filtered_numt =
      config.filter_numt && feature.numt_score > config.numt_threshold;
  feature.mapping_quality = read.mapping_quality;
  feature.flags = read.flags;
  feature.reference_name = read.reference_name;
  feature.aux_tags = read.aux_tags;

  if (!feature.filtered_numt) {
    feature.snps = call_snps_from_alignment(read, reference, config);
    if (feature.snps.empty() && read.cigar.empty()) {
      feature.snps = call_snps_from_reference_span(read, reference, config);
    }
    if (config.allow_development_tags) {
      merge_snps(feature.snps, snp_tags_from_header(read));
    }
    if (std::any_of(
            feature.snps.begin(), feature.snps.end(), [&](const auto &snp) {
              return snp.position == 0 || snp.position > reference.size();
            })) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "SNP coordinate is outside the reference for read '" +
                              read.id + "'");
    }
    apply_clinical_annotations(feature.snps, clinical_annotations);
    feature.haplogroup_markers.reserve(feature.snps.size());
    for (const auto &snp : feature.snps) {
      ObservedPhyloMutation mutation;
      mutation.position = snp.position;
      mutation.kind = PhyloMutationKind::substitution;
      mutation.alternate = snp.alternate;
      mutation.encoded =
          std::to_string(snp.position) + std::string(1, snp.alternate);
      feature.haplogroup_markers.push_back(std::move(mutation));
    }
    auto aligned_phylo_indels =
        call_phylo_indels_from_alignment(read, config, reference);
    std::vector<std::size_t> aligned_phylo_indel_positions;
    aligned_phylo_indel_positions.reserve(aligned_phylo_indels.size());
    for (const auto &mutation : aligned_phylo_indels) {
      aligned_phylo_indel_positions.push_back(mutation.position);
    }
    merge_haplogroup_markers(feature.haplogroup_markers,
                             std::move(aligned_phylo_indels));
    if (config.allow_development_tags) {
      merge_haplogroup_markers(feature.haplogroup_markers,
                               phylo_tags_from_header(read));
    }
    const auto callable_ranges =
        callable_phylo_ranges_from_alignment(read, reference, config);
    feature.haplogroup_range_known = callable_ranges.known;
    feature.haplogroup_ranges = callable_ranges.ranges;
    if (feature.haplogroup_range_known) {
      for (const auto position : aligned_phylo_indel_positions) {
        feature.haplogroup_ranges.emplace_back(position, position);
      }
      feature.haplogroup_ranges =
          normalize_ranges(std::move(feature.haplogroup_ranges));
    }
    if (config.allow_development_tags) {
      if (auto development_ranges = phylo_ranges_from_header(read)) {
        feature.haplogroup_range_known = true;
        feature.haplogroup_ranges = std::move(*development_ranges);
      }
    }
    if (config.allow_development_tags) {
      result.svs = sv_tags_from_header(read);
    }
    auto cigar_svs =
        parse_cigar_svs(read, config.sv_min_length, reference.size());
    auto split_svs =
        split_alignment_svs(read, reference.size(), config.sv_min_length);
    if (auto complex_event = coalesce_complex_sv(read.id, split_svs)) {
      feature.complex_event_ids.push_back(complex_event->id);
      result.complex_events.push_back(std::move(*complex_event));
    }
    if (!split_svs.empty()) {
      cigar_svs.erase(std::remove_if(cigar_svs.begin(), cigar_svs.end(),
                                     [](const auto &sv) {
                                       return sv.type == "soft_clip_left" ||
                                              sv.type == "soft_clip_right";
                                     }),
                      cigar_svs.end());
    }
    result.svs.insert(result.svs.end(),
                      std::make_move_iterator(cigar_svs.begin()),
                      std::make_move_iterator(cigar_svs.end()));
    result.svs.insert(result.svs.end(),
                      std::make_move_iterator(split_svs.begin()),
                      std::make_move_iterator(split_svs.end()));
    if (std::any_of(result.svs.begin(), result.svs.end(), [&](const auto &sv) {
          return sv.start == 0 || sv.end > reference.size();
        })) {
      throw AnalysisError(AnalysisErrorCode::internal_error,
                          "SV coordinate is outside the reference for read '" +
                              read.id + "'");
    }

    feature.sv_ids.reserve(result.svs.size());
    for (const auto &sv : result.svs) {
      feature.sv_ids.push_back(sv.id);
    }
    std::sort(feature.sv_ids.begin(), feature.sv_ids.end());
    feature.sv_ids.erase(
        std::unique(feature.sv_ids.begin(), feature.sv_ids.end()),
        feature.sv_ids.end());
    std::sort(feature.complex_event_ids.begin(),
              feature.complex_event_ids.end());
    feature.complex_event_ids.erase(
        std::unique(feature.complex_event_ids.begin(),
                    feature.complex_event_ids.end()),
        feature.complex_event_ids.end());
  }

  return result;
}

[[nodiscard]] std::vector<FeatureExtractionResult> extract_features_parallel(
    const std::vector<ReadRecord> &reads, const std::string &reference,
    const AnalysisConfig &config,
    const std::map<std::string, SnpCall> &clinical_annotations) {
  std::vector<FeatureExtractionResult> results(reads.size());
  if (reads.empty()) {
    return results;
  }

  const std::size_t worker_count =
      effective_thread_count(config.threads, reads.size());
  if (worker_count == 1) {
    for (std::size_t i = 0; i < reads.size(); ++i) {
      throw_if_cancelled(config);
      results[i] = extract_read_feature(reads[i], reference, config,
                                        clinical_annotations);
    }
    return results;
  }

  constexpr std::size_t kChunkSize = 32;
  std::atomic<std::size_t> next_index{0};
  std::atomic<bool> failed{false};
  std::vector<std::exception_ptr> errors(worker_count);
  std::vector<std::jthread> workers;
  workers.reserve(worker_count);

  for (std::size_t worker_id = 0; worker_id < worker_count; ++worker_id) {
    workers.emplace_back([&, worker_id] {
      try {
        while (!failed.load(std::memory_order_relaxed)) {
          throw_if_cancelled(config);
          const std::size_t begin =
              next_index.fetch_add(kChunkSize, std::memory_order_relaxed);
          if (begin >= reads.size()) {
            break;
          }
          const std::size_t end = std::min(reads.size(), begin + kChunkSize);
          for (std::size_t i = begin; i < end; ++i) {
            throw_if_cancelled(config);
            results[i] = extract_read_feature(reads[i], reference, config,
                                              clinical_annotations);
          }
        }
      } catch (...) {
        errors[worker_id] = std::current_exception();
        failed.store(true, std::memory_order_relaxed);
      }
    });
  }

  for (auto &worker : workers) {
    worker.join();
  }
  for (const auto &error : errors) {
    if (error) {
      std::rethrow_exception(error);
    }
  }

  return results;
}

} // namespace mito::detail
