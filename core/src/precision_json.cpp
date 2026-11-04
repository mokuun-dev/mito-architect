#include "detail/pipeline.hpp"
#include "detail/precision.hpp"

#include <charconv>

namespace mito::detail {
namespace {

constexpr std::size_t kMaximumListedMolecules = 200U;

class ObjectWriter {
public:
  explicit ObjectWriter(std::string &out) : out_(out) { out_.push_back('{'); }

  void key(const std::string_view name) {
    if (!first_) {
      out_.push_back(',');
    }
    first_ = false;
    out_.push_back('"');
    out_.append(name);
    out_.append("\":");
  }

  void number(const std::string_view name, const double value) {
    key(name);
    append_double(out_, value);
  }

  void integer(const std::string_view name, const std::size_t value) {
    key(name);
    out_.append(std::to_string(value));
  }

  void text(const std::string_view name, const std::string_view value) {
    key(name);
    out_.append(quoted(value));
  }

  void boolean(const std::string_view name, const bool value) {
    key(name);
    out_.append(value ? "true" : "false");
  }

  void close() { out_.push_back('}'); }

  static void append_double(std::string &out, const double value) {
    if (!std::isfinite(value)) {
      out.append("null");
      return;
    }
    std::array<char, 64> buffer{};
    const auto result = std::to_chars(buffer.data(), buffer.data() + buffer.size(),
                                      value, std::chars_format::general, 6);
    out.append(buffer.data(), result.ptr);
  }

private:
  std::string &out_;
  bool first_ = true;
};

void write_strings(std::string &out, const std::vector<std::string> &values,
                   const std::size_t limit = std::numeric_limits<std::size_t>::max()) {
  out.push_back('[');
  for (std::size_t index = 0U; index < values.size() && index < limit; ++index) {
    if (index != 0U) {
      out.push_back(',');
    }
    out.append(quoted(values[index]));
  }
  out.push_back(']');
}

void write_repeat(ObjectWriter &writer, std::string &out, const PrecisionRepeat &repeat) {
  writer.key("repeat");
  ObjectWriter object(out);
  object.integer("homopolymer_run", repeat.homopolymer_run);
  object.integer("unit", repeat.unit);
  object.integer("span", repeat.span);
  object.close();
}

void write_strands(ObjectWriter &writer, std::string &out,
                   const PrecisionStrandSupport &strands) {
  writer.key("strands");
  ObjectWriter object(out);
  object.integer("alt_forward", strands.alternate_forward);
  object.integer("alt_reverse", strands.alternate_reverse);
  object.integer("depth_forward", strands.depth_forward);
  object.integer("depth_reverse", strands.depth_reverse);
  object.close();
}

void write_effects(ObjectWriter &writer, std::string &out,
                   const std::vector<VariantEffect> &effects) {
  writer.key("effects");
  out.push_back('[');
  for (std::size_t index = 0U; index < effects.size(); ++index) {
    if (index != 0U) {
      out.push_back(',');
    }
    const auto &effect = effects[index];
    ObjectWriter object(out);
    if (!effect.gene.empty()) {
      object.text("gene", effect.gene);
    }
    object.text("biotype", effect.biotype);
    object.text("consequence", effect.consequence);
    if (!effect.hgvs_protein.empty()) {
      object.text("hgvs_p", effect.hgvs_protein);
    }
    if (!effect.codon_change.empty()) {
      object.text("codons", effect.codon_change);
    }
    if (effect.amino_acid_position != 0U) {
      object.integer("residue", effect.amino_acid_position);
    }
    if (effect.grantham) {
      object.integer("grantham", static_cast<std::size_t>(*effect.grantham));
    }
    if (effect.structure) {
      object.key("structure");
      ObjectWriter structure(out);
      structure.text("structure_id", effect.structure->pdb_id);
      structure.text("chain", effect.structure->chain);
      structure.integer("residue_index", effect.structure->residue);
      structure.text("complex", effect.structure->complex_name);
      structure.close();
    }
    object.close();
  }
  out.push_back(']');
}

void write_status(ObjectWriter &writer, std::string &out, const bool pass,
                  const std::vector<std::string> &filters,
                  const std::vector<std::string> &flags) {
  writer.text("status", pass ? "PASS" : "FILTERED");
  writer.key("filters");
  write_strings(out, filters);
  writer.key("flags");
  write_strings(out, flags);
}

void write_rates(ObjectWriter &writer, std::string &out, const std::string_view name,
                 const std::vector<std::pair<std::size_t, double>> &rates) {
  writer.key(name);
  out.push_back('[');
  for (std::size_t index = 0U; index < rates.size(); ++index) {
    if (index != 0U) {
      out.push_back(',');
    }
    ObjectWriter object(out);
    object.integer("run", rates[index].first);
    object.number("rate", rates[index].second);
    object.close();
  }
  out.push_back(']');
}

} // namespace

[[nodiscard]] std::string render_precision_calls_json(const PrecisionCallResult &result,
                                                      const AnalysisConfig &config) {
  static_cast<void>(config);
  std::string out;
  out.reserve(1024U + 512U * (result.snvs.size() + result.indels.size() + result.deletions.size()));
  ObjectWriter root(out);
  root.text("schema_version", "1.0");
  root.boolean("assessed", result.assessed);
  if (!result.assessed) {
    root.text("reason", "no eligible aligned molecules");
    root.close();
    return out;
  }
  root.key("method");
  {
    ObjectWriter method(out);
    method.number("false_discovery_rate", 0.01);
    method.text("snv_model", "context_beta_binomial_two_strand");
    method.text("indel_model", "repeat_length_error_model_pair_hmm_em");
    method.text("deletion_model", "junction_clustering_breakpoint_ml");
    method.close();
  }

  std::size_t pass_snvs = 0U;
  std::size_t pass_indels = 0U;
  std::size_t pass_deletions = 0U;
  for (const auto &call : result.snvs) {
    pass_snvs += call.pass ? 1U : 0U;
  }
  for (const auto &call : result.indels) {
    pass_indels += call.pass ? 1U : 0U;
  }
  for (const auto &call : result.deletions) {
    pass_deletions += call.pass ? 1U : 0U;
  }
  root.key("summary");
  {
    ObjectWriter summary(out);
    summary.integer("snvs_pass", pass_snvs);
    summary.integer("indels_pass", pass_indels);
    summary.integer("deletions_pass", pass_deletions);
    summary.close();
  }

  const auto &diagnostics = result.diagnostics;
  root.key("diagnostics");
  {
    ObjectWriter object(out);
    object.integer("molecules_used", diagnostics.molecules_used);
    object.integer("molecules_excluded_divergent", diagnostics.molecules_excluded_divergent);
    object.integer("reference_positions_covered", diagnostics.reference_positions_covered);
    object.number("median_depth", diagnostics.median_depth);
    object.number("mean_substitution_rate", diagnostics.mean_substitution_rate);
    object.number("median_molecule_divergence", diagnostics.median_molecule_divergence);
    object.number("divergence_limit", diagnostics.divergence_limit);
    object.integer("snv_hypotheses", diagnostics.snv_hypotheses);
    object.integer("indel_candidates", diagnostics.indel_candidates);
    object.integer("indels_rescored", diagnostics.indels_rescored);
    object.number("indel_noise_p99", diagnostics.indel_noise_p99);
    object.integer("junction_molecules", diagnostics.junction_molecules);
    object.integer("singleton_junctions", diagnostics.singleton_junctions);
    write_rates(object, out, "homopolymer_deletion_rate", diagnostics.homopolymer_deletion_rate);
    write_rates(object, out, "homopolymer_insertion_rate", diagnostics.homopolymer_insertion_rate);
    object.close();
  }

  if (result.consensus_haplogroup) {
    const auto &haplogroup = *result.consensus_haplogroup;
    root.key("consensus_haplogroup");
    ObjectWriter object(out);
    object.text("best", haplogroup.best);
    object.number("quality", haplogroup.quality);
    object.integer("observed_markers", haplogroup.observed_markers.size());
    object.key("candidates");
    out.push_back('[');
    const std::size_t listed = std::min<std::size_t>(haplogroup.candidates.size(), 5U);
    for (std::size_t index = 0U; index < listed; ++index) {
      if (index != 0U) {
        out.push_back(',');
      }
      ObjectWriter candidate(out);
      candidate.text("name", haplogroup.candidates[index].name);
      candidate.number("score", haplogroup.candidates[index].score);
      candidate.close();
    }
    out.push_back(']');
    object.close();
  }

  root.key("haplotypes");
  out.push_back('[');
  for (std::size_t index = 0U; index < result.haplotypes.size(); ++index) {
    if (index != 0U) {
      out.push_back(',');
    }
    const auto &haplotype = result.haplotypes[index];
    ObjectWriter object(out);
    object.integer("id", static_cast<std::size_t>(haplotype.id));
    object.number("fraction", haplotype.fraction);
    object.integer("molecules", haplotype.molecules);
    object.key("defining_variants");
    write_strings(out, haplotype.defining_variants);
    if (haplotype.haplogroup) {
      object.text("haplogroup", haplotype.haplogroup->best);
      object.number("haplogroup_quality", haplotype.haplogroup->quality);
    }
    object.close();
  }
  out.push_back(']');

  root.key("callable_ranges");
  out.push_back('[');
  for (std::size_t index = 0U; index < result.callable_ranges.size(); ++index) {
    if (index != 0U) {
      out.push_back(',');
    }
    out.append("[" + std::to_string(result.callable_ranges[index].first) + "," +
               std::to_string(result.callable_ranges[index].second) + "]");
  }
  out.push_back(']');

  root.key("snvs");
  out.push_back('[');
  for (std::size_t index = 0U; index < result.snvs.size(); ++index) {
    if (index != 0U) {
      out.push_back(',');
    }
    const auto &call = result.snvs[index];
    ObjectWriter object(out);
    object.text("id", "snv:" + std::to_string(call.position) + ":" +
                          std::string(1U, call.reference) + ":" +
                          std::string(1U, call.alternate));
    object.integer("position", call.position);
    object.text("ref", std::string(1U, call.reference));
    object.text("alt", std::string(1U, call.alternate));
    write_status(object, out, call.pass, call.filters, call.flags);
    object.number("fraction", call.fraction);
    object.number("ci95_low", call.ci95_low);
    object.number("ci95_high", call.ci95_high);
    object.integer("depth", call.depth);
    object.integer("alt_depth", call.alternate_depth);
    write_strands(object, out, call.strands);
    object.number("p_value", call.p_value);
    object.number("q_value", call.q_value);
    object.number("background_rate", call.background_rate);
    object.number("alt_mean_quality", call.alternate_mean_quality);
    write_repeat(object, out, call.repeat);
    write_effects(object, out, call.effects);
    object.close();
  }
  out.push_back(']');

  root.key("indels");
  out.push_back('[');
  for (std::size_t index = 0U; index < result.indels.size(); ++index) {
    if (index != 0U) {
      out.push_back(',');
    }
    const auto &call = result.indels[index];
    ObjectWriter object(out);
    object.text("id", call.id);
    object.text("type", call.deletion ? "deletion" : "insertion");
    object.integer("position", call.left);
    object.integer("position_3prime", call.right);
    object.integer("length", call.length);
    object.text("allele", call.allele);
    object.integer("vcf_position", call.vcf_position);
    object.text("vcf_ref", call.vcf_reference);
    object.text("vcf_alt", call.vcf_alternate);
    object.text("hgvs", call.hgvs);
    write_status(object, out, call.pass, call.filters, call.flags);
    object.number("fraction", call.fraction);
    object.number("ci95_low", call.ci95_low);
    object.number("ci95_high", call.ci95_high);
    object.integer("depth", call.depth);
    object.integer("raw_alt_depth", call.raw_alternate_depth);
    object.number("effective_alt_depth", call.effective_alternate_depth);
    object.number("effective_depth", call.effective_depth);
    write_strands(object, out, call.strands);
    object.number("lod", call.lod);
    object.boolean("rescored", call.rescored);
    object.number("p_value", call.p_value);
    object.number("q_value", call.q_value);
    object.number("background_rate", call.background_rate);
    write_repeat(object, out, call.repeat);
    write_effects(object, out, call.effects);
    object.integer("supporting_molecule_count", call.supporting_molecules.size());
    object.key("supporting_molecules");
    write_strings(out, call.supporting_molecules, kMaximumListedMolecules);
    object.close();
  }
  out.push_back(']');

  root.key("deletions");
  out.push_back('[');
  for (std::size_t index = 0U; index < result.deletions.size(); ++index) {
    if (index != 0U) {
      out.push_back(',');
    }
    const auto &call = result.deletions[index];
    ObjectWriter object(out);
    object.text("id", call.id);
    object.integer("start", call.start);
    object.integer("end", call.end);
    object.integer("left_aligned_start", call.left_start);
    object.integer("left_aligned_end", call.left_end);
    object.integer("length", call.length);
    object.text("microhomology", call.microhomology);
    write_status(object, out, call.pass, call.filters, call.flags);
    object.key("fraction");
    if (call.fraction) {
      ObjectWriter::append_double(out, *call.fraction);
    } else {
      out.append("null");
    }
    object.number("ci95_low", call.ci95_low);
    object.number("ci95_high", call.ci95_high);
    object.integer("junction_molecules", call.junction_molecules);
    object.integer("junction_forward", call.junction_forward);
    object.integer("junction_reverse", call.junction_reverse);
    object.integer("wild_type_left", call.wild_type_left);
    object.integer("wild_type_right", call.wild_type_right);
    object.integer("breakpoint_variants", call.breakpoint_variants);
    object.key("depth_ratio");
    if (call.depth_ratio) {
      ObjectWriter::append_double(out, *call.depth_ratio);
    } else {
      out.append("null");
    }
    object.number("p_value", call.p_value);
    object.number("q_value", call.q_value);
    object.boolean("known_common_deletion", call.known_common_deletion);
    object.key("impact");
    {
      ObjectWriter impact(out);
      impact.key("genes_lost");
      write_strings(out, call.impact.genes_lost);
      impact.key("genes_partial");
      write_strings(out, call.impact.genes_partial);
      impact.integer("trnas_affected", call.impact.transfer_rnas_affected);
      if (!call.impact.fusion.empty()) {
        impact.text("fusion", call.impact.fusion);
      }
      impact.boolean("heavy_strand_origin_lost", call.impact.heavy_strand_origin_lost);
      impact.boolean("light_strand_origin_lost", call.impact.light_strand_origin_lost);
      impact.boolean("major_arc", call.impact.major_arc);
      impact.boolean("contact_zone", call.impact.contact_zone);
      impact.close();
    }
    object.integer("supporting_molecule_count", call.supporting_molecules.size());
    object.key("supporting_molecules");
    write_strings(out, call.supporting_molecules, kMaximumListedMolecules);
    object.close();
  }
  out.push_back(']');
  root.close();
  return out;
}

} // namespace mito::detail
