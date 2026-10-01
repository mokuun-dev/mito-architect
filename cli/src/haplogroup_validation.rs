use crate::args::ValidateHaplogroupManifestArgs;
use crate::contract::{require_array, validate_result_contract};
use crate::storage::sha256_file;
use anyhow::Context;
use anyhow::Result;
use mito_ffi::AnalyzeOptions;
use mito_ffi::MitoEngine;
use serde_json::Value;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::fs;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

#[derive(Debug)]
pub(super) struct HsdProfile {
    pub(super) ranges: BTreeSet<String>,
    pub(super) markers: BTreeSet<String>,
}

pub(super) fn canonical_hsd_marker(sample: &str, raw: &str) -> Result<String> {
    let marker = raw.trim();
    if marker.is_empty() || marker.chars().any(char::is_whitespace) {
        anyhow::bail!("HSD profile '{sample}' has invalid marker {marker:?}");
    }
    let bytes = marker.as_bytes();
    if bytes.len() > 3 && bytes[bytes.len() - 3..].eq_ignore_ascii_case(b"DEL") {
        let position = marker
            .get(..bytes.len() - 3)
            .with_context(|| format!("HSD profile '{sample}' has invalid deletion {marker:?}"))?;
        let parsed = position
            .parse::<u64>()
            .with_context(|| format!("HSD profile '{sample}' has invalid deletion {marker:?}"))?;
        if parsed == 0 || parsed > 16_569 {
            anyhow::bail!("HSD profile '{sample}' has invalid deletion {marker:?}");
        }
        return Ok(format!("{parsed}d"));
    }
    Ok(marker.to_owned())
}

pub(super) fn canonical_hsd_ranges(sample: &str, raw: &str) -> Result<BTreeSet<String>> {
    let mut ranges = BTreeSet::new();
    let unquoted = raw.trim().trim_matches('"');
    for raw_range in unquoted.split(';') {
        let compact: String = raw_range
            .chars()
            .filter(|value| !value.is_whitespace())
            .collect();
        if compact.is_empty() {
            continue;
        }
        let (start, end) = compact.split_once('-').unwrap_or((&compact, &compact));
        let start = start
            .parse::<u64>()
            .with_context(|| format!("HSD profile '{sample}' has invalid range {compact:?}"))?;
        let end = end
            .parse::<u64>()
            .with_context(|| format!("HSD profile '{sample}' has invalid range {compact:?}"))?;
        if start == 0 || end == 0 || start > 16_569 || end > 16_569 {
            anyhow::bail!("HSD profile '{sample}' has invalid range {compact:?}");
        }
        if start <= end {
            ranges.insert(format!("{start}-{end}"));
        } else {
            ranges.insert(format!("1-{end}"));
            ranges.insert(format!("{start}-16569"));
        }
    }
    if ranges.is_empty() {
        anyhow::bail!("HSD profile '{sample}' has no tested range");
    }
    Ok(ranges)
}

pub(super) fn load_hsd_profiles(path: &Path) -> Result<BTreeMap<String, HsdProfile>> {
    let input = fs::File::open(path)
        .with_context(|| format!("failed to open pinned HSD input {}", path.display()))?;
    let mut profiles = BTreeMap::new();
    for (line_index, line) in BufReader::new(input).lines().enumerate() {
        let line = line.with_context(|| {
            format!(
                "failed to read pinned HSD input {} at line {}",
                path.display(),
                line_index + 1
            )
        })?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            anyhow::bail!(
                "pinned HSD input {} line {} requires sample, range, metadata, and markers",
                path.display(),
                line_index + 1
            );
        }
        let sample = fields[0].trim();
        if sample.is_empty() || sample.chars().any(char::is_whitespace) {
            anyhow::bail!(
                "pinned HSD input {} line {} has invalid sample identifier",
                path.display(),
                line_index + 1
            );
        }
        let mut markers = BTreeSet::new();
        for raw_marker in &fields[3..] {
            if raw_marker.trim().is_empty() {
                continue;
            }
            let marker = canonical_hsd_marker(sample, raw_marker)?;
            if !markers.insert(marker.clone()) {
                anyhow::bail!("HSD profile '{sample}' repeats marker {marker}");
            }
        }
        if markers.is_empty() {
            anyhow::bail!("HSD profile '{sample}' has no markers");
        }
        let profile = HsdProfile {
            ranges: canonical_hsd_ranges(sample, fields[1])?,
            markers,
        };
        if profiles.insert(sample.to_owned(), profile).is_some() {
            anyhow::bail!("duplicate HSD profile: {sample}");
        }
    }
    if profiles.is_empty() {
        anyhow::bail!("pinned HSD input is empty: {}", path.display());
    }
    Ok(profiles)
}

pub(super) fn validate_haplogroup_manifest(args: ValidateHaplogroupManifestArgs) -> Result<()> {
    if !args.manifest.is_file() {
        anyhow::bail!(
            "haplogroup manifest does not exist: {}",
            args.manifest.display()
        );
    }
    let manifest: Value = serde_json::from_str(
        &fs::read_to_string(&args.manifest)
            .with_context(|| format!("failed to read {}", args.manifest.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.manifest.display()))?;
    if manifest.get("schema_version").and_then(Value::as_str) != Some("1.0") {
        anyhow::bail!("haplogroup manifest requires schema_version 1.0");
    }
    let reference = manifest
        .get("reference")
        .and_then(Value::as_object)
        .context("haplogroup manifest requires a reference object")?;
    for field in [
        "tool",
        "version",
        "tree",
        "asset_sha256",
        "tree_sha256",
        "weights_sha256",
        "alignment_rules_sha256",
        "source_hsd",
        "source_hsd_sha256",
    ] {
        if reference.get(field).and_then(Value::as_str).is_none() {
            anyhow::bail!("haplogroup manifest reference requires string {field}");
        }
    }
    if reference.get("tool").and_then(Value::as_str) != Some("HaploGrep 3")
        || reference.get("version").and_then(Value::as_str) != Some("3.3.2")
        || reference.get("tree").and_then(Value::as_str) != Some("phylotree-rcrs@17.3")
    {
        anyhow::bail!("haplogroup manifest must pin HaploGrep 3.3.2 with phylotree-rcrs@17.3");
    }
    for field in [
        "asset_sha256",
        "tree_sha256",
        "weights_sha256",
        "alignment_rules_sha256",
        "source_hsd_sha256",
    ] {
        let digest = reference
            .get(field)
            .and_then(Value::as_str)
            .context("validated reference digest must be present")?;
        if digest.len() != 64 || !digest.chars().all(|value| value.is_ascii_hexdigit()) {
            anyhow::bail!("haplogroup manifest reference {field} is not a SHA-256 digest");
        }
    }
    let base = args.manifest.parent().unwrap_or_else(|| Path::new("."));
    let source_hsd = base.join(
        reference
            .get("source_hsd")
            .and_then(Value::as_str)
            .context("source_hsd must be present")?,
    );
    if !source_hsd.is_file() {
        anyhow::bail!(
            "pinned HaploGrep input is missing: {}",
            source_hsd.display()
        );
    }
    let source_digest = sha256_file(&source_hsd)?;
    let expected_source_digest = reference
        .get("source_hsd_sha256")
        .and_then(Value::as_str)
        .context("source_hsd_sha256 must be present")?;
    if source_digest != expected_source_digest {
        anyhow::bail!(
            "pinned HaploGrep input checksum mismatch: expected {expected_source_digest}, got {source_digest}"
        );
    }
    let hsd_profiles = load_hsd_profiles(&source_hsd)?;

    let cases = manifest
        .get("cases")
        .and_then(Value::as_array)
        .context("haplogroup manifest requires a cases array")?;
    if cases.is_empty() {
        anyhow::bail!("haplogroup manifest cases array must not be empty");
    }
    let coverage = manifest
        .get("coverage")
        .and_then(Value::as_object)
        .context("haplogroup manifest requires a coverage object")?;
    let expected_differential_profiles = coverage
        .get("differential_profiles")
        .and_then(Value::as_u64)
        .context("haplogroup manifest coverage requires differential_profiles")?;
    let expected_alignment_cases = coverage
        .get("alignment_projection_cases")
        .and_then(Value::as_u64)
        .context("haplogroup manifest coverage requires alignment_projection_cases")?;
    let expected_compound_rule_cases = coverage
        .get("compound_rule_cases")
        .and_then(Value::as_u64)
        .context("haplogroup manifest coverage requires compound_rule_cases")?;
    let expected_backmutation_cases = coverage
        .get("backmutation_cases")
        .and_then(Value::as_u64)
        .context("haplogroup manifest coverage requires backmutation_cases")?;
    let expected_lineage_groups = coverage
        .get("lineage_groups")
        .and_then(Value::as_array)
        .context("haplogroup manifest coverage requires lineage_groups")?
        .iter()
        .map(|value| {
            value
                .as_str()
                .map(str::to_owned)
                .context("haplogroup manifest coverage has a non-string lineage group")
        })
        .collect::<Result<BTreeSet<_>>>()?;
    if expected_lineage_groups.is_empty() {
        anyhow::bail!("haplogroup manifest coverage has no lineage groups");
    }

    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    let mut case_names = BTreeSet::new();
    let mut referenced_hsd_profiles = BTreeSet::new();
    let mut differential_cases = 0usize;
    let mut alignment_cases = 0usize;
    let mut compound_rule_cases = 0usize;
    let mut backmutation_cases = 0usize;
    let mut observed_lineage_groups = BTreeSet::new();
    for (index, case) in cases.iter().enumerate() {
        let object = case
            .as_object()
            .with_context(|| format!("haplogroup case {index} must be an object"))?;
        let name = object
            .get("name")
            .and_then(Value::as_str)
            .with_context(|| format!("haplogroup case {index} requires string name"))?;
        if !case_names.insert(name.to_owned()) {
            anyhow::bail!("duplicate haplogroup case name: {name}");
        }
        let source_profile = object
            .get("source_profile")
            .map(|value| {
                value.as_str().with_context(|| {
                    format!("haplogroup case '{name}' has a non-string source_profile")
                })
            })
            .transpose()?;
        let (expected_markers, expected_ranges) = if let Some(source_profile) = source_profile {
            if object.contains_key("markers") || object.contains_key("range") {
                anyhow::bail!(
                    "haplogroup case '{name}' must not duplicate markers/range from source_profile"
                );
            }
            if !referenced_hsd_profiles.insert(source_profile.to_owned()) {
                anyhow::bail!("HSD profile '{source_profile}' is referenced more than once");
            }
            let profile = hsd_profiles.get(source_profile).with_context(|| {
                format!(
                    "haplogroup case '{name}' references missing HSD profile '{source_profile}'"
                )
            })?;
            (profile.markers.clone(), profile.ranges.clone())
        } else {
            let marker_values = object
                .get("markers")
                .and_then(Value::as_array)
                .with_context(|| format!("haplogroup case '{name}' requires markers"))?;
            if marker_values.is_empty() {
                anyhow::bail!("haplogroup case '{name}' has no markers");
            }
            let mut markers = BTreeSet::new();
            for marker in marker_values {
                let marker = marker
                    .as_str()
                    .with_context(|| format!("haplogroup case '{name}' has a non-string marker"))?;
                if marker.is_empty() || marker.chars().any(char::is_whitespace) {
                    anyhow::bail!("haplogroup case '{name}' has an invalid marker: {marker:?}");
                }
                if !markers.insert(marker.to_owned()) {
                    anyhow::bail!("haplogroup case '{name}' repeats marker {marker}");
                }
            }
            let range_values = object
                .get("range")
                .and_then(Value::as_array)
                .with_context(|| format!("haplogroup case '{name}' requires range"))?;
            if range_values.is_empty() {
                anyhow::bail!("haplogroup case '{name}' has no tested range");
            }
            let mut ranges = BTreeSet::new();
            for range in range_values {
                let range = range
                    .as_str()
                    .with_context(|| format!("haplogroup case '{name}' has a non-string range"))?;
                let (start, end) = range.split_once('-').unwrap_or((range, range));
                let start = start.parse::<u64>().with_context(|| {
                    format!("haplogroup case '{name}' has invalid range {range:?}")
                })?;
                let end = end.parse::<u64>().with_context(|| {
                    format!("haplogroup case '{name}' has invalid range {range:?}")
                })?;
                if start == 0 || end < start || end > 16_569 {
                    anyhow::bail!("haplogroup case '{name}' has invalid range {range:?}");
                }
                ranges.insert(format!("{start}-{end}"));
            }
            (markers, ranges)
        };
        let expected_candidates = object
            .get("haplogrep_candidates")
            .map(|value| {
                let values = value.as_array().with_context(|| {
                    format!("haplogroup case '{name}' has invalid haplogrep_candidates")
                })?;
                values
                    .iter()
                    .map(|candidate| {
                        candidate.as_str().map(str::to_owned).with_context(|| {
                            format!("haplogroup case '{name}' has a non-string reference candidate")
                        })
                    })
                    .collect::<Result<Vec<_>>>()
            })
            .transpose()?;
        let backmutation_absent_markers = object
            .get("backmutation_absent_markers")
            .map(|value| {
                value
                    .as_array()
                    .with_context(|| {
                        format!("haplogroup case '{name}' has invalid backmutation_absent_markers")
                    })?
                    .iter()
                    .map(|marker| {
                        marker.as_str().map(str::to_owned).with_context(|| {
                            format!("haplogroup case '{name}' has a non-string backmutation marker")
                        })
                    })
                    .collect::<Result<BTreeSet<_>>>()
            })
            .transpose()?
            .unwrap_or_default();
        if !backmutation_absent_markers.is_empty() {
            if source_profile.is_none() {
                anyhow::bail!("haplogroup backmutation case '{name}' requires a source_profile");
            }
            for marker in &backmutation_absent_markers {
                if expected_markers.contains(marker) {
                    anyhow::bail!(
                        "haplogroup backmutation case '{name}' still contains removed marker {marker}"
                    );
                }
            }
            backmutation_cases += 1;
        }
        if expected_candidates.as_ref().is_some_and(Vec::is_empty) {
            anyhow::bail!("haplogroup case '{name}' has no reference top hit");
        }
        if let Some(candidates) = &expected_candidates {
            if candidates.len() != 3 {
                anyhow::bail!(
                    "haplogroup case '{name}' must pin exactly three HaploGrep candidates"
                );
            }
            let quality = object
                .get("haplogrep_quality")
                .and_then(Value::as_f64)
                .with_context(|| {
                    format!("haplogroup case '{name}' requires numeric haplogrep_quality")
                })?;
            if !quality.is_finite() || quality <= 0.0 {
                anyhow::bail!("haplogroup case '{name}' has invalid haplogrep_quality");
            }
            let lineage_group = object
                .get("lineage_group")
                .and_then(Value::as_str)
                .with_context(|| {
                    format!("haplogroup case '{name}' requires string lineage_group")
                })?;
            if lineage_group.is_empty() || lineage_group.chars().any(char::is_whitespace) {
                anyhow::bail!("haplogroup case '{name}' has invalid lineage_group");
            }
            observed_lineage_groups.insert(lineage_group.to_owned());
        }

        let input = object.get("input").and_then(Value::as_str);
        if source_profile.is_some() && input.is_some() {
            anyhow::bail!("haplogroup case '{name}' cannot combine source_profile with input");
        }
        if source_profile.is_some() && expected_candidates.is_none() {
            anyhow::bail!("haplogroup case '{name}' source_profile requires reference candidates");
        }
        if input.is_none() && expected_candidates.is_none() {
            anyhow::bail!("haplogroup case '{name}' requires either input or haplogrep_candidates");
        }
        let rule_class = object
            .get("rule_class")
            .map(|value| {
                value
                    .as_str()
                    .with_context(|| format!("haplogroup case '{name}' has invalid rule_class"))
            })
            .transpose()?;
        if input.is_some() {
            if !matches!(rule_class, Some("cigar_3prime" | "compound")) {
                anyhow::bail!(
                    "haplogroup alignment case '{name}' requires rule_class cigar_3prime or compound"
                );
            }
            let normalization_rules = object
                .get("normalization_rules")
                .and_then(Value::as_array)
                .with_context(|| {
                    format!("haplogroup alignment case '{name}' requires normalization_rules")
                })?;
            if normalization_rules.is_empty()
                || normalization_rules.iter().any(|value| {
                    value
                        .as_str()
                        .is_none_or(|rule| rule.trim().is_empty() || !rule.contains("->"))
                })
            {
                anyhow::bail!("haplogroup alignment case '{name}' has invalid normalization_rules");
            }
        } else if rule_class.is_some() {
            anyhow::bail!("haplogroup differential case '{name}' must not set rule_class");
        }
        let temporary_fixture = input.is_none();
        let fixture_path = if let Some(input) = input {
            base.join(input)
        } else {
            std::env::temp_dir().join(format!(
                "mito-haplogroup-{}-{index}.fastq",
                std::process::id()
            ))
        };
        if temporary_fixture {
            let mut fastq = format!(
                "@{name} phylo_range={}",
                expected_ranges
                    .iter()
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(",")
            );
            for marker in &expected_markers {
                fastq.push_str(" phylo=");
                fastq.push_str(marker);
            }
            fastq.push_str("\nA\n+\nI\n");
            fs::write(&fixture_path, fastq)
                .with_context(|| format!("failed to write {}", fixture_path.display()))?;
        } else if !fixture_path.is_file() {
            anyhow::bail!(
                "haplogroup case '{name}' input does not exist: {}",
                fixture_path.display()
            );
        }
        let analysis = engine.analyze_with_options(
            &fixture_path,
            None,
            AnalyzeOptions {
                filter_numt: false,
                threads: 1,
                allow_development_tags: true,
                ..AnalyzeOptions::default()
            },
        );
        if temporary_fixture {
            let _ = fs::remove_file(&fixture_path);
        }
        let json = analysis.with_context(|| format!("haplogroup case '{name}' failed"))?;
        let data: Value = serde_json::from_str(&json)
            .with_context(|| format!("haplogroup case '{name}' returned invalid JSON"))?;
        validate_result_contract(&data)?;
        let clusters = require_array(&data, "/clusters")?;
        if clusters.len() != 1 {
            anyhow::bail!(
                "haplogroup case '{name}' expected one cluster, got {}",
                clusters.len()
            );
        }
        let assignment = clusters[0]
            .get("haplogroup_assignment")
            .and_then(Value::as_object)
            .with_context(|| format!("haplogroup case '{name}' has no assignment object"))?;
        let actual_best = clusters[0]
            .get("haplogroup")
            .and_then(Value::as_str)
            .with_context(|| format!("haplogroup case '{name}' has no best assignment"))?;
        if let Some(expected_best) = expected_candidates
            .as_ref()
            .and_then(|candidates| candidates.first())
        {
            if actual_best != expected_best {
                anyhow::bail!(
                    "HaploGrep differential mismatch for '{name}': expected {expected_best}, got {actual_best}"
                );
            }
        }
        if let Some(expected_candidates) = &expected_candidates {
            let candidate_values = assignment
                .get("candidates")
                .and_then(Value::as_array)
                .with_context(|| format!("haplogroup case '{name}' has no candidates array"))?;
            if !backmutation_absent_markers.is_empty() {
                let winning_candidate = candidate_values.first().with_context(|| {
                    format!("haplogroup backmutation case '{name}' has no winning candidate")
                })?;
                for field in ["matched", "missing"] {
                    let values = winning_candidate
                        .get(field)
                        .and_then(Value::as_array)
                        .with_context(|| {
                            format!(
                                "haplogroup backmutation case '{name}' winner has no {field} array"
                            )
                        })?;
                    for marker in values.iter().filter_map(Value::as_str) {
                        if backmutation_absent_markers.contains(marker) {
                            anyhow::bail!(
                                "haplogroup backmutation case '{name}' emitted removed marker {marker} in {field}"
                            );
                        }
                    }
                }
            }
            let actual_candidates = candidate_values
                .iter()
                .map(|candidate| {
                    candidate
                        .get("name")
                        .and_then(Value::as_str)
                        .map(str::to_owned)
                        .with_context(|| {
                            format!("haplogroup case '{name}' emitted a candidate without name")
                        })
                })
                .collect::<Result<Vec<_>>>()?;
            if &actual_candidates != expected_candidates {
                anyhow::bail!(
                    "HaploGrep top-3 mismatch for '{name}': expected {:?}, got {:?}",
                    expected_candidates,
                    actual_candidates
                );
            }
        }
        let actual_markers = assignment
            .get("observed_markers")
            .and_then(Value::as_array)
            .with_context(|| format!("haplogroup case '{name}' has no observed_markers"))?
            .iter()
            .map(|marker| {
                marker.as_str().map(str::to_owned).with_context(|| {
                    format!("haplogroup case '{name}' emitted a non-string marker")
                })
            })
            .collect::<Result<BTreeSet<_>>>()?;
        if actual_markers != expected_markers {
            anyhow::bail!(
                "haplogroup marker projection mismatch for '{name}': expected {:?}, got {:?}",
                expected_markers,
                actual_markers
            );
        }
        let actual_ranges = assignment
            .get("callable_ranges")
            .and_then(Value::as_array)
            .with_context(|| format!("haplogroup case '{name}' has no callable_ranges"))?
            .iter()
            .map(|range| {
                let start = range
                    .get("start")
                    .and_then(Value::as_u64)
                    .with_context(|| {
                        format!("haplogroup case '{name}' emitted a range without start")
                    })?;
                let end = range.get("end").and_then(Value::as_u64).with_context(|| {
                    format!("haplogroup case '{name}' emitted a range without end")
                })?;
                Ok(format!("{start}-{end}"))
            })
            .collect::<Result<BTreeSet<_>>>()?;
        if actual_ranges != expected_ranges {
            anyhow::bail!(
                "haplogroup callable-range mismatch for '{name}': expected {:?}, got {:?}",
                expected_ranges,
                actual_ranges
            );
        }
        if expected_candidates.is_some() {
            differential_cases += 1;
            println!("haplogroup differential passed: {name} -> {actual_best}");
        } else {
            alignment_cases += 1;
            if rule_class == Some("compound") {
                compound_rule_cases += 1;
            }
            println!("haplogroup alignment projection passed: {name}");
        }
    }

    let available_hsd_profiles = hsd_profiles.keys().cloned().collect::<BTreeSet<_>>();
    if referenced_hsd_profiles != available_hsd_profiles {
        let missing = available_hsd_profiles
            .difference(&referenced_hsd_profiles)
            .cloned()
            .collect::<Vec<_>>();
        let unknown = referenced_hsd_profiles
            .difference(&available_hsd_profiles)
            .cloned()
            .collect::<Vec<_>>();
        anyhow::bail!(
            "haplogroup manifest/HSD profile mismatch: unreferenced {:?}, unknown {:?}",
            missing,
            unknown
        );
    }
    if u64::try_from(differential_cases).context("differential case count overflow")?
        != expected_differential_profiles
        || u64::try_from(alignment_cases).context("alignment case count overflow")?
            != expected_alignment_cases
        || u64::try_from(compound_rule_cases).context("compound rule case count overflow")?
            != expected_compound_rule_cases
        || u64::try_from(backmutation_cases).context("backmutation case count overflow")?
            != expected_backmutation_cases
    {
        anyhow::bail!(
            "haplogroup coverage count mismatch: expected {expected_differential_profiles} differential/{expected_alignment_cases} alignment/{expected_compound_rule_cases} compound/{expected_backmutation_cases} backmutation, got {differential_cases}/{alignment_cases}/{compound_rule_cases}/{backmutation_cases}"
        );
    }
    if observed_lineage_groups != expected_lineage_groups {
        anyhow::bail!(
            "haplogroup lineage coverage mismatch: expected {:?}, got {:?}",
            expected_lineage_groups,
            observed_lineage_groups
        );
    }

    println!(
        "haplogroup manifest validation passed: {differential_cases} HaploGrep differential cases, {alignment_cases} alignment-projection cases ({compound_rule_cases} compound), and {backmutation_cases} explicit backmutation cases"
    );
    Ok(())
}
