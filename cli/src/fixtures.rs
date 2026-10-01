use crate::args::{
    ValidateClinicalFixtureArgs, ValidateClinicalManifestArgs, ValidateErrorManifestArgs,
    ValidateEvidenceFixtureArgs, ValidateEvidenceGraphFixtureArgs, ValidateFixtureArgs,
    ValidateSvFixtureArgs,
};
use crate::clinical::validate_clinical_tsv_schema;
use crate::contract::{require_array, validate_result_contract};
use crate::exports::render_vcf;
use crate::projection::{project_evidence, project_evidence_graph};
use crate::storage::{normalize_text, write_json_atomic};
use anyhow::Context;
use anyhow::Result;
use mito_ffi::AnalyzeOptions;
use mito_ffi::MitoEngine;
use serde_json::Value;
use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

pub(super) fn validate_fixture(args: ValidateFixtureArgs) -> Result<()> {
    if !args.input.exists() {
        anyhow::bail!("fixture input does not exist: {}", args.input.display());
    }
    if !args.expected_vcf.exists() {
        anyhow::bail!(
            "expected VCF does not exist: {}",
            args.expected_vcf.display()
        );
    }

    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    let json = engine
        .analyze_with_options(
            &args.input,
            None,
            AnalyzeOptions {
                filter_numt: true,
                threads: 1,
                ..AnalyzeOptions::default()
            },
        )
        .context("analysis failed")?;
    let data: Value = serde_json::from_str(&json).context("analysis returned invalid JSON")?;
    assert_json_u64(
        &data,
        "/filter_stats/passed_reads",
        args.expected_passed,
        "passed read count",
    )?;
    assert_json_u64(
        &data,
        "/filter_stats/numt_filtered_reads",
        args.expected_numt,
        "NUMT filtered read count",
    )?;
    for expected_snp in &args.expected_snp {
        assert_snp_present(&data, expected_snp)?;
    }
    for expected_sv in &args.expected_sv {
        assert_sv_present(&data, expected_sv)?;
    }
    for expected_mapq in args.expected_mapq {
        assert_mapq_present(&data, expected_mapq)?;
    }
    for expected_aux in &args.expected_aux {
        assert_aux_present(&data, expected_aux)?;
    }

    let actual_vcf = normalize_text(&render_vcf(&json)?);
    let expected_vcf = normalize_text(
        &fs::read_to_string(&args.expected_vcf)
            .with_context(|| format!("failed to read {}", args.expected_vcf.display()))?,
    );
    if actual_vcf != expected_vcf {
        anyhow::bail!(
            "VCF mismatch for {}\nexpected:\n{}\nactual:\n{}",
            args.input.display(),
            expected_vcf,
            actual_vcf
        );
    }

    println!("fixture validation passed: {}", args.input.display());
    Ok(())
}

pub(super) fn validate_sv_fixture(args: ValidateSvFixtureArgs) -> Result<()> {
    if !args.input.exists() {
        anyhow::bail!("fixture input does not exist: {}", args.input.display());
    }
    if !args.expected_json.exists() {
        anyhow::bail!(
            "expected SV JSON does not exist: {}",
            args.expected_json.display()
        );
    }

    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    let json = engine
        .analyze_with_options(
            &args.input,
            None,
            AnalyzeOptions {
                filter_numt: true,
                threads: 1,
                ..AnalyzeOptions::default()
            },
        )
        .context("analysis failed")?;
    let data: Value = serde_json::from_str(&json).context("analysis returned invalid JSON")?;
    let actual = serde_json::json!({
        "schema_version": data.pointer("/metadata/sv_event_schema_version"),
        "input_alignment_records": data.pointer("/filter_stats/input_alignment_records"),
        "input_molecules": data.pointer("/filter_stats/input_molecules"),
        "passed_reads": data.pointer("/filter_stats/passed_reads"),
        "numt_filtered_reads": data.pointer("/filter_stats/numt_filtered_reads"),
        "svs": data.get("svs"),
        "complex_events": data.get("complex_events"),
    });
    let expected: Value = serde_json::from_str(
        &fs::read_to_string(&args.expected_json)
            .with_context(|| format!("failed to read {}", args.expected_json.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.expected_json.display()))?;
    if actual != expected {
        anyhow::bail!(
            "SV golden mismatch for {}\nexpected:\n{}\nactual:\n{}",
            args.input.display(),
            serde_json::to_string_pretty(&expected)?,
            serde_json::to_string_pretty(&actual)?
        );
    }

    println!("SV fixture validation passed: {}", args.input.display());
    Ok(())
}

pub(super) fn validate_evidence_fixture(args: ValidateEvidenceFixtureArgs) -> Result<()> {
    if !args.input.exists() {
        anyhow::bail!("fixture input does not exist: {}", args.input.display());
    }
    if !args.expected_json.exists() {
        anyhow::bail!(
            "expected evidence JSON does not exist: {}",
            args.expected_json.display()
        );
    }

    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    let json = engine
        .analyze_with_options(
            &args.input,
            None,
            AnalyzeOptions {
                filter_numt: true,
                threads: 1,
                ..AnalyzeOptions::default()
            },
        )
        .context("analysis failed")?;
    let data: Value = serde_json::from_str(&json).context("analysis returned invalid JSON")?;
    validate_result_contract(&data)?;
    let actual = project_evidence(&data)?;
    if let Some(path) = &args.write_projection {
        write_json_atomic(path, &actual)
            .with_context(|| format!("failed to write projection {}", path.display()))?;
    }
    let expected: Value = serde_json::from_str(
        &fs::read_to_string(&args.expected_json)
            .with_context(|| format!("failed to read {}", args.expected_json.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.expected_json.display()))?;
    if actual != expected {
        anyhow::bail!(
            "scientific evidence golden mismatch for {}\nexpected:\n{}\nactual:\n{}",
            args.input.display(),
            serde_json::to_string_pretty(&expected)?,
            serde_json::to_string_pretty(&actual)?
        );
    }

    println!(
        "scientific evidence fixture validation passed: {}",
        args.input.display()
    );
    Ok(())
}

pub(super) fn validate_evidence_graph_fixture(
    args: ValidateEvidenceGraphFixtureArgs,
) -> Result<()> {
    if !args.input.exists() {
        anyhow::bail!("fixture input does not exist: {}", args.input.display());
    }
    if !args.expected_json.exists() {
        anyhow::bail!(
            "expected evidence-graph JSON does not exist: {}",
            args.expected_json.display()
        );
    }

    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    let json = engine
        .analyze_with_options(
            &args.input,
            None,
            AnalyzeOptions {
                filter_numt: true,
                threads: 4,
                emit_evidence_graph: true,
                evidence_page_size: args.evidence_page_size,
                molecule_id_tag: args.molecule_id_tag,
                umi_tag: args.umi_tag,
                duplex_tag: args.duplex_tag,
                ..AnalyzeOptions::default()
            },
        )
        .context("schema 0.6 analysis failed")?;
    let data: Value = serde_json::from_str(&json).context("analysis returned invalid JSON")?;
    validate_result_contract(&data)?;
    let actual = project_evidence_graph(&data)?;
    if let Some(path) = &args.write_projection {
        write_json_atomic(path, &actual)
            .with_context(|| format!("failed to write projection {}", path.display()))?;
    }
    let expected: Value = serde_json::from_str(
        &fs::read_to_string(&args.expected_json)
            .with_context(|| format!("failed to read {}", args.expected_json.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.expected_json.display()))?;
    if actual != expected {
        anyhow::bail!(
            "schema 0.6 evidence-graph golden mismatch for {}\nexpected:\n{}\nactual:\n{}",
            args.input.display(),
            serde_json::to_string_pretty(&expected)?,
            serde_json::to_string_pretty(&actual)?
        );
    }

    println!(
        "schema 0.6 evidence-graph fixture validation passed: {}",
        args.input.display()
    );
    Ok(())
}

pub(super) struct ScopedEnvironmentVariable {
    pub(super) name: &'static str,
    pub(super) previous: Option<std::ffi::OsString>,
}

impl ScopedEnvironmentVariable {
    pub(super) fn set(name: &'static str, value: &Path) -> Self {
        let previous = std::env::var_os(name);
        std::env::set_var(name, value);
        Self { name, previous }
    }
}

impl Drop for ScopedEnvironmentVariable {
    fn drop(&mut self) {
        if let Some(previous) = &self.previous {
            std::env::set_var(self.name, previous);
        } else {
            std::env::remove_var(self.name);
        }
    }
}

pub(super) fn validate_clinical_fixture(args: ValidateClinicalFixtureArgs) -> Result<()> {
    if !args.input.exists() {
        anyhow::bail!("fixture input does not exist: {}", args.input.display());
    }
    if !args.annotations.exists() {
        anyhow::bail!(
            "clinical fixture does not exist: {}",
            args.annotations.display()
        );
    }
    if !args.expected_json.exists() {
        anyhow::bail!(
            "expected clinical JSON does not exist: {}",
            args.expected_json.display()
        );
    }
    validate_clinical_tsv_schema(&args.annotations)?;
    let _environment =
        ScopedEnvironmentVariable::set("MITO_CLINICAL_ANNOTATIONS", &args.annotations);
    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    let json = engine
        .analyze_with_options(
            &args.input,
            None,
            AnalyzeOptions {
                filter_numt: true,
                threads: 1,
                ..AnalyzeOptions::default()
            },
        )
        .context("analysis failed")?;
    let data: Value = serde_json::from_str(&json).context("analysis returned invalid JSON")?;
    validate_result_contract(&data)?;
    let variants = require_array(&data, "/variants")?
        .iter()
        .filter_map(|variant| {
            variant.get("annotation").map(|annotation| {
                serde_json::json!({
                    "position": variant.get("position"),
                    "ref": variant.get("ref"),
                    "alt": variant.get("alt"),
                    "annotation": annotation,
                })
            })
        })
        .collect::<Vec<_>>();
    let actual = serde_json::json!({
        "schema_version": data.pointer("/metadata/clinical_annotation_schema_version"),
        "variants": variants,
    });
    let expected: Value = serde_json::from_str(
        &fs::read_to_string(&args.expected_json)
            .with_context(|| format!("failed to read {}", args.expected_json.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.expected_json.display()))?;
    if actual != expected {
        anyhow::bail!(
            "clinical assertion golden mismatch for {}\nexpected:\n{}\nactual:\n{}",
            args.input.display(),
            serde_json::to_string_pretty(&expected)?,
            serde_json::to_string_pretty(&actual)?
        );
    }
    println!(
        "clinical assertion fixture validation passed: {}",
        args.input.display()
    );
    Ok(())
}

pub(super) fn validate_clinical_manifest(args: ValidateClinicalManifestArgs) -> Result<()> {
    if !args.manifest.exists() {
        anyhow::bail!(
            "clinical manifest does not exist: {}",
            args.manifest.display()
        );
    }
    let manifest: Value = serde_json::from_str(
        &fs::read_to_string(&args.manifest)
            .with_context(|| format!("failed to read {}", args.manifest.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.manifest.display()))?;
    if manifest.get("schema_version").and_then(Value::as_str) != Some("1.0") {
        anyhow::bail!("clinical manifest requires schema_version 1.0");
    }
    let base = args.manifest.parent().unwrap_or_else(|| Path::new("."));
    let input = manifest
        .get("input")
        .and_then(Value::as_str)
        .context("clinical manifest requires string input")?;
    let input = base.join(input);
    if !input.exists() {
        anyhow::bail!(
            "clinical manifest input does not exist: {}",
            input.display()
        );
    }
    let cases = manifest
        .get("cases")
        .and_then(Value::as_array)
        .context("clinical manifest requires cases array")?;
    if cases.is_empty() {
        anyhow::bail!("clinical manifest contains no cases");
    }
    for case in cases {
        let name = case
            .get("name")
            .and_then(Value::as_str)
            .context("clinical case requires string name")?;
        let annotations = case
            .get("annotations")
            .and_then(Value::as_str)
            .context("clinical case requires string annotations")?;
        let expected_code = case
            .get("expected_error_code")
            .and_then(Value::as_str)
            .context("clinical case requires string expected_error_code")?;
        let annotations = base.join(annotations);
        if !annotations.exists() {
            anyhow::bail!(
                "clinical case '{name}' annotations do not exist: {}",
                annotations.display()
            );
        }
        let result = {
            let _environment =
                ScopedEnvironmentVariable::set("MITO_CLINICAL_ANNOTATIONS", &annotations);
            let engine = MitoEngine::new().context("failed to create analysis engine")?;
            engine.analyze_with_options(
                &input,
                None,
                AnalyzeOptions {
                    filter_numt: true,
                    threads: 1,
                    ..AnalyzeOptions::default()
                },
            )
        };
        match result {
            Ok(_) => anyhow::bail!("clinical negative case '{name}' unexpectedly succeeded"),
            Err(error) if error.code == expected_code => {}
            Err(error) => anyhow::bail!(
                "clinical negative case '{name}' returned {}, expected {}: {}",
                error.code,
                expected_code,
                error.message
            ),
        }
    }
    println!(
        "clinical negative manifest validation passed: {} cases",
        cases.len()
    );
    Ok(())
}

pub(super) fn validate_error_manifest(args: ValidateErrorManifestArgs) -> Result<()> {
    if !args.manifest.exists() {
        anyhow::bail!("error manifest does not exist: {}", args.manifest.display());
    }
    let manifest: Value = serde_json::from_str(
        &fs::read_to_string(&args.manifest)
            .with_context(|| format!("failed to read {}", args.manifest.display()))?,
    )
    .with_context(|| format!("invalid JSON in {}", args.manifest.display()))?;
    let expected_schema = manifest
        .get("error_schema_version")
        .and_then(Value::as_str)
        .context("error manifest requires string error_schema_version")?;
    let capabilities = MitoEngine::capabilities();
    if expected_schema != capabilities.error_schema_version {
        anyhow::bail!(
            "error-schema mismatch: manifest {}, engine {}",
            expected_schema,
            capabilities.error_schema_version
        );
    }
    let cases = manifest
        .get("cases")
        .and_then(Value::as_array)
        .context("error manifest requires a cases array")?;
    if cases.is_empty() {
        anyhow::bail!("error manifest cases array must not be empty");
    }

    let base = args.manifest.parent().unwrap_or_else(|| Path::new("."));
    let engine = MitoEngine::new().context("failed to create analysis engine")?;
    let mut case_names = BTreeSet::new();
    for (index, case) in cases.iter().enumerate() {
        let object = case
            .as_object()
            .with_context(|| format!("error manifest case {index} must be an object"))?;
        let name = object
            .get("name")
            .and_then(Value::as_str)
            .with_context(|| format!("error manifest case {index} requires string name"))?;
        if !case_names.insert(name.to_owned()) {
            anyhow::bail!("duplicate error manifest case name: {name}");
        }
        let input = object
            .get("input")
            .and_then(Value::as_str)
            .with_context(|| format!("error manifest case '{name}' requires string input"))?;
        let expected_code = object
            .get("code")
            .and_then(Value::as_str)
            .with_context(|| format!("error manifest case '{name}' requires string code"))?;
        let input_path = {
            let path = Path::new(input);
            if path.is_absolute() {
                path.to_path_buf()
            } else {
                base.join(path)
            }
        };
        if !input_path.is_file() {
            anyhow::bail!(
                "error manifest case '{name}' input does not exist: {}",
                input_path.display()
            );
        }
        let error = match engine.analyze_with_options(&input_path, None, AnalyzeOptions::default())
        {
            Ok(_) => {
                anyhow::bail!("negative golden '{name}' unexpectedly produced a scientific result")
            }
            Err(error) => error,
        };
        if error.code != expected_code {
            anyhow::bail!(
                "negative golden '{name}' expected {expected_code}, got {}: {}",
                error.code,
                error.message
            );
        }
        println!("negative fixture passed: {name} [{expected_code}]");
    }

    println!(
        "error manifest validation passed: {} cases from {}",
        cases.len(),
        args.manifest.display()
    );
    Ok(())
}

pub(super) fn assert_json_u64(
    data: &Value,
    pointer: &str,
    expected: u64,
    label: &str,
) -> Result<()> {
    let actual = data
        .pointer(pointer)
        .and_then(Value::as_u64)
        .with_context(|| format!("missing numeric field {pointer}"))?;
    if actual != expected {
        anyhow::bail!("{label} mismatch: expected {expected}, got {actual}");
    }
    Ok(())
}

pub(super) fn assert_snp_present(data: &Value, expected: &str) -> Result<()> {
    let mut parts = expected.split(':');
    let position = parts
        .next()
        .context("expected SNP must be formatted as position:ref:alt")?
        .parse::<u64>()
        .context("expected SNP position is not numeric")?;
    let reference = parts
        .next()
        .context("expected SNP must include reference allele")?;
    let alternate = parts
        .next()
        .context("expected SNP must include alternate allele")?;
    if parts.next().is_some() {
        anyhow::bail!("expected SNP must be formatted as position:ref:alt");
    }

    for read in data
        .get("reads")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
    {
        for snp in read
            .get("snps")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
        {
            if snp.get("position").and_then(Value::as_u64) == Some(position)
                && snp.get("ref").and_then(Value::as_str) == Some(reference)
                && snp.get("alt").and_then(Value::as_str) == Some(alternate)
            {
                return Ok(());
            }
        }
    }
    anyhow::bail!("expected SNP not found: {expected}");
}

pub(super) fn assert_sv_present(data: &Value, expected_id: &str) -> Result<()> {
    let found = data
        .get("svs")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .any(|sv| sv.get("id").and_then(Value::as_str) == Some(expected_id));
    if found {
        Ok(())
    } else {
        anyhow::bail!("expected SV not found: {expected_id}");
    }
}

pub(super) fn assert_mapq_present(data: &Value, expected_mapq: u8) -> Result<()> {
    let found = data
        .pointer("/coverage_metrics/mapping_quality_histogram")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .any(|entry| {
            entry.get("mapq").and_then(Value::as_u64) == Some(u64::from(expected_mapq))
                && entry.get("count").and_then(Value::as_u64).unwrap_or(0) > 0
        });
    if found {
        Ok(())
    } else {
        anyhow::bail!("expected MAPQ not found in histogram: {expected_mapq}");
    }
}

pub(super) fn assert_aux_present(data: &Value, expected: &str) -> Result<()> {
    let (key, value) = expected
        .split_once('=')
        .context("expected aux tag must be formatted as TAG=value")?;
    if key.len() != 2 || !key.chars().all(|c| c.is_ascii_alphanumeric()) {
        anyhow::bail!("expected aux tag key must be a two-character SAM tag");
    }

    for read in data
        .get("reads")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
    {
        if read
            .get("aux_tags")
            .and_then(Value::as_object)
            .and_then(|tags| tags.get(key))
            .and_then(Value::as_str)
            == Some(value)
        {
            return Ok(());
        }
    }
    anyhow::bail!("expected aux tag not found: {expected}");
}
