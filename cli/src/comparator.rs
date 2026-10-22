use crate::args::{ComparatorFilterPolicy, CompareMtDnaServer2Args};
use crate::contract::{require_array, validate_result_contract};
use crate::storage::{sha256_bytes, sha256_file, write_json_atomic};
use anyhow::Context;
use anyhow::Result;
use serde_json::Value;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::fs;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

#[derive(Clone, Debug)]
pub(super) struct ComparatorVariant {
    pub(super) sample_id: String,
    pub(super) filter: String,
    pub(super) position: u64,
    pub(super) reference: String,
    pub(super) alternate: String,
    pub(super) heteroplasmy: f64,
    pub(super) coverage: u64,
    pub(super) coverage_forward: u64,
    pub(super) coverage_reverse: u64,
    pub(super) mutation: String,
}

#[derive(Clone, Debug)]
pub(super) struct MitoVariantProjection {
    pub(super) event_id: String,
    pub(super) position: u64,
    pub(super) reference: String,
    pub(super) alternate: String,
    pub(super) heteroplasmy: f64,
    pub(super) callable_depth: u64,
    pub(super) alt_depth: u64,
    pub(super) filter_status: String,
}

pub(super) fn compare_mtdna_server2(args: CompareMtDnaServer2Args) -> Result<()> {
    if args.comparator_version != "2.1.16" {
        anyhow::bail!(
            "unsupported mtDNA-Server 2 format version {}; only 2.1.16 is pinned",
            args.comparator_version
        );
    }
    validate_unit_interval(args.min_call_concordance, "--min-call-concordance")?;
    validate_unit_interval(args.max_mean_hf_delta, "--max-mean-hf-delta")?;
    if args
        .sample
        .as_deref()
        .is_some_and(|sample| sample.is_empty() || sample.trim() != sample)
    {
        anyhow::bail!("--sample must be a non-empty exact ID without surrounding whitespace");
    }

    let result_bytes = fs::read(&args.result)
        .with_context(|| format!("failed to read {}", args.result.display()))?;
    let result: Value = serde_json::from_slice(&result_bytes)
        .with_context(|| format!("{} is not valid JSON", args.result.display()))?;
    validate_result_contract(&result)
        .with_context(|| format!("{} is not a valid analysis result", args.result.display()))?;
    if result
        .pointer("/metadata/schema_version")
        .and_then(Value::as_str)
        != Some("0.6")
    {
        anyhow::bail!("mtDNA-Server 2 comparison requires a schema 0.6 result");
    }
    if result
        .pointer("/metadata/reference_length")
        .and_then(Value::as_u64)
        != Some(16_569)
    {
        anyhow::bail!(
            "mtDNA-Server 2 comparison requires the 16,569 bp human rCRS coordinate system"
        );
    }
    if let Some(accession) = result
        .pointer("/metadata/reference_accession")
        .and_then(Value::as_str)
    {
        if accession != "NC_012920.1" {
            anyhow::bail!(
                "mtDNA-Server 2 comparison requires rCRS accession NC_012920.1, found {accession}"
            );
        }
    }

    let comparator_rows = parse_mtdna_server2_variants(&args.comparator)?;
    let available_samples = comparator_rows
        .iter()
        .map(|row| row.sample_id.as_str())
        .collect::<BTreeSet<_>>();
    let selected_sample =
        match args.sample.as_deref() {
            Some(sample) if available_samples.is_empty() => sample.to_owned(),
            Some(sample) if available_samples.contains(sample) => sample.to_owned(),
            Some(sample) => anyhow::bail!(
                "sample {sample:?} is absent from {}; available samples: {}",
                args.comparator.display(),
                available_samples
                    .iter()
                    .copied()
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
            None if available_samples.len() == 1 => available_samples
                .iter()
                .next()
                .context("comparator sample set is unexpectedly empty")?
                .to_string(),
            None if available_samples.is_empty() => anyhow::bail!(
                "--sample is required when the comparator file contains only a header"
            ),
            None => {
                anyhow::bail!(
            "--sample is required for a multi-sample comparator file; available samples: {}",
            available_samples.iter().copied().collect::<Vec<_>>().join(", ")
        )
            }
        };

    let include_comparator = |row: &&ComparatorVariant| {
        row.sample_id == selected_sample
            && match args.filter_policy {
                ComparatorFilterPolicy::PassOnly => row.filter == "PASS",
                ComparatorFilterPolicy::All => true,
            }
    };
    let mut comparator_by_key = BTreeMap::<String, &ComparatorVariant>::new();
    for row in comparator_rows.iter().filter(include_comparator) {
        let key = variant_comparison_key(row.position, &row.reference, &row.alternate);
        if comparator_by_key.insert(key.clone(), row).is_some() {
            anyhow::bail!("duplicate mtDNA-Server 2 variant {key} for sample {selected_sample}");
        }
    }
    let mito_variants = parse_mito_variant_projection(&result)?;
    let mut mito_by_key = BTreeMap::<String, &MitoVariantProjection>::new();
    for variant in &mito_variants {
        let key = variant_comparison_key(variant.position, &variant.reference, &variant.alternate);
        if mito_by_key.insert(key.clone(), variant).is_some() {
            anyhow::bail!("duplicate Mito-Architect variant comparison key {key}");
        }
    }

    let mut matched = Vec::new();
    let mut mito_only = Vec::new();
    let mut comparator_only = Vec::new();
    let mut hf_delta_sum = 0.0;
    let mut max_hf_delta = 0.0_f64;
    for (key, mito) in &mito_by_key {
        if let Some(comparator) = comparator_by_key.get(key) {
            let hf_delta = (mito.heteroplasmy - comparator.heteroplasmy).abs();
            hf_delta_sum += hf_delta;
            max_hf_delta = max_hf_delta.max(hf_delta);
            matched.push(serde_json::json!({
                "key": key,
                "event_id": mito.event_id,
                "position": mito.position,
                "ref": mito.reference,
                "alt": mito.alternate,
                "mito_architect": {
                    "heteroplasmy": mito.heteroplasmy,
                    "callable_depth": mito.callable_depth,
                    "alt_depth": mito.alt_depth,
                    "filter_status": mito.filter_status,
                },
                "mtdna_server_2": {
                    "heteroplasmy": comparator.heteroplasmy,
                    "coverage": comparator.coverage,
                    "coverage_forward": comparator.coverage_forward,
                    "coverage_reverse": comparator.coverage_reverse,
                    "filter": comparator.filter,
                    "mutation": comparator.mutation,
                },
                "absolute_hf_delta": hf_delta,
            }));
        } else {
            mito_only.push(serde_json::json!({
                "key": key,
                "event_id": mito.event_id,
                "position": mito.position,
                "ref": mito.reference,
                "alt": mito.alternate,
                "heteroplasmy": mito.heteroplasmy,
                "callable_depth": mito.callable_depth,
                "alt_depth": mito.alt_depth,
                "filter_status": mito.filter_status,
            }));
        }
    }
    for (key, comparator) in &comparator_by_key {
        if !mito_by_key.contains_key(key) {
            comparator_only.push(serde_json::json!({
                "key": key,
                "position": comparator.position,
                "ref": comparator.reference,
                "alt": comparator.alternate,
                "heteroplasmy": comparator.heteroplasmy,
                "coverage": comparator.coverage,
                "coverage_forward": comparator.coverage_forward,
                "coverage_reverse": comparator.coverage_reverse,
                "filter": comparator.filter,
                "mutation": comparator.mutation,
            }));
        }
    }

    let union_count = matched.len() + mito_only.len() + comparator_only.len();
    let call_concordance = if union_count == 0 {
        1.0
    } else {
        matched.len() as f64 / union_count as f64
    };
    let mean_hf_delta = if matched.is_empty() {
        None
    } else {
        Some(hf_delta_sum / matched.len() as f64)
    };
    let filter_policy = match args.filter_policy {
        ComparatorFilterPolicy::PassOnly => "pass_only",
        ComparatorFilterPolicy::All => "all",
    };
    let report = serde_json::json!({
        "schema_version": "1.0",
        "report_type": "variant_callset_differential",
        "interpretation": "Comparator concordance is not analytical truth, sensitivity, specificity, or clinical validation.",
        "coordinate_system": {
            "name": "rCRS",
            "accession": "NC_012920.1",
            "length": 16569,
        },
        "mtdna_server_2": {
            "version": args.comparator_version,
            "format": "variants.annotated.txt",
            "format_contract": "mtdna-server-2-v2.1.16-variants-annotated",
            "upstream_release": "https://github.com/genepi/mtdna-server-2/releases/tag/v2.1.16",
            "sample_id": selected_sample,
            "filter_policy": filter_policy,
        },
        "provenance": {
            "mito_architect_result": {
                "path": args.result.to_string_lossy(),
                "sha256": sha256_bytes(&result_bytes)?,
            },
            "comparator_result": {
                "path": args.comparator.to_string_lossy(),
                "sha256": sha256_file(&args.comparator)?,
            },
        },
        "metrics": {
            "matched": matched.len(),
            "mito_architect_only": mito_only.len(),
            "mtdna_server_2_only": comparator_only.len(),
            "union": union_count,
            "call_concordance": call_concordance,
            "mean_absolute_hf_delta": mean_hf_delta,
            "max_absolute_hf_delta": if matched.is_empty() { None } else { Some(max_hf_delta) },
        },
        "matched": matched,
        "mito_architect_only": mito_only,
        "mtdna_server_2_only": comparator_only,
    });
    write_json_atomic(&args.output, &report)?;

    if let Some(minimum) = args.min_call_concordance {
        if call_concordance + f64::EPSILON < minimum {
            anyhow::bail!(
                "call concordance {call_concordance:.9} is below the development gate {minimum:.9}; report written to {}",
                args.output.display()
            );
        }
    }
    if let Some(maximum) = args.max_mean_hf_delta {
        match mean_hf_delta {
            Some(value) if value <= maximum + f64::EPSILON => {}
            Some(value) => anyhow::bail!(
                "mean absolute HF delta {value:.9} exceeds the development gate {maximum:.9}; report written to {}",
                args.output.display()
            ),
            None => anyhow::bail!(
                "mean absolute HF delta is undefined because no variants matched; report written to {}",
                args.output.display()
            ),
        }
    }
    println!(
        "mtDNA-Server 2 differential: matched={} mito_only={} comparator_only={} concordance={:.6}",
        report
            .pointer("/metrics/matched")
            .and_then(Value::as_u64)
            .unwrap_or(0),
        report
            .pointer("/metrics/mito_architect_only")
            .and_then(Value::as_u64)
            .unwrap_or(0),
        report
            .pointer("/metrics/mtdna_server_2_only")
            .and_then(Value::as_u64)
            .unwrap_or(0),
        call_concordance
    );
    Ok(())
}

pub(super) fn validate_unit_interval(value: Option<f64>, name: &str) -> Result<()> {
    if value.is_some_and(|value| !value.is_finite() || !(0.0..=1.0).contains(&value)) {
        anyhow::bail!("{name} must be a finite number between 0 and 1");
    }
    Ok(())
}

pub(super) fn parse_mtdna_server2_variants(path: &Path) -> Result<Vec<ComparatorVariant>> {
    const REQUIRED: [&str; 15] = [
        "ID",
        "Filter",
        "Pos",
        "Ref",
        "Variant",
        "VariantLevel",
        "MajorBase",
        "MajorLevel",
        "MinorBase",
        "MinorLevel",
        "Coverage",
        "CoverageFWD",
        "CoverageREV",
        "Type",
        "Mutation",
    ];
    let file =
        fs::File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    let mut lines = BufReader::new(file).lines();
    let header_line = lines
        .next()
        .transpose()
        .with_context(|| format!("failed to read {}", path.display()))?
        .context("mtDNA-Server 2 comparator file is empty")?;
    let header = header_line
        .trim_end_matches('\r')
        .split('\t')
        .collect::<Vec<_>>();
    let mut indices = BTreeMap::<&str, usize>::new();
    for (index, name) in header.iter().copied().enumerate() {
        if name.is_empty() || indices.insert(name, index).is_some() {
            anyhow::bail!("mtDNA-Server 2 header contains an empty or duplicate column");
        }
    }
    for required in REQUIRED {
        if !indices.contains_key(required) {
            anyhow::bail!(
                "mtDNA-Server 2 v2.1.16 variants.annotated.txt is missing required column {required}"
            );
        }
    }

    let field = |fields: &[&str], name: &str| -> Result<String> {
        let index = *indices
            .get(name)
            .with_context(|| format!("required comparator column {name} is unresolved"))?;
        fields
            .get(index)
            .map(|value| value.trim().to_owned())
            .with_context(|| format!("comparator row has no {name} column"))
    };
    let mut rows = Vec::new();
    for (offset, line) in lines.enumerate() {
        let line_number = offset + 2;
        let line =
            line.with_context(|| format!("failed to read {}:{line_number}", path.display()))?;
        let line = line.trim_end_matches('\r');
        if line.is_empty() {
            anyhow::bail!("empty comparator row at {}:{line_number}", path.display());
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != header.len() {
            anyhow::bail!(
                "{}:{line_number} has {} columns; expected {}",
                path.display(),
                fields.len(),
                header.len()
            );
        }
        let sample_id = field(&fields, "ID")?;
        let filter = field(&fields, "Filter")?;
        let position = parse_comparator_u64(&field(&fields, "Pos")?, "Pos", line_number)?;
        if !(1..=16_569).contains(&position) {
            anyhow::bail!("{}:{line_number} Pos is outside rCRS", path.display());
        }
        let reference = field(&fields, "Ref")?.to_ascii_uppercase();
        let alternate = field(&fields, "Variant")?.to_ascii_uppercase();
        if sample_id.is_empty()
            || filter.is_empty()
            || !valid_comparator_allele(&reference)
            || !valid_comparator_allele(&alternate)
        {
            anyhow::bail!(
                "{}:{line_number} has an invalid identity, filter, REF, or ALT",
                path.display()
            );
        }
        let heteroplasmy = parse_comparator_fraction(
            &field(&fields, "VariantLevel")?,
            "VariantLevel",
            line_number,
        )?;
        let major_base = field(&fields, "MajorBase")?.to_ascii_uppercase();
        let minor_base = field(&fields, "MinorBase")?.to_ascii_uppercase();
        if !valid_comparator_allele(&major_base) || !valid_comparator_allele(&minor_base) {
            anyhow::bail!(
                "{}:{line_number} has an invalid MajorBase or MinorBase",
                path.display()
            );
        }
        let _major_level =
            parse_comparator_fraction(&field(&fields, "MajorLevel")?, "MajorLevel", line_number)?;
        let _minor_level =
            parse_comparator_fraction(&field(&fields, "MinorLevel")?, "MinorLevel", line_number)?;
        let coverage = parse_comparator_u64(&field(&fields, "Coverage")?, "Coverage", line_number)?;
        let coverage_forward =
            parse_comparator_u64(&field(&fields, "CoverageFWD")?, "CoverageFWD", line_number)?;
        let coverage_reverse =
            parse_comparator_u64(&field(&fields, "CoverageREV")?, "CoverageREV", line_number)?;
        if coverage_forward.checked_add(coverage_reverse) != Some(coverage) {
            anyhow::bail!(
                "{}:{line_number} CoverageFWD + CoverageREV does not equal Coverage",
                path.display()
            );
        }
        if coverage == 0 {
            anyhow::bail!("{}:{line_number} variant Coverage is zero", path.display());
        }
        let variant_type = field(&fields, "Type")?;
        let mutation = field(&fields, "Mutation")?;
        if variant_type.is_empty() || mutation.is_empty() {
            anyhow::bail!(
                "{}:{line_number} has an empty Type or Mutation",
                path.display()
            );
        }
        rows.push(ComparatorVariant {
            sample_id,
            filter,
            position,
            reference,
            alternate,
            heteroplasmy,
            coverage,
            coverage_forward,
            coverage_reverse,
            mutation,
        });
    }
    Ok(rows)
}

pub(super) fn parse_comparator_u64(value: &str, field: &str, line: usize) -> Result<u64> {
    value
        .parse::<u64>()
        .with_context(|| format!("comparator line {line} {field} is not an unsigned integer"))
}

pub(super) fn parse_comparator_fraction(value: &str, field: &str, line: usize) -> Result<f64> {
    let value = value
        .parse::<f64>()
        .with_context(|| format!("comparator line {line} {field} is not numeric"))?;
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        anyhow::bail!("comparator line {line} {field} is outside [0,1]");
    }
    Ok(value)
}

pub(super) fn valid_comparator_allele(value: &str) -> bool {
    !value.is_empty()
        && value
            .bytes()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N' | b'-'))
}

pub(super) fn parse_mito_variant_projection(result: &Value) -> Result<Vec<MitoVariantProjection>> {
    require_array(result, "/variants")?
        .iter()
        .enumerate()
        .map(|(index, variant)| {
            let string = |name: &str| {
                variant
                    .get(name)
                    .and_then(Value::as_str)
                    .map(str::to_owned)
                    .with_context(|| format!("/variants/{index}/{name} must be a string"))
            };
            let unsigned = |name: &str| {
                variant
                    .get(name)
                    .and_then(Value::as_u64)
                    .with_context(|| format!("/variants/{index}/{name} must be an integer"))
            };
            let heteroplasmy = variant
                .get("heteroplasmy")
                .and_then(Value::as_f64)
                .with_context(|| format!("/variants/{index}/heteroplasmy must be a number"))?;
            if !heteroplasmy.is_finite() || !(0.0..=1.0).contains(&heteroplasmy) {
                anyhow::bail!("/variants/{index}/heteroplasmy is outside [0,1]");
            }
            Ok(MitoVariantProjection {
                event_id: string("event_id")?,
                position: unsigned("position")?,
                reference: string("ref")?.to_ascii_uppercase(),
                alternate: string("alt")?.to_ascii_uppercase(),
                heteroplasmy,
                callable_depth: unsigned("callable_depth")?,
                alt_depth: unsigned("alt_depth")?,
                filter_status: string("filter_status")?,
            })
        })
        .collect()
}

pub(super) fn variant_comparison_key(position: u64, reference: &str, alternate: &str) -> String {
    format!("{position}:{reference}:{alternate}")
}
