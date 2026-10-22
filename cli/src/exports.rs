use crate::args::AnalyzeArgs;
use crate::contract::{require_array, require_string, validate_result_contract};
use crate::storage::{sha256_bytes, sha256_file, write_json_atomic};
use anyhow::Context;
use anyhow::Result;
use serde_json::Value;
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command as ProcessCommand;

#[derive(Debug)]
pub(super) struct ExportArtifact {
    pub(super) kind: &'static str,
    pub(super) path: PathBuf,
    pub(super) bytes: u64,
    pub(super) sha256: String,
}

impl ExportArtifact {
    pub(super) fn from_path(kind: &'static str, path: &Path) -> Result<Self> {
        Ok(Self {
            kind,
            path: path.to_path_buf(),
            bytes: fs::metadata(path)
                .with_context(|| format!("failed to stat {}", path.display()))?
                .len(),
            sha256: sha256_file(path)?,
        })
    }
}

pub(super) fn current_git_commit() -> Option<String> {
    let output = ProcessCommand::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let value = String::from_utf8(output.stdout).ok()?.trim().to_owned();
    (!value.is_empty()).then_some(value)
}

pub(super) fn build_provenance_manifest(
    result: &str,
    args: &AnalyzeArgs,
    exported: &[ExportArtifact],
) -> Result<Value> {
    let data: Value = serde_json::from_str(result).context("analysis returned invalid JSON")?;
    validate_result_contract(&data)?;
    let reference_path = data
        .pointer("/metadata/reference_path")
        .and_then(Value::as_str)
        .context("result contract violation: /metadata/reference_path must be a string")?;
    let reference_path = Path::new(reference_path);
    let output_records = exported
        .iter()
        .map(|artifact| {
            serde_json::json!({
                "kind": artifact.kind,
                "path": artifact.path.display().to_string(),
                "bytes": artifact.bytes,
                "sha256": artifact.sha256,
            })
        })
        .collect::<Vec<_>>();
    Ok(serde_json::json!({
        "schema_version": "1.0",
        "timestamp_policy": "omitted_for_determinism",
        "software": {
            "name": "Mito-Architect",
            "engine_version": data.pointer("/metadata/engine_version"),
            "git_commit": current_git_commit(),
            "result_schema_version": data.pointer("/metadata/schema_version"),
            "sv_event_schema_version": data.pointer("/metadata/sv_event_schema_version"),
            "complex_sv_event_schema_version": data.pointer("/metadata/complex_sv_event_schema_version"),
            "clinical_annotation_schema_version": data.pointer("/metadata/clinical_annotation_schema_version"),
        },
        "determinism": {
            "deterministic_algorithms": true,
            "random_seed": Value::Null,
            "thread_count": data.pointer("/metadata/threads"),
        },
        "input": {
            "path": args.input.display().to_string(),
            "bytes": fs::metadata(&args.input)?.len(),
            "sha256": sha256_file(&args.input)?,
        },
        "reference": {
            "path": reference_path.display().to_string(),
            "accession": data.pointer("/metadata/reference_accession"),
            "length": data.pointer("/metadata/reference_length"),
            "sha256": sha256_file(reference_path)?,
        },
        "calling_parameters": data.pointer("/metadata/calling_parameters"),
        "resources": data.pointer("/metadata/resources"),
        "command_line": std::env::args().collect::<Vec<_>>(),
        "authoritative_result": {
            "format": "application/json",
            "bytes": result.len(),
            "sha256": sha256_bytes(result.as_bytes())?,
        },
        "exports": output_records,
    }))
}

pub(super) fn export_evidence_pages(result: &str, directory: &Path) -> Result<()> {
    let data: Value = serde_json::from_str(result).context("analysis returned invalid JSON")?;
    if data
        .pointer("/metadata/schema_version")
        .and_then(Value::as_str)
        != Some("0.6")
    {
        anyhow::bail!("--evidence-pages-dir requires --evidence-graph");
    }
    validate_result_contract(&data)?;
    fs::create_dir_all(directory)
        .with_context(|| format!("failed to create {}", directory.display()))?;
    let pages = require_array(&data, "/observation_pages")?;
    let mut manifest_pages = Vec::with_capacity(pages.len());
    for (page_index, page) in pages.iter().enumerate() {
        let file_name = format!("observations-{page_index:06}.json");
        let path = directory.join(&file_name);
        write_json_atomic(&path, page)?;
        let bytes = fs::metadata(&path)
            .with_context(|| format!("failed to stat {}", path.display()))?
            .len();
        manifest_pages.push(serde_json::json!({
            "index": page_index,
            "offset": page.get("offset"),
            "count": page.get("count"),
            "path": file_name,
            "bytes": bytes,
            "sha256": sha256_file(&path)?,
        }));
    }
    let manifest = serde_json::json!({
        "schema_version": "1.0",
        "result_schema_version": "0.6",
        "sample": data.pointer("/metadata/sample"),
        "reference_accession": data.pointer("/metadata/reference_accession"),
        "reference_length": data.pointer("/metadata/reference_length"),
        "evidence_encoding": data.get("evidence_encoding"),
        "pages": manifest_pages,
    });
    write_json_atomic(&directory.join("manifest.json"), &manifest)
}

pub(super) fn clean_export_field(value: &str) -> String {
    value
        .chars()
        .map(|character| match character {
            '\t' | '\n' | '\r' => ' ',
            other => other,
        })
        .collect::<String>()
        .trim()
        .to_owned()
}

pub(super) fn render_variant_tsv(json: &str) -> Result<String> {
    let data: Value = serde_json::from_str(json).context("analysis returned invalid JSON")?;
    validate_result_contract(&data)?;
    let schema = require_string(&data, "/metadata/schema_version")?;
    let variants = require_array(&data, "/variants")?;
    let mut out = String::from(
        "event_id\ttype\tposition\tstart\tend\tref\talt\talt_molecules\tref_molecules\tother_molecules\tcallable_molecules\theteroplasmy\tci95_low\tci95_high\talt_forward\talt_reverse\tref_forward\tref_reverse\talt_mean_mapq\talt_mean_baseq\tmulti_allelic\thomopolymer_run\tnumt_assessability\tfilter_status\tqc_flags\tsupporting_molecule_ids\tgene\tconsequence\tclinical_significance\n",
    );
    for (index, variant) in variants.iter().enumerate() {
        let string = |field: &str| {
            variant
                .get(field)
                .and_then(Value::as_str)
                .map(clean_export_field)
                .unwrap_or_default()
        };
        let integer = |field: &str| {
            variant
                .get(field)
                .and_then(Value::as_u64)
                .map(|value| value.to_string())
                .unwrap_or_default()
        };
        let number = |field: &str| {
            variant
                .get(field)
                .and_then(Value::as_f64)
                .map(|value| format!("{value:.8}"))
                .unwrap_or_default()
        };
        let event_id = variant
            .get("event_id")
            .and_then(Value::as_str)
            .map(clean_export_field)
            .unwrap_or_else(|| {
                format!(
                    "legacy:{}:{}:{}",
                    integer("position"),
                    string("ref"),
                    string("alt")
                )
            });
        let supporting = variant
            .get("supporting_molecule_ids")
            .or_else(|| variant.get("supporting_reads"))
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(Value::as_str)
            .map(clean_export_field)
            .collect::<Vec<_>>()
            .join(";");
        let qc_flags = variant
            .get("qc_flags")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(Value::as_str)
            .map(clean_export_field)
            .collect::<Vec<_>>()
            .join(";");
        let clinical = variant
            .pointer("/annotation/consensus_significance")
            .or_else(|| variant.pointer("/annotation/pathogenicity"))
            .and_then(Value::as_str)
            .map(clean_export_field)
            .unwrap_or_default();
        let fields = [
            event_id,
            if schema == "0.6" {
                string("type")
            } else {
                "SNV".to_owned()
            },
            integer("position"),
            integer("start"),
            integer("end"),
            string("ref"),
            string("alt"),
            integer("alt_depth"),
            integer("ref_depth"),
            integer("other_depth"),
            integer("callable_depth"),
            number("heteroplasmy"),
            number("ci95_low"),
            number("ci95_high"),
            variant
                .pointer("/strand_support/alt_forward")
                .and_then(Value::as_u64)
                .map(|value| value.to_string())
                .unwrap_or_default(),
            variant
                .pointer("/strand_support/alt_reverse")
                .and_then(Value::as_u64)
                .map(|value| value.to_string())
                .unwrap_or_default(),
            variant
                .pointer("/strand_support/ref_forward")
                .and_then(Value::as_u64)
                .map(|value| value.to_string())
                .unwrap_or_default(),
            variant
                .pointer("/strand_support/ref_reverse")
                .and_then(Value::as_u64)
                .map(|value| value.to_string())
                .unwrap_or_default(),
            variant
                .pointer("/mapping_quality/alternate/mean")
                .and_then(Value::as_f64)
                .map(|value| format!("{value:.4}"))
                .unwrap_or_default(),
            variant
                .pointer("/allele_quality/alternate/mean_phred")
                .and_then(Value::as_f64)
                .map(|value| format!("{value:.4}"))
                .unwrap_or_default(),
            variant
                .get("multi_allelic")
                .and_then(Value::as_bool)
                .map(|value| value.to_string())
                .unwrap_or_default(),
            variant
                .pointer("/homopolymer_context/run_length")
                .and_then(Value::as_u64)
                .map(|value| value.to_string())
                .unwrap_or_default(),
            string("numt_assessability"),
            string("filter_status"),
            qc_flags,
            supporting,
            string("gene"),
            string("consequence"),
            clinical,
        ];
        if fields
            .iter()
            .any(|field| field.contains('\t') || field.contains('\n') || field.contains('\r'))
        {
            anyhow::bail!("variant {index} contains an unsafe TSV field after sanitization");
        }
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }
    Ok(out)
}

#[derive(Debug)]
pub(super) struct VcfVariant {
    pub(super) chrom: String,
    pub(super) position: u64,
    pub(super) reference: String,
    pub(super) alternate: String,
    pub(super) support: u64,
    pub(super) depth: u64,
    pub(super) gene: Option<String>,
    pub(super) clinical_significance: Option<String>,
    pub(super) clinical_conflict: bool,
    pub(super) clinical_sources: Vec<String>,
    pub(super) ci95_low: Option<f64>,
    pub(super) ci95_high: Option<f64>,
    pub(super) event_id: Option<String>,
    pub(super) event_type: Option<String>,
    pub(super) reference_support: u64,
    pub(super) other_support: u64,
    pub(super) low_quality: u64,
    pub(super) conflict: u64,
    pub(super) qc_flags: Vec<String>,
    pub(super) numt_assessability: Option<String>,
    pub(super) normalization: Option<String>,
    pub(super) alt_mapping_quality: Option<f64>,
    pub(super) alt_base_quality: Option<f64>,
    pub(super) strand_support: Option<[u64; 4]>,
    pub(super) schema_0_6: bool,
}

impl VcfVariant {
    pub(super) fn heteroplasmy(&self) -> f64 {
        if self.depth == 0 {
            0.0
        } else {
            self.support as f64 / self.depth as f64
        }
    }
}

pub(super) fn render_vcf(json: &str) -> Result<String> {
    let data: Value = serde_json::from_str(json).context("analysis returned invalid JSON")?;
    validate_result_contract(&data)?;
    let sample = data
        .pointer("/metadata/sample")
        .and_then(Value::as_str)
        .unwrap_or("sample");
    let chrom = data
        .pointer("/metadata/reference_accession")
        .and_then(Value::as_str)
        .filter(|value| !value.is_empty() && *value != "custom")
        .unwrap_or("MT");
    let reference_length = data
        .pointer("/metadata/reference_length")
        .and_then(Value::as_u64)
        .context("result contract violation: /metadata/reference_length must be an integer")?;
    let schema_0_6 = data
        .pointer("/metadata/schema_version")
        .and_then(Value::as_str)
        == Some("0.6");
    let mut variants = BTreeMap::<(u64, String, String), VcfVariant>::new();
    let aggregates = data
        .get("variants")
        .and_then(Value::as_array)
        .context("result contract violation: /variants must be an array")?;
    for variant in aggregates {
        let Some(position) = variant
            .get(if schema_0_6 {
                "vcf_position"
            } else {
                "position"
            })
            .and_then(Value::as_u64)
        else {
            continue;
        };
        let Some(reference) = variant.get("ref").and_then(Value::as_str) else {
            continue;
        };
        let Some(alternate) = variant.get("alt").and_then(Value::as_str) else {
            continue;
        };
        let support = variant
            .get("alt_depth")
            .and_then(Value::as_u64)
            .unwrap_or(0);
        let depth = variant
            .get("callable_depth")
            .and_then(Value::as_u64)
            .unwrap_or(0);
        if schema_0_6 && variant.get("vcf_representable").and_then(Value::as_bool) == Some(false) {
            let event_id = variant
                .get("event_id")
                .and_then(Value::as_str)
                .unwrap_or("unknown");
            anyhow::bail!(
                "variant {event_id} crosses the circular origin and has no lossless linear VCF representation; use JSON or TSV"
            );
        }
        if reference.is_empty()
            || alternate.is_empty()
            || (!schema_0_6 && (reference.len() != 1 || alternate.len() != 1))
            || depth == 0
        {
            continue;
        }
        variants.insert(
            (position, reference.to_string(), alternate.to_string()),
            VcfVariant {
                chrom: chrom.to_string(),
                position,
                reference: reference.to_string(),
                alternate: alternate.to_string(),
                support,
                depth,
                gene: variant
                    .get("gene")
                    .and_then(Value::as_str)
                    .map(str::to_string),
                clinical_significance: variant
                    .pointer("/annotation/consensus_significance")
                    .and_then(Value::as_str)
                    .or_else(|| {
                        variant
                            .pointer("/annotation/pathogenicity")
                            .and_then(Value::as_str)
                    })
                    .map(str::to_string),
                clinical_conflict: variant
                    .pointer("/annotation/conflict_status")
                    .and_then(Value::as_str)
                    == Some("conflicting"),
                clinical_sources: variant
                    .pointer("/annotation/sources")
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten()
                    .filter_map(Value::as_str)
                    .map(str::to_string)
                    .collect(),
                ci95_low: variant.get("ci95_low").and_then(Value::as_f64),
                ci95_high: variant.get("ci95_high").and_then(Value::as_f64),
                event_id: variant
                    .get("event_id")
                    .and_then(Value::as_str)
                    .map(str::to_string),
                event_type: variant
                    .get("type")
                    .and_then(Value::as_str)
                    .map(str::to_string)
                    .or_else(|| {
                        variant
                            .get("event_id")
                            .and_then(Value::as_str)
                            .map(|_| "SNV".to_string())
                    }),
                reference_support: variant
                    .get("ref_depth")
                    .and_then(Value::as_u64)
                    .unwrap_or(0),
                other_support: variant
                    .get("other_depth")
                    .and_then(Value::as_u64)
                    .unwrap_or(0),
                low_quality: variant
                    .get("low_quality_depth")
                    .and_then(Value::as_u64)
                    .unwrap_or(0),
                conflict: variant
                    .get("conflict_depth")
                    .and_then(Value::as_u64)
                    .unwrap_or(0),
                qc_flags: variant
                    .get("qc_flags")
                    .and_then(Value::as_array)
                    .into_iter()
                    .flatten()
                    .filter_map(Value::as_str)
                    .map(str::to_owned)
                    .collect(),
                numt_assessability: variant
                    .get("numt_assessability")
                    .and_then(Value::as_str)
                    .map(str::to_owned),
                normalization: variant
                    .get("normalization")
                    .and_then(Value::as_str)
                    .map(str::to_owned),
                alt_mapping_quality: variant
                    .pointer("/mapping_quality/alternate/mean")
                    .and_then(Value::as_f64),
                alt_base_quality: variant
                    .pointer("/allele_quality/alternate/mean_phred")
                    .and_then(Value::as_f64),
                strand_support: [
                    "/strand_support/alt_forward",
                    "/strand_support/alt_reverse",
                    "/strand_support/ref_forward",
                    "/strand_support/ref_reverse",
                ]
                .map(|pointer| variant.pointer(pointer).and_then(Value::as_u64))
                .into_iter()
                .collect::<Option<Vec<_>>>()
                .and_then(|values| values.try_into().ok()),
                schema_0_6,
            },
        );
    }

    let mut out = String::new();
    out.push_str("##fileformat=VCFv4.3\n");
    out.push_str("##source=MitoArchitect\n");
    out.push_str(&format!(
        "##contig=<ID={},length={}>\n",
        sanitize_vcf_token(chrom),
        reference_length
    ));
    if schema_0_6 {
        out.push_str(
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate molecule support\">\n",
        );
        out.push_str("##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Reference and alternate molecule support\">\n");
        out.push_str("##INFO=<ID=ODC,Number=1,Type=Integer,Description=\"Callable molecules supporting another allele at the normalized event\">\n");
        out.push_str(
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Event-callable molecule depth\">\n",
        );
        out.push_str("##INFO=<ID=LOWQ,Number=1,Type=Integer,Description=\"Low-quality molecule observations excluded from DP\">\n");
        out.push_str("##INFO=<ID=CONFLICT,Number=1,Type=Integer,Description=\"Conflicting molecule observations excluded from DP\">\n");
        out.push_str("##INFO=<ID=MOLECULE_SUPPORT,Number=1,Type=Integer,Description=\"Physical molecules supporting the alternate event\">\n");
        out.push_str("##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Mean mapping quality of alternate-supporting molecules\">\n");
        out.push_str("##INFO=<ID=BQ,Number=1,Type=Float,Description=\"Mean base quality of alternate-supporting observations when defined\">\n");
        out.push_str("##INFO=<ID=STRAND_SUPPORT,Number=4,Type=Integer,Description=\"ALT forward, ALT reverse, REF forward, REF reverse molecule counts\">\n");
        out.push_str("##INFO=<ID=NUMT_ASSESSABLE,Number=1,Type=String,Description=\"Whether competitive nuclear-plus-mitochondrial evidence makes NUMT specificity assessable\">\n");
        out.push_str("##INFO=<ID=QC_FLAGS,Number=.,Type=String,Description=\"Observed reason-coded QC facts; thresholds are not calibrated filters\">\n");
        out.push_str("##INFO=<ID=NORMALIZATION,Number=1,Type=String,Description=\"Normalized mitochondrial event representation rule\">\n");
    } else {
        out.push_str(
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate read support\">\n",
        );
        out.push_str("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Locus-callable A/C/G/T molecule depth after quality filters\">\n");
    }
    out.push_str("##INFO=<ID=HF,Number=A,Type=Float,Description=\"Heteroplasmy fraction estimated from alternate depth / locus-callable depth\">\n");
    out.push_str("##INFO=<ID=HF_CI95,Number=2,Type=Float,Description=\"Wilson 95% confidence interval for heteroplasmy fraction\">\n");
    out.push_str(
        "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Annotated mitochondrial gene\">\n",
    );
    out.push_str("##INFO=<ID=CLNSIG,Number=1,Type=String,Description=\"Deterministic clinical significance summary; inspect source assertions before interpretation\">\n");
    out.push_str("##INFO=<ID=CLNCONFLICT,Number=0,Type=Flag,Description=\"Source assertions contain incompatible or explicitly conflicting significance groups\">\n");
    out.push_str("##INFO=<ID=CLNSRC,Number=.,Type=String,Description=\"Clinical assertion sources preserved in the JSON result\">\n");
    const EVENT_TYPE_HEADER: &str = "##INFO=<ID=EVENT_TYPE,Number=1,Type=String,Description=\"Schema 0.6 normalized event class\">\n";
    if variants
        .values()
        .any(|variant| variant.event_type.is_some())
    {
        out.push_str(EVENT_TYPE_HEADER);
    }
    out.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype placeholder for haploid mtDNA export\">\n");
    if schema_0_6 {
        out.push_str("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Event-callable molecule depth\">\n");
        out.push_str("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Reference and alternate molecule support\">\n");
    }
    out.push_str("##FORMAT=<ID=HF,Number=1,Type=Float,Description=\"Heteroplasmy fraction\">\n");
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
    out.push_str(&sanitize_vcf_token(sample));
    out.push('\n');

    for variant in variants.values() {
        let mut info = vec![
            format!("AC={}", variant.support),
            format!("DP={}", variant.depth),
            format!("HF={:.6}", variant.heteroplasmy()),
        ];
        if variant.schema_0_6 {
            info.push(format!(
                "AD={},{}",
                variant.reference_support, variant.support
            ));
            info.push(format!("ODC={}", variant.other_support));
            info.push(format!("LOWQ={}", variant.low_quality));
            info.push(format!("CONFLICT={}", variant.conflict));
            info.push(format!("MOLECULE_SUPPORT={}", variant.support));
            if let Some(value) = variant.alt_mapping_quality {
                info.push(format!("MQ={value:.4}"));
            }
            if let Some(value) = variant.alt_base_quality {
                info.push(format!("BQ={value:.4}"));
            }
            if let Some(strands) = variant.strand_support {
                info.push(format!(
                    "STRAND_SUPPORT={},{},{},{}",
                    strands[0], strands[1], strands[2], strands[3]
                ));
            }
            if let Some(value) = &variant.numt_assessability {
                info.push(format!("NUMT_ASSESSABLE={}", sanitize_vcf_token(value)));
            }
            if !variant.qc_flags.is_empty() {
                info.push(format!(
                    "QC_FLAGS={}",
                    variant
                        .qc_flags
                        .iter()
                        .map(|value| sanitize_vcf_token(value))
                        .collect::<Vec<_>>()
                        .join(",")
                ));
            }
            if let Some(value) = &variant.normalization {
                info.push(format!("NORMALIZATION={}", sanitize_vcf_token(value)));
            }
        }
        if let Some(gene) = &variant.gene {
            info.push(format!("GENE={}", sanitize_vcf_token(gene)));
        }
        if let Some(clinical_significance) = &variant.clinical_significance {
            info.push(format!(
                "CLNSIG={}",
                sanitize_vcf_token(clinical_significance)
            ));
        }
        if variant.clinical_conflict {
            info.push("CLNCONFLICT".to_string());
        }
        if !variant.clinical_sources.is_empty() {
            info.push(format!(
                "CLNSRC={}",
                variant
                    .clinical_sources
                    .iter()
                    .map(|value| sanitize_vcf_token(value))
                    .collect::<Vec<_>>()
                    .join(",")
            ));
        }
        if let (Some(low), Some(high)) = (variant.ci95_low, variant.ci95_high) {
            info.push(format!("HF_CI95={low:.6},{high:.6}"));
        }
        if let Some(event_type) = &variant.event_type {
            info.push(format!("EVENT_TYPE={}", sanitize_vcf_token(event_type)));
        }
        let event_id = variant
            .event_id
            .as_deref()
            .map(sanitize_vcf_token)
            .unwrap_or_else(|| ".".to_string());
        if variant.schema_0_6 {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t.\t.\t{}\tGT:DP:AD:HF\t.:{}:{},{}:{:.6}\n",
                sanitize_vcf_token(&variant.chrom),
                variant.position,
                event_id,
                sanitize_vcf_token(&variant.reference),
                sanitize_vcf_token(&variant.alternate),
                info.join(";"),
                variant.depth,
                variant.reference_support,
                variant.support,
                variant.heteroplasmy()
            ));
        } else {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t.\tPASS\t{}\tGT:HF\t0/1:{:.6}\n",
                sanitize_vcf_token(&variant.chrom),
                variant.position,
                event_id,
                sanitize_vcf_token(&variant.reference),
                sanitize_vcf_token(&variant.alternate),
                info.join(";"),
                variant.heteroplasmy()
            ));
        }
    }

    Ok(out)
}

pub(super) fn sanitize_vcf_token(value: &str) -> String {
    value
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '_' | '-' | '.' | ':' | '>') {
                c
            } else {
                '_'
            }
        })
        .collect()
}
