use crate::contract::require_array;
use anyhow::Context;
use anyhow::Result;
use serde_json::Value;

pub(super) fn project_evidence(data: &Value) -> Result<Value> {
    let variants = require_array(data, "/variants")?
        .iter()
        .enumerate()
        .map(|(index, variant)| {
            project_object_fields(
                variant,
                &[
                    "position",
                    "ref",
                    "alt",
                    "alt_depth",
                    "ref_depth",
                    "other_depth",
                    "callable_depth",
                    "heteroplasmy",
                    "ci95_low",
                    "ci95_high",
                    "molecule_support",
                    "strand_support",
                    "strand_bias_delta",
                    "allele_quality",
                    "read_position",
                    "supporting_reads",
                ],
                &format!("/variants/{index}"),
            )
        })
        .collect::<Result<Vec<_>>>()?;

    let reads = require_array(data, "/reads")?
        .iter()
        .enumerate()
        .map(|(index, read)| {
            let mut projected = project_object_fields(
                read,
                &["id", "flags", "mapping_quality"],
                &format!("/reads/{index}"),
            )?;
            let source_snps = read
                .get("snps")
                .and_then(Value::as_array)
                .with_context(|| format!("result contract violation: /reads/{index}/snps"))?;
            let snps = source_snps
                .iter()
                .enumerate()
                .map(|(snp_index, snp)| {
                    project_object_fields(
                        snp,
                        &["position", "ref", "alt"],
                        &format!("/reads/{index}/snps/{snp_index}"),
                    )
                })
                .collect::<Result<Vec<_>>>()?;
            projected
                .as_object_mut()
                .context("internal evidence projection must be an object")?
                .insert("snps".to_owned(), Value::Array(snps));
            Ok(projected)
        })
        .collect::<Result<Vec<_>>>()?;

    Ok(serde_json::json!({
        "schema_version": data.pointer("/metadata/schema_version"),
        "input_alignment_records": data.pointer("/filter_stats/input_alignment_records"),
        "input_molecules": data.pointer("/filter_stats/input_molecules"),
        "passed_reads": data.pointer("/filter_stats/passed_reads"),
        "numt_filtered_reads": data.pointer("/filter_stats/numt_filtered_reads"),
        "variants": variants,
        "reads": reads,
    }))
}

pub(super) fn materialize_observations(data: &Value) -> Result<Vec<Value>> {
    if let Some(observations) = data.get("observations").and_then(Value::as_array) {
        return Ok(observations.clone());
    }
    let pages = require_array(data, "/observation_pages")?;
    let column_names = [
        "molecule_id",
        "event_id",
        "alignment_id",
        "state",
        "observed_allele",
        "base_quality",
        "mapping_quality",
        "strand",
        "evidence_source",
        "read_position",
    ];
    let expected_count = data
        .pointer("/evidence_encoding/observation_count")
        .and_then(Value::as_u64)
        .context("schema 0.6 observation_count must be an integer")?;
    let expected_page_count = data
        .pointer("/evidence_encoding/observation_page_count")
        .and_then(Value::as_u64)
        .context("schema 0.6 observation_page_count must be an integer")?;
    let page_size = data
        .pointer("/evidence_encoding/observation_page_size")
        .and_then(Value::as_u64)
        .filter(|value| (1..=1_000_000).contains(value))
        .context("schema 0.6 observation_page_size is invalid")?;
    if pages.len() as u64 != expected_page_count {
        anyhow::bail!("schema 0.6 observation page cardinality does not match metadata");
    }
    let capacity = usize::try_from(expected_count)
        .context("schema 0.6 observation_count exceeds platform size")?;
    let mut observations = Vec::with_capacity(capacity);
    for (page_index, page) in pages.iter().enumerate() {
        let index = page
            .get("index")
            .and_then(Value::as_u64)
            .with_context(|| format!("/observation_pages/{page_index}/index must be an integer"))?;
        let offset = page
            .get("offset")
            .and_then(Value::as_u64)
            .with_context(|| {
                format!("/observation_pages/{page_index}/offset must be an integer")
            })?;
        let count = page
            .get("count")
            .and_then(Value::as_u64)
            .with_context(|| format!("/observation_pages/{page_index}/count must be an integer"))?;
        if index != page_index as u64 || offset != observations.len() as u64 {
            anyhow::bail!("schema 0.6 observation pages are not contiguous");
        }
        let count =
            usize::try_from(count).context("observation page count exceeds platform size")?;
        if count as u64 > page_size {
            anyhow::bail!("schema 0.6 observation page {page_index} exceeds the page-size limit");
        }
        let columns = page
            .get("columns")
            .and_then(Value::as_object)
            .with_context(|| {
                format!("/observation_pages/{page_index}/columns must be an object")
            })?;
        for name in column_names {
            if columns.get(name).and_then(Value::as_array).map(Vec::len) != Some(count) {
                anyhow::bail!(
                    "schema 0.6 observation page {page_index} column {name} has the wrong length"
                );
            }
        }
        for row_index in 0..count {
            let mut row = serde_json::Map::new();
            row.insert(
                "id".to_owned(),
                Value::String(format!("observation:{}", observations.len())),
            );
            for name in column_names {
                row.insert(name.to_owned(), columns[name][row_index].clone());
            }
            observations.push(Value::Object(row));
        }
    }
    if observations.len() != capacity {
        anyhow::bail!("schema 0.6 observation page count does not match observation_count");
    }
    Ok(observations)
}

pub(super) fn project_evidence_graph(data: &Value) -> Result<Value> {
    let project_array = |pointer: &str, fields: &[&str]| -> Result<Vec<Value>> {
        require_array(data, pointer)?
            .iter()
            .enumerate()
            .map(|(index, value)| {
                project_object_fields(value, fields, &format!("{pointer}/{index}"))
            })
            .collect()
    };

    let alignments = project_array(
        "/alignments",
        &[
            "id",
            "source_record_index",
            "molecule_id",
            "role",
            "selected_representative",
            "flags",
        ],
    )?;
    let molecules = project_array(
        "/molecules",
        &[
            "id",
            "index",
            "identity_policy",
            "assembly_status",
            "primary_candidate_count",
            "ambiguous",
            "analysis_eligible",
            "evidence_eligible",
            "callability_status",
            "callable_bases",
            "callable_fraction",
            "query_length",
            "mean_base_quality",
            "mapping_quality",
            "numt_score",
            "numt_evidence",
            "cluster_id",
            "architecture_assignment",
            "alternate_event_ids",
            "evidence_state_counts",
            "representative_alignment_id",
            "source_qnames",
            "protocol_metadata",
            "protocol_flags",
            "exclusion_reasons",
            "alignment_ids",
            "warnings",
        ],
    )?;
    let events = project_array(
        "/events",
        &[
            "id",
            "index",
            "type",
            "start",
            "end",
            "length",
            "ref",
            "alt",
            "normalization",
            "source_projection",
            "negative_evidence_rule",
            "assessability",
            "component_event_ids",
            "supporting_molecule_ids",
            "evidence_counts",
        ],
    )?;
    let observations = materialize_observations(data)?
        .iter()
        .enumerate()
        .map(|(index, value)| {
            project_object_fields(
                value,
                &[
                    "id",
                    "molecule_id",
                    "event_id",
                    "alignment_id",
                    "state",
                    "observed_allele",
                    "base_quality",
                    "mapping_quality",
                    "strand",
                    "evidence_source",
                ],
                &format!("/observation_pages/materialized/{index}"),
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let variants = project_array(
        "/variants",
        &[
            "event_id",
            "type",
            "position",
            "start",
            "end",
            "length",
            "ref",
            "alt",
            "normalization",
            "negative_evidence_rule",
            "assessability",
            "alt_depth",
            "ref_depth",
            "other_depth",
            "event_absent_depth",
            "low_quality_depth",
            "conflict_depth",
            "callable_depth",
            "heteroplasmy",
            "ci95_low",
            "ci95_high",
            "filter_status",
            "qc_flags",
            "numt_assessability",
            "multi_allelic",
            "homopolymer_context",
            "supporting_molecule_ids",
            "supporting_reads",
        ],
    )?;

    Ok(serde_json::json!({
        "schema_version": data.pointer("/metadata/schema_version"),
        "input_alignment_records": data.pointer("/filter_stats/input_alignment_records"),
        "input_molecules": data.pointer("/filter_stats/input_molecules"),
        "evidence_encoding": data.get("evidence_encoding"),
        "alignments": alignments,
        "molecules": molecules,
        "callability": data.get("callability"),
        "events": events,
        "observations": observations,
        "variants": variants,
        "phase_links": data.get("phase_links"),
        "architecture_inference": data.get("architecture_inference"),
        "architectures": data.get("architectures"),
    }))
}

pub(super) fn project_object_fields(value: &Value, fields: &[&str], path: &str) -> Result<Value> {
    let source = value
        .as_object()
        .with_context(|| format!("result contract violation: {path} must be an object"))?;
    let mut projected = serde_json::Map::new();
    for field in fields {
        projected.insert(
            (*field).to_owned(),
            source
                .get(*field)
                .cloned()
                .with_context(|| format!("result contract violation: {path}/{field} is missing"))?,
        );
    }
    Ok(Value::Object(projected))
}
