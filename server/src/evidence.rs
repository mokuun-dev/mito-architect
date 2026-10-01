use axum::body::Bytes;
use serde::Deserialize;
use serde::Serialize;
use serde_json::Value;
use std::collections::HashMap;
use std::sync::Arc;
use uuid::Uuid;

#[derive(Debug)]
pub(super) struct CompletedArtifacts {
    pub(super) result: Bytes,
    pub(super) result_summary: Bytes,
    pub(super) evidence_pages: Vec<Bytes>,
    pub(super) evidence_search_index: Option<Arc<EvidenceSearchIndex>>,
    pub(super) html_report: Bytes,
}

#[derive(Debug)]
pub(super) struct EvidenceSearchIndex {
    pub(super) strings: Vec<Arc<str>>,
    pub(super) string_ids: HashMap<Arc<str>, u32>,
    pub(super) rows: Vec<EvidenceSearchRecord>,
    pub(super) by_molecule: HashMap<u32, Vec<usize>>,
    pub(super) by_event: HashMap<u32, Vec<usize>>,
    pub(super) by_state: HashMap<u32, Vec<usize>>,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct EvidenceSearchRecord {
    pub(super) global_index: usize,
    pub(super) page_index: usize,
    pub(super) row_index: usize,
    pub(super) molecule_id: u32,
    pub(super) event_id: u32,
    pub(super) alignment_id: u32,
    pub(super) state: u32,
    pub(super) observed_allele: Option<u32>,
    pub(super) base_quality: Option<u64>,
    pub(super) mapping_quality: u64,
    pub(super) strand: u32,
    pub(super) evidence_source: u32,
    pub(super) read_position: Option<f64>,
}

#[derive(Debug, Deserialize)]
pub(super) struct EvidenceSearchQuery {
    pub(super) molecule_id: Option<String>,
    pub(super) event_id: Option<String>,
    pub(super) state: Option<String>,
    pub(super) cursor: Option<usize>,
    pub(super) limit: Option<usize>,
}

#[derive(Debug, Serialize)]
pub(super) struct EvidenceSearchFilters {
    pub(super) molecule_id: Option<String>,
    pub(super) event_id: Option<String>,
    pub(super) state: Option<String>,
}

#[derive(Debug, Serialize)]
pub(super) struct EvidenceSearchResponse {
    pub(super) schema_version: &'static str,
    pub(super) filters: EvidenceSearchFilters,
    pub(super) cursor: usize,
    pub(super) next_cursor: Option<usize>,
    pub(super) total_matches: usize,
    pub(super) rows: Vec<EvidenceSearchResult>,
}

#[derive(Debug, Serialize)]
pub(super) struct EvidenceSearchResult {
    pub(super) id: String,
    pub(super) page_index: usize,
    pub(super) row_index: usize,
    pub(super) molecule_id: String,
    pub(super) event_id: String,
    pub(super) alignment_id: String,
    pub(super) state: String,
    pub(super) observed_allele: Option<String>,
    pub(super) base_quality: Option<u64>,
    pub(super) mapping_quality: u64,
    pub(super) strand: String,
    pub(super) evidence_source: String,
    pub(super) read_position: Option<f64>,
}

pub(super) fn normalized_search_filter(value: Option<String>) -> Option<String> {
    value.and_then(|value| {
        let trimmed = value.trim();
        (!trimmed.is_empty()).then(|| trimmed.to_string())
    })
}

impl EvidenceSearchIndex {
    pub(super) fn search(
        &self,
        filters: EvidenceSearchFilters,
        cursor: usize,
        limit: usize,
    ) -> EvidenceSearchResponse {
        let molecule_id = filters
            .molecule_id
            .as_ref()
            .and_then(|value| self.string_ids.get(value.as_str()).copied());
        let event_id = filters
            .event_id
            .as_ref()
            .and_then(|value| self.string_ids.get(value.as_str()).copied());
        let state = filters
            .state
            .as_ref()
            .and_then(|value| self.string_ids.get(value.as_str()).copied());
        let unknown_filter = (filters.molecule_id.is_some() && molecule_id.is_none())
            || (filters.event_id.is_some() && event_id.is_none())
            || (filters.state.is_some() && state.is_none());
        if unknown_filter {
            return EvidenceSearchResponse {
                schema_version: "1.0",
                filters,
                cursor,
                next_cursor: None,
                total_matches: 0,
                rows: Vec::new(),
            };
        }

        let mut candidate_postings: Option<&[usize]> = None;
        for postings in [
            molecule_id.and_then(|id| self.by_molecule.get(&id)),
            event_id.and_then(|id| self.by_event.get(&id)),
            state.and_then(|id| self.by_state.get(&id)),
        ]
        .into_iter()
        .flatten()
        {
            if match candidate_postings {
                Some(current) => postings.len() < current.len(),
                None => true,
            } {
                candidate_postings = Some(postings);
            }
        }

        let mut total_matches = 0usize;
        let mut response_rows = Vec::with_capacity(limit);
        let mut visit = |record_index: usize| {
            let record = &self.rows[record_index];
            if molecule_id.is_some_and(|id| record.molecule_id != id)
                || event_id.is_some_and(|id| record.event_id != id)
                || state.is_some_and(|id| record.state != id)
            {
                return;
            }
            if total_matches >= cursor && response_rows.len() < limit {
                response_rows.push(self.response_row(record));
            }
            total_matches += 1;
        };
        if let Some(postings) = candidate_postings {
            for &record_index in postings {
                visit(record_index);
            }
        } else {
            for record_index in 0..self.rows.len() {
                visit(record_index);
            }
        }
        let returned_end = cursor.saturating_add(response_rows.len());
        EvidenceSearchResponse {
            schema_version: "1.0",
            filters,
            cursor,
            next_cursor: (returned_end < total_matches).then_some(returned_end),
            total_matches,
            rows: response_rows,
        }
    }

    pub(super) fn response_row(&self, record: &EvidenceSearchRecord) -> EvidenceSearchResult {
        let resolve = |id: u32| self.strings[id as usize].to_string();
        EvidenceSearchResult {
            id: format!("observation:{}", record.global_index),
            page_index: record.page_index,
            row_index: record.row_index,
            molecule_id: resolve(record.molecule_id),
            event_id: resolve(record.event_id),
            alignment_id: resolve(record.alignment_id),
            state: resolve(record.state),
            observed_allele: record.observed_allele.map(resolve),
            base_quality: record.base_quality,
            mapping_quality: record.mapping_quality,
            strand: resolve(record.strand),
            evidence_source: resolve(record.evidence_source),
            read_position: record.read_position,
        }
    }
}

pub(super) fn prepare_completed_artifacts(
    job_id: Uuid,
    value: &Value,
    json_text: String,
    html: String,
) -> anyhow::Result<CompletedArtifacts> {
    let root = value
        .as_object()
        .ok_or_else(|| anyhow::anyhow!("analysis result root is not an object"))?;
    let evidence_pages = value
        .pointer("/observation_pages")
        .and_then(Value::as_array)
        .map(|pages| {
            pages
                .iter()
                .map(serde_json::to_vec)
                .map(|page| page.map(Bytes::from))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()?
        .unwrap_or_default();

    let is_schema_0_6 = value
        .pointer("/metadata/schema_version")
        .and_then(Value::as_str)
        == Some("0.6");
    let evidence_search_index = if is_schema_0_6 {
        Some(Arc::new(build_evidence_search_index(value)?))
    } else {
        None
    };
    let mut summary_root = serde_json::Map::with_capacity(root.len());
    for (name, field) in root {
        if is_schema_0_6 && name == "observation_pages" {
            summary_root.insert(name.clone(), Value::Array(Vec::new()));
        } else {
            summary_root.insert(name.clone(), field.clone());
        }
    }
    let mut summary = Value::Object(summary_root);
    if is_schema_0_6 {
        if !root.contains_key("observation_pages") {
            anyhow::bail!("schema 0.6 summary has no observation_pages");
        }
        let encoding = summary
            .pointer_mut("/evidence_encoding")
            .and_then(Value::as_object_mut)
            .ok_or_else(|| anyhow::anyhow!("schema 0.6 summary has no evidence_encoding"))?;
        encoding.insert(
            "observation_storage".to_string(),
            Value::String("remote_http_pages".to_string()),
        );
        encoding.insert(
            "observation_page_endpoint".to_string(),
            Value::String(format!("/result/{job_id}/evidence/{{page_index}}")),
        );
        encoding.insert(
            "observation_search_endpoint".to_string(),
            Value::String(format!("/result/{job_id}/evidence")),
        );
    }
    Ok(CompletedArtifacts {
        result: Bytes::from(json_text),
        result_summary: Bytes::from(serde_json::to_vec(&summary)?),
        evidence_pages,
        evidence_search_index,
        html_report: Bytes::from(html),
    })
}

pub(super) fn build_evidence_search_index(value: &Value) -> anyhow::Result<EvidenceSearchIndex> {
    let pages = value
        .pointer("/observation_pages")
        .and_then(Value::as_array)
        .ok_or_else(|| anyhow::anyhow!("schema 0.6 result has no observation_pages"))?;
    let expected_count = value
        .pointer("/evidence_encoding/observation_count")
        .and_then(Value::as_u64)
        .and_then(|count| usize::try_from(count).ok())
        .ok_or_else(|| anyhow::anyhow!("schema 0.6 result has invalid observation_count"))?;
    let mut strings = Vec::<Arc<str>>::new();
    let mut string_ids = HashMap::<Arc<str>, u32>::new();
    let mut rows = Vec::<EvidenceSearchRecord>::with_capacity(expected_count);
    let mut by_molecule = HashMap::<u32, Vec<usize>>::new();
    let mut by_event = HashMap::<u32, Vec<usize>>::new();
    let mut by_state = HashMap::<u32, Vec<usize>>::new();

    for page in pages {
        let page_index = required_page_usize(page, "index")?;
        let offset = required_page_usize(page, "offset")?;
        let count = required_page_usize(page, "count")?;
        let columns = page
            .get("columns")
            .and_then(Value::as_object)
            .ok_or_else(|| anyhow::anyhow!("evidence page {page_index} has no columns"))?;
        let molecule_ids = required_page_column(columns, "molecule_id", page_index, count)?;
        let event_ids = required_page_column(columns, "event_id", page_index, count)?;
        let alignment_ids = required_page_column(columns, "alignment_id", page_index, count)?;
        let states = required_page_column(columns, "state", page_index, count)?;
        let observed_alleles = required_page_column(columns, "observed_allele", page_index, count)?;
        let base_qualities = required_page_column(columns, "base_quality", page_index, count)?;
        let mapping_qualities =
            required_page_column(columns, "mapping_quality", page_index, count)?;
        let strands = required_page_column(columns, "strand", page_index, count)?;
        let evidence_sources = required_page_column(columns, "evidence_source", page_index, count)?;
        let read_positions = required_page_column(columns, "read_position", page_index, count)?;

        for row_index in 0..count {
            let molecule_id = intern_required_column_string(
                molecule_ids,
                row_index,
                "molecule_id",
                page_index,
                &mut strings,
                &mut string_ids,
            )?;
            let event_id = intern_required_column_string(
                event_ids,
                row_index,
                "event_id",
                page_index,
                &mut strings,
                &mut string_ids,
            )?;
            let alignment_id = intern_required_column_string(
                alignment_ids,
                row_index,
                "alignment_id",
                page_index,
                &mut strings,
                &mut string_ids,
            )?;
            let state = intern_required_column_string(
                states,
                row_index,
                "state",
                page_index,
                &mut strings,
                &mut string_ids,
            )?;
            let strand = intern_required_column_string(
                strands,
                row_index,
                "strand",
                page_index,
                &mut strings,
                &mut string_ids,
            )?;
            let evidence_source = intern_required_column_string(
                evidence_sources,
                row_index,
                "evidence_source",
                page_index,
                &mut strings,
                &mut string_ids,
            )?;
            let observed_allele = match &observed_alleles[row_index] {
                Value::Null => None,
                Value::String(value) => Some(intern_evidence_string(
                    value,
                    &mut strings,
                    &mut string_ids,
                )?),
                _ => anyhow::bail!(
                    "evidence page {page_index} observed_allele[{row_index}] is invalid"
                ),
            };
            let base_quality =
                optional_u64_column(base_qualities, row_index, "base_quality", page_index)?;
            let mapping_quality = mapping_qualities[row_index].as_u64().ok_or_else(|| {
                anyhow::anyhow!(
                    "evidence page {page_index} mapping_quality[{row_index}] is invalid"
                )
            })?;
            let read_position =
                optional_f64_column(read_positions, row_index, "read_position", page_index)?;
            let record_index = rows.len();
            rows.push(EvidenceSearchRecord {
                global_index: offset + row_index,
                page_index,
                row_index,
                molecule_id,
                event_id,
                alignment_id,
                state,
                observed_allele,
                base_quality,
                mapping_quality,
                strand,
                evidence_source,
                read_position,
            });
            by_molecule
                .entry(molecule_id)
                .or_default()
                .push(record_index);
            by_event.entry(event_id).or_default().push(record_index);
            by_state.entry(state).or_default().push(record_index);
        }
    }
    if rows.len() != expected_count {
        anyhow::bail!(
            "evidence search index row count {} differs from declared observation_count {expected_count}",
            rows.len()
        );
    }
    Ok(EvidenceSearchIndex {
        strings,
        string_ids,
        rows,
        by_molecule,
        by_event,
        by_state,
    })
}

pub(super) fn required_page_usize(page: &Value, name: &str) -> anyhow::Result<usize> {
    page.get(name)
        .and_then(Value::as_u64)
        .and_then(|value| usize::try_from(value).ok())
        .ok_or_else(|| anyhow::anyhow!("evidence page has invalid {name}"))
}

pub(super) fn required_page_column<'a>(
    columns: &'a serde_json::Map<String, Value>,
    name: &str,
    page_index: usize,
    count: usize,
) -> anyhow::Result<&'a Vec<Value>> {
    let values = columns
        .get(name)
        .and_then(Value::as_array)
        .ok_or_else(|| anyhow::anyhow!("evidence page {page_index} has invalid {name} column"))?;
    if values.len() != count {
        anyhow::bail!(
            "evidence page {page_index} {name} column has {} rows, expected {count}",
            values.len()
        );
    }
    Ok(values)
}

pub(super) fn intern_required_column_string(
    values: &[Value],
    row_index: usize,
    name: &str,
    page_index: usize,
    strings: &mut Vec<Arc<str>>,
    string_ids: &mut HashMap<Arc<str>, u32>,
) -> anyhow::Result<u32> {
    let value = values[row_index].as_str().ok_or_else(|| {
        anyhow::anyhow!("evidence page {page_index} {name}[{row_index}] is invalid")
    })?;
    intern_evidence_string(value, strings, string_ids)
}

pub(super) fn intern_evidence_string(
    value: &str,
    strings: &mut Vec<Arc<str>>,
    string_ids: &mut HashMap<Arc<str>, u32>,
) -> anyhow::Result<u32> {
    if let Some(id) = string_ids.get(value) {
        return Ok(*id);
    }
    let id = u32::try_from(strings.len())
        .map_err(|_| anyhow::anyhow!("evidence string dictionary exceeds u32"))?;
    let interned: Arc<str> = Arc::from(value);
    strings.push(interned.clone());
    string_ids.insert(interned, id);
    Ok(id)
}

pub(super) fn optional_u64_column(
    values: &[Value],
    row_index: usize,
    name: &str,
    page_index: usize,
) -> anyhow::Result<Option<u64>> {
    match &values[row_index] {
        Value::Null => Ok(None),
        value => value.as_u64().map(Some).ok_or_else(|| {
            anyhow::anyhow!("evidence page {page_index} {name}[{row_index}] is invalid")
        }),
    }
}

pub(super) fn optional_f64_column(
    values: &[Value],
    row_index: usize,
    name: &str,
    page_index: usize,
) -> anyhow::Result<Option<f64>> {
    match &values[row_index] {
        Value::Null => Ok(None),
        value => value
            .as_f64()
            .filter(|number| number.is_finite())
            .map(Some)
            .ok_or_else(|| {
                anyhow::anyhow!("evidence page {page_index} {name}[{row_index}] is invalid")
            }),
    }
}
