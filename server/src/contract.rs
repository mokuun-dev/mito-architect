use crate::architecture_contract::validate_architecture_projection;
use serde_json::Value;
use std::collections::HashMap;
use std::collections::HashSet;

pub(super) fn validate_result_contract(value: &Value) -> Result<(), String> {
    require_object(value, "/metadata")?;
    require_object(value, "/filter_stats")?;
    require_string(value, "/metadata/schema_version")?;
    let schema_version = value
        .pointer("/metadata/schema_version")
        .and_then(Value::as_str)
        .ok_or_else(|| "/metadata/schema_version must be a string".to_string())?;
    if schema_version != "0.5" && schema_version != "0.6" {
        return Err(format!("unsupported result schema {schema_version}"));
    }
    require_string(value, "/metadata/sv_event_schema_version")?;
    require_string(value, "/metadata/complex_sv_event_schema_version")?;
    require_string(value, "/metadata/clinical_annotation_schema_version")?;
    require_string(value, "/metadata/engine_version")?;
    require_u64(value, "/metadata/reference_length")?;
    require_array(value, "/reads")?;
    require_array(value, "/variants")?;
    require_array(value, "/svs")?;
    require_array(value, "/complex_events")?;
    require_array(value, "/clusters")?;
    require_array(value, "/metadata/resources")?;
    require_object(value, "/filter_stats/numt_assessment")?;
    require_u64(value, "/filter_stats/passed_reads")?;
    if schema_version == "0.6" {
        require_object(value, "/evidence_encoding")?;
        require_array(value, "/alignments")?;
        require_array(value, "/molecules")?;
        require_array(value, "/callability")?;
        require_array(value, "/events")?;
        require_array(value, "/observation_pages")?;
        require_array(value, "/phase_links")?;
        require_object(value, "/architecture_inference")?;
        require_array(value, "/architectures")?;
        validate_observation_pages(value)?;
        validate_unified_variant_projection(value)?;
        validate_architecture_projection(value)?;
    }
    Ok(())
}

pub(super) fn validate_observation_pages(value: &Value) -> Result<(), String> {
    if value.get("observations").is_some() {
        return Err("schema 0.6 row observations are forbidden; use observation_pages".to_string());
    }
    let encoding = value
        .pointer("/evidence_encoding")
        .and_then(Value::as_object)
        .ok_or_else(|| "/evidence_encoding must be an object".to_string())?;
    if encoding.get("layout").and_then(Value::as_str) != Some("paged_columnar_molecule_event")
        || encoding.get("observation_storage").and_then(Value::as_str)
            != Some("embedded_columnar_pages")
        || encoding
            .get("phase_molecule_policy")
            .and_then(Value::as_str)
            != Some("evidence_eligible_only")
        || encoding
            .get("phase_molecule_reference")
            .and_then(Value::as_str)
            != Some("molecules[].index")
        || encoding.get("phase_null_model").and_then(Value::as_str)
            != Some("independent_marginals_within_jointly_callable")
    {
        return Err("unsupported schema 0.6 evidence page encoding".to_string());
    }
    let expected_count = encoding
        .get("observation_count")
        .and_then(Value::as_u64)
        .ok_or_else(|| "/evidence_encoding/observation_count must be an integer".to_string())?;
    let expected_page_count = encoding
        .get("observation_page_count")
        .and_then(Value::as_u64)
        .ok_or_else(|| {
            "/evidence_encoding/observation_page_count must be an integer".to_string()
        })?;
    let page_size = encoding
        .get("observation_page_size")
        .and_then(Value::as_u64)
        .filter(|value| (1..=1_000_000).contains(value))
        .ok_or_else(|| "/evidence_encoding/observation_page_size is invalid".to_string())?;
    let pages = value
        .pointer("/observation_pages")
        .and_then(Value::as_array)
        .ok_or_else(|| "/observation_pages must be an array".to_string())?;
    if pages.len() as u64 != expected_page_count {
        return Err("observation page cardinality does not match metadata".to_string());
    }
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
    let mut offset = 0_u64;
    for (page_index, page) in pages.iter().enumerate() {
        if page.get("index").and_then(Value::as_u64) != Some(page_index as u64)
            || page.get("offset").and_then(Value::as_u64) != Some(offset)
        {
            return Err(format!("observation page {page_index} is not contiguous"));
        }
        let count = page
            .get("count")
            .and_then(Value::as_u64)
            .ok_or_else(|| format!("observation page {page_index} count must be an integer"))?;
        if count > page_size {
            return Err(format!(
                "observation page {page_index} exceeds the page-size limit"
            ));
        }
        let count_usize = usize::try_from(count)
            .map_err(|_| format!("observation page {page_index} count is too large"))?;
        let columns = page
            .get("columns")
            .and_then(Value::as_object)
            .ok_or_else(|| format!("observation page {page_index} columns must be an object"))?;
        for name in column_names {
            if columns.get(name).and_then(Value::as_array).map(Vec::len) != Some(count_usize) {
                return Err(format!(
                    "observation page {page_index} column {name} has the wrong length"
                ));
            }
        }
        offset = offset
            .checked_add(count)
            .ok_or_else(|| "observation page count overflow".to_string())?;
    }
    if offset != expected_count {
        return Err("observation page count does not match observation_count".to_string());
    }
    Ok(())
}

pub(super) fn validate_unified_variant_projection(value: &Value) -> Result<(), String> {
    #[derive(Clone)]
    struct EventDefinition {
        event_type: String,
        start: u64,
        end: u64,
        reference: String,
        alternate: String,
        support: HashSet<String>,
    }

    let events = value
        .pointer("/events")
        .and_then(Value::as_array)
        .ok_or_else(|| "/events must be an array".to_string())?;
    let mut definitions = HashMap::<String, EventDefinition>::new();
    let mut event_order = HashMap::<String, usize>::new();
    let mut event_complete_callability = HashMap::<String, bool>::new();
    for (index, event) in events.iter().enumerate() {
        let event_id = event
            .get("id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("/events/{index}/id must be a string"))?;
        let event_type = event
            .get("type")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("/events/{index}/type must be a string"))?;
        if event_order.insert(event_id.to_owned(), index).is_some() {
            return Err(format!("duplicate evidence event {event_id}"));
        }
        let assessability = event
            .get("assessability")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("/events/{index}/assessability must be a string"))?;
        event_complete_callability.insert(
            event_id.to_owned(),
            assessability == "REFERENCE_AND_ALTERNATE",
        );
        if !matches!(event_type, "SNV" | "SMALL_INSERTION" | "SMALL_DELETION") {
            continue;
        }
        let support = event
            .get("supporting_molecule_ids")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("/events/{index}/supporting_molecule_ids must be an array"))?
            .iter()
            .map(|item| {
                item.as_str()
                    .map(str::to_owned)
                    .ok_or_else(|| "supporting molecule ID must be a string".to_string())
            })
            .collect::<Result<HashSet<_>, _>>()?;
        let definition = EventDefinition {
            event_type: event_type.to_owned(),
            start: event
                .get("start")
                .and_then(Value::as_u64)
                .ok_or_else(|| format!("/events/{index}/start must be an integer"))?,
            end: event
                .get("end")
                .and_then(Value::as_u64)
                .ok_or_else(|| format!("/events/{index}/end must be an integer"))?,
            reference: event
                .get("ref")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("/events/{index}/ref must be a string"))?
                .to_owned(),
            alternate: event
                .get("alt")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("/events/{index}/alt must be a string"))?
                .to_owned(),
            support,
        };
        if definitions
            .insert(event_id.to_owned(), definition)
            .is_some()
        {
            return Err(format!("duplicate unified variant event {event_id}"));
        }
    }

    let molecules = value
        .pointer("/molecules")
        .and_then(Value::as_array)
        .ok_or_else(|| "/molecules must be an array".to_string())?;
    let mut evidence_eligible = HashSet::new();
    let mut evidence_eligible_molecules = Vec::<(u64, String)>::new();
    for (index, molecule) in molecules.iter().enumerate() {
        let id = molecule
            .get("id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("/molecules/{index}/id must be a string"))?;
        if molecule.get("index").and_then(Value::as_u64) != Some(index as u64) {
            return Err(format!("/molecules/{index}/index is not contiguous"));
        }
        let eligible = molecule
            .get("evidence_eligible")
            .and_then(Value::as_bool)
            .ok_or_else(|| format!("/molecules/{index}/evidence_eligible must be boolean"))?;
        if eligible {
            if !evidence_eligible.insert(id.to_owned()) {
                return Err(format!("duplicate evidence-eligible molecule {id}"));
            }
            evidence_eligible_molecules.push((index as u64, id.to_owned()));
        }
    }
    let mut counts = HashMap::<String, [u64; 6]>::new();
    let mut phase_observations = HashMap::<String, HashMap<String, String>>::new();
    for (page_index, page) in value
        .pointer("/observation_pages")
        .and_then(Value::as_array)
        .ok_or_else(|| "/observation_pages must be an array".to_string())?
        .iter()
        .enumerate()
    {
        let columns = page
            .get("columns")
            .and_then(Value::as_object)
            .ok_or_else(|| format!("observation page {page_index} columns must be an object"))?;
        let event_ids = columns
            .get("event_id")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("observation page {page_index} event_id must be an array"))?;
        let molecule_ids = columns
            .get("molecule_id")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("observation page {page_index} molecule_id must be an array"))?;
        let states = columns
            .get("state")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("observation page {page_index} state must be an array"))?;
        let sources = columns
            .get("evidence_source")
            .and_then(Value::as_array)
            .ok_or_else(|| {
                format!("observation page {page_index} evidence_source must be an array")
            })?;
        for row in 0..event_ids.len() {
            let molecule_id = molecule_ids[row]
                .as_str()
                .ok_or_else(|| "observation molecule ID must be a string".to_string())?;
            let event_id = event_ids[row]
                .as_str()
                .ok_or_else(|| "observation event ID must be a string".to_string())?;
            let Some(definition) = definitions.get(event_id) else {
                continue;
            };
            let state = states[row]
                .as_str()
                .ok_or_else(|| "observation state must be a string".to_string())?;
            let source = sources[row]
                .as_str()
                .ok_or_else(|| "observation evidence source must be a string".to_string())?;
            if !evidence_eligible.contains(molecule_id) {
                return Err(format!(
                    "observation for event {event_id} belongs to an evidence-ineligible molecule"
                ));
            }
            if phase_observations
                .entry(molecule_id.to_owned())
                .or_default()
                .insert(event_id.to_owned(), state.to_owned())
                .is_some()
            {
                return Err(format!(
                    "duplicate molecule/event observation {molecule_id}/{event_id}"
                ));
            }
            let projected = counts.entry(event_id.to_owned()).or_default();
            match state {
                "ALTERNATE" => projected[0] += 1,
                "REFERENCE" => projected[1] += 1,
                "EVENT_ABSENT" => {
                    projected[3] += 1;
                    if definition.event_type == "SNV" || source == "cigar_alternative_small_indel" {
                        projected[2] += 1;
                    } else {
                        projected[1] += 1;
                    }
                }
                "LOW_QUALITY" => projected[4] += 1,
                "CONFLICT" => projected[5] += 1,
                _ => return Err(format!("unsupported observation state {state}")),
            }
        }
    }

    let variants = value
        .pointer("/variants")
        .and_then(Value::as_array)
        .ok_or_else(|| "/variants must be an array".to_string())?;
    if variants.len() != definitions.len() {
        return Err("unified variant/event cardinality mismatch".to_string());
    }
    let mut seen = HashSet::new();
    for (index, variant) in variants.iter().enumerate() {
        let event_id = variant
            .get("event_id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("/variants/{index}/event_id must be a string"))?;
        if !seen.insert(event_id.to_owned()) {
            return Err(format!("duplicate unified variant {event_id}"));
        }
        let definition = definitions
            .get(event_id)
            .ok_or_else(|| format!("variant {event_id} has no source event"))?;
        if variant.get("type").and_then(Value::as_str) != Some(definition.event_type.as_str())
            || variant.get("position").and_then(Value::as_u64) != Some(definition.start)
            || variant.get("start").and_then(Value::as_u64) != Some(definition.start)
            || variant.get("end").and_then(Value::as_u64) != Some(definition.end)
            || variant.get("ref").and_then(Value::as_str) != Some(definition.reference.as_str())
            || variant.get("alt").and_then(Value::as_str) != Some(definition.alternate.as_str())
        {
            return Err(format!("variant/event identity mismatch for {event_id}"));
        }
        let projected = counts.get(event_id).copied().unwrap_or_default();
        let callable = projected[0]
            .checked_add(projected[1])
            .and_then(|total| total.checked_add(projected[2]))
            .ok_or_else(|| "variant count overflow".to_string())?;
        for (field, expected) in [
            ("alt_depth", projected[0]),
            ("ref_depth", projected[1]),
            ("other_depth", projected[2]),
            ("event_absent_depth", projected[3]),
            ("low_quality_depth", projected[4]),
            ("conflict_depth", projected[5]),
            ("callable_depth", callable),
        ] {
            if variant.get(field).and_then(Value::as_u64) != Some(expected) {
                return Err(format!("variant count mismatch for {event_id}/{field}"));
            }
        }
        let support = variant
            .get("supporting_molecule_ids")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("/variants/{index}/supporting_molecule_ids must be an array"))?
            .iter()
            .map(|item| {
                item.as_str()
                    .map(str::to_owned)
                    .ok_or_else(|| "supporting molecule ID must be a string".to_string())
            })
            .collect::<Result<HashSet<_>, _>>()?;
        if support != definition.support {
            return Err(format!("variant support mismatch for {event_id}"));
        }
    }
    if seen.len() != definitions.len() {
        return Err("variant projection does not cover every SNV/small-indel event".to_string());
    }
    let mut expected_phase_pairs = HashSet::<(String, String)>::new();
    for observations in phase_observations.values() {
        for (alternate_event, state) in observations {
            if state != "ALTERNATE" {
                continue;
            }
            for other_event in observations.keys() {
                if alternate_event == other_event {
                    continue;
                }
                let alternate_order = event_order
                    .get(alternate_event)
                    .ok_or_else(|| "phase alternate event is unresolved".to_string())?;
                let other_order = event_order
                    .get(other_event)
                    .ok_or_else(|| "phase neighbor event is unresolved".to_string())?;
                expected_phase_pairs.insert(if alternate_order < other_order {
                    (alternate_event.clone(), other_event.clone())
                } else {
                    (other_event.clone(), alternate_event.clone())
                });
            }
        }
    }
    for (index, link) in value
        .pointer("/phase_links")
        .and_then(Value::as_array)
        .ok_or_else(|| "/phase_links must be an array".to_string())?
        .iter()
        .enumerate()
    {
        let event_a = link
            .get("event_a_id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("phase link {index} event_a_id must be a string"))?;
        let event_b = link
            .get("event_b_id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("phase link {index} event_b_id must be a string"))?;
        let phase_id = link
            .get("id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("phase link {index} ID must be a string"))?;
        if event_order.get(event_a) >= event_order.get(event_b)
            || phase_id != format!("phase:{event_a}|{event_b}")
            || !expected_phase_pairs.remove(&(event_a.to_owned(), event_b.to_owned()))
        {
            return Err(format!(
                "phase link {index} is not a canonical evidence-derived pair"
            ));
        }
        let complete = event_complete_callability.get(event_a) == Some(&true)
            && event_complete_callability.get(event_b) == Some(&true);
        let expected_assessability = if complete {
            "COMPLETE_FOR_BOTH_EVENTS"
        } else {
            "SUPPORT_CONDITIONED"
        };
        if link.get("assessability").and_then(Value::as_str) != Some(expected_assessability) {
            return Err(format!("phase link {index} assessability is inconsistent"));
        }

        let mut projected = [0_u64; 5];
        let mut support = Vec::new();
        let mut uncertain = Vec::new();
        for (molecule_index, molecule_id) in &evidence_eligible_molecules {
            let Some(observations) = phase_observations.get(molecule_id) else {
                continue;
            };
            let (Some(state_a), Some(state_b)) =
                (observations.get(event_a), observations.get(event_b))
            else {
                continue;
            };
            let callable_a = matches!(state_a.as_str(), "REFERENCE" | "ALTERNATE" | "EVENT_ABSENT");
            let callable_b = matches!(state_b.as_str(), "REFERENCE" | "ALTERNATE" | "EVENT_ABSENT");
            if !callable_a || !callable_b {
                projected[4] += 1;
                uncertain.push(*molecule_index);
                continue;
            }
            match (state_a == "ALTERNATE", state_b == "ALTERNATE") {
                (true, true) => {
                    projected[0] += 1;
                    support.push(*molecule_index);
                }
                (true, false) => projected[1] += 1,
                (false, true) => projected[2] += 1,
                (false, false) => projected[3] += 1,
            }
        }
        let jointly_callable = projected[..4].iter().try_fold(0_u64, |total, value| {
            total
                .checked_add(*value)
                .ok_or_else(|| "phase count overflow".to_string())
        })?;
        for (field, expected_count) in [
            ("both_alternate", projected[0]),
            ("a_alternate_b_absent", projected[1]),
            ("a_absent_b_alternate", projected[2]),
            ("neither_alternate", projected[3]),
            ("jointly_uncertain", projected[4]),
            ("jointly_callable", jointly_callable),
        ] {
            if link.get(field).and_then(Value::as_u64) != Some(expected_count) {
                return Err(format!(
                    "phase link {index}/{field} is not evidence-derived"
                ));
            }
        }
        if server_phase_u64_array(link, "supporting_molecule_indices", index)? != support
            || server_phase_u64_array(link, "uncertain_molecule_indices", index)? != uncertain
        {
            return Err(format!(
                "phase link {index} molecule traceability is inconsistent"
            ));
        }
        let expected_qc = [
            (!complete).then_some("SUPPORT_CONDITIONED"),
            (projected[4] != 0).then_some("UNCERTAIN_COOCCURRENCE_EXCLUDED"),
        ]
        .into_iter()
        .flatten()
        .map(str::to_owned)
        .collect::<Vec<_>>();
        if server_phase_string_array(link, "qc_flags", index)? != expected_qc {
            return Err(format!("phase link {index} QC facts are inconsistent"));
        }
        let observed = if jointly_callable == 0 {
            0.0
        } else {
            projected[0] as f64 / jointly_callable as f64
        };
        let expected = if jointly_callable == 0 {
            0.0
        } else {
            let n = jointly_callable as f64;
            ((projected[0] + projected[1]) as f64 / n) * ((projected[0] + projected[2]) as f64 / n)
        };
        let (ci_low, ci_high) = server_phase_wilson_interval(projected[0], jointly_callable);
        for (field, expected_value) in [
            ("co_alternate_fraction", observed),
            ("co_alternate_ci95_low", ci_low),
            ("co_alternate_ci95_high", ci_high),
            ("expected_co_alternate_fraction", expected),
            ("linkage_delta", observed - expected),
        ] {
            let actual = link
                .get(field)
                .and_then(Value::as_f64)
                .ok_or_else(|| format!("phase link {index}/{field} must be a number"))?;
            if (actual - expected_value).abs() > 2e-9 {
                return Err(format!("phase link {index}/{field} is inconsistent"));
            }
        }
    }
    if !expected_phase_pairs.is_empty() {
        return Err("phase projection omits evidence-derived candidate pairs".to_string());
    }
    Ok(())
}

pub(super) fn server_phase_string_array(
    link: &Value,
    field: &str,
    index: usize,
) -> Result<Vec<String>, String> {
    link.get(field)
        .and_then(Value::as_array)
        .ok_or_else(|| format!("phase link {index}/{field} must be an array"))?
        .iter()
        .map(|value| {
            value
                .as_str()
                .map(str::to_owned)
                .ok_or_else(|| format!("phase link {index}/{field} must contain strings"))
        })
        .collect()
}

pub(super) fn server_phase_u64_array(
    link: &Value,
    field: &str,
    index: usize,
) -> Result<Vec<u64>, String> {
    link.get(field)
        .and_then(Value::as_array)
        .ok_or_else(|| format!("phase link {index}/{field} must be an array"))?
        .iter()
        .map(|value| {
            value
                .as_u64()
                .ok_or_else(|| format!("phase link {index}/{field} must contain integers"))
        })
        .collect()
}

pub(super) fn server_phase_wilson_interval(successes: u64, total: u64) -> (f64, f64) {
    if total == 0 {
        return (0.0, 0.0);
    }
    const Z: f64 = 1.959_963_984_540_054;
    let n = total as f64;
    let p = successes as f64 / n;
    let z2 = Z * Z;
    let denominator = 1.0 + z2 / n;
    let center = (p + z2 / (2.0 * n)) / denominator;
    let margin = Z * ((p * (1.0 - p) / n) + z2 / (4.0 * n * n)).sqrt() / denominator;
    ((center - margin).max(0.0), (center + margin).min(1.0))
}

pub(super) fn require_object(value: &Value, pointer: &str) -> Result<(), String> {
    value
        .pointer(pointer)
        .and_then(Value::as_object)
        .map(|_| ())
        .ok_or_else(|| format!("{pointer} must be an object"))
}

pub(super) fn require_array(value: &Value, pointer: &str) -> Result<(), String> {
    value
        .pointer(pointer)
        .and_then(Value::as_array)
        .map(|_| ())
        .ok_or_else(|| format!("{pointer} must be an array"))
}

pub(super) fn require_string(value: &Value, pointer: &str) -> Result<(), String> {
    value
        .pointer(pointer)
        .and_then(Value::as_str)
        .map(|_| ())
        .ok_or_else(|| format!("{pointer} must be a string"))
}

pub(super) fn require_u64(value: &Value, pointer: &str) -> Result<(), String> {
    value
        .pointer(pointer)
        .and_then(Value::as_u64)
        .map(|_| ())
        .ok_or_else(|| format!("{pointer} must be an unsigned integer"))
}
