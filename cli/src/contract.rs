use crate::architecture_contract::validate_architecture_contract;
use crate::projection::materialize_observations;
use anyhow::Context;
use anyhow::Result;
use serde_json::Value;
use std::collections::BTreeMap;
use std::collections::BTreeSet;

pub(super) fn validate_result_contract(data: &Value) -> Result<()> {
    require_object(data, "/metadata")?;
    require_object(data, "/filter_stats")?;
    let schema_version = require_string(data, "/metadata/schema_version")?;
    if schema_version != "0.5" && schema_version != "0.6" {
        anyhow::bail!("result contract violation: unsupported schema {schema_version}");
    }
    require_string(data, "/metadata/sv_event_schema_version")?;
    require_string(data, "/metadata/complex_sv_event_schema_version")?;
    require_string(data, "/metadata/clinical_annotation_schema_version")?;
    require_string(data, "/metadata/engine_version")?;
    require_u64(data, "/metadata/reference_length")?;
    require_array(data, "/reads")?;
    require_array(data, "/variants")?;
    require_array(data, "/svs")?;
    require_array(data, "/complex_events")?;
    require_array(data, "/clusters")?;
    require_array(data, "/metadata/resources")?;
    require_object(data, "/filter_stats/numt_assessment")?;
    require_u64(data, "/filter_stats/passed_reads")?;
    if schema_version == "0.6" {
        require_object(data, "/evidence_encoding")?;
        require_array(data, "/alignments")?;
        require_array(data, "/molecules")?;
        require_array(data, "/callability")?;
        require_array(data, "/events")?;
        require_array(data, "/observation_pages")?;
        require_array(data, "/phase_links")?;
        require_object(data, "/architecture_inference")?;
        require_array(data, "/architectures")?;
        validate_schema_0_6_contract(data)?;
    }
    Ok(())
}

pub(super) fn validate_schema_0_6_contract(data: &Value) -> Result<()> {
    if data.get("observations").is_some() {
        anyhow::bail!("schema 0.6 row observations are forbidden; use observation_pages");
    }
    let reference_length = require_u64(data, "/metadata/reference_length")?;
    let encoding = require_object(data, "/evidence_encoding")?;
    if encoding.get("layout").and_then(Value::as_str) != Some("paged_columnar_molecule_event")
        || encoding.get("observation_storage").and_then(Value::as_str)
            != Some("embedded_columnar_pages")
        || encoding.get("missing_pair_state").and_then(Value::as_str) != Some("NOT_CALLABLE")
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
        anyhow::bail!("schema 0.6 evidence encoding is unsupported");
    }
    let observation_limit = encoding
        .get("observation_limit")
        .and_then(Value::as_u64)
        .context("schema 0.6 observation_limit must be an integer")?;
    let phase_link_limit = encoding
        .get("phase_link_limit")
        .and_then(Value::as_u64)
        .context("schema 0.6 phase_link_limit must be an integer")?;
    let page_size = encoding
        .get("observation_page_size")
        .and_then(Value::as_u64)
        .context("schema 0.6 observation_page_size must be an integer")?;
    if page_size == 0 || page_size > 1_000_000 {
        anyhow::bail!("schema 0.6 observation_page_size is outside the supported envelope");
    }

    let alignments = require_array(data, "/alignments")?;
    let molecules = require_array(data, "/molecules")?;
    let callability = require_array(data, "/callability")?;
    let events = require_array(data, "/events")?;
    let observations = materialize_observations(data)?;
    let phase_links = require_array(data, "/phase_links")?;
    if observations.len() as u64 > observation_limit || phase_links.len() as u64 > phase_link_limit
    {
        anyhow::bail!("schema 0.6 evidence payload exceeds its declared resource limit");
    }
    if molecules.len() != callability.len() {
        anyhow::bail!("schema 0.6 callability/molecule cardinality mismatch");
    }

    let alignment_ids = alignments
        .iter()
        .enumerate()
        .map(|(index, alignment)| {
            alignment
                .get("id")
                .and_then(Value::as_str)
                .map(str::to_owned)
                .with_context(|| format!("/alignments/{index}/id must be a string"))
        })
        .collect::<Result<BTreeSet<_>>>()?;
    if alignment_ids.len() != alignments.len() {
        anyhow::bail!("schema 0.6 alignment IDs are not unique");
    }
    let mut alignment_owners = BTreeMap::new();
    for (index, alignment) in alignments.iter().enumerate() {
        let alignment_id = alignment
            .get("id")
            .and_then(Value::as_str)
            .with_context(|| format!("/alignments/{index}/id must be a string"))?;
        let molecule_id = alignment
            .get("molecule_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/alignments/{index}/molecule_id must be a string"))?;
        alignment_owners.insert(alignment_id, molecule_id);
    }
    let molecule_ids = molecules
        .iter()
        .enumerate()
        .map(|(index, molecule)| {
            molecule
                .get("id")
                .and_then(Value::as_str)
                .map(str::to_owned)
                .with_context(|| format!("/molecules/{index}/id must be a string"))
        })
        .collect::<Result<BTreeSet<_>>>()?;
    if molecule_ids.len() != molecules.len() {
        anyhow::bail!("schema 0.6 molecule IDs are not unique");
    }
    let mut projected_alignments = BTreeSet::new();
    let mut phase_eligible_molecule_ids = BTreeSet::new();
    let mut phase_eligible_molecules = BTreeMap::<usize, String>::new();
    for (index, molecule) in molecules.iter().enumerate() {
        let molecule_id = molecule
            .get("id")
            .and_then(Value::as_str)
            .with_context(|| format!("/molecules/{index}/id must be a string"))?;
        if molecule.get("index").and_then(Value::as_u64) != Some(index as u64) {
            anyhow::bail!("schema 0.6 molecule index is not contiguous at {molecule_id}");
        }
        let representative = molecule
            .get("representative_alignment_id")
            .and_then(Value::as_str)
            .with_context(|| {
                format!("/molecules/{index}/representative_alignment_id must be a string")
            })?;
        let identity_policy = molecule
            .get("identity_policy")
            .and_then(Value::as_str)
            .filter(|value| !value.is_empty())
            .with_context(|| format!("/molecules/{index}/identity_policy must be non-empty"))?;
        let source_qnames = molecule
            .get("source_qnames")
            .and_then(Value::as_array)
            .with_context(|| format!("/molecules/{index}/source_qnames must be an array"))?;
        let source_qname_set = source_qnames
            .iter()
            .map(|value| value.as_str().context("source QNAME must be a string"))
            .collect::<Result<BTreeSet<_>>>()?;
        if source_qname_set.is_empty() || source_qname_set.len() != source_qnames.len() {
            anyhow::bail!("schema 0.6 molecule {molecule_id} has invalid source QNAMEs");
        }
        let protocol_metadata = molecule
            .get("protocol_metadata")
            .and_then(Value::as_object)
            .with_context(|| format!("/molecules/{index}/protocol_metadata must be an object"))?;
        if protocol_metadata
            .values()
            .any(|value| value.as_str().is_none())
        {
            anyhow::bail!("schema 0.6 molecule {molecule_id} protocol metadata must be strings");
        }
        if identity_policy.starts_with("sam_tag:")
            && protocol_metadata
                .get("molecule_id_tag")
                .and_then(Value::as_str)
                .is_none()
        {
            anyhow::bail!("schema 0.6 tagged molecule {molecule_id} lacks tag provenance");
        }
        let analysis_eligible = molecule
            .get("analysis_eligible")
            .and_then(Value::as_bool)
            .with_context(|| format!("/molecules/{index}/analysis_eligible must be boolean"))?;
        let evidence_eligible = molecule
            .get("evidence_eligible")
            .and_then(Value::as_bool)
            .with_context(|| format!("/molecules/{index}/evidence_eligible must be boolean"))?;
        if evidence_eligible && !analysis_eligible {
            anyhow::bail!(
                "schema 0.6 molecule {molecule_id} cannot be evidence-eligible while analysis-ineligible"
            );
        }
        if evidence_eligible {
            phase_eligible_molecule_ids.insert(molecule_id.to_owned());
            phase_eligible_molecules.insert(index, molecule_id.to_owned());
        }
        let exclusion_reasons = molecule
            .get("exclusion_reasons")
            .and_then(Value::as_array)
            .with_context(|| format!("/molecules/{index}/exclusion_reasons must be an array"))?;
        let exclusion_set = exclusion_reasons
            .iter()
            .map(|value| value.as_str().context("exclusion reason must be a string"))
            .collect::<Result<BTreeSet<_>>>()?;
        if exclusion_set.len() != exclusion_reasons.len()
            || (analysis_eligible && !exclusion_set.is_empty())
            || (!analysis_eligible && exclusion_set.is_empty())
        {
            anyhow::bail!("schema 0.6 molecule {molecule_id} exclusion contract is inconsistent");
        }
        let mut molecule_alignments = BTreeSet::new();
        for alignment_id in molecule
            .get("alignment_ids")
            .and_then(Value::as_array)
            .with_context(|| format!("/molecules/{index}/alignment_ids must be an array"))?
        {
            let alignment_id = alignment_id
                .as_str()
                .context("molecule alignment ID must be a string")?;
            if alignment_owners.get(alignment_id).copied() != Some(molecule_id)
                || !molecule_alignments.insert(alignment_id)
                || !projected_alignments.insert(alignment_id)
            {
                anyhow::bail!(
                    "schema 0.6 molecule/alignment ownership mismatch for {alignment_id}"
                );
            }
        }
        if !molecule_alignments.contains(representative) {
            anyhow::bail!("schema 0.6 representative alignment is unresolved for {molecule_id}");
        }
    }
    if projected_alignments != alignment_ids.iter().map(String::as_str).collect() {
        anyhow::bail!("schema 0.6 not every alignment resolves to exactly one molecule");
    }
    let mut callability_molecules = BTreeSet::new();
    for (index, summary) in callability.iter().enumerate() {
        let molecule_id = summary
            .get("molecule_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/callability/{index}/molecule_id must be a string"))?;
        if !molecule_ids.contains(molecule_id) {
            anyhow::bail!("schema 0.6 callability references unknown molecule {molecule_id}");
        }
        if !callability_molecules.insert(molecule_id) {
            anyhow::bail!("schema 0.6 callability duplicates molecule {molecule_id}");
        }
        let ranges = summary
            .get("ranges")
            .and_then(Value::as_array)
            .with_context(|| format!("/callability/{index}/ranges must be an array"))?;
        let mut previous_end = 0_u64;
        for (range_index, range) in ranges.iter().enumerate() {
            let start = range
                .get("start")
                .and_then(Value::as_u64)
                .with_context(|| {
                    format!("/callability/{index}/ranges/{range_index}/start must be an integer")
                })?;
            let end = range.get("end").and_then(Value::as_u64).with_context(|| {
                format!("/callability/{index}/ranges/{range_index}/end must be an integer")
            })?;
            if start == 0 || end < start || end > reference_length || start <= previous_end {
                anyhow::bail!(
                    "schema 0.6 callability ranges are invalid at molecule {molecule_id}"
                );
            }
            previous_end = end;
        }
        for (alignment_index, alignment) in summary
            .get("alignments")
            .and_then(Value::as_array)
            .with_context(|| format!("/callability/{index}/alignments must be an array"))?
            .iter()
            .enumerate()
        {
            let alignment_id = alignment
                .get("alignment_id")
                .and_then(Value::as_str)
                .with_context(|| format!("/callability/{index}/alignments/{alignment_index}/alignment_id must be a string"))?;
            if alignment_owners.get(alignment_id).copied() != Some(molecule_id) {
                anyhow::bail!(
                    "schema 0.6 callability alignment {alignment_id} belongs to another molecule"
                );
            }
        }
    }

    let mut event_ids = BTreeSet::new();
    let mut event_order = BTreeMap::<String, usize>::new();
    let mut event_complete_callability = BTreeMap::<String, bool>::new();
    let mut variant_events = BTreeMap::<String, (String, u64, u64, String, String)>::new();
    for (index, event) in events.iter().enumerate() {
        let event_id = event
            .get("id")
            .and_then(Value::as_str)
            .with_context(|| format!("/events/{index}/id must be a string"))?;
        if event.get("index").and_then(Value::as_u64) != Some(index as u64) {
            anyhow::bail!("schema 0.6 event index is not contiguous at {event_id}");
        }
        if !event_ids.insert(event_id.to_owned()) {
            anyhow::bail!("schema 0.6 event ID is duplicated: {event_id}");
        }
        event_order.insert(event_id.to_owned(), index);
        let event_type = event
            .get("type")
            .and_then(Value::as_str)
            .with_context(|| format!("/events/{index}/type must be a string"))?;
        if matches!(event_type, "SNV" | "SMALL_INSERTION" | "SMALL_DELETION") {
            variant_events.insert(
                event_id.to_owned(),
                (
                    event_type.to_owned(),
                    event
                        .get("start")
                        .and_then(Value::as_u64)
                        .with_context(|| {
                            format!("/events/{index}/start must be an integer for a variant event")
                        })?,
                    event.get("end").and_then(Value::as_u64).with_context(|| {
                        format!("/events/{index}/end must be an integer for a variant event")
                    })?,
                    event
                        .get("ref")
                        .and_then(Value::as_str)
                        .with_context(|| format!("/events/{index}/ref must be a string"))?
                        .to_owned(),
                    event
                        .get("alt")
                        .and_then(Value::as_str)
                        .with_context(|| format!("/events/{index}/alt must be a string"))?
                        .to_owned(),
                ),
            );
        }
        let assessability = event
            .get("assessability")
            .and_then(Value::as_str)
            .with_context(|| format!("/events/{index}/assessability must be a string"))?;
        let negative_rule = event
            .get("negative_evidence_rule")
            .and_then(Value::as_str)
            .filter(|value| !value.is_empty())
            .with_context(|| format!("/events/{index}/negative_evidence_rule must be non-empty"))?;
        event_complete_callability.insert(
            event_id.to_owned(),
            assessability == "REFERENCE_AND_ALTERNATE",
        );
        if (assessability == "ALTERNATE_SUPPORT_ONLY")
            != (negative_rule == "support_only_no_negative_inference")
        {
            anyhow::bail!("schema 0.6 event {event_id} has an inconsistent negative-evidence rule");
        }
        for molecule_id in event
            .get("supporting_molecule_ids")
            .and_then(Value::as_array)
            .with_context(|| format!("/events/{index}/supporting_molecule_ids must be an array"))?
        {
            let molecule_id = molecule_id
                .as_str()
                .context("supporting molecule ID must be a string")?;
            if !molecule_ids.contains(molecule_id) {
                anyhow::bail!(
                    "schema 0.6 event {event_id} references unknown molecule {molecule_id}"
                );
            }
        }
    }

    let mut molecule_event_pairs = BTreeSet::new();
    let mut projected_support: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut projected_counts: BTreeMap<String, [u64; 5]> = BTreeMap::new();
    let mut projected_variant_counts: BTreeMap<String, [u64; 6]> = BTreeMap::new();
    let mut projected_molecule_counts: BTreeMap<String, [u64; 5]> = BTreeMap::new();
    let mut projected_molecule_alternates: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let mut phase_observations = BTreeMap::<String, BTreeMap<String, String>>::new();
    for (index, observation) in observations.iter().enumerate() {
        let expected_id = format!("observation:{index}");
        if observation.get("id").and_then(Value::as_str) != Some(expected_id.as_str()) {
            anyhow::bail!("schema 0.6 observation ID is not contiguous at index {index}");
        }
        let molecule_id = observation
            .get("molecule_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/observations/{index}/molecule_id must be a string"))?;
        let event_id = observation
            .get("event_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/observations/{index}/event_id must be a string"))?;
        let alignment_id = observation
            .get("alignment_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/observations/{index}/alignment_id must be a string"))?;
        let state = observation
            .get("state")
            .and_then(Value::as_str)
            .with_context(|| format!("/observations/{index}/state must be a string"))?;
        if !molecule_ids.contains(molecule_id)
            || !event_ids.contains(event_id)
            || !alignment_ids.contains(alignment_id)
        {
            anyhow::bail!("schema 0.6 observation {index} has an unresolved reference");
        }
        if !phase_eligible_molecule_ids.contains(molecule_id) {
            anyhow::bail!(
                "schema 0.6 observation {index} belongs to a molecule excluded from evidence"
            );
        }
        if !matches!(
            state,
            "REFERENCE" | "ALTERNATE" | "EVENT_ABSENT" | "LOW_QUALITY" | "CONFLICT"
        ) {
            anyhow::bail!("schema 0.6 observation {index} has unsupported state {state}");
        }
        let evidence_source = observation
            .get("evidence_source")
            .and_then(Value::as_str)
            .filter(|value| !value.is_empty())
            .with_context(|| format!("/observations/{index}/evidence_source must be non-empty"))?;
        if !molecule_event_pairs.insert((molecule_id.to_owned(), event_id.to_owned())) {
            anyhow::bail!("schema 0.6 molecule/event observation is duplicated");
        }
        phase_observations
            .entry(molecule_id.to_owned())
            .or_default()
            .insert(event_id.to_owned(), state.to_owned());
        if state == "ALTERNATE" {
            projected_support
                .entry(event_id.to_owned())
                .or_default()
                .insert(molecule_id.to_owned());
        }
        let count_index = match state {
            "ALTERNATE" => 0,
            "REFERENCE" => 1,
            "EVENT_ABSENT" => 2,
            "LOW_QUALITY" => 3,
            "CONFLICT" => 4,
            _ => unreachable!("observation state was validated above"),
        };
        let count = &mut projected_counts.entry(event_id.to_owned()).or_default()[count_index];
        *count = count
            .checked_add(1)
            .context("schema 0.6 observation count overflow")?;
        let molecule_count = &mut projected_molecule_counts
            .entry(molecule_id.to_owned())
            .or_default()[count_index];
        *molecule_count = molecule_count
            .checked_add(1)
            .context("schema 0.6 molecule observation count overflow")?;
        if state == "ALTERNATE" {
            projected_molecule_alternates
                .entry(molecule_id.to_owned())
                .or_default()
                .insert(event_id.to_owned());
        }
        if let Some((event_type, _, _, _, _)) = variant_events.get(event_id) {
            let counts = projected_variant_counts
                .entry(event_id.to_owned())
                .or_default();
            match state {
                "ALTERNATE" => counts[0] += 1,
                "REFERENCE" => counts[1] += 1,
                "EVENT_ABSENT" => {
                    counts[3] += 1;
                    if event_type == "SNV" || evidence_source == "cigar_alternative_small_indel" {
                        counts[2] += 1;
                    } else {
                        counts[1] += 1;
                    }
                }
                "LOW_QUALITY" => counts[4] += 1,
                "CONFLICT" => counts[5] += 1,
                _ => unreachable!("observation state was validated above"),
            }
        }
    }
    let mut verified_event_support = BTreeMap::<String, BTreeSet<String>>::new();
    for (index, event) in events.iter().enumerate() {
        let event_id = event.get("id").and_then(Value::as_str).unwrap_or_default();
        let declared = event
            .get("supporting_molecule_ids")
            .and_then(Value::as_array)
            .with_context(|| format!("/events/{index}/supporting_molecule_ids must be an array"))?
            .iter()
            .map(|value| {
                value
                    .as_str()
                    .map(str::to_owned)
                    .context("supporting molecule ID must be a string")
            })
            .collect::<Result<BTreeSet<_>>>()?;
        if declared != projected_support.remove(event_id).unwrap_or_default() {
            anyhow::bail!("schema 0.6 support projection mismatch for {event_id}");
        }
        verified_event_support.insert(event_id.to_owned(), declared);
        let counts = projected_counts.remove(event_id).unwrap_or_default();
        let declared_counts = event
            .get("evidence_counts")
            .and_then(Value::as_object)
            .with_context(|| format!("/events/{index}/evidence_counts must be an object"))?;
        let expected = [
            ("alternate", counts[0]),
            ("reference", counts[1]),
            ("event_absent", counts[2]),
            ("low_quality", counts[3]),
            ("conflict", counts[4]),
            ("callable", counts[0] + counts[1] + counts[2]),
        ];
        for (field, expected_count) in expected {
            if declared_counts.get(field).and_then(Value::as_u64) != Some(expected_count) {
                anyhow::bail!("schema 0.6 evidence count mismatch for {event_id}/{field}");
            }
        }
    }

    for (index, molecule) in molecules.iter().enumerate() {
        let molecule_id = molecule
            .get("id")
            .and_then(Value::as_str)
            .with_context(|| format!("/molecules/{index}/id must be a string"))?;
        let projected = projected_molecule_counts
            .remove(molecule_id)
            .unwrap_or_default();
        let declared = molecule
            .get("evidence_state_counts")
            .and_then(Value::as_object)
            .with_context(|| {
                format!("/molecules/{index}/evidence_state_counts must be an object")
            })?;
        for (field, expected) in [
            ("alternate", projected[0]),
            ("reference", projected[1]),
            ("event_absent", projected[2]),
            ("low_quality", projected[3]),
            ("conflict", projected[4]),
        ] {
            if declared.get(field).and_then(Value::as_u64) != Some(expected) {
                anyhow::bail!(
                    "schema 0.6 molecule evidence count mismatch for {molecule_id}/{field}"
                );
            }
        }
        let alternate_event_ids = molecule
            .get("alternate_event_ids")
            .and_then(Value::as_array)
            .with_context(|| format!("/molecules/{index}/alternate_event_ids must be an array"))?
            .iter()
            .map(|value| {
                value
                    .as_str()
                    .map(str::to_owned)
                    .context("molecule alternate event ID must be a string")
            })
            .collect::<Result<BTreeSet<_>>>()?;
        if alternate_event_ids
            != projected_molecule_alternates
                .remove(molecule_id)
                .unwrap_or_default()
        {
            anyhow::bail!(
                "schema 0.6 molecule alternate-event projection mismatch for {molecule_id}"
            );
        }
        for field in ["query_length", "mapping_quality"] {
            if molecule.get(field).and_then(Value::as_u64).is_none() {
                anyhow::bail!("schema 0.6 molecule {molecule_id}/{field} must be an integer");
            }
        }
        for field in ["mean_base_quality", "numt_score"] {
            if molecule.get(field).and_then(Value::as_f64).is_none() {
                anyhow::bail!("schema 0.6 molecule {molecule_id}/{field} must be a number");
            }
        }
        require_array(molecule, "/numt_evidence")?;
    }

    let variants = require_array(data, "/variants")?;
    if variants.len() != variant_events.len() {
        anyhow::bail!("schema 0.6 unified variant/event cardinality mismatch");
    }
    let mut projected_variant_ids = BTreeSet::new();
    for (index, variant) in variants.iter().enumerate() {
        let event_id = variant
            .get("event_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/variants/{index}/event_id must be a string"))?;
        if !projected_variant_ids.insert(event_id.to_owned()) {
            anyhow::bail!("schema 0.6 variant event ID is duplicated: {event_id}");
        }
        let (event_type, start, end, reference, alternate) = variant_events
            .get(event_id)
            .with_context(|| format!("variant {event_id} has no normalized source event"))?;
        if variant.get("type").and_then(Value::as_str) != Some(event_type.as_str())
            || variant.get("position").and_then(Value::as_u64) != Some(*start)
            || variant.get("start").and_then(Value::as_u64) != Some(*start)
            || variant.get("end").and_then(Value::as_u64) != Some(*end)
            || variant.get("ref").and_then(Value::as_str) != Some(reference.as_str())
            || variant.get("alt").and_then(Value::as_str) != Some(alternate.as_str())
        {
            anyhow::bail!("schema 0.6 variant/event identity mismatch for {event_id}");
        }
        let counts = projected_variant_counts
            .get(event_id)
            .copied()
            .unwrap_or_default();
        let callable = counts[0]
            .checked_add(counts[1])
            .and_then(|value| value.checked_add(counts[2]))
            .context("schema 0.6 variant callable count overflow")?;
        for (field, expected) in [
            ("alt_depth", counts[0]),
            ("ref_depth", counts[1]),
            ("other_depth", counts[2]),
            ("event_absent_depth", counts[3]),
            ("low_quality_depth", counts[4]),
            ("conflict_depth", counts[5]),
            ("callable_depth", callable),
        ] {
            if variant.get(field).and_then(Value::as_u64) != Some(expected) {
                anyhow::bail!("schema 0.6 variant count mismatch for {event_id}/{field}");
            }
        }
        let declared_support = variant
            .get("supporting_molecule_ids")
            .and_then(Value::as_array)
            .with_context(|| format!("/variants/{index}/supporting_molecule_ids must be an array"))?
            .iter()
            .map(|value| {
                value
                    .as_str()
                    .map(str::to_owned)
                    .context("variant supporting molecule ID must be a string")
            })
            .collect::<Result<BTreeSet<_>>>()?;
        if verified_event_support.get(event_id) != Some(&declared_support) {
            anyhow::bail!("schema 0.6 variant support mismatch for {event_id}");
        }
        let heteroplasmy = variant
            .get("heteroplasmy")
            .and_then(Value::as_f64)
            .with_context(|| format!("/variants/{index}/heteroplasmy must be a number"))?;
        let expected_hf = if callable == 0 {
            0.0
        } else {
            counts[0] as f64 / callable as f64
        };
        if (heteroplasmy - expected_hf).abs() > 2e-9 {
            anyhow::bail!("schema 0.6 variant HF mismatch for {event_id}");
        }
        require_string(variant, "/filter_status")?;
        require_array(variant, "/qc_flags")?;
        require_string(variant, "/numt_assessability")?;
    }
    if projected_variant_ids != variant_events.keys().cloned().collect() {
        anyhow::bail!("schema 0.6 variant projection does not cover every SNV/small-indel event");
    }
    let mut expected_phase_pairs = BTreeSet::<(String, String)>::new();
    for observations in phase_observations.values() {
        for (alternate_event, state) in observations {
            if state != "ALTERNATE" {
                continue;
            }
            for other_event in observations.keys() {
                if other_event == alternate_event {
                    continue;
                }
                let alternate_order = *event_order
                    .get(alternate_event)
                    .context("phase alternate event order is unresolved")?;
                let other_order = *event_order
                    .get(other_event)
                    .context("phase neighbor event order is unresolved")?;
                let pair = if alternate_order < other_order {
                    (alternate_event.clone(), other_event.clone())
                } else {
                    (other_event.clone(), alternate_event.clone())
                };
                expected_phase_pairs.insert(pair);
            }
        }
    }
    let mut phase_ids = BTreeSet::new();
    for (index, link) in phase_links.iter().enumerate() {
        let phase_id = link
            .get("id")
            .and_then(Value::as_str)
            .with_context(|| format!("/phase_links/{index}/id must be a string"))?;
        if !phase_ids.insert(phase_id) {
            anyhow::bail!("schema 0.6 phase link ID is duplicated: {phase_id}");
        }
        let event_a = link
            .get("event_a_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/phase_links/{index}/event_a_id must be a string"))?;
        let event_b = link
            .get("event_b_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/phase_links/{index}/event_b_id must be a string"))?;
        if !event_ids.contains(event_a) || !event_ids.contains(event_b) || event_a == event_b {
            anyhow::bail!("schema 0.6 phase link {index} has invalid event references");
        }
        if event_order.get(event_a) >= event_order.get(event_b)
            || phase_id != format!("phase:{event_a}|{event_b}")
        {
            anyhow::bail!("schema 0.6 phase link {index} is not in canonical event order");
        }
        if !expected_phase_pairs.remove(&(event_a.to_owned(), event_b.to_owned())) {
            anyhow::bail!(
                "schema 0.6 phase link {index} is not derived from an alternate observation"
            );
        }
        let complete = event_complete_callability.get(event_a) == Some(&true)
            && event_complete_callability.get(event_b) == Some(&true);
        let expected_assessability = if complete {
            "COMPLETE_FOR_BOTH_EVENTS"
        } else {
            "SUPPORT_CONDITIONED"
        };
        if link.get("assessability").and_then(Value::as_str) != Some(expected_assessability) {
            anyhow::bail!("schema 0.6 phase assessability mismatch at link {index}");
        }

        let mut counts = [0_u64; 5];
        let mut supporting_molecule_indices = Vec::new();
        let mut uncertain_molecule_indices = Vec::new();
        for (molecule_index, molecule_id) in &phase_eligible_molecules {
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
                counts[4] = counts[4]
                    .checked_add(1)
                    .context("phase uncertain count overflow")?;
                uncertain_molecule_indices.push(*molecule_index as u64);
                continue;
            }
            match (
                state_a.as_str() == "ALTERNATE",
                state_b.as_str() == "ALTERNATE",
            ) {
                (true, true) => {
                    counts[0] += 1;
                    supporting_molecule_indices.push(*molecule_index as u64);
                }
                (true, false) => counts[1] += 1,
                (false, true) => counts[2] += 1,
                (false, false) => counts[3] += 1,
            }
        }
        let jointly_callable = counts[0]
            .checked_add(counts[1])
            .and_then(|value| value.checked_add(counts[2]))
            .and_then(|value| value.checked_add(counts[3]))
            .context("phase jointly-callable count overflow")?;
        for (field, expected_count) in [
            ("both_alternate", counts[0]),
            ("a_alternate_b_absent", counts[1]),
            ("a_absent_b_alternate", counts[2]),
            ("neither_alternate", counts[3]),
            ("jointly_uncertain", counts[4]),
            ("jointly_callable", jointly_callable),
        ] {
            if link.get(field).and_then(Value::as_u64) != Some(expected_count) {
                anyhow::bail!(
                    "schema 0.6 phase observation projection mismatch at link {index}/{field}"
                );
            }
        }
        let declared_support = phase_u64_array(link, "supporting_molecule_indices", index)?;
        let declared_uncertain = phase_u64_array(link, "uncertain_molecule_indices", index)?;
        if declared_support != supporting_molecule_indices
            || declared_uncertain != uncertain_molecule_indices
        {
            anyhow::bail!("schema 0.6 phase molecule traceability mismatch at link {index}");
        }
        let expected_qc_flags = [
            (!complete).then_some("SUPPORT_CONDITIONED"),
            (counts[4] != 0).then_some("UNCERTAIN_COOCCURRENCE_EXCLUDED"),
        ]
        .into_iter()
        .flatten()
        .map(str::to_owned)
        .collect::<Vec<_>>();
        if phase_string_array(link, "qc_flags", index)? != expected_qc_flags {
            anyhow::bail!("schema 0.6 phase QC facts mismatch at link {index}");
        }

        let observed = if jointly_callable == 0 {
            0.0
        } else {
            counts[0] as f64 / jointly_callable as f64
        };
        let expected = if jointly_callable == 0 {
            0.0
        } else {
            let n = jointly_callable as f64;
            ((counts[0] + counts[1]) as f64 / n) * ((counts[0] + counts[2]) as f64 / n)
        };
        let (ci_low, ci_high) = phase_wilson_interval(counts[0], jointly_callable);
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
                .with_context(|| format!("phase {field} must be a number"))?;
            if (actual - expected_value).abs() > 2e-9 {
                anyhow::bail!("schema 0.6 phase statistic mismatch at link {index}/{field}");
            }
        }
    }
    if !expected_phase_pairs.is_empty() {
        anyhow::bail!("schema 0.6 phase projection omits candidate event pairs");
    }
    validate_architecture_contract(data, &event_ids, &phase_eligible_molecules)?;
    Ok(())
}

pub(super) fn phase_string_array(link: &Value, field: &str, index: usize) -> Result<Vec<String>> {
    link.get(field)
        .and_then(Value::as_array)
        .with_context(|| format!("/phase_links/{index}/{field} must be an array"))?
        .iter()
        .map(|value| {
            value
                .as_str()
                .map(str::to_owned)
                .with_context(|| format!("/phase_links/{index}/{field} values must be strings"))
        })
        .collect()
}

pub(super) fn phase_u64_array(link: &Value, field: &str, index: usize) -> Result<Vec<u64>> {
    link.get(field)
        .and_then(Value::as_array)
        .with_context(|| format!("/phase_links/{index}/{field} must be an array"))?
        .iter()
        .map(|value| {
            value
                .as_u64()
                .with_context(|| format!("/phase_links/{index}/{field} values must be integers"))
        })
        .collect()
}

pub(super) fn phase_wilson_interval(successes: u64, total: u64) -> (f64, f64) {
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

pub(super) fn require_object<'a>(
    data: &'a Value,
    pointer: &str,
) -> Result<&'a serde_json::Map<String, Value>> {
    data.pointer(pointer)
        .and_then(Value::as_object)
        .with_context(|| format!("result contract violation: {pointer} must be an object"))
}

pub(super) fn require_array<'a>(data: &'a Value, pointer: &str) -> Result<&'a Vec<Value>> {
    data.pointer(pointer)
        .and_then(Value::as_array)
        .with_context(|| format!("result contract violation: {pointer} must be an array"))
}

pub(super) fn require_string<'a>(data: &'a Value, pointer: &str) -> Result<&'a str> {
    data.pointer(pointer)
        .and_then(Value::as_str)
        .with_context(|| format!("result contract violation: {pointer} must be a string"))
}

pub(super) fn require_u64(data: &Value, pointer: &str) -> Result<u64> {
    data.pointer(pointer)
        .and_then(Value::as_u64)
        .with_context(|| {
            format!("result contract violation: {pointer} must be an unsigned integer")
        })
}
