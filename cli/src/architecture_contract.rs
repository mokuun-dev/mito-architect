use crate::contract::{require_array, require_object};
use anyhow::Context;
use anyhow::Result;
use serde_json::Value;
use std::collections::BTreeMap;
use std::collections::BTreeSet;

pub(super) fn validate_architecture_contract(
    data: &Value,
    event_ids: &BTreeSet<String>,
    eligible_molecules: &BTreeMap<usize, String>,
) -> Result<()> {
    let architectures = require_array(data, "/architectures")?;
    let molecules = require_array(data, "/molecules")?;
    let callability = require_array(data, "/callability")?;
    let mut callable_molecule_ids = BTreeSet::new();
    let detailed_callability = callability.iter().all(|entry| {
        entry.get("known").and_then(Value::as_bool).is_some()
            && entry.get("callable_bases").and_then(Value::as_u64).is_some()
    });
    if detailed_callability {
        let positive_molecule_ids = molecules
            .iter()
            .enumerate()
            .map(|(index, molecule)| {
                let id = molecule
                    .get("id")
                    .and_then(Value::as_str)
                    .with_context(|| format!("/molecules/{index}/id must be a string"))?;
                let alternates = molecule
                    .get("alternate_event_ids")
                    .and_then(Value::as_array)
                    .with_context(|| format!("/molecules/{index}/alternate_event_ids must be an array"))?;
                Ok((id, !alternates.is_empty()))
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .filter_map(|(id, positive)| positive.then_some(id))
            .collect::<BTreeSet<_>>();
        for (index, entry) in callability.iter().enumerate() {
            let molecule_id = entry
                .get("molecule_id")
                .and_then(Value::as_str)
                .with_context(|| format!("/callability/{index}/molecule_id must be a string"))?;
            let known = entry
                .get("known")
                .and_then(Value::as_bool)
                .with_context(|| format!("/callability/{index}/known must be boolean"))?;
            let callable_bases = entry
                .get("callable_bases")
                .and_then(Value::as_u64)
                .with_context(|| format!("/callability/{index}/callable_bases must be an integer"))?;
            // Positive event evidence, including a deletion, can be informative
            // even when no aligned reference base passes the quality threshold.
            if known && (callable_bases != 0 || positive_molecule_ids.contains(molecule_id)) {
                callable_molecule_ids.insert(molecule_id);
            }
        }
    }
    let architecture_eligible_molecules = if detailed_callability {
        eligible_molecules
            .iter()
            .filter(|(_, id)| callable_molecule_ids.contains(id.as_str()))
            .map(|(index, id)| (*index, id.clone()))
            .collect::<BTreeMap<_, _>>()
    } else {
        eligible_molecules.clone()
    };
    let mut members_by_architecture = BTreeMap::<String, BTreeSet<u64>>::new();
    let mut globally_assigned = BTreeSet::new();
    for (index, architecture) in architectures.iter().enumerate() {
        let id = architecture
            .get("architecture_id")
            .and_then(Value::as_str)
            .with_context(|| format!("/architectures/{index}/architecture_id must be a string"))?;
        let members = architecture
            .get("member_molecule_indices")
            .and_then(Value::as_array)
            .with_context(|| format!("architecture {id} members must be an array"))?
            .iter()
            .map(|value| {
                value
                    .as_u64()
                    .context("architecture member index must be an integer")
            })
            .collect::<Result<BTreeSet<_>>>()?;
        if architecture.get("molecule_count").and_then(Value::as_u64) != Some(members.len() as u64)
        {
            anyhow::bail!("architecture {id} molecule_count mismatch");
        }
        for member in &members {
            let member_index =
                usize::try_from(*member).context("architecture member index overflow")?;
            if !architecture_eligible_molecules.contains_key(&member_index) || !globally_assigned.insert(*member)
            {
                anyhow::bail!("architecture {id} has an ineligible or duplicate member {member}");
            }
        }
        let expected_fraction = if architecture_eligible_molecules.is_empty() {
            0.0
        } else {
            members.len() as f64 / architecture_eligible_molecules.len() as f64
        };
        let fraction = architecture
            .get("estimated_fraction")
            .and_then(Value::as_f64)
            .with_context(|| format!("architecture {id} estimated_fraction must be a number"))?;
        if (fraction - expected_fraction).abs() > 2e-9 {
            anyhow::bail!("architecture {id} estimated_fraction mismatch");
        }
        for field in [
            "defining_event_signature",
            "seed_event_signature",
            "optional_event_ids",
        ] {
            let mut unique = BTreeSet::new();
            for event_id in architecture
                .get(field)
                .and_then(Value::as_array)
                .with_context(|| format!("architecture {id}/{field} must be an array"))?
            {
                let event_id = event_id
                    .as_str()
                    .with_context(|| format!("architecture {id}/{field} values must be strings"))?;
                if !event_ids.contains(event_id) || !unique.insert(event_id) {
                    anyhow::bail!("architecture {id}/{field} is not a unique event projection");
                }
            }
        }
        let mut profile_event_ids = BTreeSet::new();
        for (profile_index, profile) in architecture
            .get("event_profile")
            .and_then(Value::as_array)
            .with_context(|| format!("architecture {id} event_profile must be an array"))?
            .iter()
            .enumerate()
        {
            let event_id = profile
                .get("event_id")
                .and_then(Value::as_str)
                .with_context(|| {
                    format!(
                        "architecture {id} event_profile/{profile_index}/event_id must be a string"
                    )
                })?;
            if !event_ids.contains(event_id) || !profile_event_ids.insert(event_id) {
                anyhow::bail!(
                    "architecture {id} event_profile contains an unknown or duplicate event"
                );
            }
            let count = |field: &str| {
                profile
                    .get(field)
                    .and_then(Value::as_u64)
                    .with_context(|| {
                        format!(
                            "architecture {id} event_profile/{profile_index}/{field} must be an integer"
                        )
                    })
            };
            let alternate = count("alternate")?;
            let absent = count("absent")?;
            let uncertain = count("uncertain")?;
            let not_callable = count("not_callable")?;
            let total = alternate
                .checked_add(absent)
                .and_then(|value| value.checked_add(uncertain))
                .and_then(|value| value.checked_add(not_callable))
                .with_context(|| format!("architecture {id} event profile count overflow"))?;
            if total != members.len() as u64 {
                anyhow::bail!(
                    "architecture {id} event_profile/{profile_index} does not partition members"
                );
            }
            let callable = alternate + absent;
            if callable == 0 {
                anyhow::bail!(
                    "architecture {id} event_profile/{profile_index} has no callable support"
                );
            }
            let actual_fraction = profile
                .get("alternate_fraction")
                .and_then(Value::as_f64)
                .with_context(|| {
                    format!(
                        "architecture {id} event_profile/{profile_index}/alternate_fraction must be a number"
                    )
                })?;
            let expected_fraction = alternate as f64 / callable as f64;
            if !actual_fraction.is_finite() || (actual_fraction - expected_fraction).abs() > 2e-9 {
                anyhow::bail!(
                    "architecture {id} event_profile/{profile_index} alternate fraction mismatch"
                );
            }
            if !matches!(
                profile.get("consensus_state").and_then(Value::as_str),
                Some("ALTERNATE" | "ABSENT" | "OPTIONAL")
            ) {
                anyhow::bail!(
                    "architecture {id} event_profile/{profile_index} consensus state is invalid"
                );
            }
        }
        if architecture
            .pointer("/method/not_callable_is_reference")
            .and_then(Value::as_bool)
            != Some(false)
        {
            anyhow::bail!("architecture {id} must preserve NOT_CALLABLE semantics");
        }
        if members_by_architecture
            .insert(id.to_owned(), members)
            .is_some()
        {
            anyhow::bail!("duplicate architecture ID {id}");
        }
    }

    let mut counts = [0_u64; 4];
    for (index, molecule) in molecules.iter().enumerate() {
        let status = molecule
            .pointer("/architecture_assignment/status")
            .and_then(Value::as_str)
            .with_context(|| format!("molecule {index} architecture status must be a string"))?;
        match status {
            "ASSIGNED" => {
                counts[0] += 1;
                let id = molecule
                    .pointer("/architecture_assignment/architecture_id")
                    .and_then(Value::as_str)
                    .with_context(|| {
                        format!("assigned molecule {index} needs an architecture ID")
                    })?;
                if !members_by_architecture
                    .get(id)
                    .is_some_and(|members| members.contains(&(index as u64)))
                {
                    anyhow::bail!("assigned molecule {index} is absent from architecture {id}");
                }
            }
            "AMBIGUOUS" => {
                counts[1] += 1;
                let candidates = molecule
                    .pointer("/architecture_assignment/candidate_architecture_ids")
                    .and_then(Value::as_array)
                    .with_context(|| {
                        format!("ambiguous molecule {index} candidates must be an array")
                    })?;
                if candidates.len() < 2
                    || candidates.iter().any(|candidate| match candidate.as_str() {
                        Some(id) => !members_by_architecture.contains_key(id),
                        None => true,
                    })
                {
                    anyhow::bail!("ambiguous molecule {index} has invalid candidates");
                }
            }
            "UNASSIGNED" => counts[2] += 1,
            "INELIGIBLE" => counts[3] += 1,
            _ => anyhow::bail!("unsupported architecture assignment status {status}"),
        }
        let eligible = architecture_eligible_molecules.contains_key(&index);
        if eligible == (status == "INELIGIBLE") {
            anyhow::bail!("molecule {index} architecture eligibility mismatch");
        }
    }
    let inference = require_object(data, "/architecture_inference")?;
    for (field, expected) in [
        ("eligible_molecules", architecture_eligible_molecules.len() as u64),
        ("assigned_molecules", counts[0]),
        ("ambiguous_molecules", counts[1]),
        ("unassigned_molecules", counts[2]),
        ("candidate_count", architectures.len() as u64),
    ] {
        if inference.get(field).and_then(Value::as_u64) != Some(expected) {
            anyhow::bail!("architecture inference {field} mismatch");
        }
    }
    if counts[0] + counts[1] + counts[2] != architecture_eligible_molecules.len() as u64 {
        anyhow::bail!("architecture assignments do not partition eligible molecules");
    }
    Ok(())
}
