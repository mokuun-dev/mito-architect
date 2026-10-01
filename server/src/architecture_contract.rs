use serde_json::Value;
use std::collections::HashMap;
use std::collections::HashSet;

pub(super) fn validate_architecture_projection(value: &Value) -> Result<(), String> {
    let events = value
        .pointer("/events")
        .and_then(Value::as_array)
        .ok_or_else(|| "/events must be an array".to_string())?;
    let event_ids = events
        .iter()
        .enumerate()
        .map(|(index, event)| {
            event
                .get("id")
                .and_then(Value::as_str)
                .map(str::to_owned)
                .ok_or_else(|| format!("/events/{index}/id must be a string"))
        })
        .collect::<Result<HashSet<_>, _>>()?;
    let molecules = value
        .pointer("/molecules")
        .and_then(Value::as_array)
        .ok_or_else(|| "/molecules must be an array".to_string())?;
    let callability = value
        .pointer("/callability")
        .and_then(Value::as_array)
        .ok_or_else(|| "/callability must be an array".to_string())?;
    let detailed_callability = callability.iter().all(|entry| {
        entry.get("known").and_then(Value::as_bool).is_some()
            && entry.get("callable_bases").and_then(Value::as_u64).is_some()
    });
    let callable_molecule_ids = callability
        .iter()
        .enumerate()
        .filter_map(|(index, entry)| {
            let id = entry.get("molecule_id")?.as_str()?;
            let known = entry.get("known")?.as_bool()?;
            let callable_bases = entry.get("callable_bases")?.as_u64()?;
            (known && callable_bases != 0).then_some((index, id.to_owned()))
        })
        .map(|(_, id)| id)
        .collect::<HashSet<_>>();
    let mut molecule_ids = Vec::with_capacity(molecules.len());
    let mut eligible_indices = HashSet::new();
    for (index, molecule) in molecules.iter().enumerate() {
        molecule_ids.push(
            molecule
                .get("id")
                .and_then(Value::as_str)
                .ok_or_else(|| format!("/molecules/{index}/id must be a string"))?
                .to_owned(),
        );
        if molecule
            .get("evidence_eligible")
            .and_then(Value::as_bool)
            .ok_or_else(|| format!("/molecules/{index}/evidence_eligible must be boolean"))?
        {
            eligible_indices.insert(index as u64);
        }
    }
    if detailed_callability {
        eligible_indices.retain(|index| callable_molecule_ids.contains(&molecule_ids[*index as usize]));
    }

    let architectures = value
        .pointer("/architectures")
        .and_then(Value::as_array)
        .ok_or_else(|| "/architectures must be an array".to_string())?;
    let mut architecture_members = HashMap::<String, HashSet<u64>>::new();
    let mut architecture_ambiguous = HashMap::<String, HashSet<String>>::new();
    let mut assigned_indices = HashSet::new();
    for (index, architecture) in architectures.iter().enumerate() {
        let id = architecture
            .get("architecture_id")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("/architectures/{index}/architecture_id must be a string"))?;
        let members = architecture
            .get("member_molecule_indices")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("architecture {id} member_molecule_indices must be an array"))?
            .iter()
            .map(|value| {
                value
                    .as_u64()
                    .ok_or_else(|| format!("architecture {id} member index must be an integer"))
            })
            .collect::<Result<HashSet<_>, _>>()?;
        if architecture.get("molecule_count").and_then(Value::as_u64) != Some(members.len() as u64)
        {
            return Err(format!("architecture {id} molecule_count is inconsistent"));
        }
        for member in &members {
            if !eligible_indices.contains(member) || !assigned_indices.insert(*member) {
                return Err(format!(
                    "architecture {id} member {member} is ineligible or assigned twice"
                ));
            }
        }
        let expected_fraction = if eligible_indices.is_empty() {
            0.0
        } else {
            members.len() as f64 / eligible_indices.len() as f64
        };
        let actual_fraction = architecture
            .get("estimated_fraction")
            .and_then(Value::as_f64)
            .ok_or_else(|| format!("architecture {id} estimated_fraction must be a number"))?;
        if (actual_fraction - expected_fraction).abs() > 2e-9 {
            return Err(format!("architecture {id} abundance is inconsistent"));
        }
        for field in [
            "defining_event_signature",
            "seed_event_signature",
            "optional_event_ids",
        ] {
            let mut seen = HashSet::new();
            for event_id in architecture
                .get(field)
                .and_then(Value::as_array)
                .ok_or_else(|| format!("architecture {id}/{field} must be an array"))?
            {
                let event_id = event_id
                    .as_str()
                    .ok_or_else(|| format!("architecture {id}/{field} must contain strings"))?;
                if !event_ids.contains(event_id) || !seen.insert(event_id) {
                    return Err(format!(
                        "architecture {id}/{field} is not a unique event projection"
                    ));
                }
            }
        }
        let mut profile_event_ids = HashSet::new();
        let event_profile = architecture
            .get("event_profile")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("architecture {id} event_profile must be an array"))?;
        for (profile_index, profile) in event_profile.iter().enumerate() {
            let event_id = profile
                .get("event_id")
                .and_then(Value::as_str)
                .ok_or_else(|| {
                    format!(
                        "architecture {id} event_profile/{profile_index}/event_id must be a string"
                    )
                })?;
            if !event_ids.contains(event_id) || !profile_event_ids.insert(event_id) {
                return Err(format!(
                    "architecture {id} event_profile contains an unknown or duplicate event"
                ));
            }
            let count = |field: &str| {
                profile.get(field).and_then(Value::as_u64).ok_or_else(|| {
                    format!("architecture {id} event_profile/{profile_index}/{field} must be an integer")
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
                .ok_or_else(|| format!("architecture {id} event_profile count overflow"))?;
            if total != members.len() as u64 {
                return Err(format!(
                    "architecture {id} event_profile/{profile_index} does not partition members"
                ));
            }
            let callable = alternate + absent;
            if callable == 0 {
                return Err(format!(
                    "architecture {id} event_profile/{profile_index} has no callable support"
                ));
            }
            let fraction = profile
                .get("alternate_fraction")
                .and_then(Value::as_f64)
                .ok_or_else(|| {
                    format!("architecture {id} event_profile/{profile_index}/alternate_fraction must be a number")
                })?;
            let expected_fraction = alternate as f64 / callable as f64;
            if !fraction.is_finite() || (fraction - expected_fraction).abs() > 2e-9 {
                return Err(format!(
                    "architecture {id} event_profile/{profile_index} alternate fraction is inconsistent"
                ));
            }
            match profile.get("consensus_state").and_then(Value::as_str) {
                Some("ALTERNATE" | "ABSENT" | "OPTIONAL") => {}
                _ => {
                    return Err(format!(
                        "architecture {id} event_profile/{profile_index} consensus state is invalid"
                    ));
                }
            }
        }
        if architecture
            .pointer("/method/not_callable_is_reference")
            .and_then(Value::as_bool)
            != Some(false)
        {
            return Err(format!(
                "architecture {id} collapses NOT_CALLABLE into reference"
            ));
        }
        let ambiguous = architecture
            .get("ambiguous_molecule_ids")
            .and_then(Value::as_array)
            .ok_or_else(|| format!("architecture {id} ambiguous_molecule_ids must be an array"))?
            .iter()
            .map(|value| {
                value
                    .as_str()
                    .map(str::to_owned)
                    .ok_or_else(|| format!("architecture {id} ambiguous molecule must be a string"))
            })
            .collect::<Result<HashSet<_>, _>>()?;
        if ambiguous
            .iter()
            .any(|molecule_id| !molecule_ids.contains(molecule_id))
        {
            return Err(format!(
                "architecture {id} references an unknown ambiguous molecule"
            ));
        }
        if architecture_members
            .insert(id.to_owned(), members)
            .is_some()
            || architecture_ambiguous
                .insert(id.to_owned(), ambiguous)
                .is_some()
        {
            return Err(format!("duplicate architecture ID {id}"));
        }
    }

    let mut projected_counts = [0_u64; 4];
    for (index, molecule) in molecules.iter().enumerate() {
        let assignment = molecule
            .get("architecture_assignment")
            .and_then(Value::as_object)
            .ok_or_else(|| format!("molecule {index} architecture_assignment must be an object"))?;
        let status = assignment
            .get("status")
            .and_then(Value::as_str)
            .ok_or_else(|| format!("molecule {index} architecture status must be a string"))?;
        match status {
            "ASSIGNED" => {
                projected_counts[0] += 1;
                let architecture_id = assignment
                    .get("architecture_id")
                    .and_then(Value::as_str)
                    .ok_or_else(|| format!("assigned molecule {index} needs an architecture ID"))?;
                if !architecture_members
                    .get(architecture_id)
                    .is_some_and(|members| members.contains(&(index as u64)))
                {
                    return Err(format!(
                        "assigned molecule {index} is absent from {architecture_id}"
                    ));
                }
            }
            "AMBIGUOUS" => {
                projected_counts[1] += 1;
                let candidates = assignment
                    .get("candidate_architecture_ids")
                    .and_then(Value::as_array)
                    .ok_or_else(|| {
                        format!("ambiguous molecule {index} candidates must be an array")
                    })?;
                if candidates.len() < 2 {
                    return Err(format!(
                        "ambiguous molecule {index} needs at least two candidates"
                    ));
                }
                for candidate in candidates {
                    let candidate = candidate.as_str().ok_or_else(|| {
                        format!("ambiguous molecule {index} candidate must be a string")
                    })?;
                    if !architecture_ambiguous
                        .get(candidate)
                        .is_some_and(|ids| ids.contains(&molecule_ids[index]))
                    {
                        return Err(format!(
                            "ambiguous molecule {index} is absent from candidate {candidate}"
                        ));
                    }
                }
            }
            "UNASSIGNED" => projected_counts[2] += 1,
            "INELIGIBLE" => projected_counts[3] += 1,
            _ => return Err(format!("unsupported molecule architecture status {status}")),
        }
        if eligible_indices.contains(&(index as u64)) == (status == "INELIGIBLE") {
            return Err(format!(
                "molecule {index} architecture eligibility is inconsistent"
            ));
        }
    }

    let inference = value
        .pointer("/architecture_inference")
        .and_then(Value::as_object)
        .ok_or_else(|| "/architecture_inference must be an object".to_string())?;
    for (field, expected) in [
        ("eligible_molecules", eligible_indices.len() as u64),
        ("assigned_molecules", projected_counts[0]),
        ("ambiguous_molecules", projected_counts[1]),
        ("unassigned_molecules", projected_counts[2]),
        ("candidate_count", architectures.len() as u64),
    ] {
        if inference.get(field).and_then(Value::as_u64) != Some(expected) {
            return Err(format!(
                "architecture inference count {field} is inconsistent"
            ));
        }
    }
    if projected_counts[0] + projected_counts[1] + projected_counts[2]
        != eligible_indices.len() as u64
    {
        return Err(
            "architecture assignment counts do not partition eligible molecules".to_string(),
        );
    }
    Ok(())
}
