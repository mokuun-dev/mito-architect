use crate::args::{ClinicalSnapshotAction, ClinicalSnapshotArgs, UpdateClinicalArgs};
use crate::storage::{
    now_unix_secs, read_json_file, sha256_file, sync_directory, sync_file, write_json_atomic,
};
use anyhow::Context;
use anyhow::Result;
use serde_json::Value;
use std::collections::BTreeSet;
use std::fs;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command as ProcessCommand;
use std::process::Stdio;
use std::time::SystemTime;
use std::time::UNIX_EPOCH;

pub(super) const CLINICAL_TSV_HEADER: &str = "position\tref\talt\tgene\tconsequence\tprotein\tresidue\tstructure_id\tstructure_chain\tstructure_residue\tstructure_complex\tsource\tassertion_id\tallele_id\tdisease\tclinical_significance\treview_status\tassertion_date\tsource_url\treferences\tresource_version\tretrieved_at";

pub(super) fn update_clinical(args: UpdateClinicalArgs) -> Result<PathBuf> {
    if args.clinvar_live {
        return update_clinvar_cache(args.output, args.clinvar_gz);
    }

    let source = args.source.unwrap_or_else(default_bundled_clinical_tsv);
    if !source.exists() {
        anyhow::bail!(
            "clinical annotation source does not exist: {}",
            source.display()
        );
    }
    validate_clinical_tsv_schema(&source)?;

    let explicit_output = args.output.is_some();
    let mut output = args.output.unwrap_or_else(|| {
        default_cache_dir()
            .join("mito-architect")
            .join("clinical_annotations.tsv")
    });
    if let Some(parent) = output.parent() {
        if let Err(error) = fs::create_dir_all(parent) {
            if explicit_output {
                return Err(error)
                    .with_context(|| format!("failed to create {}", parent.display()));
            }
            output = std::env::temp_dir()
                .join("mito-architect")
                .join("clinical_annotations.tsv");
            if let Some(fallback_parent) = output.parent() {
                fs::create_dir_all(fallback_parent)
                    .with_context(|| format!("failed to create {}", fallback_parent.display()))?;
            }
        }
    }
    let output_name = output
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("clinical_annotations.tsv");
    let temporary_output =
        output.with_file_name(format!(".{output_name}.tmp-{}", std::process::id()));
    fs::copy(&source, &temporary_output).with_context(|| {
        format!(
            "failed to copy clinical annotations from {} to {}",
            source.display(),
            temporary_output.display()
        )
    })?;
    if let Err(error) = fs::rename(&temporary_output, &output) {
        let _ = fs::remove_file(&temporary_output);
        return Err(error).with_context(|| {
            format!(
                "failed to atomically replace clinical cache {}",
                output.display()
            )
        });
    }
    write_clinical_cache_metadata(&output, &source, false)?;
    Ok(output)
}

pub(super) fn validate_clinical_tsv_schema(path: &Path) -> Result<()> {
    let input = fs::File::open(path)
        .with_context(|| format!("failed to open clinical TSV {}", path.display()))?;
    let mut lines = BufReader::new(input).lines();
    let header = lines.next().transpose()?.context("clinical TSV is empty")?;
    if header.trim_end_matches('\r') != CLINICAL_TSV_HEADER {
        anyhow::bail!("clinical TSV header does not match assertion schema 1.0");
    }

    let mut row_count = 0usize;
    for (index, line) in lines.enumerate() {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 22 {
            anyhow::bail!(
                "clinical TSV row {} has {} fields; expected 22",
                index + 2,
                fields.len()
            );
        }
        let position = fields[0]
            .parse::<u64>()
            .with_context(|| format!("invalid clinical position at row {}", index + 2))?;
        let valid_allele = |value: &str| {
            value.len() == 1
                && matches!(
                    value.as_bytes()[0].to_ascii_uppercase(),
                    b'A' | b'C' | b'G' | b'T'
                )
        };
        if !(1..=16_569).contains(&position)
            || !valid_allele(fields[1])
            || !valid_allele(fields[2])
            || fields[1].eq_ignore_ascii_case(fields[2])
            || fields[11].trim().is_empty()
        {
            anyhow::bail!(
                "invalid clinical variant/assertion key at row {}",
                index + 2
            );
        }
        if !fields[18].is_empty()
            && !fields[18].starts_with("https://")
            && !fields[18].starts_with("http://")
        {
            anyhow::bail!("clinical source URL must be HTTP(S) at row {}", index + 2);
        }
        row_count += 1;
    }
    if row_count == 0 {
        anyhow::bail!("clinical TSV contains no assertions");
    }
    Ok(())
}

pub(super) fn manage_clinical_snapshot(args: ClinicalSnapshotArgs) -> Result<()> {
    let store = args.store.unwrap_or_else(default_clinical_snapshot_store);
    let result = match args.action {
        ClinicalSnapshotAction::Stage => stage_clinical_snapshot(
            &store,
            args.source.as_deref(),
            args.snapshot_id.as_deref(),
            args.license_id.as_deref(),
            args.source_policy.as_deref(),
            args.activate,
        )?,
        ClinicalSnapshotAction::Activate => {
            let snapshot_id = args
                .snapshot_id
                .as_deref()
                .context("clinical snapshot activate requires --snapshot-id")?;
            activate_clinical_snapshot(&store, snapshot_id, args.max_age_days)?
        }
        ClinicalSnapshotAction::Rollback => rollback_clinical_snapshot(&store, args.max_age_days)?,
        ClinicalSnapshotAction::Verify => {
            let snapshot_id = match args.snapshot_id.as_deref() {
                Some(snapshot_id) => snapshot_id.to_owned(),
                None => active_snapshot_id(&store)?.context(
                    "clinical snapshot verify requires --snapshot-id or an active snapshot",
                )?,
            };
            verify_clinical_snapshot(&store, &snapshot_id, args.max_age_days)?
        }
        ClinicalSnapshotAction::Status => clinical_snapshot_status(&store, args.max_age_days)?,
    };
    println!("{}", serde_json::to_string_pretty(&result)?);
    Ok(())
}

pub(super) fn default_clinical_snapshot_store() -> PathBuf {
    default_cache_dir()
        .join("mito-architect")
        .join("clinical-snapshots")
}

pub(super) fn validate_snapshot_id(snapshot_id: &str) -> Result<()> {
    if snapshot_id.is_empty()
        || snapshot_id.len() > 128
        || snapshot_id.starts_with('.')
        || snapshot_id.ends_with('.')
        || !snapshot_id
            .chars()
            .all(|value| value.is_ascii_alphanumeric() || matches!(value, '.' | '_' | '-'))
    {
        anyhow::bail!(
            "snapshot ID must be 1-128 ASCII letters, digits, '.', '_', or '-', without edge dots"
        );
    }
    Ok(())
}

pub(super) fn require_governance_label<'a>(name: &str, value: Option<&'a str>) -> Result<&'a str> {
    let value = value
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .with_context(|| format!("clinical snapshot stage requires --{name}"))?;
    if value.len() > 256
        || matches!(
            value.to_ascii_lowercase().as_str(),
            "unknown" | "unrecorded" | "not-recorded" | "none"
        )
        || value.chars().any(char::is_control)
    {
        anyhow::bail!("clinical snapshot --{name} is not a valid governed value");
    }
    Ok(value)
}

pub(super) fn snapshot_dir(store: &Path, snapshot_id: &str) -> PathBuf {
    store.join("snapshots").join(snapshot_id)
}

pub(super) fn snapshot_data_path(store: &Path, snapshot_id: &str) -> PathBuf {
    snapshot_dir(store, snapshot_id).join("clinical_annotations.tsv")
}

pub(super) fn snapshot_manifest_path(store: &Path, snapshot_id: &str) -> PathBuf {
    snapshot_dir(store, snapshot_id).join("manifest.json")
}

pub(super) fn clinical_tsv_summary(path: &Path) -> Result<Value> {
    validate_clinical_tsv_schema(path)?;
    let input = fs::File::open(path)
        .with_context(|| format!("failed to open clinical TSV {}", path.display()))?;
    let mut rows = 0u64;
    let mut variants = BTreeSet::new();
    let mut sources = BTreeSet::new();
    for (index, line) in BufReader::new(input).lines().enumerate() {
        let line = line?;
        if index == 0 || line.trim().is_empty() {
            continue;
        }
        let fields = line.trim_end_matches('\r').split('\t').collect::<Vec<_>>();
        rows = rows
            .checked_add(1)
            .context("clinical assertion count overflow")?;
        variants.insert(format!("{}:{}:{}", fields[0], fields[1], fields[2]));
        sources.insert(fields[11].to_owned());
    }
    Ok(serde_json::json!({
        "assertion_count": rows,
        "variant_count": variants.len(),
        "sources": sources,
    }))
}

pub(super) fn stage_clinical_snapshot(
    store: &Path,
    source: Option<&Path>,
    snapshot_id: Option<&str>,
    license_id: Option<&str>,
    source_policy: Option<&str>,
    activate: bool,
) -> Result<Value> {
    let source = source.context("clinical snapshot stage requires --source")?;
    if !source.is_file() {
        anyhow::bail!(
            "clinical snapshot source does not exist: {}",
            source.display()
        );
    }
    let snapshot_id = snapshot_id.context("clinical snapshot stage requires --snapshot-id")?;
    validate_snapshot_id(snapshot_id)?;
    let license_id = require_governance_label("license-id", license_id)?;
    let source_policy = require_governance_label("source-policy", source_policy)?;
    let summary = clinical_tsv_summary(source)?;

    let snapshots = store.join("snapshots");
    fs::create_dir_all(&snapshots)
        .with_context(|| format!("failed to create {}", snapshots.display()))?;
    let destination = snapshot_dir(store, snapshot_id);
    if destination.exists() {
        anyhow::bail!("clinical snapshot already exists: {snapshot_id}");
    }
    let temporary = snapshots.join(format!(
        ".{snapshot_id}.tmp-{}-{}",
        std::process::id(),
        now_unix_secs()?
    ));
    if temporary.exists() {
        anyhow::bail!(
            "clinical snapshot staging path already exists: {}",
            temporary.display()
        );
    }
    fs::create_dir(&temporary)
        .with_context(|| format!("failed to create {}", temporary.display()))?;

    let stage_result = (|| -> Result<Value> {
        let staged_data = temporary.join("clinical_annotations.tsv");
        fs::copy(source, &staged_data).with_context(|| {
            format!(
                "failed to stage clinical snapshot from {} to {}",
                source.display(),
                staged_data.display()
            )
        })?;
        let staged_summary = clinical_tsv_summary(&staged_data)?;
        if staged_summary != summary {
            anyhow::bail!("clinical snapshot summary changed during staging");
        }
        let source_sha256 = sha256_file(source)?;
        let data_sha256 = sha256_file(&staged_data)?;
        if source_sha256 != data_sha256 {
            anyhow::bail!("clinical snapshot bytes changed during staging");
        }
        let staged_at_unix = now_unix_secs()?;
        let manifest = serde_json::json!({
            "schema_version": "1.0",
            "clinical_annotation_schema_version": "1.0",
            "snapshot_id": snapshot_id,
            "status": "staged",
            "source": source.display().to_string(),
            "source_sha256": source_sha256,
            "data_sha256": data_sha256,
            "staged_at_unix": staged_at_unix,
            "license_id": license_id,
            "source_policy": source_policy,
            "summary": summary,
        });
        fs::write(
            temporary.join("manifest.json"),
            serde_json::to_vec_pretty(&manifest)?,
        )
        .with_context(|| format!("failed to write manifest for snapshot {snapshot_id}"))?;
        sync_file(&staged_data)?;
        sync_file(&temporary.join("manifest.json"))?;
        sync_directory(&temporary)?;
        fs::rename(&temporary, &destination).with_context(|| {
            format!(
                "failed to atomically publish clinical snapshot {}",
                destination.display()
            )
        })?;
        sync_directory(&snapshots)?;
        Ok(manifest)
    })();
    if stage_result.is_err() {
        let _ = fs::remove_dir_all(&temporary);
    }
    let manifest = stage_result?;

    if activate {
        return activate_clinical_snapshot(store, snapshot_id, None);
    }
    Ok(serde_json::json!({
        "operation": "stage",
        "store": store,
        "snapshot": manifest,
        "data_path": snapshot_data_path(store, snapshot_id),
    }))
}

pub(super) fn read_snapshot_manifest(store: &Path, snapshot_id: &str) -> Result<Value> {
    validate_snapshot_id(snapshot_id)?;
    let path = snapshot_manifest_path(store, snapshot_id);
    let manifest = read_json_file(&path, "clinical snapshot manifest")?;
    if manifest.get("schema_version").and_then(Value::as_str) != Some("1.0")
        || manifest
            .get("clinical_annotation_schema_version")
            .and_then(Value::as_str)
            != Some("1.0")
        || manifest.get("snapshot_id").and_then(Value::as_str) != Some(snapshot_id)
        || manifest.get("status").and_then(Value::as_str) != Some("staged")
    {
        anyhow::bail!("clinical snapshot manifest contract mismatch: {snapshot_id}");
    }
    for field in ["source_sha256", "data_sha256"] {
        let digest = manifest
            .get(field)
            .and_then(Value::as_str)
            .with_context(|| format!("clinical snapshot manifest requires {field}"))?;
        if digest.len() != 64 || !digest.chars().all(|value| value.is_ascii_hexdigit()) {
            anyhow::bail!("clinical snapshot manifest has invalid {field}");
        }
    }
    for field in ["license_id", "source_policy"] {
        require_governance_label(
            &field.replace('_', "-"),
            manifest.get(field).and_then(Value::as_str),
        )?;
    }
    manifest
        .get("staged_at_unix")
        .and_then(Value::as_u64)
        .context("clinical snapshot manifest requires staged_at_unix")?;
    manifest
        .get("summary")
        .and_then(Value::as_object)
        .context("clinical snapshot manifest requires summary")?;
    Ok(manifest)
}

pub(super) fn verify_clinical_snapshot(
    store: &Path,
    snapshot_id: &str,
    max_age_days: Option<u64>,
) -> Result<Value> {
    let manifest = read_snapshot_manifest(store, snapshot_id)?;
    let data_path = snapshot_data_path(store, snapshot_id);
    let summary = clinical_tsv_summary(&data_path)?;
    if manifest.get("summary") != Some(&summary) {
        anyhow::bail!("clinical snapshot summary mismatch: {snapshot_id}");
    }
    let actual_sha256 = sha256_file(&data_path)?;
    let expected_sha256 = manifest
        .get("data_sha256")
        .and_then(Value::as_str)
        .context("clinical snapshot manifest requires data_sha256")?;
    if actual_sha256 != expected_sha256 {
        anyhow::bail!(
            "clinical snapshot checksum mismatch for {snapshot_id}: expected {expected_sha256}, got {actual_sha256}"
        );
    }
    let staged_at = manifest
        .get("staged_at_unix")
        .and_then(Value::as_u64)
        .context("clinical snapshot manifest requires staged_at_unix")?;
    let now = now_unix_secs()?;
    if staged_at > now.saturating_add(300) {
        anyhow::bail!("clinical snapshot staging time is in the future: {snapshot_id}");
    }
    let age_seconds = now.saturating_sub(staged_at);
    if let Some(max_age_days) = max_age_days {
        let max_age_seconds = max_age_days
            .checked_mul(24 * 60 * 60)
            .context("clinical snapshot maximum age overflow")?;
        if age_seconds > max_age_seconds {
            anyhow::bail!(
                "clinical snapshot {snapshot_id} is stale: age {age_seconds}s exceeds {max_age_seconds}s"
            );
        }
    }
    Ok(serde_json::json!({
        "operation": "verify",
        "verified": true,
        "snapshot_id": snapshot_id,
        "data_path": data_path,
        "data_sha256": actual_sha256,
        "age_seconds": age_seconds,
        "max_age_days": max_age_days,
        "manifest": manifest,
    }))
}

pub(super) fn snapshot_state_path(store: &Path) -> PathBuf {
    store.join("state.json")
}

pub(super) fn read_snapshot_state(store: &Path) -> Result<Option<Value>> {
    let path = snapshot_state_path(store);
    if !path.exists() {
        return Ok(None);
    }
    if !path.is_file() {
        anyhow::bail!("clinical snapshot state is not a file: {}", path.display());
    }
    let state = read_json_file(&path, "clinical snapshot state")?;
    if state.get("schema_version").and_then(Value::as_str) != Some("1.0") {
        anyhow::bail!("clinical snapshot state requires schema_version 1.0");
    }
    if let Some(active) = state.get("active_snapshot").and_then(Value::as_str) {
        validate_snapshot_id(active)?;
    }
    if let Some(previous) = state.get("previous_snapshot").and_then(Value::as_str) {
        validate_snapshot_id(previous)?;
    }
    Ok(Some(state))
}

pub(super) fn active_snapshot_id(store: &Path) -> Result<Option<String>> {
    Ok(read_snapshot_state(store)?.and_then(|state| {
        state
            .get("active_snapshot")
            .and_then(Value::as_str)
            .map(str::to_owned)
    }))
}

pub(super) fn activate_clinical_snapshot(
    store: &Path,
    snapshot_id: &str,
    max_age_days: Option<u64>,
) -> Result<Value> {
    let verified = verify_clinical_snapshot(store, snapshot_id, max_age_days)?;
    let current_state = read_snapshot_state(store)?;
    let current = current_state
        .as_ref()
        .and_then(|state| state.get("active_snapshot"))
        .and_then(Value::as_str);
    let previous = current
        .filter(|current| *current != snapshot_id)
        .map(str::to_owned)
        .or_else(|| {
            current_state.as_ref().and_then(|state| {
                state
                    .get("previous_snapshot")
                    .and_then(Value::as_str)
                    .map(str::to_owned)
            })
        });
    let state = serde_json::json!({
        "schema_version": "1.0",
        "active_snapshot": snapshot_id,
        "previous_snapshot": previous,
        "updated_at_unix": now_unix_secs()?,
    });
    write_json_atomic(&snapshot_state_path(store), &state)?;
    Ok(serde_json::json!({
        "operation": "activate",
        "state": state,
        "active_path": snapshot_data_path(store, snapshot_id),
        "verification": verified,
    }))
}

pub(super) fn rollback_clinical_snapshot(store: &Path, max_age_days: Option<u64>) -> Result<Value> {
    let state = read_snapshot_state(store)?.context("no clinical snapshot state to roll back")?;
    let active = state
        .get("active_snapshot")
        .and_then(Value::as_str)
        .context("clinical snapshot state has no active snapshot")?;
    let previous = state
        .get("previous_snapshot")
        .and_then(Value::as_str)
        .context("clinical snapshot state has no previous snapshot")?;
    let verified = verify_clinical_snapshot(store, previous, max_age_days)?;
    let next_state = serde_json::json!({
        "schema_version": "1.0",
        "active_snapshot": previous,
        "previous_snapshot": active,
        "updated_at_unix": now_unix_secs()?,
    });
    write_json_atomic(&snapshot_state_path(store), &next_state)?;
    Ok(serde_json::json!({
        "operation": "rollback",
        "state": next_state,
        "active_path": snapshot_data_path(store, previous),
        "verification": verified,
    }))
}

pub(super) fn clinical_snapshot_status(store: &Path, max_age_days: Option<u64>) -> Result<Value> {
    let Some(state) = read_snapshot_state(store)? else {
        return Ok(serde_json::json!({
            "operation": "status",
            "store": store,
            "active_snapshot": null,
            "active_path": null,
            "verified": false,
        }));
    };
    let active = state
        .get("active_snapshot")
        .and_then(Value::as_str)
        .context("clinical snapshot state has no active snapshot")?;
    let verification = verify_clinical_snapshot(store, active, max_age_days)?;
    Ok(serde_json::json!({
        "operation": "status",
        "store": store,
        "active_snapshot": active,
        "active_path": snapshot_data_path(store, active),
        "state": state,
        "verified": true,
        "verification": verification,
    }))
}

pub(super) fn update_clinvar_cache(
    output: Option<PathBuf>,
    clinvar_gz: Option<PathBuf>,
) -> Result<PathBuf> {
    let explicit_output = output.is_some();
    let mut output = output.unwrap_or_else(|| {
        default_cache_dir()
            .join("mito-architect")
            .join("clinical_annotations.tsv")
    });

    if let Err(error) = try_prepare_parent(&output) {
        if explicit_output {
            return Err(error);
        }
        output = std::env::temp_dir()
            .join("mito-architect")
            .join("clinical_annotations.tsv");
        try_prepare_parent(&output)?;
    }

    let (source, remove_source) = match clinvar_gz {
        Some(path) => (path, false),
        None => (download_clinvar_summary()?, true),
    };
    let output_name = output
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("clinical_annotations.tsv");
    let temporary_output =
        output.with_file_name(format!(".{output_name}.tmp-{}", std::process::id()));
    let retrieved_at_unix = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .context("system time is before the Unix epoch")?
        .as_secs();
    let write_result = write_clinvar_mtdna_tsv(&source, &temporary_output, retrieved_at_unix)
        .and_then(|()| validate_clinical_tsv_schema(&temporary_output))
        .and_then(|()| {
            fs::rename(&temporary_output, &output).with_context(|| {
                format!(
                    "failed to atomically replace {} with refreshed clinical cache",
                    output.display()
                )
            })
        });
    if write_result.is_err() {
        let _ = fs::remove_file(&temporary_output);
    }
    let metadata_result = if write_result.is_ok() {
        write_clinical_cache_metadata(&output, &source, remove_source)
    } else {
        Ok(())
    };
    if remove_source {
        let _ = fs::remove_file(&source);
    }
    write_result?;
    metadata_result?;
    Ok(output)
}

pub(super) fn write_clinical_cache_metadata(
    output: &Path,
    source: &Path,
    downloaded: bool,
) -> Result<()> {
    let retrieved_at_unix = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .context("system time is before the Unix epoch")?
        .as_secs();
    let metadata = serde_json::json!({
        "schema_version": 2,
        "clinical_annotation_schema_version": "1.0",
        "resource": "ClinVar variant_summary mitochondrial SNV assertions",
        "source": if downloaded {
            "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz".to_string()
        } else {
            source.display().to_string()
        },
        "retrieved_at_unix": retrieved_at_unix,
        "source_sha256": sha256_file(source)?,
        "cache_sha256": sha256_file(output)?,
        "normalization": "GRCh38 MT/M single-nucleotide records; PositionVCF and VCF alleles; one source record per assertion row; bundled curated assertions appended without field-wise conflict loss",
    });
    let file_name = output
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("clinical_annotations.tsv");
    let metadata_path = output.with_file_name(format!("{file_name}.metadata.json"));
    let temporary_path = metadata_path.with_extension(format!("json.tmp-{}", std::process::id()));
    fs::write(&temporary_path, serde_json::to_vec_pretty(&metadata)?)
        .with_context(|| format!("failed to write {}", temporary_path.display()))?;
    fs::rename(&temporary_path, &metadata_path).with_context(|| {
        format!(
            "failed to atomically replace clinical metadata {}",
            metadata_path.display()
        )
    })?;
    Ok(())
}

pub(super) fn try_prepare_parent(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    Ok(())
}

pub(super) fn download_clinvar_summary() -> Result<PathBuf> {
    let output = std::env::temp_dir().join(format!(
        "clinvar_variant_summary_{}.txt.gz",
        std::process::id()
    ));
    let status = ProcessCommand::new("curl")
        .args([
            "--fail",
            "--silent",
            "--show-error",
            "-L",
            "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
            "-o",
        ])
        .arg(&output)
        .status()
        .context("failed to start curl for ClinVar download")?;
    if !status.success() {
        anyhow::bail!("curl failed while downloading ClinVar variant_summary.txt.gz");
    }
    Ok(output)
}

pub(super) fn write_clinvar_mtdna_tsv(
    source_gz: &Path,
    output: &Path,
    retrieved_at_unix: u64,
) -> Result<()> {
    let mut gzip = ProcessCommand::new("gzip")
        .arg("-dc")
        .arg(source_gz)
        .stdout(Stdio::piped())
        .spawn()
        .with_context(|| format!("failed to decompress {}", source_gz.display()))?;
    let stdout = gzip.stdout.take().context("gzip did not provide stdout")?;
    let mut reader = BufReader::new(stdout);
    let mut header = String::new();
    reader
        .read_line(&mut header)
        .context("ClinVar summary is empty")?;
    let columns = header
        .trim_start_matches('#')
        .trim_end()
        .split('\t')
        .enumerate()
        .map(|(index, name)| (name.to_string(), index))
        .collect::<std::collections::HashMap<_, _>>();
    let mut out = fs::File::create(output)
        .with_context(|| format!("failed to create {}", output.display()))?;
    writeln!(out, "{CLINICAL_TSV_HEADER}")?;

    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        {
            let fields = line.trim_end().split('\t').collect::<Vec<_>>();
            let get = |name: &str| -> &str {
                columns
                    .get(name)
                    .and_then(|index| fields.get(*index))
                    .copied()
                    .unwrap_or_default()
            };
            if get("Assembly") == "GRCh38" && matches!(get("Chromosome"), "MT" | "M") {
                let position = get("PositionVCF");
                let reference = get("ReferenceAlleleVCF");
                let alternate = get("AlternateAlleleVCF");
                if !position.is_empty()
                    && reference.len() == 1
                    && alternate.len() == 1
                    && reference != "na"
                    && alternate != "na"
                {
                    let variation_id = clean_tsv(get("VariationID"));
                    let source_url = if !variation_id.is_empty()
                        && variation_id.chars().all(|value| value.is_ascii_digit())
                    {
                        format!("https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/")
                    } else {
                        String::new()
                    };
                    let row = [
                        clean_tsv(position),
                        clean_tsv(reference),
                        clean_tsv(alternate),
                        clean_tsv(get("GeneSymbol")),
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                        String::new(),
                        "ClinVar".to_string(),
                        clean_tsv(&get("RCVaccession").replace('|', ";")),
                        clean_tsv(get("AlleleID")),
                        clean_tsv(get("PhenotypeList")),
                        clean_tsv(get("ClinicalSignificance")),
                        clean_tsv(get("ReviewStatus")),
                        clean_tsv(get("LastEvaluated")),
                        source_url,
                        clean_tsv(&get("RCVaccession").replace('|', ";")),
                        "clinvar-variant-summary-grch38".to_string(),
                        format!("unix:{retrieved_at_unix}"),
                    ];
                    writeln!(out, "{}", row.join("\t"))?;
                }
            }
        }
        line.clear();
    }

    let curated = fs::read_to_string(default_bundled_clinical_tsv())
        .context("failed to read bundled curated clinical annotations")?;
    for row in curated.lines().skip(1) {
        if !row.trim().is_empty() {
            writeln!(out, "{row}")?;
        }
    }

    let status = gzip.wait().context("failed to wait for gzip")?;
    if !status.success() {
        anyhow::bail!("gzip failed while reading {}", source_gz.display());
    }
    Ok(())
}

pub(super) fn clean_tsv(value: &str) -> String {
    value.replace(['\t', '\n', '\r'], " ").trim().to_string()
}

pub(super) fn default_bundled_clinical_tsv() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../core/data/clinical_annotations.tsv")
        .components()
        .collect::<PathBuf>()
}

pub(super) fn default_cache_dir() -> PathBuf {
    if let Some(path) = std::env::var_os("XDG_CACHE_HOME") {
        return PathBuf::from(path);
    }
    std::env::var_os("HOME")
        .map(PathBuf::from)
        .map(|home| home.join(".cache"))
        .unwrap_or_else(|| Path::new(".").join(".cache"))
}
