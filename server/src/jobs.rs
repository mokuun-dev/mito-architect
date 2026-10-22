use crate::contract::validate_result_contract;
use crate::error::FailureBody;
use crate::evidence::{prepare_completed_artifacts, CompletedArtifacts};
use crate::report::render_html_report;
use crate::state::{AppState, JobAnalysisOptions, JobRecord, JobStatus};
use mito_ffi::AnalyzeOptions;
use mito_ffi::MitoEngine;
use serde_json::Value;
use std::collections::HashMap;
use std::path::Path;
use std::path::PathBuf;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering;
use std::sync::Arc;
use std::time::Duration;
use std::time::SystemTime;
use std::time::UNIX_EPOCH;
use tokio::fs;
use tokio::sync::OwnedSemaphorePermit;
use uuid::Uuid;

pub(super) const JOB_TTL_SECS: u64 = 24 * 60 * 60;

pub(super) async fn run_job(
    state: AppState,
    job_id: Uuid,
    input_path: PathBuf,
    cancel_flag: Arc<AtomicBool>,
    analysis_options: JobAnalysisOptions,
    _job_permit: OwnedSemaphorePermit,
) {
    let Ok(_permit) = state.analysis_slots.clone().acquire_owned().await else {
        remove_job_input(&input_path).await;
        update_job(
            &state,
            job_id,
            JobStatus::Error,
            100,
            Some(FailureBody::new(
                "MITO-API-E2001",
                "analysis scheduler is unavailable",
                true,
            )),
        )
        .await;
        return;
    };
    if is_cancelled(&state, job_id).await {
        remove_job_input(&input_path).await;
        state.cancel_flags.write().await.remove(&job_id);
        return;
    }
    update_job(&state, job_id, JobStatus::Processing, 12, None).await;

    let options = AnalyzeOptions {
        filter_numt: true,
        threads: state.worker_threads,
        min_mapping_quality: state.min_mapping_quality,
        min_base_quality: state.min_base_quality,
        excluded_snp_flags: state.excluded_snp_flags,
        numt_threshold: state.numt_threshold,
        allow_development_tags: false,
        emit_evidence_graph: analysis_options.emit_evidence_graph,
        max_evidence_observations: state.max_evidence_observations,
        max_phase_links: state.max_phase_links,
        max_phase_work: state.max_phase_work,
        max_phase_molecule_references: state.max_phase_molecule_references,
        max_result_bytes: state.max_result_bytes,
        max_input_records: state.max_input_records,
        max_input_bases: state.max_input_bases,
        allow_external_cram_references: false,
        evidence_page_size: analysis_options.evidence_page_size,
        min_architecture_molecules: state.min_architecture_molecules,
        max_candidate_architectures: state.max_candidate_architectures,
        architecture_max_distance: state.architecture_max_distance,
        architecture_ambiguity_margin: state.architecture_ambiguity_margin,
        architecture_min_overlap_fraction: state.architecture_min_overlap_fraction,
        architecture_consensus_fraction: state.architecture_consensus_fraction,
        architecture_optional_fraction: state.architecture_optional_fraction,
        architecture_stability_replicates: state.architecture_stability_replicates,
        architecture_seed: state.architecture_seed,
        molecule_id_tag: analysis_options.molecule_id_tag,
        umi_tag: analysis_options.umi_tag,
        duplex_tag: analysis_options.duplex_tag,
    };
    let input_for_worker = input_path.clone();
    let worker_cancel_flag = cancel_flag.clone();
    // Analysis and every CPU-heavy post-processing step (JSON parsing,
    // contract validation, report rendering, paging) run on the blocking pool
    // so a large result never stalls the async runtime that serves requests.
    let outcome = tokio::task::spawn_blocking(move || {
        complete_analysis(job_id, &input_for_worker, options, &worker_cancel_flag)
    })
    .await;
    // The upload is no longer needed once the engine has returned; do not keep
    // genomic input on disk until the TTL sweep.
    remove_job_input(&input_path).await;

    if is_cancelled(&state, job_id).await {
        state.cancel_flags.write().await.remove(&job_id);
        return;
    }
    match outcome {
        Ok(Ok(artifacts)) => {
            store_completed_job(&state, job_id, artifacts).await;
        }
        Ok(Err(failure)) => {
            update_job(&state, job_id, JobStatus::Error, 100, Some(failure)).await;
        }
        Err(error) => {
            eprintln!("job {job_id} analysis task failed: {error}");
            update_job(
                &state,
                job_id,
                JobStatus::Error,
                100,
                Some(internal_analysis_failure()),
            )
            .await;
        }
    }
    state.cancel_flags.write().await.remove(&job_id);
}

fn internal_analysis_failure() -> FailureBody {
    FailureBody::new("MITO-E9001", "internal analysis error", false)
}

/// Runs the engine and prepares transport artifacts on a blocking thread.
fn complete_analysis(
    job_id: Uuid,
    input_path: &Path,
    options: AnalyzeOptions,
    cancel_flag: &AtomicBool,
) -> Result<CompletedArtifacts, FailureBody> {
    let cancelled = || FailureBody::new("MITO-E1501", "cancelled by user", false);
    let engine =
        MitoEngine::new().map_err(|error| public_engine_failure(error, job_id, input_path))?;
    let json_text = engine
        .analyze_with_cancel_flag(input_path, None, options, cancel_flag)
        .map_err(|error| public_engine_failure(error, job_id, input_path))?;
    drop(engine);
    let value = serde_json::from_str::<Value>(&json_text).map_err(|error| {
        eprintln!("job {job_id} returned invalid JSON: {error}");
        internal_analysis_failure()
    })?;
    if cancel_flag.load(Ordering::Relaxed) {
        return Err(cancelled());
    }
    validate_result_contract(&value).map_err(|error| {
        eprintln!("job {job_id} result contract violation: {error}");
        internal_analysis_failure()
    })?;
    if cancel_flag.load(Ordering::Relaxed) {
        return Err(cancelled());
    }
    let html = render_html_report(&value);
    prepare_completed_artifacts(job_id, value, json_text, html).map_err(|error| {
        eprintln!("job {job_id} transport preparation failed: {error}");
        internal_analysis_failure()
    })
}

async fn remove_job_input(input_path: &Path) {
    if let Some(job_dir) = input_path.parent() {
        let _ = fs::remove_dir_all(job_dir).await;
    }
}

/// Stores a finished result, evicting the oldest other finished jobs while the
/// retained bytes exceed max_stored_bytes. A result that alone exceeds the
/// budget is not retained and the job fails with a resource error.
async fn store_completed_job(state: &AppState, job_id: Uuid, artifacts: CompletedArtifacts) {
    let size = artifacts.retained_bytes();
    if size > state.max_stored_bytes {
        update_job(
            state,
            job_id,
            JobStatus::Error,
            100,
            Some(FailureBody::new(
                "MITO-E1601",
                "analysis result exceeds the server retention limit (MITO_MAX_STORED_BYTES)",
                false,
            )),
        )
        .await;
        return;
    }
    let evicted = {
        let mut jobs = state.jobs.write().await;
        let Some(job) = jobs.get_mut(&job_id) else {
            return;
        };
        if job.status == JobStatus::Cancelled {
            return;
        }
        job.status = JobStatus::Done;
        job.progress = 100;
        job.result = Some(artifacts.result);
        job.result_summary = Some(artifacts.result_summary);
        job.evidence_pages = artifacts.evidence_pages;
        job.evidence_search_index = artifacts.evidence_search_index;
        job.html_report = Some(artifacts.html_report);
        job.stored_bytes = size;
        evict_finished_jobs(&mut jobs, Some(job_id), |jobs| {
            jobs.values().map(|job| job.stored_bytes).sum::<usize>() <= state.max_stored_bytes
        })
    };
    discard_evicted_jobs(state, evicted).await;
}

/// Removes finished jobs, oldest first and never `keep`, until `fits` holds.
/// Returns the removed job IDs and input paths for cleanup outside the lock.
pub(super) fn evict_finished_jobs(
    jobs: &mut HashMap<Uuid, JobRecord>,
    keep: Option<Uuid>,
    fits: impl Fn(&HashMap<Uuid, JobRecord>) -> bool,
) -> Vec<(Uuid, PathBuf)> {
    let mut evicted = Vec::new();
    while !fits(jobs) {
        let oldest = jobs
            .iter()
            .filter(|(id, job)| Some(**id) != keep && job.status.is_finished())
            .min_by_key(|(id, job)| (job.created_at, **id))
            .map(|(id, _)| *id);
        let Some(oldest) = oldest else {
            break;
        };
        if let Some(job) = jobs.remove(&oldest) {
            evicted.push((oldest, job.input_path));
        }
    }
    evicted
}

pub(super) async fn discard_evicted_jobs(state: &AppState, evicted: Vec<(Uuid, PathBuf)>) {
    if evicted.is_empty() {
        return;
    }
    {
        let mut cancel_flags = state.cancel_flags.write().await;
        for (job_id, _) in &evicted {
            cancel_flags.remove(job_id);
        }
    }
    for (_, input_path) in evicted {
        remove_job_input(&input_path).await;
    }
}

pub(super) fn public_engine_failure(
    error: mito_ffi::MitoError,
    job_id: Uuid,
    input_path: &Path,
) -> FailureBody {
    let message = match error.code.as_str() {
        "MITO-E1101" => "input could not be opened".to_string(),
        "MITO-E1201" => "reference could not be opened".to_string(),
        "MITO-E1301" => "required analysis resource could not be opened".to_string(),
        "MITO-E9001" => {
            eprintln!("job {job_id} internal engine error: {}", error.message);
            "internal analysis error".to_string()
        }
        // Diagnostics may name the stored upload; show only its file name.
        _ => {
            let stored = input_path.to_string_lossy();
            let file_name = input_path
                .file_name()
                .map(|name| name.to_string_lossy().into_owned())
                .unwrap_or_default();
            if stored.is_empty() {
                error.message
            } else {
                error.message.replace(stored.as_ref(), &file_name)
            }
        }
    };
    FailureBody::new(error.code, message, false)
}

pub(super) async fn is_cancelled(state: &AppState, job_id: Uuid) -> bool {
    state
        .jobs
        .read()
        .await
        .get(&job_id)
        .is_some_and(|job| job.cancel_requested || job.status == JobStatus::Cancelled)
}

pub(super) async fn update_job(
    state: &AppState,
    job_id: Uuid,
    status: JobStatus,
    progress: u8,
    error: Option<FailureBody>,
) {
    if let Some(job) = state.jobs.write().await.get_mut(&job_id) {
        if job.status == JobStatus::Cancelled && status != JobStatus::Cancelled {
            return;
        }
        job.status = status;
        job.progress = progress.min(100);
        if error.is_some() {
            job.error = error;
        }
    }
}

pub(super) fn spawn_cleanup(state: AppState) {
    tokio::spawn(async move {
        loop {
            tokio::time::sleep(Duration::from_secs(60 * 60)).await;
            let expired = {
                let jobs = state.jobs.read().await;
                jobs.iter()
                    .filter_map(|(job_id, job)| {
                        (now_epoch_secs().saturating_sub(job.created_at) > JOB_TTL_SECS)
                            .then_some((*job_id, job.input_path.clone()))
                    })
                    .collect::<Vec<_>>()
            };

            if expired.is_empty() {
                continue;
            }

            {
                let cancel_flags = state.cancel_flags.read().await;
                for (job_id, _) in &expired {
                    if let Some(flag) = cancel_flags.get(job_id) {
                        flag.store(true, Ordering::Relaxed);
                    }
                }
            }
            {
                let mut jobs = state.jobs.write().await;
                for (job_id, _) in &expired {
                    jobs.remove(job_id);
                }
            }
            {
                let mut cancel_flags = state.cancel_flags.write().await;
                for (job_id, _) in &expired {
                    cancel_flags.remove(job_id);
                }
            }
            for (_, input_path) in expired {
                if let Some(parent) = input_path.parent() {
                    let _ = fs::remove_dir_all(parent).await;
                }
            }
        }
    });
}

pub(super) fn now_epoch_secs() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |duration| duration.as_secs())
}
