use crate::contract::validate_result_contract;
use crate::error::FailureBody;
use crate::evidence::{prepare_completed_artifacts, CompletedArtifacts};
use crate::report::render_html_report;
use crate::state::{AppState, JobAnalysisOptions, JobStatus};
use mito_ffi::AnalyzeOptions;
use mito_ffi::MitoEngine;
use serde_json::Value;
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
        update_job(
            &state,
            job_id,
            JobStatus::Error,
            100,
            None,
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
        state.cancel_flags.write().await.remove(&job_id);
        return;
    }
    update_job(&state, job_id, JobStatus::Processing, 12, None, None).await;

    let input_for_worker = input_path.clone();
    let threads = state.worker_threads;
    let min_mapping_quality = state.min_mapping_quality;
    let min_base_quality = state.min_base_quality;
    let excluded_snp_flags = state.excluded_snp_flags;
    let numt_threshold = state.numt_threshold;
    let max_evidence_observations = state.max_evidence_observations;
    let max_phase_links = state.max_phase_links;
    let max_phase_work = state.max_phase_work;
    let max_phase_molecule_references = state.max_phase_molecule_references;
    let max_result_bytes = state.max_result_bytes;
    let min_architecture_molecules = state.min_architecture_molecules;
    let max_candidate_architectures = state.max_candidate_architectures;
    let architecture_max_distance = state.architecture_max_distance;
    let architecture_ambiguity_margin = state.architecture_ambiguity_margin;
    let architecture_min_overlap_fraction = state.architecture_min_overlap_fraction;
    let architecture_consensus_fraction = state.architecture_consensus_fraction;
    let architecture_optional_fraction = state.architecture_optional_fraction;
    let architecture_stability_replicates = state.architecture_stability_replicates;
    let architecture_seed = state.architecture_seed;
    let analysis = tokio::task::spawn_blocking(move || {
        let engine = MitoEngine::new()?;
        engine.analyze_with_cancel_flag(
            &input_for_worker,
            None,
            AnalyzeOptions {
                filter_numt: true,
                threads,
                min_mapping_quality,
                min_base_quality,
                excluded_snp_flags,
                numt_threshold,
                allow_development_tags: false,
                emit_evidence_graph: analysis_options.emit_evidence_graph,
                max_evidence_observations,
                max_phase_links,
                max_phase_work,
                max_phase_molecule_references,
                max_result_bytes,
                evidence_page_size: analysis_options.evidence_page_size,
                min_architecture_molecules,
                max_candidate_architectures,
                architecture_max_distance,
                architecture_ambiguity_margin,
                architecture_min_overlap_fraction,
                architecture_consensus_fraction,
                architecture_optional_fraction,
                architecture_stability_replicates,
                architecture_seed,
                molecule_id_tag: analysis_options.molecule_id_tag,
                umi_tag: analysis_options.umi_tag,
                duplex_tag: analysis_options.duplex_tag,
            },
            &cancel_flag,
        )
    })
    .await;

    if is_cancelled(&state, job_id).await {
        state.cancel_flags.write().await.remove(&job_id);
        return;
    }

    match analysis {
        Ok(Ok(json_text)) => match serde_json::from_str::<Value>(&json_text) {
            Ok(value) => {
                if is_cancelled(&state, job_id).await {
                    state.cancel_flags.write().await.remove(&job_id);
                    return;
                }
                if let Err(error) = validate_result_contract(&value) {
                    eprintln!("job {job_id} result contract violation: {error}");
                    update_job(
                        &state,
                        job_id,
                        JobStatus::Error,
                        100,
                        None,
                        Some(FailureBody::new(
                            "MITO-E9001",
                            "internal analysis error",
                            false,
                        )),
                    )
                    .await;
                } else {
                    if is_cancelled(&state, job_id).await {
                        state.cancel_flags.write().await.remove(&job_id);
                        return;
                    }
                    let html = render_html_report(&value);
                    if is_cancelled(&state, job_id).await {
                        state.cancel_flags.write().await.remove(&job_id);
                        return;
                    }
                    match prepare_completed_artifacts(job_id, &value, json_text, html) {
                        Ok(artifacts) => {
                            if is_cancelled(&state, job_id).await {
                                state.cancel_flags.write().await.remove(&job_id);
                                return;
                            }
                            update_job(&state, job_id, JobStatus::Done, 100, Some(artifacts), None)
                                .await;
                        }
                        Err(error) => {
                            eprintln!("job {job_id} transport preparation failed: {error}");
                            update_job(
                                &state,
                                job_id,
                                JobStatus::Error,
                                100,
                                None,
                                Some(FailureBody::new(
                                    "MITO-E9001",
                                    "internal analysis error",
                                    false,
                                )),
                            )
                            .await;
                        }
                    }
                }
            }
            Err(error) => {
                eprintln!("job {job_id} returned invalid JSON: {error}");
                update_job(
                    &state,
                    job_id,
                    JobStatus::Error,
                    100,
                    None,
                    Some(FailureBody::new(
                        "MITO-E9001",
                        "internal analysis error",
                        false,
                    )),
                )
                .await;
            }
        },
        Ok(Err(error)) => {
            let failure = public_engine_failure(error, job_id);
            update_job(&state, job_id, JobStatus::Error, 100, None, Some(failure)).await;
        }
        Err(error) => {
            eprintln!("job {job_id} analysis task failed: {error}");
            update_job(
                &state,
                job_id,
                JobStatus::Error,
                100,
                None,
                Some(FailureBody::new(
                    "MITO-E9001",
                    "internal analysis error",
                    false,
                )),
            )
            .await;
        }
    }

    state.cancel_flags.write().await.remove(&job_id);
}

pub(super) fn public_engine_failure(error: mito_ffi::MitoError, job_id: Uuid) -> FailureBody {
    let message = match error.code.as_str() {
        "MITO-E1101" => "input could not be opened".to_string(),
        "MITO-E1201" => "reference could not be opened".to_string(),
        "MITO-E1301" => "required analysis resource could not be opened".to_string(),
        "MITO-E9001" => {
            eprintln!("job {job_id} internal engine error: {}", error.message);
            "internal analysis error".to_string()
        }
        _ => error.message,
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
    artifacts: Option<CompletedArtifacts>,
    error: Option<FailureBody>,
) {
    if let Some(job) = state.jobs.write().await.get_mut(&job_id) {
        if job.status == JobStatus::Cancelled && status != JobStatus::Cancelled {
            return;
        }
        job.status = status;
        job.progress = progress.min(100);
        if let Some(artifacts) = artifacts {
            job.result = Some(artifacts.result);
            job.result_summary = Some(artifacts.result_summary);
            job.evidence_pages = artifacts.evidence_pages;
            job.evidence_search_index = artifacts.evidence_search_index;
            job.html_report = Some(artifacts.html_report);
        }
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
