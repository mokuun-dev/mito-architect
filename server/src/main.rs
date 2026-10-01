mod architecture_contract;
mod auth;
mod config;
mod contract;
mod error;
mod evidence;
mod jobs;
mod report;
mod routes;
mod state;
#[cfg(test)]
mod tests;
mod upload;

use crate::auth::{health, ready, require_api_key};
use crate::config::{
    cors_layer, env_u16, env_u64, env_u8, nonnegative_env_usize, optional_env_number,
    positive_env_usize, unit_interval_env, DEFAULT_MAX_UPLOAD_BYTES,
};
use crate::jobs::spawn_cleanup;
use crate::routes::{
    cancel, download, evidence_page, result, result_summary, search_evidence, status,
};
use crate::state::AppState;
use crate::upload::upload;
use axum::extract::DefaultBodyLimit;
use axum::middleware;
use axum::routing::get;
use axum::routing::post;
use axum::Router;
use mito_ffi::MitoEngine;
use std::collections::HashMap;
use std::net::SocketAddr;
use std::path::PathBuf;
use std::sync::Arc;
use tokio::fs;
use tokio::sync::RwLock;
use tokio::sync::Semaphore;
use tower_http::cors::CorsLayer;

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    let tmp_dir = PathBuf::from("./tmp");
    fs::create_dir_all(&tmp_dir).await?;

    let available_threads = std::thread::available_parallelism().map_or(1, usize::from);
    let worker_threads = positive_env_usize("MITO_JOB_THREADS")?.unwrap_or(available_threads);
    let max_concurrent_jobs = positive_env_usize("MITO_MAX_CONCURRENT_JOBS")?.unwrap_or(1);
    let max_active_jobs = positive_env_usize("MITO_MAX_ACTIVE_JOBS")?.unwrap_or(16);
    let max_upload_bytes =
        positive_env_usize("MITO_MAX_UPLOAD_BYTES")?.unwrap_or(DEFAULT_MAX_UPLOAD_BYTES);
    let max_evidence_observations =
        positive_env_usize("MITO_MAX_EVIDENCE_OBSERVATIONS")?.unwrap_or(5_000_000);
    let max_phase_links = positive_env_usize("MITO_MAX_PHASE_LINKS")?.unwrap_or(1_000_000);
    let max_phase_work = positive_env_usize("MITO_MAX_PHASE_WORK")?.unwrap_or(20_000_000);
    let max_phase_molecule_references =
        positive_env_usize("MITO_MAX_PHASE_MOLECULE_REFERENCES")?.unwrap_or(5_000_000);
    let max_result_bytes =
        positive_env_usize("MITO_MAX_RESULT_BYTES")?.unwrap_or(128 * 1024 * 1024);
    let evidence_page_size = positive_env_usize("MITO_EVIDENCE_PAGE_SIZE")?.unwrap_or(4096);
    if evidence_page_size > 1_000_000 {
        anyhow::bail!("MITO_EVIDENCE_PAGE_SIZE must not exceed 1000000");
    }
    let min_mapping_quality = env_u8("MITO_MIN_MAPQ")?.unwrap_or(20);
    let min_base_quality = env_u8("MITO_MIN_BASE_QUALITY")?.unwrap_or(10);
    let excluded_snp_flags = env_u16("MITO_EXCLUDED_SNP_FLAGS")?.unwrap_or(0xF00);
    let numt_threshold = optional_env_number("MITO_NUMT_THRESHOLD")?.unwrap_or(0.30);
    if !(0.0..=1.0).contains(&numt_threshold) {
        anyhow::bail!("MITO_NUMT_THRESHOLD must be between 0 and 1");
    }
    let min_architecture_molecules =
        positive_env_usize("MITO_MIN_ARCHITECTURE_MOLECULES")?.unwrap_or(2);
    if min_architecture_molecules < 2 {
        anyhow::bail!("MITO_MIN_ARCHITECTURE_MOLECULES must be at least 2");
    }
    let max_candidate_architectures =
        positive_env_usize("MITO_MAX_CANDIDATE_ARCHITECTURES")?.unwrap_or(256);
    if max_candidate_architectures > 4096 {
        anyhow::bail!("MITO_MAX_CANDIDATE_ARCHITECTURES must not exceed 4096");
    }
    let architecture_max_distance = unit_interval_env("MITO_ARCHITECTURE_MAX_DISTANCE", 0.20)?;
    let architecture_ambiguity_margin =
        unit_interval_env("MITO_ARCHITECTURE_AMBIGUITY_MARGIN", 0.05)?;
    let architecture_min_overlap_fraction =
        unit_interval_env("MITO_ARCHITECTURE_MIN_OVERLAP", 0.50)?;
    let architecture_consensus_fraction = unit_interval_env("MITO_ARCHITECTURE_CONSENSUS", 0.80)?;
    if architecture_consensus_fraction <= 0.5 {
        anyhow::bail!("MITO_ARCHITECTURE_CONSENSUS must be greater than 0.5");
    }
    let architecture_optional_fraction = unit_interval_env("MITO_ARCHITECTURE_OPTIONAL", 0.20)?;
    if architecture_optional_fraction > architecture_consensus_fraction {
        anyhow::bail!("MITO_ARCHITECTURE_OPTIONAL must not exceed consensus");
    }
    let architecture_stability_replicates =
        nonnegative_env_usize("MITO_ARCHITECTURE_STABILITY_REPLICATES")?.unwrap_or(32);
    if architecture_stability_replicates > 10_000 {
        anyhow::bail!("MITO_ARCHITECTURE_STABILITY_REPLICATES must not exceed 10000");
    }
    let architecture_seed = env_u64("MITO_ARCHITECTURE_SEED")?.unwrap_or(0x4d49_544f_4152_4348);
    let capabilities = MitoEngine::capabilities();
    if !capabilities.htslib {
        anyhow::bail!("mito-server requires a native core built with htslib for BAM/CRAM support");
    }
    let api_key = std::env::var("MITO_API_KEY")
        .ok()
        .filter(|value| !value.is_empty())
        .map(Arc::<str>::from);
    let addr = std::env::var("MITO_SERVER_ADDR")
        .ok()
        .and_then(|value| value.parse::<SocketAddr>().ok())
        .unwrap_or_else(|| SocketAddr::from(([127, 0, 0, 1], 8080)));
    if !addr.ip().is_loopback() && api_key.is_none() {
        anyhow::bail!("MITO_API_KEY is required when MITO_SERVER_ADDR is not loopback");
    }
    let state = AppState {
        jobs: Arc::new(RwLock::new(HashMap::new())),
        cancel_flags: Arc::new(RwLock::new(HashMap::new())),
        analysis_slots: Arc::new(Semaphore::new(max_concurrent_jobs)),
        job_slots: Arc::new(Semaphore::new(max_active_jobs)),
        worker_threads,
        min_mapping_quality,
        min_base_quality,
        excluded_snp_flags,
        numt_threshold,
        max_evidence_observations,
        max_phase_links,
        max_phase_work,
        max_phase_molecule_references,
        max_result_bytes,
        evidence_page_size,
        min_architecture_molecules,
        max_candidate_architectures,
        architecture_max_distance,
        architecture_ambiguity_margin,
        architecture_min_overlap_fraction,
        architecture_consensus_fraction,
        architecture_optional_fraction,
        architecture_stability_replicates,
        architecture_seed,
        engine_version: Arc::from(capabilities.engine_version),
        schema_version: Arc::from(capabilities.schema_version),
        error_schema_version: Arc::from(capabilities.error_schema_version),
        htslib_enabled: capabilities.htslib,
        api_key,
        max_upload_bytes,
        tmp_dir,
    };

    spawn_cleanup(state.clone());
    let cors = cors_layer()?;
    let app = build_app(state, cors);

    let listener = tokio::net::TcpListener::bind(addr).await?;
    println!("mito-server listening on http://{addr}");
    axum::serve(listener, app).await?;
    Ok(())
}

fn build_app(state: AppState, cors: CorsLayer) -> Router {
    let max_upload_bytes = state.max_upload_bytes;
    let protected = Router::new()
        .route("/upload", post(upload))
        .route("/status/:job_id", get(status))
        .route("/result/:job_id", get(result))
        .route("/result/:job_id/summary", get(result_summary))
        .route("/result/:job_id/evidence", get(search_evidence))
        .route("/result/:job_id/evidence/:page_index", get(evidence_page))
        .route("/download/:job_id", get(download))
        .route("/cancel/:job_id", post(cancel))
        .route_layer(middleware::from_fn_with_state(
            state.clone(),
            require_api_key,
        ));
    Router::new()
        .route("/healthz", get(health))
        .route("/readyz", get(ready))
        .merge(protected)
        .layer(DefaultBodyLimit::max(max_upload_bytes))
        .layer(cors)
        .with_state(state)
}
