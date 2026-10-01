use crate::error::FailureBody;
use crate::evidence::EvidenceSearchIndex;
use axum::body::Bytes;
use serde::Serialize;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::atomic::AtomicBool;
use std::sync::Arc;
use tokio::sync::RwLock;
use tokio::sync::Semaphore;
use uuid::Uuid;

#[derive(Clone)]
pub(super) struct AppState {
    pub(super) jobs: Arc<RwLock<HashMap<Uuid, JobRecord>>>,
    pub(super) cancel_flags: Arc<RwLock<HashMap<Uuid, Arc<AtomicBool>>>>,
    pub(super) analysis_slots: Arc<Semaphore>,
    pub(super) job_slots: Arc<Semaphore>,
    pub(super) worker_threads: usize,
    pub(super) min_mapping_quality: u8,
    pub(super) min_base_quality: u8,
    pub(super) excluded_snp_flags: u16,
    pub(super) numt_threshold: f64,
    pub(super) max_evidence_observations: usize,
    pub(super) max_phase_links: usize,
    pub(super) max_phase_work: usize,
    pub(super) max_phase_molecule_references: usize,
    pub(super) max_result_bytes: usize,
    pub(super) evidence_page_size: usize,
    pub(super) min_architecture_molecules: usize,
    pub(super) max_candidate_architectures: usize,
    pub(super) architecture_max_distance: f64,
    pub(super) architecture_ambiguity_margin: f64,
    pub(super) architecture_min_overlap_fraction: f64,
    pub(super) architecture_consensus_fraction: f64,
    pub(super) architecture_optional_fraction: f64,
    pub(super) architecture_stability_replicates: usize,
    pub(super) architecture_seed: u64,
    pub(super) engine_version: Arc<str>,
    pub(super) schema_version: Arc<str>,
    pub(super) error_schema_version: Arc<str>,
    pub(super) htslib_enabled: bool,
    pub(super) api_key: Option<Arc<str>>,
    pub(super) max_upload_bytes: usize,
    pub(super) tmp_dir: PathBuf,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub(super) enum JobStatus {
    Queued,
    Processing,
    Done,
    Error,
    Cancelled,
}

#[derive(Debug)]
pub(super) struct JobRecord {
    pub(super) status: JobStatus,
    pub(super) progress: u8,
    pub(super) input_path: PathBuf,
    pub(super) result: Option<Bytes>,
    pub(super) result_summary: Option<Bytes>,
    pub(super) evidence_pages: Vec<Bytes>,
    pub(super) evidence_search_index: Option<Arc<EvidenceSearchIndex>>,
    pub(super) html_report: Option<Bytes>,
    pub(super) error: Option<FailureBody>,
    pub(super) cancel_requested: bool,
    pub(super) created_at: u64,
}

#[derive(Debug, Serialize)]
pub(super) struct UploadResponse {
    pub(super) job_id: Uuid,
    pub(super) status: JobStatus,
    pub(super) result_schema: &'static str,
}

#[derive(Clone, Debug)]
pub(super) struct JobAnalysisOptions {
    pub(super) emit_evidence_graph: bool,
    pub(super) evidence_page_size: usize,
    pub(super) molecule_id_tag: String,
    pub(super) umi_tag: String,
    pub(super) duplex_tag: String,
}

#[derive(Debug, Serialize)]
pub(super) struct StatusResponse {
    pub(super) job_id: Uuid,
    pub(super) status: JobStatus,
    pub(super) progress: u8,
    pub(super) error: Option<FailureBody>,
}
