use crate::error::{internal_error, ApiError};
use crate::jobs::{now_epoch_secs, run_job};
use crate::state::{AppState, JobAnalysisOptions, JobRecord, JobStatus, UploadResponse};
use axum::extract::Multipart;
use axum::extract::State;
use axum::http::StatusCode;
use axum::Json;
use std::path::Path;
use std::path::PathBuf;
use std::sync::atomic::AtomicBool;
use std::sync::Arc;
use tokio::fs;
use tokio::io::AsyncWriteExt;
use uuid::Uuid;

pub(super) async fn upload(
    State(state): State<AppState>,
    multipart: Result<Multipart, axum::extract::multipart::MultipartRejection>,
) -> Result<Json<UploadResponse>, ApiError> {
    let mut multipart = multipart.map_err(|error| {
        eprintln!("multipart extraction failed: {error}");
        ApiError::new(StatusCode::BAD_REQUEST, "invalid multipart upload")
    })?;
    let job_permit = state
        .job_slots
        .clone()
        .try_acquire_owned()
        .map_err(|_| ApiError::new(StatusCode::TOO_MANY_REQUESTS, "analysis queue is full"))?;
    let mut analysis_options = JobAnalysisOptions {
        emit_evidence_graph: false,
        evidence_page_size: state.evidence_page_size,
        molecule_id_tag: String::new(),
        umi_tag: String::new(),
        duplex_tag: String::new(),
    };
    let mut uploaded: Option<(Uuid, PathBuf, Arc<AtomicBool>)> = None;
    loop {
        let next_field = match multipart.next_field().await {
            Ok(field) => field,
            Err(error) => {
                if let Some((job_id, _, _)) = &uploaded {
                    let _ = fs::remove_dir_all(state.tmp_dir.join(job_id.to_string())).await;
                }
                return Err(ApiError::new(StatusCode::BAD_REQUEST, error.to_string()));
            }
        };
        let Some(mut field) = next_field else {
            break;
        };
        let field_name = field.name().unwrap_or_default().to_owned();
        if matches!(
            field_name.as_str(),
            "evidence_graph" | "evidence_page_size" | "molecule_id_tag" | "umi_tag" | "duplex_tag"
        ) {
            let value = match field.text().await {
                Ok(value) => value,
                Err(error) => {
                    if let Some((job_id, _, _)) = &uploaded {
                        let _ = fs::remove_dir_all(state.tmp_dir.join(job_id.to_string())).await;
                    }
                    return Err(ApiError::new(StatusCode::BAD_REQUEST, error.to_string()));
                }
            };
            match field_name.as_str() {
                "evidence_graph" => {
                    analysis_options.emit_evidence_graph = match value.trim() {
                        "true" | "1" => true,
                        "false" | "0" => false,
                        _ => {
                            if let Some((job_id, _, _)) = &uploaded {
                                let _ = fs::remove_dir_all(state.tmp_dir.join(job_id.to_string()))
                                    .await;
                            }
                            return Err(ApiError::new(
                                StatusCode::BAD_REQUEST,
                                "evidence_graph must be true, false, 1, or 0",
                            ));
                        }
                    };
                }
                "evidence_page_size" => {
                    let parsed = value
                        .trim()
                        .parse::<usize>()
                        .ok()
                        .filter(|value| (1..=1_000_000).contains(value));
                    let Some(parsed) = parsed else {
                        if let Some((job_id, _, _)) = &uploaded {
                            let _ =
                                fs::remove_dir_all(state.tmp_dir.join(job_id.to_string())).await;
                        }
                        return Err(ApiError::new(
                            StatusCode::BAD_REQUEST,
                            "evidence_page_size must be between 1 and 1000000",
                        ));
                    };
                    analysis_options.evidence_page_size = parsed;
                }
                "molecule_id_tag" => analysis_options.molecule_id_tag = value.trim().to_owned(),
                "umi_tag" => analysis_options.umi_tag = value.trim().to_owned(),
                "duplex_tag" => analysis_options.duplex_tag = value.trim().to_owned(),
                _ => unreachable!("multipart protocol field was matched above"),
            }
            for tag in [
                &analysis_options.molecule_id_tag,
                &analysis_options.umi_tag,
                &analysis_options.duplex_tag,
            ] {
                if !valid_optional_sam_tag(tag) {
                    if let Some((job_id, _, _)) = &uploaded {
                        let _ = fs::remove_dir_all(state.tmp_dir.join(job_id.to_string())).await;
                    }
                    return Err(ApiError::new(
                        StatusCode::BAD_REQUEST,
                        "protocol SAM tags must match [A-Za-z][A-Za-z0-9]",
                    ));
                }
            }
            continue;
        }
        if field.name() != Some("file") {
            continue;
        }
        if let Some((job_id, _, _)) = &uploaded {
            let _ = fs::remove_dir_all(state.tmp_dir.join(job_id.to_string())).await;
            return Err(ApiError::new(
                StatusCode::BAD_REQUEST,
                "multipart upload must include exactly one file field",
            ));
        }

        let original_name = field.file_name().unwrap_or("input.fastq").to_string();
        validate_upload_name(&original_name)?;
        let job_id = Uuid::new_v4();
        let cancel_flag = Arc::new(AtomicBool::new(false));
        let job_dir = state.tmp_dir.join(job_id.to_string());
        fs::create_dir_all(&job_dir).await.map_err(internal_error)?;
        let input_path = job_dir.join(sanitize_file_name(&original_name));
        let mut output = match fs::File::create(&input_path).await {
            Ok(output) => output,
            Err(error) => {
                let _ = fs::remove_dir_all(&job_dir).await;
                return Err(internal_error(error));
            }
        };
        let mut written = 0usize;

        loop {
            let chunk = match field.chunk().await {
                Ok(Some(chunk)) => chunk,
                Ok(None) => break,
                Err(error) => {
                    drop(output);
                    let _ = fs::remove_dir_all(&job_dir).await;
                    return Err(ApiError::new(StatusCode::BAD_REQUEST, error.to_string()));
                }
            };
            written = match written.checked_add(chunk.len()) {
                Some(total) => total,
                None => {
                    drop(output);
                    let _ = fs::remove_dir_all(&job_dir).await;
                    return Err(ApiError::new(
                        StatusCode::PAYLOAD_TOO_LARGE,
                        "upload size overflow",
                    ));
                }
            };
            if written > state.max_upload_bytes {
                drop(output);
                let _ = fs::remove_dir_all(&job_dir).await;
                return Err(ApiError::new(
                    StatusCode::PAYLOAD_TOO_LARGE,
                    format!("upload exceeds the {} byte limit", state.max_upload_bytes),
                ));
            }
            if let Err(error) = output.write_all(&chunk).await {
                drop(output);
                let _ = fs::remove_dir_all(&job_dir).await;
                return Err(internal_error(error));
            }
        }

        if written == 0 {
            drop(output);
            let _ = fs::remove_dir_all(&job_dir).await;
            return Err(ApiError::new(
                StatusCode::BAD_REQUEST,
                "uploaded file is empty",
            ));
        }
        if let Err(error) = output.flush().await {
            drop(output);
            let _ = fs::remove_dir_all(&job_dir).await;
            return Err(internal_error(error));
        }
        drop(output);
        uploaded = Some((job_id, input_path, cancel_flag));
    }

    let Some((job_id, input_path, cancel_flag)) = uploaded else {
        return Err(ApiError::new(
            StatusCode::BAD_REQUEST,
            "multipart upload must include a file field named 'file'",
        ));
    };
    if !analysis_options.emit_evidence_graph
        && (!analysis_options.molecule_id_tag.is_empty()
            || !analysis_options.umi_tag.is_empty()
            || !analysis_options.duplex_tag.is_empty())
    {
        let _ = fs::remove_dir_all(state.tmp_dir.join(job_id.to_string())).await;
        return Err(ApiError::new(
            StatusCode::BAD_REQUEST,
            "protocol SAM tags require evidence_graph=true",
        ));
    }
    let record = JobRecord {
        status: JobStatus::Queued,
        progress: 0,
        input_path: input_path.clone(),
        result: None,
        result_summary: None,
        evidence_pages: Vec::new(),
        evidence_search_index: None,
        html_report: None,
        error: None,
        cancel_requested: false,
        created_at: now_epoch_secs(),
    };
    state.jobs.write().await.insert(job_id, record);
    state
        .cancel_flags
        .write()
        .await
        .insert(job_id, cancel_flag.clone());

    let job_state = state.clone();
    let job_options = analysis_options.clone();
    tokio::spawn(async move {
        run_job(
            job_state,
            job_id,
            input_path,
            cancel_flag,
            job_options,
            job_permit,
        )
        .await;
    });

    Ok(Json(UploadResponse {
        job_id,
        status: JobStatus::Queued,
        result_schema: if analysis_options.emit_evidence_graph {
            "0.6"
        } else {
            "0.5"
        },
    }))
}

pub(super) fn validate_upload_name(file_name: &str) -> Result<(), ApiError> {
    let path = Path::new(file_name);
    let extension = path
        .extension()
        .and_then(|value| value.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    match extension.as_str() {
        "fastq" | "fq" | "sam" | "bam" | "cram" => Ok(()),
        _ => Err(ApiError::new(
            StatusCode::BAD_REQUEST,
            "expected .fastq, .fq, .sam, .bam, or .cram input",
        )),
    }
}

pub(super) fn valid_optional_sam_tag(tag: &str) -> bool {
    if tag.is_empty() {
        return true;
    }
    let bytes = tag.as_bytes();
    bytes.len() == 2 && bytes[0].is_ascii_alphabetic() && bytes[1].is_ascii_alphanumeric()
}

pub(super) fn sanitize_file_name(file_name: &str) -> String {
    let sanitized: String = file_name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-') {
                c
            } else {
                '_'
            }
        })
        .collect();
    if sanitized.is_empty() {
        "input.fastq".to_string()
    } else {
        sanitized
    }
}
