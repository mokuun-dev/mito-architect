use crate::error::{internal_error, ApiError, FailureBody};
use crate::evidence::{normalized_search_filter, EvidenceSearchFilters, EvidenceSearchQuery};
use crate::state::{AppState, JobStatus, StatusResponse};
use axum::body::Body;
use axum::extract::Path as AxumPath;
use axum::extract::Query;
use axum::extract::State;
use axum::http::header;
use axum::http::HeaderValue;
use axum::http::StatusCode;
use axum::response::IntoResponse;
use axum::response::Response;
use axum::Json;
use std::sync::atomic::Ordering;
use uuid::Uuid;

pub(super) async fn status(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<String>,
) -> Result<Json<StatusResponse>, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let jobs = state.jobs.read().await;
    let job = jobs
        .get(&job_id)
        .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;
    Ok(Json(StatusResponse {
        job_id,
        status: job.status.clone(),
        progress: job.progress,
        error: job.error.clone(),
    }))
}

pub(super) async fn cancel(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<String>,
) -> Result<Json<StatusResponse>, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let response = {
        let mut jobs = state.jobs.write().await;
        let job = jobs
            .get_mut(&job_id)
            .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;

        match job.status {
            JobStatus::Done => {
                return Err(ApiError::new(
                    StatusCode::CONFLICT,
                    "completed jobs cannot be cancelled",
                ));
            }
            JobStatus::Error | JobStatus::Cancelled => {}
            JobStatus::Queued | JobStatus::Processing => {
                job.cancel_requested = true;
                job.status = JobStatus::Cancelled;
                job.progress = 100;
                job.error = Some(FailureBody::new("MITO-E1501", "cancelled by user", false));
            }
        }

        StatusResponse {
            job_id,
            status: job.status.clone(),
            progress: job.progress,
            error: job.error.clone(),
        }
    };
    if let Some(flag) = state.cancel_flags.read().await.get(&job_id) {
        flag.store(true, Ordering::Relaxed);
    }

    Ok(Json(response))
}

pub(super) async fn result(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<String>,
) -> Result<Response, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let jobs = state.jobs.read().await;
    let job = jobs
        .get(&job_id)
        .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;

    match (&job.status, &job.result) {
        (JobStatus::Done, Some(value)) => {
            let mut response = Response::new(Body::from(value.clone()));
            response.headers_mut().insert(
                header::CONTENT_TYPE,
                HeaderValue::from_static("application/json; charset=utf-8"),
            );
            response.headers_mut().insert(
                header::CONTENT_DISPOSITION,
                HeaderValue::from_str(&format!("attachment; filename=\"mito-{job_id}.json\""))
                    .map_err(internal_error)?,
            );
            Ok(response)
        }
        (JobStatus::Error, _) => {
            let failure = job
                .error
                .clone()
                .unwrap_or_else(|| FailureBody::new("MITO-E9001", "analysis failed", false));
            Err(ApiError::coded(
                StatusCode::CONFLICT,
                failure.code,
                failure.message,
                failure.retryable,
            ))
        }
        (JobStatus::Cancelled, _) => Err(ApiError::coded(
            StatusCode::CONFLICT,
            "MITO-E1501",
            "job was cancelled",
            false,
        )),
        _ => Err(ApiError::new(
            StatusCode::ACCEPTED,
            "analysis result is not ready yet",
        )),
    }
}

pub(super) async fn result_summary(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<String>,
) -> Result<Response, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let jobs = state.jobs.read().await;
    let job = jobs
        .get(&job_id)
        .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;

    match (&job.status, &job.result_summary) {
        (JobStatus::Done, Some(value)) => {
            let mut response = Response::new(Body::from(value.clone()));
            response.headers_mut().insert(
                header::CONTENT_TYPE,
                HeaderValue::from_static("application/json; charset=utf-8"),
            );
            response.headers_mut().insert(
                header::CACHE_CONTROL,
                HeaderValue::from_static("private, immutable, max-age=86400"),
            );
            Ok(response)
        }
        (JobStatus::Error, _) => {
            let failure = job
                .error
                .clone()
                .unwrap_or_else(|| FailureBody::new("MITO-E9001", "analysis failed", false));
            Err(ApiError::coded(
                StatusCode::CONFLICT,
                failure.code,
                failure.message,
                failure.retryable,
            ))
        }
        (JobStatus::Cancelled, _) => Err(ApiError::coded(
            StatusCode::CONFLICT,
            "MITO-E1501",
            "job was cancelled",
            false,
        )),
        _ => Err(ApiError::new(
            StatusCode::ACCEPTED,
            "analysis result is not ready yet",
        )),
    }
}

pub(super) async fn search_evidence(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<String>,
    Query(query): Query<EvidenceSearchQuery>,
) -> Result<Response, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let limit = query.limit.unwrap_or(100);
    if !(1..=500).contains(&limit) {
        return Err(ApiError::new(
            StatusCode::BAD_REQUEST,
            "evidence search limit must be between 1 and 500",
        ));
    }
    let filters = EvidenceSearchFilters {
        molecule_id: normalized_search_filter(query.molecule_id),
        event_id: normalized_search_filter(query.event_id),
        state: normalized_search_filter(query.state).map(|state| state.to_ascii_uppercase()),
    };
    if let Some(state) = &filters.state {
        if !matches!(
            state.as_str(),
            "REFERENCE" | "ALTERNATE" | "EVENT_ABSENT" | "LOW_QUALITY" | "CONFLICT"
        ) {
            return Err(ApiError::new(
                StatusCode::BAD_REQUEST,
                "evidence search state is not a stored observation state",
            ));
        }
    }
    let search_index = {
        let jobs = state.jobs.read().await;
        let job = jobs
            .get(&job_id)
            .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;
        match job.status {
            JobStatus::Done => job.evidence_search_index.clone().ok_or_else(|| {
                ApiError::new(
                    StatusCode::CONFLICT,
                    "analysis result has no searchable evidence graph",
                )
            })?,
            JobStatus::Error | JobStatus::Cancelled => {
                return Err(ApiError::new(
                    StatusCode::CONFLICT,
                    "analysis has no searchable evidence graph",
                ));
            }
            JobStatus::Queued | JobStatus::Processing => {
                return Err(ApiError::new(
                    StatusCode::ACCEPTED,
                    "analysis result is not ready yet",
                ));
            }
        }
    };
    let mut response =
        Json(search_index.search(filters, query.cursor.unwrap_or(0), limit)).into_response();
    response.headers_mut().insert(
        header::CACHE_CONTROL,
        HeaderValue::from_static("private, immutable, max-age=86400"),
    );
    Ok(response)
}

pub(super) async fn evidence_page(
    State(state): State<AppState>,
    AxumPath((job_id, page_index)): AxumPath<(String, String)>,
) -> Result<Response, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let page_index = page_index.parse::<usize>().map_err(|_| {
        ApiError::new(
            StatusCode::BAD_REQUEST,
            "evidence page index must be a non-negative integer",
        )
    })?;
    let page = {
        let jobs = state.jobs.read().await;
        let job = jobs
            .get(&job_id)
            .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;
        match job.status {
            JobStatus::Done => {
                job.evidence_pages.get(page_index).cloned().ok_or_else(|| {
                    ApiError::new(StatusCode::NOT_FOUND, "evidence page not found")
                })?
            }
            JobStatus::Error | JobStatus::Cancelled => {
                return Err(ApiError::new(
                    StatusCode::CONFLICT,
                    "analysis has no evidence pages",
                ));
            }
            JobStatus::Queued | JobStatus::Processing => {
                return Err(ApiError::new(
                    StatusCode::ACCEPTED,
                    "analysis result is not ready yet",
                ));
            }
        }
    };
    let mut response = Response::new(Body::from(page));
    response.headers_mut().insert(
        header::CONTENT_TYPE,
        HeaderValue::from_static("application/json; charset=utf-8"),
    );
    response.headers_mut().insert(
        header::CACHE_CONTROL,
        HeaderValue::from_static("private, immutable, max-age=86400"),
    );
    Ok(response)
}

pub(super) async fn download(
    State(state): State<AppState>,
    AxumPath(job_id): AxumPath<String>,
) -> Result<Response, ApiError> {
    let job_id = parse_job_id(&job_id)?;
    let jobs = state.jobs.read().await;
    let job = jobs
        .get(&job_id)
        .ok_or_else(|| ApiError::new(StatusCode::NOT_FOUND, "job not found"))?;
    let html = job
        .html_report
        .clone()
        .ok_or_else(|| ApiError::new(StatusCode::ACCEPTED, "HTML report is not ready yet"))?;

    let mut response = Response::new(Body::from(html));
    response.headers_mut().insert(
        header::CONTENT_TYPE,
        HeaderValue::from_static("text/html; charset=utf-8"),
    );
    response.headers_mut().insert(
        header::CONTENT_DISPOSITION,
        HeaderValue::from_str(&format!("attachment; filename=\"mito-{job_id}.html\""))
            .map_err(internal_error)?,
    );
    Ok(response)
}

pub(super) fn parse_job_id(value: &str) -> Result<Uuid, ApiError> {
    value
        .parse::<Uuid>()
        .map_err(|_| ApiError::new(StatusCode::BAD_REQUEST, "job ID must be a UUID"))
}
