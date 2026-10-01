use crate::error::ApiError;
use crate::state::AppState;
use axum::extract::Request;
use axum::extract::State;
use axum::http::header;
use axum::http::StatusCode;
use axum::middleware::Next;
use axum::response::Response;
use axum::Json;
use serde_json::Value;
use tokio::fs;

pub(super) async fn health() -> Json<Value> {
    Json(serde_json::json!({ "status": "ok" }))
}

pub(super) async fn ready(State(state): State<AppState>) -> Result<Json<Value>, ApiError> {
    fs::metadata(&state.tmp_dir)
        .await
        .map_err(|error| ApiError::new(StatusCode::SERVICE_UNAVAILABLE, error.to_string()))?;
    Ok(Json(serde_json::json!({
        "status": "ready",
        "engine_version": state.engine_version.as_ref(),
        "schema_version": state.schema_version.as_ref(),
        "error_schema_version": state.error_schema_version.as_ref(),
        "htslib_enabled": state.htslib_enabled,
        "available_analysis_slots": state.analysis_slots.available_permits(),
        "available_job_slots": state.job_slots.available_permits()
    })))
}

pub(super) async fn require_api_key(
    State(state): State<AppState>,
    request: Request,
    next: Next,
) -> Result<Response, ApiError> {
    let Some(expected) = state.api_key.as_deref() else {
        return Ok(next.run(request).await);
    };
    let supplied = request
        .headers()
        .get(header::AUTHORIZATION)
        .and_then(|value| value.to_str().ok())
        .and_then(|value| value.strip_prefix("Bearer "));
    if !supplied.is_some_and(|value| constant_time_eq(value.as_bytes(), expected.as_bytes())) {
        return Err(ApiError::new(
            StatusCode::UNAUTHORIZED,
            "invalid or missing API key",
        ));
    }
    Ok(next.run(request).await)
}

pub(super) fn constant_time_eq(left: &[u8], right: &[u8]) -> bool {
    let mut difference = left.len() ^ right.len();
    let maximum = left.len().max(right.len());
    for index in 0..maximum {
        let lhs = left.get(index).copied().unwrap_or(0);
        let rhs = right.get(index).copied().unwrap_or(0);
        difference |= usize::from(lhs ^ rhs);
    }
    difference == 0
}
