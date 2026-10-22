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

/// Rejects browser requests from origins outside MITO_CORS_ORIGINS and, when
/// no API key protects the API, Host names outside the loopback allowlist.
/// CORS alone does not stop a cross-site multipart POST (a "simple" request)
/// from reaching the unauthenticated loopback API, and DNS rebinding presents
/// the attacker's own host name. Clients that send neither header (CLI tools)
/// are unaffected.
pub(super) async fn require_trusted_origin(
    State(state): State<AppState>,
    request: Request,
    next: Next,
) -> Result<Response, ApiError> {
    if let Some(origin) = request.headers().get(header::ORIGIN) {
        let allowed = state
            .allowed_origins
            .iter()
            .any(|candidate| candidate.as_bytes() == origin.as_bytes());
        if !allowed {
            return Err(ApiError::new(
                StatusCode::FORBIDDEN,
                "request origin is not allowed",
            ));
        }
    }
    if state.api_key.is_none() {
        let host = request
            .headers()
            .get(header::HOST)
            .and_then(|value| value.to_str().ok())
            .or_else(|| {
                request
                    .uri()
                    .authority()
                    .map(|authority| authority.as_str())
            });
        if let Some(host) = host {
            if !host_is_allowed(host, &state.allowed_hosts) {
                return Err(ApiError::new(
                    StatusCode::FORBIDDEN,
                    "request host is not allowed",
                ));
            }
        }
    }
    Ok(next.run(request).await)
}

/// Compares the host name of a Host header (port removed) with the allowlist.
pub(super) fn host_is_allowed(host: &str, allowed: &[String]) -> bool {
    let host = host.trim().to_ascii_lowercase();
    let name = if host.starts_with('[') {
        match host.find(']') {
            Some(end) => &host[..=end],
            None => return false,
        }
    } else {
        match host.rsplit_once(':') {
            Some((name, port)) if port.bytes().all(|byte| byte.is_ascii_digit()) => name,
            Some(_) => return false,
            None => host.as_str(),
        }
    };
    allowed.iter().any(|candidate| candidate == name)
}
