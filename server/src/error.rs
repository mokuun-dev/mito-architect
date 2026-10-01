use axum::http::StatusCode;
use axum::response::IntoResponse;
use axum::response::Response;
use axum::Json;
use serde::Serialize;

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub(super) struct FailureBody {
    pub(super) schema_version: &'static str,
    pub(super) code: String,
    pub(super) message: String,
    pub(super) retryable: bool,
}

impl FailureBody {
    pub(super) fn new(
        code: impl Into<String>,
        message: impl Into<String>,
        retryable: bool,
    ) -> Self {
        Self {
            schema_version: "1.0",
            code: code.into(),
            message: message.into(),
            retryable,
        }
    }
}

#[derive(Debug)]
pub(super) struct ApiError {
    pub(super) status: StatusCode,
    pub(super) body: FailureBody,
}

impl ApiError {
    pub(super) fn new(status: StatusCode, message: impl Into<String>) -> Self {
        let (code, retryable) = match status {
            StatusCode::BAD_REQUEST => ("MITO-API-E1001", false),
            StatusCode::UNAUTHORIZED => ("MITO-API-E1002", false),
            StatusCode::NOT_FOUND => ("MITO-API-E1003", false),
            StatusCode::CONFLICT => ("MITO-API-E1004", false),
            StatusCode::PAYLOAD_TOO_LARGE => ("MITO-API-E1005", false),
            StatusCode::TOO_MANY_REQUESTS => ("MITO-API-E1006", true),
            StatusCode::ACCEPTED => ("MITO-API-E1007", true),
            StatusCode::SERVICE_UNAVAILABLE => ("MITO-API-E2001", true),
            _ => ("MITO-API-E9001", false),
        };
        Self::coded(status, code, message, retryable)
    }

    pub(super) fn coded(
        status: StatusCode,
        code: impl Into<String>,
        message: impl Into<String>,
        retryable: bool,
    ) -> Self {
        Self {
            status,
            body: FailureBody::new(code, message, retryable),
        }
    }
}

impl IntoResponse for ApiError {
    fn into_response(self) -> Response {
        (self.status, Json(serde_json::json!({ "error": self.body }))).into_response()
    }
}

pub(super) fn internal_error(error: impl std::fmt::Display) -> ApiError {
    eprintln!("internal server error: {error}");
    ApiError::coded(
        StatusCode::INTERNAL_SERVER_ERROR,
        "MITO-API-E9001",
        "internal server error",
        false,
    )
}
