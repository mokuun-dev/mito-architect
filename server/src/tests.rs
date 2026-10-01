use crate::architecture_contract::validate_architecture_projection;
use crate::auth::constant_time_eq;
use crate::build_app;
use crate::contract::validate_result_contract;
use crate::error::ApiError;
use crate::jobs::public_engine_failure;
use crate::state::{AppState, JobRecord, JobStatus};
use crate::upload::{sanitize_file_name, validate_upload_name};

use axum::body::{to_bytes, Body, Bytes};
use axum::http::{header, Request, StatusCode};
use serde_json::{json, Value};
use std::collections::HashMap;
use std::sync::Arc;
use tokio::sync::{RwLock, Semaphore};
use tower::ServiceExt;
use tower_http::cors::CorsLayer;
use uuid::Uuid;

fn test_state(api_key: Option<&str>) -> AppState {
    let tmp_dir = std::env::temp_dir().join(format!("mito-server-http-test-{}", Uuid::new_v4()));
    std::fs::create_dir_all(&tmp_dir).expect("test temporary directory should be created");
    AppState {
        jobs: Arc::new(RwLock::new(HashMap::new())),
        cancel_flags: Arc::new(RwLock::new(HashMap::new())),
        analysis_slots: Arc::new(Semaphore::new(2)),
        job_slots: Arc::new(Semaphore::new(4)),
        worker_threads: 1,
        min_mapping_quality: 20,
        min_base_quality: 10,
        excluded_snp_flags: 0xF00,
        numt_threshold: 0.30,
        max_evidence_observations: 5_000_000,
        max_phase_links: 1_000_000,
        max_phase_work: 20_000_000,
        max_phase_molecule_references: 5_000_000,
        max_result_bytes: 128 * 1024 * 1024,
        evidence_page_size: 4096,
        min_architecture_molecules: 2,
        max_candidate_architectures: 256,
        architecture_max_distance: 0.20,
        architecture_ambiguity_margin: 0.05,
        architecture_min_overlap_fraction: 0.50,
        architecture_consensus_fraction: 0.80,
        architecture_optional_fraction: 0.20,
        architecture_stability_replicates: 32,
        architecture_seed: 0x4d49_544f_4152_4348,
        engine_version: Arc::from("0.5.0-dev"),
        schema_version: Arc::from("0.5"),
        error_schema_version: Arc::from("1.0"),
        htslib_enabled: true,
        api_key: api_key.map(Arc::<str>::from),
        max_upload_bytes: 1024 * 1024,
        tmp_dir,
    }
}

async fn response_json(response: axum::response::Response) -> serde_json::Value {
    let bytes = to_bytes(response.into_body(), 1024 * 1024)
        .await
        .expect("response body should be readable");
    serde_json::from_slice(&bytes).expect("response body should be JSON")
}

fn get_request(uri: &str, token: Option<&str>) -> Request<Body> {
    let mut builder = Request::builder().method("GET").uri(uri);
    if let Some(token) = token {
        builder = builder.header(header::AUTHORIZATION, format!("Bearer {token}"));
    }
    builder
        .body(Body::empty())
        .expect("test request should be valid")
}

fn post_request(uri: &str, token: Option<&str>) -> Request<Body> {
    let mut builder = Request::builder().method("POST").uri(uri);
    if let Some(token) = token {
        builder = builder.header(header::AUTHORIZATION, format!("Bearer {token}"));
    }
    builder
        .body(Body::empty())
        .expect("test request should be valid")
}

fn multipart_request(file_name: &str, content: &str) -> Request<Body> {
    let boundary = "mito-architect-test-boundary";
    let body = format!(
        "--{boundary}\r\nContent-Disposition: form-data; name=\"file\"; filename=\"{file_name}\"\r\nContent-Type: application/octet-stream\r\n\r\n{content}\r\n--{boundary}--\r\n"
    );
    Request::builder()
        .method("POST")
        .uri("/upload")
        .header(header::AUTHORIZATION, "Bearer laboratory-secret")
        .header(
            header::CONTENT_TYPE,
            format!("multipart/form-data; boundary={boundary}"),
        )
        .body(Body::from(body))
        .expect("multipart test request should be valid")
}

fn multipart_evidence_request(file_first: bool) -> Request<Body> {
    let boundary = "mito-architect-evidence-boundary";
    let file = format!(
        "--{boundary}\r\nContent-Disposition: form-data; name=\"file\"; filename=\"phase.sam\"\r\nContent-Type: application/octet-stream\r\n\r\n@SQ\tSN:MT\tLN:16569\nphase\t0\tMT\t1\t60\t12M\t*\t0\t0\tGAACGCAGGTCT\tIIIIIIIIIIII\tMI:Z:M1\tRX:Z:AAA\tDX:Z:duplex\n\r\n"
    );
    let options = format!(
        "--{boundary}\r\nContent-Disposition: form-data; name=\"evidence_graph\"\r\n\r\ntrue\r\n--{boundary}\r\nContent-Disposition: form-data; name=\"evidence_page_size\"\r\n\r\n1\r\n--{boundary}\r\nContent-Disposition: form-data; name=\"molecule_id_tag\"\r\n\r\nMI\r\n--{boundary}\r\nContent-Disposition: form-data; name=\"umi_tag\"\r\n\r\nRX\r\n--{boundary}\r\nContent-Disposition: form-data; name=\"duplex_tag\"\r\n\r\nDX\r\n"
    );
    let body = if file_first {
        format!("{file}{options}--{boundary}--\r\n")
    } else {
        format!("{options}{file}--{boundary}--\r\n")
    };
    Request::builder()
        .method("POST")
        .uri("/upload")
        .header(header::AUTHORIZATION, "Bearer laboratory-secret")
        .header(
            header::CONTENT_TYPE,
            format!("multipart/form-data; boundary={boundary}"),
        )
        .body(Body::from(body))
        .expect("evidence multipart request should be valid")
}

#[test]
fn result_contract_requires_runtime_metadata() {
    let valid = json!({
        "metadata": {
            "schema_version": "0.5",
            "sv_event_schema_version": "1.0",
            "complex_sv_event_schema_version": "1.0",
            "clinical_annotation_schema_version": "1.0",
            "engine_version": "0.5.0-dev",
            "reference_length": 16569,
            "resources": []
        },
        "filter_stats": {
            "passed_reads": 2,
            "numt_assessment": { "mode": "competitive_alignment" }
        },
        "reads": [],
        "variants": [],
        "svs": [],
        "complex_events": [],
        "clusters": []
    });
    validate_result_contract(&valid).expect("valid result should satisfy contract");

    let missing_schema = json!({
        "metadata": {
            "engine_version": "0.5.0-dev",
            "reference_length": 16569
        },
        "filter_stats": {
            "passed_reads": 2,
            "numt_assessment": { "mode": "competitive_alignment" }
        },
        "reads": [],
        "variants": [],
        "svs": [],
        "clusters": []
    });
    let error =
        validate_result_contract(&missing_schema).expect_err("schema_version must be required");
    assert!(error.contains("/metadata/schema_version"));
}

#[test]
fn upload_names_are_sanitized_and_extension_checked() {
    assert_eq!(
        sanitize_file_name("../patient sample.sam"),
        ".._patient_sample.sam"
    );
    validate_upload_name("sample.fastq").expect("FASTQ should be accepted");
    validate_upload_name("sample.bam").expect("BAM should be accepted");

    let error = validate_upload_name("sample.exe").expect_err("unknown inputs rejected");
    assert_eq!(error.status, StatusCode::BAD_REQUEST);
}

#[test]
fn api_key_comparison_rejects_length_and_content_mismatches() {
    assert!(constant_time_eq(b"laboratory-secret", b"laboratory-secret"));
    assert!(!constant_time_eq(
        b"laboratory-secret",
        b"laboratory-secreu"
    ));
    assert!(!constant_time_eq(
        b"laboratory-secret",
        b"laboratory-secret-extra"
    ));
}

#[test]
fn api_errors_have_a_stable_machine_readable_envelope() {
    let error = ApiError::new(StatusCode::PAYLOAD_TOO_LARGE, "too large");
    assert_eq!(error.body.schema_version, "1.0");
    assert_eq!(error.body.code, "MITO-API-E1005");
    assert_eq!(error.body.message, "too large");
    assert!(!error.body.retryable);
}

#[test]
fn engine_open_failures_do_not_expose_server_paths() {
    let failure = public_engine_failure(
        mito_ffi::MitoError::new(
            "MITO-E1101",
            "could not open input file: /tmp/private/sample.bam",
        ),
        uuid::Uuid::nil(),
    );
    assert_eq!(failure.code, "MITO-E1101");
    assert_eq!(failure.message, "input could not be opened");
    assert!(!failure.message.contains("/tmp/private"));
}

#[tokio::test]
async fn router_preserves_http_contract_across_auth_and_job_states() {
    let state = test_state(Some("laboratory-secret"));
    let tmp_dir = state.tmp_dir.clone();
    let app = build_app(state.clone(), CorsLayer::new());

    let response = app
        .clone()
        .oneshot(get_request("/healthz", None))
        .await
        .expect("health request should complete");
    assert_eq!(response.status(), StatusCode::OK);
    assert_eq!(response_json(response).await, json!({ "status": "ok" }));

    let response = app
        .clone()
        .oneshot(get_request("/readyz", None))
        .await
        .expect("readiness request should complete");
    assert_eq!(response.status(), StatusCode::OK);
    assert_eq!(
        response_json(response).await,
        json!({
            "status": "ready",
            "engine_version": "0.5.0-dev",
            "schema_version": "0.5",
            "error_schema_version": "1.0",
            "htslib_enabled": true,
            "available_analysis_slots": 2,
            "available_job_slots": 4
        })
    );

    let unknown_id = Uuid::nil();
    let response = app
        .clone()
        .oneshot(get_request(&format!("/status/{unknown_id}"), None))
        .await
        .expect("unauthorized request should complete");
    assert_eq!(response.status(), StatusCode::UNAUTHORIZED);
    assert_eq!(
        response_json(response).await,
        json!({
            "error": {
                "schema_version": "1.0",
                "code": "MITO-API-E1002",
                "message": "invalid or missing API key",
                "retryable": false
            }
        })
    );

    let response = app
        .clone()
        .oneshot(get_request(
            &format!("/status/{unknown_id}"),
            Some("laboratory-secret"),
        ))
        .await
        .expect("not-found request should complete");
    assert_eq!(response.status(), StatusCode::NOT_FOUND);
    assert_eq!(
        response_json(response).await,
        json!({
            "error": {
                "schema_version": "1.0",
                "code": "MITO-API-E1003",
                "message": "job not found",
                "retryable": false
            }
        })
    );

    let response = app
        .clone()
        .oneshot(get_request("/status/not-a-uuid", Some("laboratory-secret")))
        .await
        .expect("malformed job ID request should complete");
    assert_eq!(response.status(), StatusCode::BAD_REQUEST);
    assert_eq!(
        response_json(response).await,
        json!({
            "error": {
                "schema_version": "1.0",
                "code": "MITO-API-E1001",
                "message": "job ID must be a UUID",
                "retryable": false
            }
        })
    );

    let response = app
        .clone()
        .oneshot(post_request("/upload", Some("laboratory-secret")))
        .await
        .expect("malformed multipart request should complete");
    assert_eq!(response.status(), StatusCode::BAD_REQUEST);
    assert_eq!(
        response_json(response).await,
        json!({
            "error": {
                "schema_version": "1.0",
                "code": "MITO-API-E1001",
                "message": "invalid multipart upload",
                "retryable": false
            }
        })
    );

    let response = app
        .clone()
        .oneshot(multipart_request("payload.exe", "x"))
        .await
        .expect("invalid extension upload should complete");
    assert_eq!(response.status(), StatusCode::BAD_REQUEST);
    assert_eq!(
        response_json(response).await,
        json!({
            "error": {
                "schema_version": "1.0",
                "code": "MITO-API-E1001",
                "message": "expected .fastq, .fq, .sam, .bam, or .cram input",
                "retryable": false
            }
        })
    );

    let job_id = Uuid::new_v4();
    state.jobs.write().await.insert(
        job_id,
        JobRecord {
            status: JobStatus::Queued,
            progress: 0,
            input_path: tmp_dir.join("input.sam"),
            result: None,
            result_summary: None,
            evidence_pages: Vec::new(),
            evidence_search_index: None,
            html_report: None,
            error: None,
            cancel_requested: false,
            created_at: 1,
        },
    );
    let response = app
        .clone()
        .oneshot(get_request(
            &format!("/status/{job_id}"),
            Some("laboratory-secret"),
        ))
        .await
        .expect("status request should complete");
    assert_eq!(response.status(), StatusCode::OK);
    assert_eq!(
        response_json(response).await,
        json!({
            "job_id": job_id,
            "status": "queued",
            "progress": 0,
            "error": null
        })
    );

    let response = app
        .clone()
        .oneshot(get_request(
            &format!("/result/{job_id}"),
            Some("laboratory-secret"),
        ))
        .await
        .expect("pending result request should complete");
    assert_eq!(response.status(), StatusCode::ACCEPTED);
    assert_eq!(
        response_json(response).await,
        json!({
            "error": {
                "schema_version": "1.0",
                "code": "MITO-API-E1007",
                "message": "analysis result is not ready yet",
                "retryable": true
            }
        })
    );

    {
        let mut jobs = state.jobs.write().await;
        let job = jobs.get_mut(&job_id).expect("job should exist");
        job.status = JobStatus::Done;
        job.progress = 100;
        job.result = Some(Bytes::from_static(
            br#"{"metadata":{"schema_version":"0.5"}}"#,
        ));
        job.result_summary = job.result.clone();
        job.html_report = Some(Bytes::from_static(b"<!doctype html><title>result</title>"));
    }
    let response = app
        .clone()
        .oneshot(get_request(
            &format!("/result/{job_id}"),
            Some("laboratory-secret"),
        ))
        .await
        .expect("completed result request should complete");
    assert_eq!(response.status(), StatusCode::OK);
    assert_eq!(
        response
            .headers()
            .get(header::CONTENT_TYPE)
            .and_then(|value| value.to_str().ok()),
        Some("application/json; charset=utf-8")
    );
    assert_eq!(
        response_json(response).await,
        json!({ "metadata": { "schema_version": "0.5" } })
    );

    let response = app
        .clone()
        .oneshot(get_request(
            &format!("/download/{job_id}"),
            Some("laboratory-secret"),
        ))
        .await
        .expect("download request should complete");
    assert_eq!(response.status(), StatusCode::OK);
    assert_eq!(
        response
            .headers()
            .get(header::CONTENT_TYPE)
            .and_then(|value| value.to_str().ok()),
        Some("text/html; charset=utf-8")
    );
    let content_disposition = response
        .headers()
        .get(header::CONTENT_DISPOSITION)
        .and_then(|value| value.to_str().ok())
        .expect("download should include content disposition");
    assert_eq!(
        content_disposition,
        format!("attachment; filename=\"mito-{job_id}.html\"")
    );
    let body = to_bytes(response.into_body(), 1024)
        .await
        .expect("download body should be readable");
    assert_eq!(
        body,
        Bytes::from_static(b"<!doctype html><title>result</title>")
    );

    let response = app
        .oneshot(post_request(
            &format!("/cancel/{job_id}"),
            Some("laboratory-secret"),
        ))
        .await
        .expect("completed cancellation request should complete");
    assert_eq!(response.status(), StatusCode::CONFLICT);
    assert_eq!(
        response_json(response).await,
        json!({
            "error": {
                "schema_version": "1.0",
                "code": "MITO-API-E1004",
                "message": "completed jobs cannot be cancelled",
                "retryable": false
            }
        })
    );

    std::fs::remove_dir_all(tmp_dir).expect("test temporary directory should be removed");
}

#[tokio::test]
async fn evidence_graph_selection_is_independent_of_multipart_field_order() {
    for file_first in [false, true] {
        let state = test_state(Some("laboratory-secret"));
        let tmp_dir = state.tmp_dir.clone();
        let app = build_app(state.clone(), CorsLayer::new());
        let response = app
            .clone()
            .oneshot(multipart_evidence_request(file_first))
            .await
            .expect("evidence upload should complete");
        assert_eq!(response.status(), StatusCode::OK);
        let uploaded = response_json(response).await;
        assert_eq!(
            uploaded
                .get("result_schema")
                .and_then(serde_json::Value::as_str),
            Some("0.6")
        );
        let job_id = uploaded
            .get("job_id")
            .and_then(serde_json::Value::as_str)
            .expect("upload should return a job ID");

        let mut completed = false;
        for _ in 0..2_000 {
            let response = app
                .clone()
                .oneshot(get_request(
                    &format!("/status/{job_id}"),
                    Some("laboratory-secret"),
                ))
                .await
                .expect("status request should complete");
            let status = response_json(response).await;
            match status.get("status").and_then(serde_json::Value::as_str) {
                Some("done") => {
                    completed = true;
                    break;
                }
                Some("error") | Some("cancelled") => {
                    panic!("schema 0.6 upload failed: {status}");
                }
                _ => tokio::time::sleep(std::time::Duration::from_millis(10)).await,
            }
        }
        assert!(completed, "schema 0.6 upload did not complete in time");
        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/result/{job_id}"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("result request should complete");
        assert_eq!(response.status(), StatusCode::OK);
        let result = response_json(response).await;
        assert_eq!(
            result
                .pointer("/metadata/schema_version")
                .and_then(serde_json::Value::as_str),
            Some("0.6")
        );
        assert_eq!(
            result
                .pointer("/molecules/0/identity_policy")
                .and_then(serde_json::Value::as_str),
            Some("sam_tag:MI")
        );
        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/result/{job_id}/summary"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("result summary request should complete");
        assert_eq!(response.status(), StatusCode::OK);
        let summary = response_json(response).await;
        assert_eq!(
            summary
                .pointer("/evidence_encoding/observation_storage")
                .and_then(Value::as_str),
            Some("remote_http_pages")
        );
        let expected_page_endpoint = format!("/result/{job_id}/evidence/{{page_index}}");
        assert_eq!(
            summary
                .pointer("/evidence_encoding/observation_page_endpoint")
                .and_then(Value::as_str),
            Some(expected_page_endpoint.as_str())
        );
        let expected_search_endpoint = format!("/result/{job_id}/evidence");
        assert_eq!(
            summary
                .pointer("/evidence_encoding/observation_search_endpoint")
                .and_then(Value::as_str),
            Some(expected_search_endpoint.as_str())
        );
        assert_eq!(
            summary
                .pointer("/observation_pages")
                .and_then(Value::as_array)
                .map(Vec::len),
            Some(0)
        );
        let expected_search_count = summary
            .pointer("/evidence_encoding/observation_count")
            .and_then(Value::as_u64)
            .expect("summary should declare the global observation count");
        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/result/{job_id}/evidence/0"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("evidence page request should complete");
        assert_eq!(response.status(), StatusCode::OK);
        let page = response_json(response).await;
        assert_eq!(page.get("index").and_then(Value::as_u64), Some(0));
        assert_eq!(page.get("count").and_then(Value::as_u64), Some(1));
        let expected_state = page
            .pointer("/columns/state/0")
            .and_then(Value::as_str)
            .expect("evidence page should contain one stored state")
            .to_string();
        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/result/{job_id}/evidence?molecule_id=MI:M1&limit=1"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("global evidence search should complete");
        assert_eq!(response.status(), StatusCode::OK);
        let search = response_json(response).await;
        assert_eq!(
            search.get("schema_version").and_then(Value::as_str),
            Some("1.0")
        );
        assert_eq!(
            search.get("total_matches").and_then(Value::as_u64),
            Some(expected_search_count)
        );
        assert_eq!(
            search
                .pointer("/rows/0/molecule_id")
                .and_then(Value::as_str),
            Some("MI:M1")
        );
        assert_eq!(
            search.pointer("/rows/0/page_index").and_then(Value::as_u64),
            Some(0)
        );
        assert_eq!(
            search.pointer("/rows/0/row_index").and_then(Value::as_u64),
            Some(0)
        );
        assert_eq!(
            search.pointer("/rows/0/state").and_then(Value::as_str),
            Some(expected_state.as_str())
        );
        assert_eq!(search.get("next_cursor").and_then(Value::as_u64), Some(1));
        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/result/{job_id}/evidence?molecule_id=MI:M1&cursor=1&limit=1"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("second evidence-search cursor should complete");
        assert_eq!(response.status(), StatusCode::OK);
        let second_search_page = response_json(response).await;
        assert_eq!(
            second_search_page
                .pointer("/rows/0/page_index")
                .and_then(Value::as_u64),
            Some(1)
        );
        assert!(second_search_page
            .get("next_cursor")
            .is_some_and(Value::is_null));
        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/result/{job_id}/evidence?event_id=missing"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("empty evidence search should complete");
        assert_eq!(response.status(), StatusCode::OK);
        assert_eq!(
            response_json(response)
                .await
                .get("total_matches")
                .and_then(Value::as_u64),
            Some(0)
        );
        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/result/{job_id}/evidence?state=NOT_CALLABLE"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("invalid sparse-state search should complete");
        assert_eq!(response.status(), StatusCode::BAD_REQUEST);
        let response = app
            .clone()
            .oneshot(get_request(
                &format!("/result/{job_id}/evidence/999"),
                Some("laboratory-secret"),
            ))
            .await
            .expect("missing evidence page request should complete");
        assert_eq!(response.status(), StatusCode::NOT_FOUND);
        std::fs::remove_dir_all(tmp_dir).expect("test temporary directory should be removed");
    }
}

#[test]
fn architecture_event_profiles_fail_closed_on_count_drift() {
    let mut result: Value = serde_json::from_str(include_str!(
        "../../fixtures/truth_architectures.expected.json"
    ))
    .expect("architecture golden should parse");
    validate_architecture_projection(&result).expect("architecture golden should validate");

    result["architectures"][0]["event_profile"][0]["not_callable"] = json!(1);
    let error = validate_architecture_projection(&result)
        .expect_err("architecture member-count drift must fail closed");
    assert!(error.contains("does not partition members"));
}
