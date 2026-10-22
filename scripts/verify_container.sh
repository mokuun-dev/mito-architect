#!/usr/bin/env bash
# Smoke-tests a built container image: bundled resources, the CLI, and the
# HTTP API with its unauthenticated-mode guards.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${1:-mito-architect:local}"
PORT="${MITO_SMOKE_PORT:-18080}"
NAME="mito-smoke-$$"
API="http://127.0.0.1:$PORT"

cleanup() {
  docker rm -f "$NAME" "$NAME-guard" >/dev/null 2>&1 || true
}
trap cleanup EXIT

docker run --rm -v "$ROOT_DIR/fixtures:/fixtures:ro" "$IMAGE" \
  mito-cli analyze --input /fixtures/truth_mixed.sam --json \
  | jq -e '
      .metadata.schema_version == "0.5"
      and .filter_stats.passed_reads == 2
    ' >/dev/null
echo "CLI analysis inside the image passed"

# The image binds 0.0.0.0, so without a key or the explicit opt-in the server
# must refuse to start (exit status 1) instead of serving.
set +e
timeout 30 docker run --rm --name "$NAME-guard" "$IMAGE" >/dev/null 2>&1
guard_status=$?
set -e
if [[ "$guard_status" != "1" ]]; then
  echo "server without MITO_API_KEY exited with $guard_status instead of refusing to start" >&2
  exit 1
fi

docker run -d --name "$NAME" -p "127.0.0.1:$PORT:8080" \
  -e MITO_ALLOW_UNAUTHENTICATED=1 "$IMAGE" >/dev/null
for _ in $(seq 1 30); do
  if curl -fsS "$API/readyz" >/dev/null 2>&1; then
    break
  fi
  sleep 1
done
curl -fsS "$API/readyz" | jq -e '.status == "ready" and .htslib_enabled' >/dev/null

JOB_ID="$(curl -fsS -F "file=@$ROOT_DIR/fixtures/truth_mixed.sam" "$API/upload" | jq -r '.job_id')"
STATUS=""
for _ in $(seq 1 60); do
  STATUS="$(curl -fsS "$API/status/$JOB_ID" | jq -r '.status')"
  if [[ "$STATUS" == "done" || "$STATUS" == "error" ]]; then
    break
  fi
  sleep 1
done
if [[ "$STATUS" != "done" ]]; then
  echo "container job finished with status '$STATUS'" >&2
  docker logs "$NAME" >&2 || true
  exit 1
fi
curl -fsS "$API/result/$JOB_ID" \
  | jq -e '.metadata.schema_version == "0.5" and .filter_stats.passed_reads == 2' >/dev/null

for header in "Host: evil.example" "Origin: https://evil.example"; do
  code="$(curl -s -o /dev/null -w '%{http_code}' -H "$header" "$API/status/$JOB_ID")"
  if [[ "$code" != "403" ]]; then
    echo "request with '$header' returned $code instead of 403" >&2
    exit 1
  fi
done
echo "Container API smoke passed"
