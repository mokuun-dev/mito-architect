#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLI="$ROOT_DIR/target/debug/mito-cli"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

if [[ ! -x "$CLI" ]]; then
  echo "Build target/debug/mito-cli before architecture verification." >&2
  exit 2
fi
if ! command -v jq >/dev/null 2>&1; then
  echo "jq is required for architecture verification." >&2
  exit 2
fi

"$CLI" validate-evidence-graph-fixture \
  --input "$ROOT_DIR/fixtures/truth_architectures.sam" \
  --expected-json "$ROOT_DIR/fixtures/truth_architectures.expected.json"

analyze() {
  local input="$1"
  local threads="$2"
  local output="$3"
  "$CLI" analyze --input "$input" --evidence-graph --threads "$threads" --json > "$output"
}

analyze "$ROOT_DIR/fixtures/truth_architectures.sam" 1 "$WORK_DIR/base-t1.json"
analyze "$ROOT_DIR/fixtures/truth_architectures.sam" 4 "$WORK_DIR/base-t4.json"
analyze "$ROOT_DIR/fixtures/truth_architectures_permuted.sam" 4 "$WORK_DIR/permuted-t4.json"
analyze "$ROOT_DIR/fixtures/truth_noisy_architecture.sam" 1 "$WORK_DIR/noisy.json"
analyze "$ROOT_DIR/fixtures/truth_evidence_graph.sam" 1 "$WORK_DIR/low-quality.json"

project() {
  jq -S '{
    architecture_inference,
    architectures,
    molecule_assignments: [.molecules[] | {id, architecture_assignment}]
  }' "$1"
}

project "$WORK_DIR/base-t1.json" > "$WORK_DIR/base-t1.projection.json"
project "$WORK_DIR/base-t4.json" > "$WORK_DIR/base-t4.projection.json"
project "$WORK_DIR/permuted-t4.json" > "$WORK_DIR/permuted-t4.projection.json"
cmp -s "$WORK_DIR/base-t1.projection.json" "$WORK_DIR/base-t4.projection.json"
cmp -s "$WORK_DIR/base-t1.projection.json" "$WORK_DIR/permuted-t4.projection.json"

jq -e '
  .architecture_inference.status == "CANDIDATES_INFERRED" and
  .architecture_inference.candidate_count == 3 and
  .architecture_inference.eligible_molecules == 12 and
  .architecture_inference.assigned_molecules == 11 and
  .architecture_inference.ambiguous_molecules == 0 and
  .architecture_inference.unassigned_molecules == 1 and
  ([.architectures[].molecule_count] | sort) == [3, 4, 4] and
  ([.architectures[].defining_event_signature] | any(. == [])) and
  ([.architectures[].defining_event_signature] | any(. == ["snv:3:T:A"])) and
  ([.architectures[].defining_event_signature] |
    any(. == ["snv:3:T:A", "sv:deletion:11-20"])) and
  ([.architectures[].method.not_callable_is_reference] | all(. == false)) and
  ([.molecules[] | select(.id == "partial-snv-only") |
    .architecture_assignment.status] == ["UNASSIGNED"])
' "$WORK_DIR/base-t1.json" >/dev/null

jq -e '
  .architecture_inference.status == "CANDIDATES_INFERRED" and
  .architecture_inference.eligible_molecules == 4 and
  .architecture_inference.assigned_molecules == 4 and
  .architecture_inference.unassigned_molecules == 0 and
  (.architectures | length) == 1 and
  .architectures[0].defining_event_signature == ["snv:3:T:A"] and
  .architectures[0].molecule_count == 4
' "$WORK_DIR/noisy.json" >/dev/null

jq -e '
  [.molecules[] | select(.id == "low-quality") |
    .architecture_assignment.status] == ["INELIGIBLE"] and
  .architecture_inference.eligible_molecules == 4
' "$WORK_DIR/low-quality.json" >/dev/null

expect_analysis_failure() {
  local expected_code="$1"
  shift
  if "$CLI" analyze --input "$ROOT_DIR/fixtures/truth_architectures.sam" \
      --evidence-graph --json "$@" \
      >"$WORK_DIR/rejected.stdout" 2>"$WORK_DIR/rejected.stderr"; then
    echo "architecture analysis unexpectedly accepted unsafe configuration: $*" >&2
    exit 1
  fi
  if ! grep -F "[$expected_code]" "$WORK_DIR/rejected.stderr" >/dev/null; then
    echo "architecture failure did not preserve $expected_code: $*" >&2
    cat "$WORK_DIR/rejected.stderr" >&2
    exit 1
  fi
}

expect_analysis_failure MITO-E1001 --min-architecture-molecules 1
expect_analysis_failure MITO-E1001 --max-candidate-architectures 4097
expect_analysis_failure MITO-E1001 --architecture-stability-replicates 10001
expect_analysis_failure MITO-E1601 --max-candidate-architectures 2

echo "architecture inference verification passed (golden, noise, threads, order, callable refusal, limits)"
