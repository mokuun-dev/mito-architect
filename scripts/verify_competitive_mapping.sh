#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

for command_name in minimap2 samtools jq sha256sum; do
  if ! command -v "$command_name" >/dev/null 2>&1; then
    echo "$command_name is required for competitive-mapping verification." >&2
    exit 2
  fi
done

printf '>chr1\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n' \
  > "$WORK_DIR/nuclear.fa"
MITO_NUCLEAR_FASTA="$WORK_DIR/nuclear.fa" \
MITO_REFERENCE_DIR="$WORK_DIR/reference" \
  bash "$ROOT_DIR/scripts/build_competitive_reference.sh" >/dev/null

printf 'field\tvalue\nschema_version\tcompetitive-reference-0.9\n' \
  > "$WORK_DIR/not-competitive-manifest.tsv"

SEQUENCE="$(awk '!/^>/{printf "%s", $0}' "$ROOT_DIR/core/data/rcrs.fasta" | cut -c1-1000)"
QUALITY="$(printf '%*s' "${#SEQUENCE}" '' | tr ' ' 'I')"
{
  printf '@competitive-mtdna-read\n%s\n+\n%s\n' "$SEQUENCE" "$QUALITY"
} > "$WORK_DIR/read.fastq"

if MITO_REFERENCE_MANIFEST="$WORK_DIR/not-competitive-manifest.tsv" \
  bash "$ROOT_DIR/scripts/map_competitive.sh" \
    "$WORK_DIR/read.fastq" "$WORK_DIR/reference/nuclear-plus-rcrs.mmi" \
    "$WORK_DIR/rejected-output" >/dev/null 2>&1; then
  echo "mapping accepted a manifest that did not attest competitive reference content" >&2
  exit 1
fi
test ! -e "$WORK_DIR/rejected-output"

MITO_MAP_THREADS=2 MITO_SAMPLE_ID=competitive-smoke \
MITO_REFERENCE_MANIFEST="$WORK_DIR/reference/reference-manifest.txt" \
  bash "$ROOT_DIR/scripts/map_competitive.sh" \
    "$WORK_DIR/read.fastq" "$WORK_DIR/reference/nuclear-plus-rcrs.mmi" \
    "$WORK_DIR/output" >/dev/null

samtools quickcheck -v "$WORK_DIR/output/competitive.bam"
test -s "$WORK_DIR/output/competitive.bam.bai"
MANIFEST="$WORK_DIR/output/competitive-mapping-manifest.tsv"
test -s "$MANIFEST"
for field in schema_version input_sha256 reference_index_sha256 \
  reference_manifest_sha256 minimap2_version samtools_version mapping_command \
  sorting_command output_bam_sha256 output_bai_sha256 numt_assessability; do
  if [[ "$(awk -F '\t' -v key="$field" '$1 == key { count += 1 } END { print count + 0 }' "$MANIFEST")" != "1" ]]; then
    echo "competitive mapping manifest field is missing or duplicated: $field" >&2
    exit 1
  fi
done
[[ "$(awk -F '\t' '$1 == "schema_version" { print $2 }' "$MANIFEST")" == "competitive-mapping-1.1" ]]
[[ "$(awk -F '\t' '$1 == "numt_assessability" { print $2 }' "$MANIFEST")" == "ASSESSABLE_BY_VERSIONED_COMPETITIVE_ALIGNMENT_INPUT" ]]

if bash "$ROOT_DIR/scripts/map_competitive.sh" \
  "$WORK_DIR/read.fastq" "$ROOT_DIR/core/data/rcrs.fasta" \
  "$WORK_DIR/mt-only-rejected" >/dev/null 2>&1; then
  echo "mtDNA-only mapping was accepted without an explicit opt-in" >&2
  exit 1
fi
MITO_ALLOW_MT_ONLY=1 MITO_SAMPLE_ID=mt-only-smoke \
  bash "$ROOT_DIR/scripts/map_competitive.sh" \
    "$WORK_DIR/read.fastq" "$ROOT_DIR/core/data/rcrs.fasta" \
    "$WORK_DIR/mt-only" >/dev/null
[[ "$(awk -F '\t' '$1 == "numt_assessability" { print $2 }' "$WORK_DIR/mt-only/competitive-mapping-manifest.tsv")" == "NOT_ASSESSABLE_MT_ONLY" ]]
"$ROOT_DIR/target/debug/mito-cli" analyze \
  --input "$WORK_DIR/mt-only/competitive.bam" --json \
  > "$WORK_DIR/mt-only-result.json"
jq -e '.filter_stats.numt_assessment.specificity_assessable == false' \
  "$WORK_DIR/mt-only-result.json" >/dev/null

"$ROOT_DIR/target/debug/mito-cli" analyze \
  --input "$WORK_DIR/output/competitive.bam" --evidence-graph --json \
  > "$WORK_DIR/result.json"
if ! jq -e '
  .filter_stats.numt_assessment.mode == "competitive_alignment" and
  .filter_stats.numt_assessment.nuclear_contigs_present == true and
  .filter_stats.numt_assessment.specificity_assessable == true
' "$WORK_DIR/result.json" >/dev/null; then
  echo "competitive BAM did not produce assessable NUMT provenance:" >&2
  jq '.filter_stats.numt_assessment' "$WORK_DIR/result.json" >&2
  exit 1
fi

echo "competitive mapping verification passed (index, BAM/BAI, provenance, NUMT assessability)"
