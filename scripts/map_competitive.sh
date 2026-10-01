#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 INPUT.fastq[.gz] COMPETITIVE_REFERENCE.mmi OUTPUT_DIR" >&2
  echo "Environment: MITO_MAP_THREADS, MITO_SAMPLE_ID, MITO_REFERENCE_MANIFEST, MITO_FORCE=1" >&2
  echo "MITO_ALLOW_MT_ONLY=1 permits only the bundled rCRS and marks NUMT specificity not assessable." >&2
}

if [[ $# -ne 3 ]]; then
  usage
  exit 2
fi

INPUT_PATH="$1"
REFERENCE_INDEX="$2"
OUTPUT_DIR="$3"
THREADS="${MITO_MAP_THREADS:-1}"
SAMPLE_ID="${MITO_SAMPLE_ID:-$(basename "$INPUT_PATH")}"
REFERENCE_MANIFEST="${MITO_REFERENCE_MANIFEST:-$(dirname "$REFERENCE_INDEX")/reference-manifest.txt}"

if [[ ! -f "$INPUT_PATH" ]]; then
  echo "Input reads do not exist: $INPUT_PATH" >&2
  exit 2
fi
if [[ ! -f "$REFERENCE_INDEX" ]]; then
  echo "Competitive reference index does not exist: $REFERENCE_INDEX" >&2
  exit 2
fi
if [[ ! "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
  echo "MITO_MAP_THREADS must be a positive integer." >&2
  exit 2
fi
for command_name in minimap2 samtools sha256sum realpath; do
  if ! command -v "$command_name" >/dev/null 2>&1; then
    echo "$command_name is required for competitive mapping." >&2
    exit 2
  fi
done

# SAM read-group fields are tab-delimited and must not accept control bytes.
SAMPLE_ID="${SAMPLE_ID%%.*}"
if [[ -z "$SAMPLE_ID" || "$SAMPLE_ID" =~ [^A-Za-z0-9._-] ]]; then
  echo "MITO_SAMPLE_ID must match [A-Za-z0-9._-]+." >&2
  exit 2
fi

OUTPUT_PARENT="$(dirname "$OUTPUT_DIR")"
OUTPUT_NAME="$(basename "$OUTPUT_DIR")"
if [[ "$OUTPUT_NAME" == "." || "$OUTPUT_NAME" == "/" || -z "$OUTPUT_NAME" ]]; then
  echo "OUTPUT_DIR must name a dedicated output directory." >&2
  exit 2
fi
mkdir -p "$OUTPUT_PARENT"
OUTPUT_DIR="$(cd "$OUTPUT_PARENT" && pwd)/$OUTPUT_NAME"
OUTPUT_BAM="$OUTPUT_DIR/competitive.bam"
OUTPUT_BAI="$OUTPUT_BAM.bai"
OUTPUT_MANIFEST="$OUTPUT_DIR/competitive-mapping-manifest.tsv"

if [[ "${MITO_FORCE:-0}" != "1" ]] &&
   [[ -e "$OUTPUT_DIR" ]]; then
  echo "Competitive mapping output directory already exists; set MITO_FORCE=1 to replace it." >&2
  exit 2
fi

REFERENCE_REAL="$(realpath "$REFERENCE_INDEX")"
MT_ONLY_MODE=0
if [[ ! -f "$REFERENCE_MANIFEST" && "${MITO_ALLOW_MT_ONLY:-0}" == "1" &&
      "$REFERENCE_REAL" == "$(realpath "$(dirname "${BASH_SOURCE[0]}")/../core/data/rcrs.fasta")" ]]; then
  MT_ONLY_MODE=1
fi
if [[ "$MT_ONLY_MODE" == "0" && ! -f "$REFERENCE_MANIFEST" ]]; then
  echo "Competitive mapping requires the reference manifest from build_competitive_reference.sh." >&2
  exit 2
fi
if [[ "$MT_ONLY_MODE" == "0" ]] && ! awk -F '\t' '
  $1 == "schema_version" && $2 == "competitive-reference-1.0" { schema = 1 }
  $1 == "nuclear_fasta_sha256" && $2 ~ /^[0-9a-f]{64}$/ { nuclear = 1 }
  $1 == "mitochondrial_fasta_sha256" && $2 ~ /^[0-9a-f]{64}$/ { mito = 1 }
  $1 == "reference_index_path" { reference_path = $2 }
  END { exit !(schema && nuclear && mito && reference_path != "") }
' "$REFERENCE_MANIFEST"; then
  echo "Reference manifest does not attest a nuclear-plus-mitochondrial reference." >&2
  exit 2
fi
if [[ "$MT_ONLY_MODE" == "0" ]]; then
  MANIFEST_REFERENCE_INDEX="$(awk -F '\t' '$1 == "reference_index_path" { print $2; exit }' "$REFERENCE_MANIFEST")"
  if [[ "$(realpath "$MANIFEST_REFERENCE_INDEX")" != "$REFERENCE_REAL" ]]; then
    echo "Reference manifest is for a different minimap2 index." >&2
    exit 2
  fi
fi

cleanup() {
  [[ -n "${STAGE_DIR:-}" && -d "$STAGE_DIR" ]] && rm -rf "$STAGE_DIR"
}
trap cleanup EXIT

STAGE_DIR="$(mktemp -d "$OUTPUT_PARENT/.${OUTPUT_NAME}.stage.XXXXXX")"
PART_BAM="$STAGE_DIR/competitive.bam"
OUTPUT_BAI="$PART_BAM.bai"
OUTPUT_MANIFEST="$STAGE_DIR/competitive-mapping-manifest.tsv"
SORT_PREFIX="$STAGE_DIR/competitive-sort"

INPUT_REAL="$(realpath "$INPUT_PATH")"
READ_GROUP="@RG\\tID:${SAMPLE_ID}\\tSM:${SAMPLE_ID}\\tPL:ONT"
MINIMAP_COMMAND=(
  minimap2 -a -x map-ont --secondary=yes -Y --MD
  -R "$READ_GROUP" -t "$THREADS" "$REFERENCE_REAL" "$INPUT_REAL"
)
SORT_COMMAND=(
  samtools sort -@ "$THREADS" -T "$SORT_PREFIX" -O BAM -o "$PART_BAM" -
)

"${MINIMAP_COMMAND[@]}" | "${SORT_COMMAND[@]}"
samtools quickcheck -v "$PART_BAM"
samtools index -@ "$THREADS" "$PART_BAM" "$OUTPUT_BAI"
samtools quickcheck -v "$PART_BAM"

render_command() {
  local rendered=""
  local token
  for token in "$@"; do
    printf -v token '%q' "$token"
    rendered+="${rendered:+ }${token}"
  done
  printf '%s' "$rendered"
}

INPUT_SHA256="$(sha256sum "$INPUT_REAL" | awk '{print $1}')"
REFERENCE_SHA256="$(sha256sum "$REFERENCE_REAL" | awk '{print $1}')"
BAM_SHA256="$(sha256sum "$PART_BAM" | awk '{print $1}')"
BAI_SHA256="$(sha256sum "$OUTPUT_BAI" | awk '{print $1}')"
MINIMAP_VERSION="$(minimap2 --version | head -n 1)"
SAMTOOLS_VERSION="$(samtools --version | head -n 1)"
MAPPING_COMMAND="$(render_command "${MINIMAP_COMMAND[@]}")"
SORTING_COMMAND="$(render_command "${SORT_COMMAND[@]}")"
if [[ "$MT_ONLY_MODE" == "1" ]]; then
  REFERENCE_MANIFEST=NOT_PROVIDED
  REFERENCE_MANIFEST_SHA256=NOT_PROVIDED
  NUMT_ASSESSABILITY=NOT_ASSESSABLE_MT_ONLY
else
  REFERENCE_MANIFEST="$(realpath "$REFERENCE_MANIFEST")"
  REFERENCE_MANIFEST_SHA256="$(sha256sum "$REFERENCE_MANIFEST" | awk '{print $1}')"
  NUMT_ASSESSABILITY=ASSESSABLE_BY_VERSIONED_COMPETITIVE_ALIGNMENT_INPUT
fi

MANIFEST_PART="$OUTPUT_MANIFEST.part"
{
  printf 'field\tvalue\n'
  printf 'schema_version\tcompetitive-mapping-1.1\n'
  printf 'generated_at_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'sample_id\t%s\n' "$SAMPLE_ID"
  printf 'input_path\t%s\n' "$INPUT_REAL"
  printf 'input_sha256\t%s\n' "$INPUT_SHA256"
  printf 'reference_index_path\t%s\n' "$REFERENCE_REAL"
  printf 'reference_index_sha256\t%s\n' "$REFERENCE_SHA256"
  printf 'reference_manifest_path\t%s\n' "$REFERENCE_MANIFEST"
  printf 'reference_manifest_sha256\t%s\n' "$REFERENCE_MANIFEST_SHA256"
  printf 'minimap2_version\t%s\n' "$MINIMAP_VERSION"
  printf 'samtools_version\t%s\n' "$SAMTOOLS_VERSION"
  printf 'mapping_command\t%s\n' "$MAPPING_COMMAND"
  printf 'sorting_command\t%s\n' "$SORTING_COMMAND"
  printf 'threads\t%s\n' "$THREADS"
  printf 'output_bam_sha256\t%s\n' "$BAM_SHA256"
  printf 'output_bai_sha256\t%s\n' "$BAI_SHA256"
  printf 'numt_assessability\t%s\n' "$NUMT_ASSESSABILITY"
} > "$MANIFEST_PART"
mv "$MANIFEST_PART" "$OUTPUT_MANIFEST"

# Publish BAM, BAI and their provenance as one directory rename. Consumers
# never observe a final BAM without its index and matching manifest.
if [[ -e "$OUTPUT_DIR" ]]; then
  rm -rf "$OUTPUT_DIR"
fi
mv "$STAGE_DIR" "$OUTPUT_DIR"
STAGE_DIR=""

trap - EXIT
echo "Competitive BAM: $OUTPUT_BAM"
echo "Mapping manifest: $OUTPUT_DIR/competitive-mapping-manifest.tsv"
