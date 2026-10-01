#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NUCLEAR_FASTA="${MITO_NUCLEAR_FASTA:-}"
OUTPUT_DIR="${MITO_REFERENCE_DIR:-$ROOT_DIR/.data/reference}"

if [[ -z "$NUCLEAR_FASTA" || ! -f "$NUCLEAR_FASTA" ]]; then
  echo "Set MITO_NUCLEAR_FASTA to a versioned human nuclear FASTA." >&2
  exit 2
fi
if ! command -v minimap2 >/dev/null 2>&1; then
  echo "minimap2 is required to build the competitive reference index." >&2
  exit 2
fi

mkdir -p "$OUTPUT_DIR"
COMBINED_FASTA="$OUTPUT_DIR/nuclear-plus-rcrs.fasta"
TEMP_FASTA="$COMBINED_FASTA.part"

# Remove any pre-existing mitochondrial contig before appending the bundled
# NC_012920.1 sequence, preventing ambiguous duplicate reference names.
awk '
  /^>/ {
    name = substr($1, 2)
    keep = !(name == "chrM" || name == "MT" || name == "M" || name == "NC_012920.1")
  }
  keep { print }
' "$NUCLEAR_FASTA" > "$TEMP_FASTA"
awk '1' "$ROOT_DIR/core/data/rcrs.fasta" >> "$TEMP_FASTA"
mv "$TEMP_FASTA" "$COMBINED_FASTA"

minimap2 -d "$OUTPUT_DIR/nuclear-plus-rcrs.mmi" "$COMBINED_FASTA"
{
  printf 'field\tvalue\n'
  printf 'schema_version\tcompetitive-reference-1.0\n'
  printf 'minimap2_version\t%s\n' "$(minimap2 --version)"
  printf 'nuclear_fasta_path\t%s\n' "$(realpath "$NUCLEAR_FASTA")"
  printf 'nuclear_fasta_sha256\t%s\n' "$(sha256sum "$NUCLEAR_FASTA" | awk '{print $1}')"
  printf 'mitochondrial_fasta_path\t%s\n' "$(realpath "$ROOT_DIR/core/data/rcrs.fasta")"
  printf 'mitochondrial_fasta_sha256\t%s\n' "$(sha256sum "$ROOT_DIR/core/data/rcrs.fasta" | awk '{print $1}')"
  printf 'combined_fasta_path\t%s\n' "$(realpath "$COMBINED_FASTA")"
  printf 'combined_fasta_sha256\t%s\n' "$(sha256sum "$COMBINED_FASTA" | awk '{print $1}')"
  printf 'reference_index_path\t%s\n' "$(realpath "$OUTPUT_DIR/nuclear-plus-rcrs.mmi")"
  printf 'reference_index_sha256\t%s\n' "$(sha256sum "$OUTPUT_DIR/nuclear-plus-rcrs.mmi" | awk '{print $1}')"
} > "$OUTPUT_DIR/reference-manifest.txt"

echo "Competitive reference: $OUTPUT_DIR/nuclear-plus-rcrs.mmi"
echo "Provenance manifest: $OUTPUT_DIR/reference-manifest.txt"
