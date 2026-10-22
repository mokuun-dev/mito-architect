#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ACCESSION="${MITO_PUBLIC_ACCESSION:-SRR18110025}"
if [[ "$ACCESSION" != SRR18110025 ]]; then
  echo 'Only the metadata-pinned SRR18110025 run is supported.' >&2
  exit 2
fi
OUT_DIR="${MITO_PUBLIC_OUT_DIR:-$ROOT_DIR/.data/public/$ACCESSION}"
ALIGNMENT_REFERENCE="${MITO_ALIGNMENT_REFERENCE:-$ROOT_DIR/core/data/rcrs.fasta}"
for tool in python3 minimap2 samtools; do
  command -v "$tool" >/dev/null || { echo "Missing tool: $tool" >&2; exit 2; }
done
[[ -s "$ALIGNMENT_REFERENCE" ]] || { echo "Missing alignment reference: $ALIGNMENT_REFERENCE" >&2; exit 2; }
MT_ONLY_MODE=0
if [[ "$(realpath "$ALIGNMENT_REFERENCE")" == "$(realpath "$ROOT_DIR/core/data/rcrs.fasta")" ]]; then
  MT_ONLY_MODE=1
fi
python3 "$ROOT_DIR/scripts/public_data/fetch_reads.py" \
  --output "$OUT_DIR/reads" \
  --max-reads "${MITO_PUBLIC_MAX_READS:-200}" \
  --max-download-bytes "${MITO_PUBLIC_MAX_DOWNLOAD_BYTES:-8388608}" \
  --max-decoded-bytes "${MITO_PUBLIC_MAX_DECODED_BYTES:-67108864}"
MITO_ALLOW_MT_ONLY="$MT_ONLY_MODE" MITO_MAP_THREADS="${MITO_PUBLIC_THREADS:-2}" MITO_SAMPLE_ID="$ACCESSION" \
  bash "$ROOT_DIR/scripts/map_competitive.sh" \
  "$OUT_DIR/reads/$ACCESSION.fastq" "$ALIGNMENT_REFERENCE" "$OUT_DIR/alignment"
echo "Analyze: cargo run -p mito-cli --offline -- analyze -i $OUT_DIR/alignment/competitive.bam --json"
echo 'Optional --evidence-graph may reach the phase-pair budget on noisy full-length reads (MITO-E1601).'
echo 'An rCRS-only reference does not establish NUMT specificity.'
