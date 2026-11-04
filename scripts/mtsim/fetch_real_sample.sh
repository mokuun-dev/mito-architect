#!/usr/bin/env bash
# Reads and alignment for the real sample of the example panel: the first 3,000
# reads of the public GM12878 nanopore run. The server can be slow to respond.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT_DIR"
OUT=".data/public/SRR18110025_3k"
THREADS="${MITO_PUBLIC_THREADS:-4}"
for tool in python3 minimap2 samtools; do
  command -v "$tool" >/dev/null || { echo "Missing tool: $tool" >&2; exit 2; }
done
mkdir -p "$OUT"
python3 scripts/public_data/fetch_reads.py --output "$OUT/reads" --max-reads 3000 \
  --max-download-bytes 120000000 --max-decoded-bytes 400000000 --timeout 280
minimap2 -a -x map-ont -Y --MD --secondary=no -t "$THREADS" core/data/rcrs.fasta \
  "$OUT/reads/SRR18110025.fastq" 2>/dev/null \
  | samtools sort -@ "$THREADS" -O BAM -o "$OUT/reads.bam" -
samtools index "$OUT/reads.bam"
echo "Wrote $OUT/reads.bam"
