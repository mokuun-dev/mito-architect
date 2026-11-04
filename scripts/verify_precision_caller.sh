#!/usr/bin/env bash
# End-to-end check of the statistical caller on freshly simulated reads: real
# alignments, known truth, and results that must not depend on the thread count.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"
for tool in minimap2 samtools jq; do
  command -v "$tool" >/dev/null || { echo "error: $tool is required" >&2; exit 2; }
done

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# Two lineages from an edit list on rCRS, so this needs no network access.
cat > "$WORK/scenario.json" <<'JSON'
{"id": "verify", "platform": "ont_r10", "coverage": 500, "seed": 7,
 "read_length": {"model": "lognormal", "mean": 8000, "sigma": 0.5, "max": 16000},
 "components": [
  {"name": "wild_type", "fraction": 0.65,
   "extra_edits": [{"position": 263, "deleted": 1, "inserted": "G"}]},
  {"name": "variant", "fraction": 0.20,
   "extra_edits": [{"position": 263, "deleted": 1, "inserted": "G"},
                   {"position": 11778, "deleted": 1, "inserted": "A"},
                   {"position": 6100, "deleted": 5}]},
  {"name": "deleted", "fraction": 0.15,
   "extra_edits": [{"position": 263, "deleted": 1, "inserted": "G"},
                   {"position": 8483, "deleted": 4977}]}
 ]}
JSON
python3 scripts/mtsim/mtsim.py --config "$WORK/scenario.json" --output "$WORK" --threads 2 >/dev/null

analyze() {
  cargo run --quiet -p mito-cli --offline -- analyze --input "$WORK/verify.bam" --json \
    --compact-min-molecules 50 --threads "$1"
}
analyze 1 > "$WORK/one.json"
analyze 4 > "$WORK/four.json"

# Thread count must not change a single call.
for file in one four; do
  jq -S 'del(.metadata.threads, .metadata.requested_threads)' "$WORK/$file.json" > "$WORK/$file.canonical.json"
done
cmp "$WORK/one.canonical.json" "$WORK/four.canonical.json" \
  || { echo "error: results differ between 1 and 4 threads" >&2; exit 1; }

jq -e '
  .precision_calls as $calls
  | ($calls.snvs | map(select(.status == "PASS"))) as $snvs
  | ($calls.indels | map(select(.status == "PASS"))) as $indels
  | ($calls.deletions | map(select(.status == "PASS"))) as $deletions
  | ($snvs | map(select(.position == 263 and .alt == "G" and .fraction > 0.95)) | length == 1)
    and ($snvs | map(select(.position == 11778 and .alt == "A" and .fraction > 0.12 and .fraction < 0.30)) | length == 1)
    and ($snvs | map(select(.position == 11778)) | first | .effects | map(select(.hgvs_p == "p.Arg340His")) | length == 1)
    and ($snvs | map(select(.position == 11778)) | first | .effects | map(select(.structure.structure_id == "9I4I" and .structure.chain == "r")) | length == 1)
    and ($indels | map(select(.hgvs == "m.6101_6105del" and .effects[0].consequence == "frameshift_variant")) | length == 1)
    and ($deletions | length == 1)
    and ($deletions[0] | .start == 8483 and .end == 13459 and .microhomology == "ACCTCCCTCACCA" and .known_common_deletion)
    and ($deletions[0].impact.fusion | startswith("MT-ATP8::MT-ND5"))
    and (.svs | map(select(.id == "deletion:8483-13459")) | length == 1)
    and ($calls.summary.snvs_pass <= 4)
' "$WORK/one.json" >/dev/null || { echo "error: statistical calls do not match the simulated truth" >&2; jq -c '.precision_calls.summary' "$WORK/one.json" >&2; exit 1; }

# A fixed variant inside the direct repeat (8473C, as in haplogroup D4) limits how
# far the junction can shift, but must not move the event or change its length.
cat > "$WORK/repeat.json" <<'JSON'
{"id": "repeat", "platform": "hifi", "coverage": 200, "seed": 11,
 "read_length": {"model": "lognormal", "mean": 8000, "sigma": 0.4, "max": 16000},
 "components": [
  {"name": "wild_type", "fraction": 0.7, "extra_edits": [{"position": 8473, "deleted": 1, "inserted": "C"}]},
  {"name": "deleted", "fraction": 0.3,
   "extra_edits": [{"position": 8473, "deleted": 1, "inserted": "C"}, {"position": 8483, "deleted": 4977}]}
 ]}
JSON
python3 scripts/mtsim/mtsim.py --config "$WORK/repeat.json" --output "$WORK" --threads 2 >/dev/null
cargo run --quiet -p mito-cli --offline -- analyze --input "$WORK/repeat.bam" --json \
  --compact-min-molecules 50 --threads 2 > "$WORK/repeat-result.json"
jq -e '
  [.precision_calls.deletions[] | select(.status == "PASS")] as $deletions
  | ($deletions | length == 1)
    and ($deletions[0] | .start == 8483 and .end == 13459 and .length == 4977
         and .known_common_deletion and .microhomology == "CCCTCACCA")
' "$WORK/repeat-result.json" >/dev/null \
  || { echo "error: a variant inside the direct repeat moved the called deletion" >&2; \
       jq -c '[.precision_calls.deletions[] | {start, end, length, microhomology, status}]' "$WORK/repeat-result.json" >&2; exit 1; }

echo "Statistical caller matches the simulated truth and is independent of thread count"
