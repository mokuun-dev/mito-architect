#!/usr/bin/env python3
"""Check the published demo panel: integrity, privacy, contract, and call quality.

The panel is generated engine output, so a regression in calling or a stale file
shows up here. Quality floors are deliberately below current values: they catch
breakage without pretending that simulated scores are accuracy on real data.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "web" / "public" / "demo"
DEFINITION = ROOT / "scripts" / "mtsim" / "demo_panel.json"

# Paths on a build machine must not be published, nor must personal data appear.
FORBIDDEN = [re.compile(pattern) for pattern in (r"/home/", r"/Users/", r"C:\\\\Users", r"/tmp/", r"@[a-z0-9-]+\.[a-z]{2,}")]
MAXIMUM_SAMPLE_BYTES = 800_000

EXPECTED_CALLS = {
    "j1c-ont-lhon": {"snv": (11778, "A")},
    "b4a-hifi-narp": {"snv": (8993, "G")},
    "u5b-ont-melas-frameshift": {"snv": (3243, "G"), "hgvs": "m.9205_9206del"},
    "t2e-hifi-common-deletion": {"deletion": (8483, 13459)},
    # This lineage has 8473C inside the direct repeat of the common deletion.
    "d4a-ont-multiple-deletions": {"deletion": (8483, 13459)},
}

# Minimum recall and precision on simulated reads, outside long C-tracts.
QUALITY_FLOOR = 0.90
MAXIMUM_FALSE_CALLS = 2


def fail(message: str) -> None:
    print(f"demo panel: {message}", file=sys.stderr)
    raise SystemExit(1)


def check_sample(entry: dict) -> None:
    path = PANEL / entry["file"]
    if not path.exists():
        fail(f"{entry['file']} is listed in the index but missing")
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != entry["sha256"]:
        fail(f"{entry['file']} does not match the checksum in the index")
    if len(payload) > MAXIMUM_SAMPLE_BYTES:
        fail(f"{entry['file']} is {len(payload)} bytes, above the {MAXIMUM_SAMPLE_BYTES} limit for the web")
    text = payload.decode()
    for pattern in FORBIDDEN:
        match = pattern.search(text)
        if match:
            fail(f"{entry['file']} contains a forbidden string: {match.group(0)!r}")
    result = json.loads(text)
    for key in ("metadata", "filter_stats", "genes", "coverage", "svs", "variants", "clusters", "reads", "precision_calls"):
        if key not in result:
            fail(f"{entry['file']} lacks {key}")
    calls = result["precision_calls"]
    if not calls.get("assessed"):
        fail(f"{entry['file']}: the statistical caller did not run")
    # A listed read that is not in the table would break molecule links in the page.
    listed = {read["id"] for read in result["reads"]}
    for cluster in result["clusters"]:
        missing = [read for read in cluster["reads"] if read not in listed]
        if missing:
            fail(f"{entry['file']}: cluster {cluster['id']} lists unknown read {missing[0]}")

    expected = EXPECTED_CALLS.get(entry["id"], {})
    passed_snvs = {(call["position"], call["alt"]) for call in calls["snvs"] if call["status"] == "PASS"}
    if "snv" in expected and expected["snv"] not in passed_snvs:
        fail(f"{entry['id']}: expected SNV {expected['snv']} was not called")
    if "hgvs" in expected and not any(
            call["hgvs"] == expected["hgvs"] and call["status"] == "PASS" for call in calls["indels"]):
        fail(f"{entry['id']}: expected indel {expected['hgvs']} was not called")
    if "deletion" in expected and not any(
            (call["start"], call["end"]) == expected["deletion"] and call["status"] == "PASS"
            for call in calls["deletions"]):
        fail(f"{entry['id']}: expected deletion {expected['deletion']} was not called")

    truth = entry.get("truth_check")
    if entry["simulated"] and truth is None:
        fail(f"{entry['id']}: a simulated sample must carry its truth check")
    if truth:
        for key in ("snv_easy", "indel_easy", "large_deletion"):
            score = truth[key]
            if score["tp"] + score["fn"] > 0 and score["recall"] < QUALITY_FLOOR:
                fail(f"{entry['id']}: {key} recall {score['recall']} is below {QUALITY_FLOOR}")
            # With a handful of true events, one false call says little about precision.
            if score["fp"] > MAXIMUM_FALSE_CALLS and score["precision"] < QUALITY_FLOOR:
                fail(f"{entry['id']}: {key} precision {score['precision']} is below {QUALITY_FLOOR}")


def main() -> int:
    index_path = PANEL / "index.json"
    if not index_path.exists():
        fail("web/public/demo/index.json is missing; run scripts/mtsim/build_demo_panel.py")
    index = json.loads(index_path.read_text())
    defined = {sample["id"] for sample in json.loads(DEFINITION.read_text())["samples"]}
    listed = {entry["id"] for entry in index["samples"]}
    if defined != listed:
        fail(f"the index and scripts/mtsim/demo_panel.json differ: {sorted(defined ^ listed)}")
    for entry in index["samples"]:
        check_sample(entry)
    extra = {path.name for path in PANEL.glob("*.json")} - {"index.json"} - {e["file"] for e in index["samples"]}
    if extra:
        fail(f"files not listed in the index: {sorted(extra)}")
    if not any(not entry["simulated"] for entry in index["samples"]):
        fail("the panel must keep at least one sample of real reads")
    lineages = {(entry.get("haplogroup") or "?")[0] for entry in index["samples"]}
    print(f"Verified {len(index['samples'])} demo samples across haplogroup lineages {''.join(sorted(lineages))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
