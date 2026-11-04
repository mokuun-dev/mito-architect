#!/usr/bin/env python3
"""Generate the web demo panel from real genomes and the actual engine.

Every demo result is real engine output on simulated or public reads, so the
panel cannot drift from what the software does. Results are reduced for the web
by listing a subset of molecules per cluster; counts are left as computed.

    build_demo_panel.py --engine "target/release/mito-cli analyze --json --input {bam} ..."
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import mtsim  # noqa: E402
import score as scoring  # noqa: E402

ROOT = HERE.parents[1]
PANEL = HERE / "demo_panel.json"
OUTPUT = ROOT / "web" / "public" / "demo"
GENOMES = ROOT / ".data" / "genbank"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="

MAX_READS_PER_CLUSTER = 80
MAX_SUPPORT_LISTED = 40


def fetch_genome(accession: str) -> Path:
    path = GENOMES / f"{accession}.fasta"
    if not path.exists():
        GENOMES.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(EFETCH + accession, timeout=120) as response:
            path.write_bytes(response.read())
    return path


def reduce_result(result: dict, sample: dict) -> dict:
    """Lists fewer molecules; every count and statistic stays as the engine wrote it."""
    keep: set[str] = set()
    for cluster in result.get("clusters", []):
        keep.update(cluster.get("reads", [])[:MAX_READS_PER_CLUSTER])
    total_reads = len(result.get("reads", []))
    for variant in result.get("variants", []):
        kept = [read for read in variant.get("supporting_reads", []) if read in keep]
        variant["supporting_reads"] = kept[:MAX_SUPPORT_LISTED]
        if "supporting_molecule_ids" in variant:
            variant["supporting_molecule_ids"] = variant["supporting_reads"]
    for sv in result.get("svs", []):
        kept = [read for read in sv.get("supporting_reads", []) if read in keep]
        sv["supporting_read_count"] = len(sv.get("supporting_reads", []))
        sv["supporting_reads"] = kept[:MAX_SUPPORT_LISTED] or sv.get("supporting_reads", [])[:MAX_SUPPORT_LISTED]
    for cluster in result.get("clusters", []):
        cluster["reads"] = [read for read in cluster.get("reads", []) if read in keep][:MAX_READS_PER_CLUSTER]
    result["reads"] = [read for read in result.get("reads", []) if read["id"] in keep]
    for read in result["reads"]:
        read.pop("aux_tags", None)
        read.pop("haplogroup_callable_ranges", None)
    calls = result.get("precision_calls", {})
    for key in ("deletions", "indels"):
        for call in calls.get(key, []):
            call["supporting_molecules"] = call.get("supporting_molecules", [])[:MAX_SUPPORT_LISTED]
    # Filtered candidates stay inspectable, but the page does not need hundreds.
    for key in ("snvs", "indels"):
        items = calls.get(key, [])
        passed = [item for item in items if item["status"] == "PASS"]
        filtered = sorted((item for item in items if item["status"] != "PASS"),
                          key=lambda item: item["p_value"])
        calls[key] = sorted(passed + filtered[:40], key=lambda item: (item["position"], item["id"]))
    metadata = result["metadata"]
    # The build machine's paths have no place in a published result.
    metadata["input_path"] = f"{sample['id']}.bam"
    metadata["reference_path"] = "rcrs.fasta"
    metadata["sample"] = sample["id"]
    metadata["demo"] = {
        "listed_reads": len(result["reads"]),
        "total_reads": total_reads,
        "note": "Molecule lists are shortened for the web; counts and statistics are unchanged.",
    }
    for resource in metadata.get("resources", []):
        resource["path"] = Path(resource.get("path", "")).name
    return result


def build_sample(sample: dict, engine: list[str], work: Path, threads: int) -> tuple[dict, dict | None]:
    if "real" in sample:
        bam = ROOT / sample["real"]["bam"]
        truth = None
    else:
        config = {"id": sample["id"], **sample["scenario"]}
        for component in config["components"]:
            if component.get("genome") not in (None, "rCRS"):
                fetch_genome(component["genome"])
        truth = mtsim.run_scenario(config, work / sample["id"], GENOMES, threads=threads)
        bam = work / sample["id"] / f"{sample['id']}.bam"
    command = [part.replace("{bam}", str(bam)) for part in engine]
    output = subprocess.run(command, check=True, capture_output=True).stdout
    return json.loads(output), truth


def summarize(sample: dict, result: dict, truth: dict | None, payload: bytes) -> dict:
    calls = result["precision_calls"]
    entry = {
        "id": sample["id"],
        "title": sample["title"],
        "summary": sample["summary"],
        "shows": sample["shows"],
        "accessions": sample["accessions"],
        "simulated": truth is not None,
        "platform": sample.get("scenario", {}).get("platform", "ont_public"),
        "file": f"{sample['id']}.json",
        "sha256": hashlib.sha256(payload).hexdigest(),
        "haplogroup": (calls.get("consensus_haplogroup") or {}).get("best"),
        "haplogroup_quality": (calls.get("consensus_haplogroup") or {}).get("quality"),
        "haplotypes": len(calls.get("haplotypes", [])),
        "median_depth": calls["diagnostics"]["median_depth"],
        "calls": calls["summary"],
    }
    if truth is not None:
        report = scoring.score(result, truth)
        entry["truth_check"] = {key: report[key] for key in ("snv", "snv_easy", "indel", "indel_easy", "large_deletion")}
    else:
        entry["source"] = sample["real"]["source"]
    return entry


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--engine", required=True, help="command printing result JSON, with {bam} marking the input")
    parser.add_argument("--work", type=Path, required=True, help="scratch directory for reads and BAM files")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--only", help="rebuild only this sample ID")
    args = parser.parse_args()
    panel = json.loads(PANEL.read_text())
    engine = args.engine.split()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    index_path = OUTPUT / "index.json"
    existing = {}
    if args.only and index_path.exists():
        existing = {entry["id"]: entry for entry in json.loads(index_path.read_text())["samples"]}
    entries = []
    for sample in panel["samples"]:
        if args.only and sample["id"] != args.only:
            if sample["id"] in existing:
                entries.append(existing[sample["id"]])
            continue
        raw, truth = build_sample(sample, engine, args.work, args.threads)
        report_input = json.loads(json.dumps(raw))
        reduced = reduce_result(raw, sample)
        payload = json.dumps(reduced, separators=(",", ":"), sort_keys=False).encode()
        (OUTPUT / f"{sample['id']}.json").write_bytes(payload)
        entry = summarize(sample, report_input, truth, payload)
        entries.append(entry)
        print(f"{sample['id']}: {len(payload) / 1e3:.0f} KB, haplogroup {entry['haplogroup']}, "
              f"calls {entry['calls']}, haplotypes {entry['haplotypes']}", file=sys.stderr)
    index = {"version": panel["version"], "sources": panel["sources"], "samples": entries}
    index_path.write_text(json.dumps(index, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
