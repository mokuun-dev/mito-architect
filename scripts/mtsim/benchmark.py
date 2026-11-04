#!/usr/bin/env python3
"""Benchmark the engine over a grid of simulated samples with known truth.

Runs mixtures of two real mitogenomes with spiked small indels and large
deletions at decreasing levels, on each simulated platform, and writes a JSON
report and a Markdown table. Results measure the caller against the simulator's
error model; they do not replace validation on real sequencing data.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import mtsim  # noqa: E402
import score as scoring  # noqa: E402

ROOT = HERE.parents[1]

# Small indels away from the long C-tracts: a mix of sizes and of contexts.
SPIKED_INDELS = [
    {"position": 1000, "deleted": 1},
    {"position": 4100, "deleted": 2},
    {"position": 6100, "deleted": 5},
    {"position": 7500, "inserted": "G"},
    {"position": 9500, "inserted": "TAC"},
    {"position": 11000, "deleted": 9},
    {"position": 12600, "deleted": 3},
    {"position": 15000, "inserted": "AG"},
]
LARGE_DELETIONS = [
    {"position": 8483, "deleted": 4977},  # The common deletion.
    {"position": 7700, "deleted": 7000},
    {"position": 10000, "deleted": 2500},
]


def mixture(platform: str, minor: float, coverage: int, seed: int, numt: float = 0.0,
            read_length: dict | None = None) -> dict:
    """Major haplotype plus a minor one carrying spiked indels and a deletion."""
    tag = f"{platform}_m{minor:g}_c{coverage}_s{seed}" + (f"_numt{numt:g}" if numt else "")
    deletion = LARGE_DELETIONS[seed % len(LARGE_DELETIONS)]
    config = {
        "id": tag.replace(".", "p"),
        "platform": platform,
        "coverage": coverage,
        "seed": seed,
        "read_length": read_length or {"model": "lognormal", "mean": 8000, "sigma": 0.5, "max": 16000},
        "components": [
            {"name": "major", "genome": "AF346963", "fraction": 1.0 - 2 * minor},
            {"name": "minor", "genome": "AF347015", "fraction": minor,
             "extra_edits": SPIKED_INDELS},
            {"name": "deleted", "genome": "AF346963", "fraction": minor,
             "extra_edits": [deletion]},
        ],
    }
    if numt:
        config["numt_fraction"] = numt
    return config


def grid(profile: str) -> list[dict]:
    configs = []
    if profile == "quick":
        for platform, minor, coverage in [("ont_r10", 0.1, 600), ("hifi", 0.05, 600)]:
            configs.append(mixture(platform, minor, coverage, seed=3))
        return configs
    for platform in ("hifi", "ont_r10", "ont_r9"):
        for minor in (0.3, 0.1, 0.05, 0.02, 0.01):
            configs.append(mixture(platform, minor, 1000, seed=1))
        for minor in (0.02, 0.01, 0.005):
            configs.append(mixture(platform, minor, 4000, seed=2))
        configs.append(mixture(platform, 0.1, 100, seed=4))
        configs.append(mixture(platform, 0.05, 1000, seed=5, numt=0.02))
    return configs


def run_engine(engine: list[str], bam: Path) -> tuple[dict, float]:
    started = time.time()
    command = [part.replace("{bam}", str(bam)) for part in engine]
    output = subprocess.run(command, check=True, capture_output=True).stdout
    return json.loads(output), time.time() - started


def raw_table(engine: list[str], bam: Path, directory: Path) -> dict:
    """Only the raw variant table, read through a filter to bound memory."""
    path = directory / "raw.json"
    command = [part.replace("{bam}", str(bam)) for part in engine]
    with path.open("wb") as handle:
        subprocess.run(command, check=True, stdout=handle, stderr=subprocess.DEVNULL)
    filtered = subprocess.run(
        ["jq", "-c", "{variants: [.variants[] | {position, alt, heteroplasmy, alt_depth}], clusters: (.clusters | length)}",
         str(path)], check=True, capture_output=True).stdout
    path.unlink()
    return json.loads(filtered)


def run_one(config: dict, directory: Path, genomes: Path, engine: list[str], raw_engine: list[str]) -> dict:
    sample_directory = directory / config["id"]
    truth = mtsim.run_scenario(config, sample_directory, genomes, threads=2)
    bam = sample_directory / f"{config['id']}.bam"
    result, seconds = run_engine(engine, bam)
    report = scoring.score(result, truth)
    # The raw table is the comparison baseline; at high depth it is too large to keep.
    raw = raw_table(raw_engine, bam, sample_directory) if config["coverage"] <= 1000 else None
    report["legacy_raw_1pct_3reads"] = scoring.legacy_score(raw, truth, 0.01, 3) if raw else None
    report["legacy_raw_2pct_5reads"] = scoring.legacy_score(raw, truth, 0.02, 5) if raw else None
    report["engine_seconds"] = round(seconds, 2)
    report["raw_variants"] = len(raw["variants"]) if raw else None
    report["raw_clusters"] = raw["clusters"] if raw else None
    report["clusters"] = len(result["clusters"])
    report["config"] = {key: config[key] for key in ("platform", "coverage", "seed")}
    report["config"]["minor"] = config["components"][1]["fraction"]
    report["config"]["numt"] = config.get("numt_fraction", 0.0)
    report["haplogroup"] = (result.get("precision_calls", {}).get("consensus_haplogroup") or {}).get("best")
    report["excluded_divergent"] = result.get("precision_calls", {}).get("diagnostics", {}).get(
        "molecules_excluded_divergent")
    shutil.rmtree(sample_directory)
    return report


def raw_f1(score: dict | None) -> str:
    return "-" if score is None else f"{score['f1']:.3f}"


def markdown(reports: list[dict]) -> str:
    lines = [
        "| Platform | Minor | Depth | NUMT-like | SNV F1 | SNV P/R | Indel F1 | Indel P/R | "
        "Deletion F1 | Deletion level error | Raw F1 (1%, 3 reads) | Raw F1 (2%, 5 reads) | Time, s |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for report in reports:
        c = report["config"]
        s = report["snv_easy"]
        i = report["indel_easy"]
        d = report["large_deletion"]
        lines.append(
            f"| {c['platform']} | {c['minor']:g} | {c['coverage']} | {c['numt']:g} | {s['f1']:.3f} | "
            f"{s['precision']:.3f}/{s['recall']:.3f} | {i['f1']:.3f} | {i['precision']:.3f}/{i['recall']:.3f} | "
            f"{d['f1']:.2f} | {d['fraction_mae'] if d['fraction_mae'] is not None else '-'} | "
            f"{raw_f1(report['legacy_raw_1pct_3reads'])} | {raw_f1(report['legacy_raw_2pct_5reads'])} | "
            f"{report['engine_seconds']} |")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--engine", required=True, help="command printing result JSON, with {bam} marking the input path")
    parser.add_argument("--raw-engine", help="command for the raw table baseline; defaults to --engine")
    parser.add_argument("--output", type=Path, required=True, help="directory for the reports")
    parser.add_argument("--work", type=Path, default=ROOT / ".data" / "benchmark-work",
                        help="scratch directory on disk for reads and alignments")
    parser.add_argument("--genomes", type=Path, default=ROOT / ".data" / "genbank")
    parser.add_argument("--profile", choices=["quick", "full"], default="quick")
    parser.add_argument("--jobs", type=int, default=3)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    engine = args.engine.split()
    raw_engine = (args.raw_engine or args.engine).split()
    configs = grid(args.profile)
    args.work.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        reports = list(pool.map(lambda c: run_one(c, args.work, args.genomes, engine, raw_engine), configs))
    (args.output / "benchmark.json").write_text(json.dumps(reports, indent=1) + "\n")
    (args.output / "benchmark.md").write_text(markdown(reports))
    print(markdown(reports))
    return 0


if __name__ == "__main__":
    sys.exit(main())
