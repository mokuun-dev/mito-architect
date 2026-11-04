#!/usr/bin/env python3
"""Score an engine result against the truth file written by mtsim.py."""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


# Long C-tracts of the control region. A haplotype here has several equally good
# alignments, so a truth edit and a called edit can be two spellings of the same
# molecule. Calls are scored with and without them, as benchmarks stratify low
# complexity regions.
HARD_REGIONS = [(300, 317), (513, 526), (956, 966), (3105, 3110), (5895, 5900),
                (8270, 8290), (12418, 12426), (16180, 16195)]


def in_hard_region(position: int, length: int = 1) -> bool:
    return any(position <= high and position + length - 1 >= low for low, high in HARD_REGIONS)


def prf(tp: int, fp: int, fn: int) -> dict:
    precision = tp / (tp + fp) if tp + fp else 1.0
    recall = tp / (tp + fn) if tp + fn else 1.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {"tp": tp, "fp": fp, "fn": fn, "precision": round(precision, 4),
            "recall": round(recall, 4), "f1": round(f1, 4)}


def fraction_error(pairs: list[tuple[float, float]]) -> float | None:
    if not pairs:
        return None
    return round(sum(abs(called - true) for called, true in pairs) / len(pairs), 4)


def score_snvs(calls: dict[tuple, float], truth: dict[tuple, float], minimum: float = 0.0) -> dict:
    wanted = {key for key, value in truth.items() if value >= minimum}
    called = set(calls)
    # A call matching a truth variant below the evaluated stratum is not a false
    # positive for that stratum.
    fp = len(called - set(truth))
    result = prf(len(called & wanted), fp, len(wanted - called))
    result["fraction_mae"] = fraction_error([(calls[k], truth[k]) for k in called & wanted])
    return result


def score_indels(calls: dict[tuple, float], truth: dict[tuple, float]) -> dict:
    result = prf(len(set(calls) & set(truth)), len(set(calls) - set(truth)),
                 len(set(truth) - set(calls)))
    result["fraction_mae"] = fraction_error([(calls[k], truth[k]) for k in set(calls) & set(truth)])
    return result


def score_deletions(calls: list[dict], truth: list[dict], placement_tolerance: int = 15) -> dict:
    """The length must match exactly. Placement may differ within a direct repeat,
    since truth is normalized on rCRS and calls on the sample's own sequence."""
    matched_truth = set()
    matched_calls = set()
    pairs = []
    for t_index, record in enumerate(truth):
        for c_index, call in enumerate(calls):
            if c_index in matched_calls:
                continue
            same_length = call["length"] == record["length"]
            same_place = abs(call["left_aligned_start"] - record["position"]) <= placement_tolerance
            if same_length and same_place:
                matched_truth.add(t_index)
                matched_calls.add(c_index)
                if call.get("fraction") is not None:
                    pairs.append((call["fraction"], record["fraction"]))
                break
    result = prf(len(matched_truth), len(calls) - len(matched_calls),
                 len(truth) - len(matched_truth))
    result["fraction_mae"] = fraction_error(pairs)
    return result


def score(result: dict, truth: dict) -> dict:
    variants = truth["variants"]
    truth_snvs = {(v["position"], v["alt"]): v["fraction"] for v in variants if v["kind"] == "snv"}
    truth_small = {}
    truth_large = []
    for variant in variants:
        if variant["kind"] in ("deletion", "insertion"):
            if variant["length"] >= 50:
                truth_large.append(variant)
            else:
                truth_small[(variant["kind"], variant["position"], variant["allele"])] = variant["fraction"]
    precision_calls = result.get("precision_calls") or {}
    # A variant in a region with no reads is not missed: it is not callable.
    ranges = precision_calls.get("callable_ranges")
    if ranges is not None:
        callable_position = lambda position: any(low <= position <= high for low, high in ranges)
        truth_snvs = {k: v for k, v in truth_snvs.items() if callable_position(k[0])}
        truth_small = {k: v for k, v in truth_small.items() if callable_position(k[1])}
        truth_large = [v for v in truth_large
                       if callable_position(v["position"]) and callable_position(v["position"] + v["length"])]
    else:
        callable_position = lambda position: True
    # Calls and truth are compared over the same callable positions.
    snvs = {(c["position"], c["alt"]): c["fraction"]
            for c in precision_calls.get("snvs", []) if c["status"] == "PASS" and callable_position(c["position"])}
    indels = {(c["type"], c["position"], c["allele"]): c["fraction"]
              for c in precision_calls.get("indels", []) if c["status"] == "PASS" and callable_position(c["position"])}
    deletions = [c for c in precision_calls.get("deletions", []) if c["status"] == "PASS"]
    report = {
        "sample": truth["sample"],
        "snv": score_snvs(snvs, truth_snvs),
        "snv_easy": score_snvs({k: v for k, v in snvs.items() if not in_hard_region(k[0])},
                               {k: v for k, v in truth_snvs.items() if not in_hard_region(k[0])}),
        "indel": score_indels(indels, truth_small),
        "indel_easy": score_indels(
            {k: v for k, v in indels.items() if not in_hard_region(k[1], len(k[2]))},
            {k: v for k, v in truth_small.items() if not in_hard_region(k[1], len(k[2]))}),
        "large_deletion": score_deletions(deletions, truth_large),
        "strata": {},
    }
    for low, high, label in [(0.5, 1.01, ">=50%"), (0.1, 0.5, "10-50%"), (0.05, 0.1, "5-10%"),
                             (0.02, 0.05, "2-5%"), (0.0, 0.02, "<2%")]:
        keys = {k for k, v in truth_snvs.items() if low <= v < high}
        if keys:
            report["strata"][label] = {"truth": len(keys), "recalled": len(keys & set(snvs))}
    report["false_snvs"] = sorted(set(snvs) - set(truth_snvs))[:30]
    report["missed_snvs"] = sorted((k[0], k[1], truth_snvs[k]) for k in set(truth_snvs) - set(snvs))[:30]
    report["false_indels"] = sorted(set(indels) - set(truth_small))[:30]
    report["missed_indels"] = sorted(set(truth_small) - set(indels))[:30]
    return report


def legacy_score(result: dict, truth: dict, minimum_fraction: float, minimum_depth: int) -> dict:
    """The raw per-position table with a fixed threshold, as a comparison."""
    truth_snvs = {(v["position"], v["alt"]): v["fraction"] for v in truth["variants"] if v["kind"] == "snv"}
    called = {(v["position"], v["alt"]): v["heteroplasmy"] for v in result["variants"]
              if v.get("heteroplasmy", 0) >= minimum_fraction and v.get("alt_depth", 0) >= minimum_depth}
    # The filtered table is compared on sites outside long C-tracts, as the caller is.
    truth_snvs = {k: v for k, v in truth_snvs.items() if not in_hard_region(k[0])}
    called = {k: v for k, v in called.items() if not in_hard_region(k[0])}
    return score_snvs(called, truth_snvs)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--legacy", action="store_true", help="also score the raw variant table")
    args = parser.parse_args()
    result = json.loads(args.result.read_text())
    truth = json.loads(args.truth.read_text())
    report = score(result, truth)
    if args.legacy:
        report["legacy_raw_1pct_3reads"] = legacy_score(result, truth, 0.01, 3)
        report["legacy_raw_2pct_5reads"] = legacy_score(result, truth, 0.02, 5)
    json.dump(report, sys.stdout, indent=1)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
