#!/usr/bin/env python3
"""Offline end-to-end truth test: 40 full rCRS molecules, 10 with two linked SNVs."""
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REFERENCE_SHA256 = "fc392cde8e63b4d2e3a870bb97cc0626dea33d46dfb8abdebffada040f42ec92"
EVENTS = {"snv:3:T:G", "snv:3243:A:G"}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def verify(result):
    require(result["metadata"]["schema_version"] == "0.6", "schema changed")
    require(result["metadata"]["reference_length"] == 16569, "full reference required")
    require(result["filter_stats"]["input_molecules"] == 40, "molecule identity lost")
    require(result["filter_stats"]["passed_reads"] == 40, "unexpected filtering")
    require(result["coverage_metrics"]["mean_depth"] == 40, "incorrect whole-genome coverage")
    require(all(bin["depth"] == 40 for bin in result["coverage"]), "coverage has gaps")
    require({v["event_id"] for v in result["variants"]} == EVENTS, "incorrect variant set")
    for variant in result["variants"]:
        require((variant["alt_depth"], variant["ref_depth"], variant["callable_depth"],
                 variant["heteroplasmy"]) == (10, 30, 40, 0.25), "incorrect allele counts")
        require(variant["numt_assessability"] == "NOT_ASSESSABLE", "mt-only NUMT overclaim")
        require(set(variant["supporting_molecule_ids"]) ==
                {f"molecule-{i:02d}" for i in range(10)}, "lost molecule traceability")
    require(len(result["phase_links"]) == 1, "incorrect phase graph")
    link = result["phase_links"][0]
    require({link["event_a_id"], link["event_b_id"]} == EVENTS, "incorrect phase endpoints")
    require(tuple(link[k] for k in ("jointly_callable", "both_alternate", "neither_alternate",
                                    "a_alternate_b_absent", "a_absent_b_alternate")) ==
            (40, 10, 30, 0, 0), "incorrect phase contingency table")
    architectures = {tuple(a["defining_event_signature"]): a for a in result["architectures"]}
    require(set(architectures) == {(), tuple(sorted(EVENTS))}, "incorrect architectures")
    for signature, count, fraction in [((), 30, 0.75), (tuple(sorted(EVENTS)), 10, 0.25)]:
        a = architectures[signature]
        require((a["molecule_count"], a["estimated_fraction"]) == (count, fraction),
                "incorrect architecture fraction")
    return {"molecules": 40, "reference_length": 16569, "depth": 40,
            "variants": sorted(EVENTS), "alternate_molecules": 10,
            "architecture_fractions": [0.75, 0.25]}


def biological_projection(result):
    """Order-independent scientific projection for the thread/order truth gate."""
    molecule_ids = {molecule["index"]: molecule["id"] for molecule in result["molecules"]}
    variants = []
    for variant in result["variants"]:
        variants.append({
            key: variant[key] for key in (
                "event_id", "ref", "alt", "alt_depth", "ref_depth", "callable_depth",
                "heteroplasmy", "numt_assessability", "supporting_molecule_ids"
            )
        })
    links = []
    for link in result["phase_links"]:
        links.append({
            key: link[key] for key in (
                "event_a_id", "event_b_id", "jointly_callable", "jointly_uncertain",
                "both_alternate", "a_alternate_b_absent", "a_absent_b_alternate",
                "neither_alternate", "co_alternate_fraction", "linkage_delta"
            )
        } | {
            "supporting_molecule_ids": sorted(
                molecule_ids[index] for index in link["supporting_molecule_indices"]),
            "uncertain_molecule_ids": sorted(
                molecule_ids[index] for index in link["uncertain_molecule_indices"]),
        })
    architectures = []
    for architecture in result["architectures"]:
        architectures.append({
            "defining_event_signature": architecture["defining_event_signature"],
            "optional_event_ids": architecture["optional_event_ids"],
            "molecule_count": architecture["molecule_count"],
            "estimated_fraction": architecture["estimated_fraction"],
            "member_molecule_ids": sorted(
                molecule_ids[index] for index in architecture["member_molecule_indices"]),
        })
    return {
        "coverage": result["coverage"],
        "coverage_metrics": result["coverage_metrics"],
        "variants": sorted(variants, key=lambda variant: variant["event_id"]),
        "phase_links": sorted(links, key=lambda link: (link["event_a_id"], link["event_b_id"])),
        "architectures": sorted(architectures,
                                key=lambda architecture: architecture["defining_event_signature"]),
    }


def main():
    reference = ROOT / "core/data/rcrs.fasta"
    require(hashlib.sha256(reference.read_bytes()).hexdigest() == REFERENCE_SHA256,
            "rCRS bytes changed; review fixture provenance")
    sequence = "".join(reference.read_text().splitlines()[1:])
    require(len(sequence) == 16569 and sequence[2] == "T" and sequence[3242] == "A",
            "unexpected rCRS coordinates")
    directory = ROOT / "fixtures/generated/full_rcrs"
    directory.mkdir(parents=True, exist_ok=True)
    records = []
    for i in range(40):
        bases = list(sequence)
        if i < 10:
            bases[2] = bases[3242] = "G"
        records.append(f"molecule-{i:02d}\t0\tchrM\t1\t60\t16569M\t*\t0\t0\t"
                       + "".join(bases) + "\t" + "I" * len(bases) + "\n")
    results = []
    projections = []
    for label, rows, threads in [("ordered", records, 1), ("reversed", records[::-1], 4)]:
        path = directory / f"{label}.sam"
        path.write_text("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrM\tLN:16569\n" + "".join(rows))
        output = directory / f"{label}.json"
        with output.open("w") as stream:
            subprocess.run([str(ROOT / "target/debug/mito-cli"), "analyze", "-i", str(path),
                            "--json", "--evidence-graph", "--threads", str(threads)],
                           cwd=ROOT, stdout=stream, check=True, timeout=120)
        result = json.loads(output.read_text())
        results.append(verify(result))
        projections.append(biological_projection(result))
    require(results[0] == results[1], "input order or worker count changed truth")
    require(projections[0] == projections[1],
            "input order or worker count changed biological projection")
    summary = {"status": "passed", "reference_sha256": REFERENCE_SHA256,
               "scope": "Synthetic engineering truth; no analytical sensitivity claim",
               "cases": ["ordered/1-thread", "reversed/4-threads"], **results[0]}
    (directory / "verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("Full rCRS mixture truth passed (40 molecules, linked 25% SNVs, 2 architectures)")


if __name__ == "__main__":
    main()
