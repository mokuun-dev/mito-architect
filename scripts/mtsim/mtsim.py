#!/usr/bin/env python3
"""Deterministic long-read simulator for human mtDNA with a known truth set.

This is engineering test infrastructure, not a sequencing model of record. The
error model reproduces the failure modes a mtDNA caller has to survive: length
errors that grow with homopolymer length, context-dependent substitutions that
differ between strands, a few recurrent systematic error motifs, per-read
quality variation, reads crossing the rCRS origin, molecules with large
deletions, and rare NUMT-like reads. Rates are order-of-magnitude settings for
each platform; benchmarks on these reads show relative behaviour and regressions
and do not establish accuracy on real sequencers.

Standard library only, so the verification gate has no Python dependencies.
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import json
import math
import random
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REFERENCE = ROOT / "core" / "data" / "rcrs.fasta"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
BASES = "ACGT"


@dataclass(frozen=True)
class Platform:
    """Error-model settings for one sequencing platform."""

    name: str
    minimap2_preset: str
    substitution: float
    insertion: float
    deletion: float
    # P(length error in a run of h) = hp_base * exp(hp_growth * (h - 3)), capped.
    hp_base: float
    hp_growth: float
    hp_cap: float
    hp_deletion_share: float
    context_sigma: float
    hotspot_fraction: float
    hotspot_multiplier: float
    read_quality_mean: float
    read_quality_sd: float
    base_quality_sd: float
    error_quality_drop: float


PLATFORMS = {
    # Nanopore R10-class duplex-free reads: roughly Q20 with dominant length
    # errors in homopolymers.
    "ont_r10": Platform(
        name="ont_r10",
        minimap2_preset="map-ont",
        substitution=0.0035,
        insertion=0.0012,
        deletion=0.0020,
        hp_base=0.018,
        hp_growth=0.42,
        hp_cap=0.45,
        hp_deletion_share=0.70,
        context_sigma=0.55,
        hotspot_fraction=0.006,
        hotspot_multiplier=7.0,
        read_quality_mean=20.0,
        read_quality_sd=2.5,
        base_quality_sd=5.0,
        error_quality_drop=9.0,
    ),
    # Nanopore R9-class reads: noisier, as in older public mtDNA datasets.
    "ont_r9": Platform(
        name="ont_r9",
        minimap2_preset="map-ont",
        substitution=0.012,
        insertion=0.006,
        deletion=0.010,
        hp_base=0.08,
        hp_growth=0.35,
        hp_cap=0.60,
        hp_deletion_share=0.75,
        context_sigma=0.70,
        hotspot_fraction=0.010,
        hotspot_multiplier=6.0,
        read_quality_mean=15.0,
        read_quality_sd=2.5,
        base_quality_sd=5.0,
        error_quality_drop=7.0,
    ),
    # PacBio HiFi: rare errors, mostly homopolymer length errors.
    "hifi": Platform(
        name="hifi",
        minimap2_preset="map-hifi",
        substitution=0.00012,
        insertion=0.00006,
        deletion=0.00008,
        hp_base=0.0012,
        hp_growth=0.55,
        hp_cap=0.20,
        hp_deletion_share=0.55,
        context_sigma=0.40,
        hotspot_fraction=0.002,
        hotspot_multiplier=5.0,
        read_quality_mean=38.0,
        read_quality_sd=3.0,
        base_quality_sd=6.0,
        error_quality_drop=16.0,
    ),
}


@dataclass(frozen=True, order=True)
class Edit:
    """One change relative to rCRS: replace deleted bases at 1-based pos."""

    position: int
    deleted: int
    inserted: str = ""

    def kind(self) -> str:
        if self.deleted == 1 and len(self.inserted) == 1:
            return "snv"
        if self.deleted and not self.inserted:
            return "deletion"
        if self.inserted and not self.deleted:
            return "insertion"
        return "complex"


@dataclass
class Component:
    """A molecular subpopulation: a sequence, its share of molecules, its edits."""

    name: str
    fraction: float
    edits: list[Edit] = field(default_factory=list)
    numt_like: bool = False


def load_fasta(path: Path) -> str:
    lines = [line.strip() for line in path.read_text().splitlines()]
    return "".join(line for line in lines if line and not line.startswith(">")).upper()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def apply_edits(reference: str, edits: list[Edit]) -> str:
    """Build a molecule sequence; rCRS placeholder N bases are always dropped."""
    pieces: list[str] = []
    cursor = 0
    for edit in sorted(edits):
        start = edit.position - 1
        if start < cursor:
            raise ValueError(f"overlapping edits at rCRS position {edit.position}")
        pieces.append(reference[cursor:start])
        pieces.append(edit.inserted)
        cursor = start + edit.deleted
    pieces.append(reference[cursor:])
    return "".join(pieces).replace("N", "")


def parse_cs_long(cs: str, reference_start: int) -> list[Edit]:
    """Convert a minimap2 --cs=long tag into edits with rCRS coordinates."""
    edits: list[Edit] = []
    position = reference_start
    for match in re.finditer(r"=([A-Za-z]+)|\*([a-z])([a-z])|\+([a-z]+)|-([a-z]+)", cs):
        same, sub_ref, sub_alt, inserted, deleted = match.groups()
        if same:
            position += len(same)
        elif sub_ref:
            if sub_ref != "n" and sub_alt != "n":
                edits.append(Edit(position, 1, sub_alt.upper()))
            position += 1
        elif inserted:
            edits.append(Edit(position, 0, inserted.upper()))
        else:
            if set(deleted) != {"n"}:
                # A placeholder N inside a longer gap stays in the reference run.
                edits.append(Edit(position, len(deleted)))
            position += len(deleted)
    return edits


def genome_edits(genome_fasta: Path, reference_fasta: Path = DEFAULT_REFERENCE) -> list[Edit]:
    """Derive an edit list for a whole mitogenome by aligning it to rCRS."""
    result = subprocess.run(
        ["minimap2", "-c", "--cs=long", "-x", "asm20", "--secondary=no", "-t", "2",
         str(reference_fasta), str(genome_fasta)],
        check=True, capture_output=True, text=True,
    )
    alignments = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        cs = next((f[5:] for f in fields if f.startswith("cs:Z:")), None)
        if cs is None:
            continue
        # PAF reference start is 0-based.
        alignments.append((int(fields[8]) - int(fields[7]), int(fields[7]) + 1, fields[4], cs))
    if not alignments:
        raise RuntimeError(f"no alignment of {genome_fasta} to rCRS")
    alignments.sort(reverse=True)
    span, start, strand, cs = alignments[0]
    if strand != "+" or span < 15000:
        raise RuntimeError(f"{genome_fasta}: expected one forward alignment, got {span} bp {strand}")
    return parse_cs_long(cs, start)


def normalize_indel(reference: str, position: int, deleted: int, inserted: str) -> tuple[str, int, str]:
    """Left-normalize an indel on rCRS (Tan et al. 2015) for truth comparison.

    Returns (kind, leftmost 1-based position, allele). A deletion's position is
    its first deleted base; an insertion's position is the base before it.
    """
    length = len(reference)
    if deleted and not inserted:
        allele = reference[position - 1 : position - 1 + deleted]
        start = position
        for _ in range(length):
            before = reference[(start - 2) % length]
            if before != allele[-1]:
                break
            start = (start - 2) % length + 1
            allele = before + allele[:-1]
        return ("deletion", start, allele)
    if inserted and not deleted:
        allele = inserted
        anchor = position - 1 if position > 1 else length
        for _ in range(length):
            before = reference[anchor - 1]
            if before != allele[-1]:
                break
            anchor = anchor - 1 if anchor > 1 else length
            allele = before + allele[:-1]
        return ("insertion", anchor, allele)
    raise ValueError("normalize_indel needs a pure insertion or deletion")


def truth_records(reference: str, components: list[Component]) -> list[dict]:
    """Collapse component edits into per-event truth with molecule fractions."""
    fractions: dict[tuple, float] = {}
    payload: dict[tuple, dict] = {}
    total = sum(c.fraction for c in components if not c.numt_like) or 1.0
    for component in components:
        if component.numt_like:
            continue
        for edit in component.edits:
            kind = edit.kind()
            if kind == "snv":
                key = ("snv", edit.position, edit.inserted)
                record = {"kind": "snv", "position": edit.position,
                          "ref": reference[edit.position - 1], "alt": edit.inserted}
            elif kind in ("deletion", "insertion"):
                normalized_kind, position, allele = normalize_indel(
                    reference, edit.position, edit.deleted, edit.inserted)
                key = (normalized_kind, position, allele)
                record = {"kind": normalized_kind, "position": position, "allele": allele,
                          "length": len(allele)}
            else:
                key = ("complex", edit.position, edit.deleted, edit.inserted)
                record = {"kind": "complex", "position": edit.position,
                          "deleted": edit.deleted, "inserted": edit.inserted}
            fractions[key] = fractions.get(key, 0.0) + component.fraction / total
            payload[key] = record
    records = []
    for key, record in payload.items():
        record["fraction"] = round(fractions[key], 6)
        records.append(record)
    records.sort(key=lambda r: (r["position"], r["kind"]))
    return records


def _motif_multiplier(motif: str, platform: Platform) -> float:
    """Deterministic per-motif substitution multiplier (systematic errors)."""
    digest = hashlib.sha256(f"{platform.name}:{motif}".encode()).digest()
    uniform_a = (int.from_bytes(digest[:4], "big") + 0.5) / 2**32
    uniform_b = (int.from_bytes(digest[4:8], "big") + 0.5) / 2**32
    gaussian = math.sqrt(-2.0 * math.log(uniform_a)) * math.cos(2.0 * math.pi * uniform_b)
    multiplier = math.exp(platform.context_sigma * gaussian)
    hotspot = (int.from_bytes(digest[8:12], "big") + 0.5) / 2**32
    if hotspot < platform.hotspot_fraction:
        multiplier *= platform.hotspot_multiplier
    return multiplier


class Template:
    """One strand of one circular molecule, with position-wise error rates."""

    def __init__(self, sequence: str, platform: Platform, motif_cache: dict[str, float]):
        length = len(sequence)
        self.length = length
        self.sequence = sequence + sequence
        doubled = self.sequence
        # Substitution rate depends on the 5-mer read through the pore, which is
        # why the two strands of one genome position have different error rates.
        rates = [0.0] * length
        for index in range(length):
            motif = (sequence[(index - 2) % length] + sequence[(index - 1) % length] +
                     sequence[index] + doubled[index + 1] + doubled[index + 2])
            multiplier = motif_cache.get(motif)
            if multiplier is None:
                multiplier = _motif_multiplier(motif, platform)
                motif_cache[motif] = multiplier
            rates[index] = min(0.45, platform.substitution * multiplier)
        self.substitution = rates + rates
        # Remaining run length from each position; a run is rolled once at its
        # first base. Positions inside a run hold 0.
        runs = [0] * (2 * length)
        index = 0
        while index < 2 * length:
            end = index + 1
            while end < 2 * length and doubled[end] == doubled[index]:
                end += 1
            runs[index] = end - index
            index = end
        self.runs = runs


def hamming_prefix(left: str, right: str) -> int:
    """Mismatches over the shared prefix; a coarse read-accuracy measure."""
    return sum(a != b for a, b in zip(left, right))


def _hp_error(platform: Platform, run: int) -> float:
    if run < 3:
        return 0.0
    return min(platform.hp_cap, platform.hp_base * math.exp(platform.hp_growth * (run - 3)))


@functools.lru_cache(maxsize=None)
def _hp_error_table(platform: Platform) -> tuple[float, ...]:
    return tuple(_hp_error(platform, run) for run in range(64))


def _read_length(rng: random.Random, model: dict, molecule_length: int) -> int:
    kind = model.get("model", "lognormal")
    maximum = min(int(model.get("max", molecule_length)), molecule_length)
    if kind == "fixed":
        return min(int(model["length"]), maximum)
    if kind == "full":
        return molecule_length
    mean = float(model.get("mean", 8000))
    sigma = float(model.get("sigma", 0.5))
    drawn = rng.lognormvariate(math.log(mean) - sigma * sigma / 2.0, sigma)
    return max(int(model.get("min", 500)), min(maximum, int(drawn)))


def _simulate_read(template: Template, start: int, length: int, platform: Platform,
                   rng: random.Random) -> tuple[str, str]:
    sequence = template.sequence
    substitution_rates = template.substitution
    runs = template.runs
    uniform = rng.random
    gauss = rng.gauss
    insertion = platform.insertion
    deletion = platform.deletion
    read_quality = gauss(platform.read_quality_mean, platform.read_quality_sd)
    base_sd = platform.base_quality_sd
    drop = platform.error_quality_drop
    hp_share = platform.hp_deletion_share
    hp_rates = _hp_error_table(platform)
    bases: list[str] = []
    qualities: list[int] = []

    def emit(base: str, penalty: float) -> None:
        bases.append(base)
        quality = int(read_quality + gauss(0.0, base_sd) - penalty)
        qualities.append(2 if quality < 2 else (60 if quality > 60 else quality))

    position = start
    end = start + length
    while position < end:
        run = runs[position]
        if run >= 3 and uniform() < hp_rates[min(run, 63)]:
            if uniform() < hp_share:
                position += 1
                if position >= end:
                    break
            else:
                emit(sequence[position], drop)
        base = sequence[position]
        draw = uniform()
        rate = substitution_rates[position]
        if draw < rate:
            alternative = BASES[int(uniform() * 4.0) & 3]
            while alternative == base:
                alternative = BASES[int(uniform() * 4.0) & 3]
            emit(alternative, drop)
        elif draw < rate + deletion:
            pass
        elif draw < rate + deletion + insertion:
            emit(BASES[int(uniform() * 4.0) & 3], drop)
            emit(base, 0.0)
        else:
            emit(base, 3.0 if run >= 4 else 0.0)
        position += 1
    return "".join(bases), "".join(chr(q + 33) for q in qualities)


def _read_start(rng: random.Random, region: tuple[int, int] | None, is_forward: bool,
                length: int, molecule_length: int) -> int | None:
    """Start of a read on its template strand; a region models an amplicon panel.

    Region bounds are molecule coordinates, close to rCRS coordinates for the
    purpose of leaving part of the genome without reads.
    """
    if region is None:
        return rng.randrange(molecule_length)
    first, last = region
    span = last - first + 1
    if length > span:
        return None
    offset = rng.randrange(span - length + 1)
    forward_start = first - 1 + offset
    return forward_start if is_forward else molecule_length - forward_start - length


def simulate_fastq(components: list[Component], reference: str, platform: Platform,
                   coverage: float, read_length: dict, seed: int, sample: str,
                   output: Path, region: tuple[int, int] | None = None) -> dict:
    """Write reads to FASTQ and return realized per-component molecule counts."""
    rng = random.Random(seed)
    motif_cache: dict[str, float] = {}
    templates: list[tuple[Component, Template, Template]] = []
    total_fraction = sum(c.fraction for c in components)
    for component in components:
        sequence = apply_edits(reference, component.edits)
        templates.append((component,
                          Template(sequence, platform, motif_cache),
                          Template(reverse_complement(sequence), platform, motif_cache)))
    target_bases = coverage * (len(reference) - 1 if region is None else region[1] - region[0] + 1)
    weights = [c.fraction / total_fraction for c in components]
    bases_written = 0
    read_index = 0
    realized = {c.name: 0 for c in components}
    with output.open("w") as handle:
        while bases_written < target_bases:
            draw = rng.random()
            cumulative = 0.0
            chosen = len(templates) - 1
            for index, weight in enumerate(weights):
                cumulative += weight
                if draw < cumulative:
                    chosen = index
                    break
            component, forward, reverse = templates[chosen]
            is_forward = rng.random() < 0.5
            template = forward if is_forward else reverse
            length = _read_length(rng, read_length, template.length)
            start = _read_start(rng, region, is_forward, length, template.length)
            if start is None:
                continue
            sequence, qualities = _simulate_read(template, start, length, platform, rng)
            if len(sequence) < 100:
                continue
            read_index += 1
            realized[component.name] += 1
            bases_written += len(sequence)
            handle.write(f"@{sample}_{read_index:07d}\n{sequence}\n+\n{qualities}\n")
    return {"reads": read_index, "bases": bases_written, "component_reads": realized}


def align(fastq: Path, bam: Path, preset: str, threads: int = 4,
          reference: Path = DEFAULT_REFERENCE) -> None:
    """Map reads with the settings used for public data; write an indexed BAM."""
    minimap2 = subprocess.Popen(
        ["minimap2", "-a", "-x", preset, "-Y", "--MD", "--secondary=no", "-t", str(threads),
         str(reference), str(fastq)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run(["samtools", "sort", "-@", str(threads), "-O", "BAM", "-o", str(bam), "-"],
                   stdin=minimap2.stdout, check=True, stderr=subprocess.DEVNULL)
    minimap2.stdout.close()
    if minimap2.wait() != 0:
        raise RuntimeError("minimap2 failed")
    subprocess.run(["samtools", "index", str(bam)], check=True)


def numt_like_edits(reference: str, seed: int, divergence: float = 0.08,
                    start: int = 5000, end: int = 9000) -> list[Edit]:
    """A fixed diverged copy of one genome segment, standing in for a NUMT.

    Every read of a NUMT carries the same differences, so at low abundance they
    look like linked low-level heteroplasmy, which is the artifact to reject.
    """
    rng = random.Random(seed)
    edits = []
    for position in range(start, end):
        base = reference[position - 1]
        if base in BASES and rng.random() < divergence:
            alternative = rng.choice([b for b in BASES if b != base])
            edits.append(Edit(position, 1, alternative))
    # Truncate to the segment by deleting the rest of the molecule's sequence.
    edits.append(Edit(end, len(reference) - end))
    edits.append(Edit(1, start - 1))
    return edits


def edit_from_json(item: dict) -> Edit:
    return Edit(int(item["position"]), int(item.get("deleted", 0)), str(item.get("inserted", "")))


def build_components(config: dict, reference: str, genome_directory: Path) -> list[Component]:
    components = []
    for item in config["components"]:
        edits: list[Edit] = []
        genome = item.get("genome")
        if genome and genome != "rCRS":
            edits.extend(genome_edits(genome_directory / f"{genome}.fasta"))
        for extra in item.get("extra_edits", []):
            edit = edit_from_json(extra)
            # A spiked change replaces any background change it overlaps.
            span_end = edit.position + max(edit.deleted, 1)
            edits = [e for e in edits
                     if e.position + max(e.deleted, 1) <= edit.position or e.position >= span_end]
            edits.append(edit)
        components.append(Component(item["name"], float(item["fraction"]), sorted(edits)))
    numt_fraction = float(config.get("numt_fraction", 0.0))
    if numt_fraction > 0.0:
        components.append(Component(
            "divergent_copy", numt_fraction,
            numt_like_edits(reference, int(config.get("seed", 1)) + 97), numt_like=True))
    return components


def run_scenario(config: dict, output_directory: Path, genome_directory: Path,
                 reference_path: Path = DEFAULT_REFERENCE, threads: int = 4) -> dict:
    """Simulate and align one scenario; write reads, BAM and truth JSON."""
    output_directory.mkdir(parents=True, exist_ok=True)
    reference = load_fasta(reference_path)
    platform = PLATFORMS[config["platform"]]
    sample = config["id"]
    components = build_components(config, reference, genome_directory)
    fastq = output_directory / f"{sample}.fastq"
    stats = simulate_fastq(components, reference, platform, float(config["coverage"]),
                           config.get("read_length", {"model": "lognormal", "mean": 8000, "sigma": 0.5}),
                           int(config.get("seed", 1)), sample, fastq,
                           tuple(config["region"]) if "region" in config else None)
    bam = output_directory / f"{sample}.bam"
    align(fastq, bam, platform.minimap2_preset, threads, reference_path)
    truth = {
        "sample": sample,
        "platform": platform.name,
        "coverage": config["coverage"],
        "seed": config.get("seed", 1),
        "components": [{"name": c.name, "fraction": c.fraction, "edits": len(c.edits),
                        "numt_like": c.numt_like} for c in components],
        "variants": truth_records(reference, components),
        "reads": stats["reads"],
        "bases": stats["bases"],
        "component_reads": stats["component_reads"],
    }
    (output_directory / f"{sample}.truth.json").write_text(json.dumps(truth, indent=1) + "\n")
    return truth


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", type=Path, required=True, help="scenario JSON file")
    parser.add_argument("--output", type=Path, required=True, help="output directory")
    parser.add_argument("--genomes", type=Path, default=ROOT / ".data" / "genbank",
                        help="directory holding <accession>.fasta background genomes")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    truth = run_scenario(json.loads(args.config.read_text()), args.output, args.genomes,
                         threads=args.threads)
    print(json.dumps({key: truth[key] for key in ("sample", "reads", "bases", "component_reads")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
