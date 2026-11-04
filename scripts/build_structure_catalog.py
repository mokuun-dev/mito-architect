#!/usr/bin/env python3
"""Build and verify the mtDNA-encoded protein to PDB chain catalog.

For each of the 13 mtDNA-encoded proteins this script translates rCRS with the
vertebrate mitochondrial code, downloads the experimental human structure as
mmCIF, and compares every modelled residue by author number. The output records
which residues a coordinate actually exists for, any residue where the model
differs from rCRS, and the numbering offset (which must be zero for a mapping
to be emitted). The engine and the 3D viewer rely on these checks; a variant
never highlights a residue that this script did not verify.

Usage:
    build_structure_catalog.py           regenerate core/data/mt_structure_catalog.tsv
    build_structure_catalog.py --check   fail if the bundled catalog is out of date
"""
from __future__ import annotations

import argparse
import shlex
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "core" / "data" / "mt_structure_catalog.tsv"
REFERENCE = ROOT / "core" / "data" / "rcrs.fasta"

# gene, start, end (rCRS 1-based, inclusive), strand, UniProt, PDB, author chain, complex
GENES = [
    ("MT-ND1", 3307, 4262, "+", "P03886", "9I4I", "s", "Complex I NADH dehydrogenase"),
    ("MT-ND2", 4470, 5511, "+", "P03891", "9I4I", "i", "Complex I NADH dehydrogenase"),
    ("MT-ND3", 10059, 10404, "+", "P03897", "9I4I", "j", "Complex I NADH dehydrogenase"),
    ("MT-ND4L", 10470, 10766, "+", "P03901", "9I4I", "k", "Complex I NADH dehydrogenase"),
    ("MT-ND4", 10760, 12137, "+", "P03905", "9I4I", "r", "Complex I NADH dehydrogenase"),
    ("MT-ND5", 12337, 14148, "+", "P03915", "9I4I", "l", "Complex I NADH dehydrogenase"),
    ("MT-ND6", 14149, 14673, "-", "P03923", "9I4I", "m", "Complex I NADH dehydrogenase"),
    ("MT-CYB", 14747, 15887, "+", "P00156", "9HZL", "J", "Complex III cytochrome bc1"),
    ("MT-CO1", 5904, 7445, "+", "P00395", "9I7U", "A", "Complex IV cytochrome c oxidase"),
    ("MT-CO2", 7586, 8269, "+", "P00403", "9I7U", "B", "Complex IV cytochrome c oxidase"),
    ("MT-CO3", 9207, 9990, "+", "P00414", "9I7U", "C", "Complex IV cytochrome c oxidase"),
    ("MT-ATP8", 8366, 8572, "+", "P03928", "8H9S", "Q", "Complex V ATP synthase"),
    ("MT-ATP6", 8527, 9207, "+", "P00846", "8H9S", "N", "Complex V ATP synthase"),
]

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "FME": "M",  # N-formylmethionine, the mitochondrial initiator residue.
}

# NCBI translation table 2, vertebrate mitochondrial code.
CODE = {}
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
for _i, _a in enumerate(_BASES):
    for _j, _b in enumerate(_BASES):
        for _k, _c in enumerate(_BASES):
            CODE[_a + _b + _c] = _AAS[16 * _i + 4 * _j + _k]


def load_reference() -> str:
    lines = REFERENCE.read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">")).upper()


def translate_gene(reference: str, start: int, end: int, strand: str) -> str:
    segment = reference[start - 1 : end]
    if strand == "-":
        segment = segment.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    protein = []
    for offset in range(0, len(segment), 3):
        # Incomplete stop codons are completed by polyadenylation.
        codon = (segment[offset : offset + 3] + "AA")[:3]
        amino = CODE[codon]
        if offset == 0:
            amino = "M"  # Mitochondrial alternative initiators encode Met.
        if amino == "*":
            break
        protein.append(amino)
    return "".join(protein)


def fetch_mmcif(pdb_id: str, cache: Path) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f"{pdb_id.lower()}.cif"
    if not path.exists() or path.stat().st_size == 0:
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        with urllib.request.urlopen(url, timeout=300) as response:
            path.write_bytes(response.read())
    return path


def read_chain_residues(path: Path, chain: str) -> dict[int, str]:
    """Author residue number to one-letter code for residues with a CA atom."""
    residues: dict[int, str] = {}
    columns: list[str] = []
    in_header = False
    with path.open() as handle:
        for line in handle:
            if line.startswith("loop_"):
                columns, in_header = [], True
                continue
            if in_header and line.startswith("_atom_site."):
                columns.append(line.strip().split(".", 1)[1])
                continue
            if in_header and columns:
                in_header = False
            if columns and line.startswith(("ATOM", "HETATM")):
                fields = shlex.split(line)
                record = dict(zip(columns, fields))
                if (record.get("auth_asym_id") == chain
                        and record.get("label_atom_id") == "CA"
                        and record.get("pdbx_PDB_model_num", "1") == "1"):
                    name = record["auth_comp_id"]
                    if name in THREE_TO_ONE:
                        residues[int(record["auth_seq_id"])] = THREE_TO_ONE[name]
            elif columns and line.startswith("#"):
                if residues or any(c == "auth_asym_id" for c in columns):
                    break
    return residues


def as_ranges(numbers: list[int]) -> str:
    ranges = []
    start = previous = None
    for number in sorted(numbers):
        if start is None:
            start = previous = number
        elif number == previous + 1:
            previous = number
        else:
            ranges.append(f"{start}-{previous}")
            start = previous = number
    if start is not None:
        ranges.append(f"{start}-{previous}")
    return ",".join(ranges)


def build(cache: Path) -> str:
    reference = load_reference()
    rows = ["\t".join(["gene", "uniprot", "pdb_id", "chain", "complex", "protein_length",
                       "modelled_ranges", "model_differs"])]
    for gene, start, end, strand, uniprot, pdb_id, chain, complex_name in GENES:
        protein = translate_gene(reference, start, end, strand)
        residues = read_chain_residues(fetch_mmcif(pdb_id, cache), chain)
        if not residues:
            raise SystemExit(f"{gene}: chain {chain} of {pdb_id} has no modelled residues")
        best_offset = max(
            range(-5, 6),
            key=lambda shift: sum(
                1 for number, code in residues.items()
                if 1 <= number + shift <= len(protein) and protein[number + shift - 1] == code),
        )
        if best_offset != 0:
            raise SystemExit(f"{gene}: {pdb_id}:{chain} numbering is offset by {best_offset}")
        differs = []
        for number, code in sorted(residues.items()):
            if number < 1 or number > len(protein):
                raise SystemExit(f"{gene}: {pdb_id}:{chain} residue {number} is outside the protein")
            if protein[number - 1] != code:
                differs.append(f"{protein[number - 1]}{number}{code}")
        identity = 1.0 - len(differs) / len(residues)
        if identity < 0.97:
            raise SystemExit(f"{gene}: only {identity:.3f} identity with {pdb_id}:{chain}")
        rows.append("\t".join([gene, uniprot, pdb_id, chain, complex_name, str(len(protein)),
                               as_ranges(list(residues)), ",".join(differs)]))
        print(f"{gene:8s} {pdb_id}:{chain:2s} length={len(protein):3d} "
              f"modelled={len(residues):3d} ({100 * len(residues) / len(protein):5.1f}%) "
              f"differs={differs}", file=sys.stderr)
    return "\n".join(rows) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="compare instead of writing")
    parser.add_argument("--cache", type=Path, default=ROOT / ".data" / "mmcif")
    args = parser.parse_args()
    catalog = build(args.cache)
    if args.check:
        if not CATALOG.exists() or CATALOG.read_text() != catalog:
            print(f"{CATALOG} is out of date; run {Path(__file__).name}", file=sys.stderr)
            return 1
        print("Structure catalog matches rCRS translation and the mmCIF models")
        return 0
    CATALOG.write_text(catalog)
    print(f"Wrote {CATALOG}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
