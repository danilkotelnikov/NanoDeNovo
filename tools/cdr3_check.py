#!/usr/bin/env python3
# version: w06 | date: 14.XI.2026
"""Measure how far the third loop (CDR3) of each folded design sits from its backbone.

A whole-chain match score is dominated by the conserved nanobody framework, so
most designs pass it. The decision that matters is the third loop, the long
flexible one that makes the antigen contact. This tool fits the framework alone
and then reports how far the loop moved.

Method, per design:
1. Read the C-alpha atoms and the residue names of the folded model and of its
   backbone. Residue i of one matches residue i of the other, because both were
   built on the same residue count.
2. Find the CDR3 window in the design sequence: the stretch between the last
   cysteine and the FR4 start (W or F, then GQG). The pattern is a heuristic.
   A design that does not fit it is written with an empty loop and counted.
3. Fit the model onto the backbone using every residue outside the window
   (Kabsch). Report `fw_rmsd` over those residues and `cdr3_rmsd` over the window.

Usage:
    py -3 -I cdr3_check.py --models DIR --scaffolds DIR --out results/cdr3_20.csv
    py -3 -I cdr3_check.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

THREE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
NAME_RE = re.compile(r"^(?P<bid>\d+)_top\.(?P<rank>\d+)\.(?P<tag>\w+)\.pdb$")
# Framework 4 is about 11 residues and ends the chain. Its start has the shape
# residue, G, any, G (WGQG, WGEG, FGQG, YGQG, RGQG, VGDG and so on).
FR4_SPAN = (14, 8)
FIELDS = ["name", "length", "cdr3_start", "cdr3_len", "fw_rmsd", "cdr3_rmsd"]


def read_ca(path: Path) -> tuple[str, np.ndarray]:
    """Return the one-letter sequence and the (n, 3) C-alpha coordinates."""
    letters, coords = [], []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                letters.append(THREE.get(line[17:20], "X"))
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return "".join(letters), np.array(coords, dtype=float)


def find_cdr3(sequence: str) -> tuple[int, int] | None:
    """Return (start, length) of the third loop, or None when it cannot be placed.

    The end of the loop is the framework 4 start: a G-any-G shape that lies
    8 to 14 residues from the chain end, the one nearest to 11. The start of the
    loop is the first cysteine within 30 residues before that, which is the
    conserved cysteine of framework 3. A heuristic, so the caller counts the misses.
    """
    n = len(sequence)
    candidates = [p for p in range(max(0, n - FR4_SPAN[0]), n - FR4_SPAN[1])
                  if p + 3 < n and sequence[p + 1] == "G" and sequence[p + 3] == "G"]
    if not candidates:
        return None
    end = min(candidates, key=lambda q: abs(q - (n - 11)))
    window = sequence[max(0, end - 31):max(0, end - 3)]
    offset = window.find("C")
    if offset < 0:
        return None
    start = max(0, end - 31) + offset + 1
    return start, end - start


def kabsch(moving: np.ndarray, fixed: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the rotation and the two centres that best map `moving` onto `fixed`."""
    mc, fc = moving.mean(axis=0), fixed.mean(axis=0)
    u, _, vt = np.linalg.svd((moving - mc).T @ (fixed - fc))
    sign = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, sign]) @ u.T
    return rot, np.stack([mc, fc])


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def measure(model_xyz: np.ndarray, scaffold_xyz: np.ndarray, start: int, length: int) -> tuple[float, float]:
    """Fit on the framework, then return (framework RMSD, CDR3 RMSD)."""
    inside = np.zeros(len(model_xyz), dtype=bool)
    inside[start:start + length] = True
    rot, centres = kabsch(model_xyz[~inside], scaffold_xyz[~inside])
    placed = (model_xyz - centres[0]) @ rot.T + centres[1]
    return rmsd(placed[~inside], scaffold_xyz[~inside]), rmsd(placed[inside], scaffold_xyz[inside])


def run(models: Path, scaffolds: Path, out: Path) -> int:
    rows, no_loop, mismatched = [], 0, []
    for model in sorted(models.glob("*.pdb")):
        match = NAME_RE.match(model.name)
        scaffold = scaffolds / f"gen_{match['bid']}.pdb" if match else None
        if match is None or not scaffold.exists():
            continue
        sequence, mxyz = read_ca(model)
        _, sxyz = read_ca(scaffold)
        if len(mxyz) != len(sxyz) or len(mxyz) == 0:
            mismatched.append(model.name)
            continue
        loop = find_cdr3(sequence)
        if loop is None:
            no_loop += 1
            rows.append({"name": model.stem, "length": len(mxyz), "cdr3_start": "",
                         "cdr3_len": "", "fw_rmsd": "", "cdr3_rmsd": ""})
            continue
        fw, cdr = measure(mxyz, sxyz, *loop)
        rows.append({"name": model.stem, "length": len(mxyz), "cdr3_start": loop[0],
                     "cdr3_len": loop[1], "fw_rmsd": round(fw, 3), "cdr3_rmsd": round(cdr, 3)})
    if not rows:
        print(f"No design with a backbone file in {models}", file=sys.stderr)
        return 1
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(out)
    print(f"Designs measured     : {len(rows) - no_loop}")
    print(f"No CDR3 motif found  : {no_loop}")
    print(f"Length mismatch      : {len(mismatched)}")
    print(f"Written              : {out}")
    return 0


def selftest() -> int:
    rng = np.random.default_rng(7)
    # The conserved cysteine opens the loop; the G-any-G of framework 4 closes it.
    seq = "QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADLRGGYSYWYWGQGTQVTVSS"
    loop = find_cdr3(seq)
    assert loop is not None and seq[loop[0]:loop[0] + loop[1]] == "AADLRGGYSYWY", loop
    assert find_cdr3("AAAAAAAAAA") is None
    # The framework 4 variants seen in the designs must all close the loop.
    for tail, loop_text in (("WGEGTLLKLATEK", "RQTKSSPWT"), ("YGQGSLITVEE", "GVLDGDGTILV"),
                            ("RGTGAEITVTE", "AVLSGPDEGKEIIV"), ("VGDGVELTVVE", "QQTATSPYTV")):
        full = "EDEATYYC" + loop_text + tail
        got = find_cdr3(full)
        assert got is not None, tail
        # The loop may include the residue that precedes G-any-G; it must start after the C.
        assert full[got[0] - 1] == "C" and abs(got[1] - len(loop_text)) <= 1, (tail, got)

    n = 40
    fixed = np.cumsum(rng.normal(size=(n, 3)) * 1.5, axis=0)
    # A copy rotated and moved, with the loop residues 20-29 pushed out by 4 angstrom.
    theta = 0.7
    rot = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    moved = fixed @ rot.T + np.array([10.0, -3.0, 5.0])
    kick = np.zeros_like(moved)
    kick[20:30] = (rot @ np.array([4.0, 0.0, 0.0]))
    moved = moved + kick
    fw, cdr = measure(moved, fixed, 20, 10)
    # The framework must fit to zero, and the loop must show the 4 angstrom kick.
    assert fw < 1e-6, fw
    assert abs(cdr - 4.0) < 1e-6, cdr
    # A pure rotation with no kick must give zero everywhere.
    fw0, cdr0 = measure(fixed @ rot.T, fixed, 20, 10)
    assert fw0 < 1e-6 and cdr0 < 1e-6, (fw0, cdr0)
    # A mirror image must not fit: the sign fix keeps the rotation proper.
    mirror = fixed * np.array([1.0, 1.0, -1.0])
    fwm, _ = measure(mirror, fixed, 20, 10)
    assert fwm > 0.5, fwm
    print("selftest: pass")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Measure the third loop of each folded design.")
    parser.add_argument("--models", type=Path, help="Folder with folded designs")
    parser.add_argument("--scaffolds", type=Path, help="Folder with gen_NNNN.pdb backbones")
    parser.add_argument("--out", type=Path, help="CSV file to write")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        return selftest()
    if not (args.models and args.scaffolds and args.out):
        parser.error("--models, --scaffolds and --out are required")
    return run(args.models, args.scaffolds, args.out)


if __name__ == "__main__":
    sys.exit(main())
