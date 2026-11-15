#!/usr/bin/env python3
# version: 0.3.0-a43ea1-4 | date: 15.XI.2026
"""Static heavy-atom interface audit for nanobody-antigen complexes.

Computes deterministic structural geometry before docking or MD relaxation:
- Heavy-atom clashes (< 2.2 A) and minimum binder-antigen heavy-atom distance.
- Interface contacts (heavy atom <= 5.0 A; CA-CA <= 8.0 A).
- Polar H-bond candidate pairs (N/O across interface <= 3.5 A).
- Buried surface area (Shrake-Rupley SASA: 0.5 * (SASA_binder + SASA_ag - SASA_complex)).
- CDR3 loop placement (reusing cdr3_check.find_cdr3) and CDR3 contact share.

Usage:
    py -3 -I interface_check.py --inputs ./rosetta_dock\\*.pdb \
        --out w06-nanodenovo\\results\\interface_audit.csv
    py -3 -I interface_check.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

# Reuse CDR3 heuristic from sibling tool rather than duplicating it.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cdr3_check import THREE, find_cdr3  # noqa: E402

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "SE": 1.90}
FIELDS = [
    "name", "path", "status", "binder_chain", "antigen_chain",
    "binder_len", "antigen_len", "clashes_lt_2_2A", "min_dist_A",
    "heavy_contacts_5A", "ca_contacts_8A", "polar_hbonds_3_5A",
    "bsa_A2", "cdr3_len", "cdr3_contact_share",
]


def sphere_points(n: int = 92) -> np.ndarray:
    """Fibonacci golden-spiral unit sphere points."""
    idx = np.arange(0, n, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * idx / n)
    theta = math.pi * (1.0 + math.sqrt(5.0)) * idx
    return np.column_stack((np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)))


SPHERE = sphere_points(92)


def read_atoms(pdb_path: Path) -> list[dict]:
    """Read non-hydrogen ATOM records from one PDB file, keeping column-27 insertion codes."""
    atoms: list[dict] = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            altloc = line[16].strip() if len(line) > 16 else ""
            if altloc not in ("", "A", "1"):
                continue
            aname = line[12:16].strip()
            elem = (line[76:78].strip() or aname[0]).upper()
            if elem == "H" or aname.startswith("H"):
                continue
            chain = line[21].strip() or "_"
            resseq = int(line[22:26])
            icode = line[26].strip() if len(line) > 26 else ""
            resname = line[17:20].strip()
            atoms.append(
                {
                    "chain": chain,
                    "resseq": resseq,
                    "icode": icode,
                    "reskey": (chain, resseq, icode),
                    "resname": resname,
                    "aname": aname,
                    "elem": elem,
                    "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                }
            )
    return atoms


def pick_chains(atoms: list[dict], binder_opt: str | None, antigen_opt: str | None) -> tuple[str | None, str | None]:
    """Select binder (default H or shortest protein chain) and antigen chain."""
    by_chain: dict[str, set[tuple[str, int, str]]] = {}
    for a in atoms:
        by_chain.setdefault(a["chain"], set()).add(a["reskey"])
    if len(by_chain) < 2:
        return None, None
    b = binder_opt if binder_opt in by_chain else ("H" if "H" in by_chain else min(by_chain, key=lambda c: len(by_chain[c])))
    others = [c for c in by_chain if c != b]
    ag = antigen_opt if antigen_opt in by_chain and antigen_opt != b else ("T" if "T" in others else ("A" if "A" in others else max(others, key=lambda c: len(by_chain[c]))))
    return b, ag


def sasa_atoms(coords: np.ndarray, elems: list[str]) -> float:
    """Compute total Shrake-Rupley SASA (1.4 A solvent probe) for one atom coordinate array."""
    if len(coords) == 0:
        return 0.0
    radii = np.array([VDW.get(e, 1.70) + 1.40 for e in elems], dtype=float)
    total = 0.0
    # Chunked distance screening keeps memory small on ~3,000-atom complexes.
    for i in range(len(coords)):
        d = np.sqrt(np.sum((coords - coords[i]) ** 2, axis=1))
        r = radii[i]
        nb = np.where((d > 1e-6) & (d < r + radii + 1e-5))[0]
        if len(nb) == 0:
            total += 4.0 * math.pi * r * r
            continue
        pts = coords[i] + r * SPHERE
        d_nb = np.sqrt(np.sum((pts[:, None, :] - coords[nb][None, :, :]) ** 2, axis=-1))
        free = np.all(d_nb >= radii[nb][None, :] - 1e-6, axis=1)
        total += 4.0 * math.pi * r * r * float(np.mean(free))
    return total


def analyze_complex(pdb_path: Path, binder_opt: str | None = None, antigen_opt: str | None = None) -> dict:
    """Compute clash, contact, H-bond, BSA, and CDR3 metrics for one complex."""
    base = {k: "" for k in FIELDS}
    base["name"] = pdb_path.stem
    base["path"] = str(pdb_path)
    try:
        atoms = read_atoms(pdb_path)
    except OSError as exc:
        base["status"] = f"read_error:{type(exc).__name__}"
        return base

    b_ch, ag_ch = pick_chains(atoms, binder_opt, antigen_opt)
    if not b_ch or not ag_ch:
        base["status"] = "single_chain_skipped"
        return base

    b_atoms = [a for a in atoms if a["chain"] == b_ch]
    ag_atoms = [a for a in atoms if a["chain"] == ag_ch]

    b_res_order: list[tuple[str, int, str]] = []
    b_res_seq: list[str] = []
    seen_res: set[tuple[str, int, str]] = set()
    for a in b_atoms:
        if a["reskey"] not in seen_res:
            seen_res.add(a["reskey"])
            b_res_order.append(a["reskey"])
            b_res_seq.append(THREE.get(a["resname"], "X"))
    ag_res_set = {a["reskey"] for a in ag_atoms}

    b_xyz = np.array([a["xyz"] for a in b_atoms], dtype=float)
    ag_xyz = np.array([a["xyz"] for a in ag_atoms], dtype=float)
    dmat = np.sqrt(np.sum((b_xyz[:, None, :] - ag_xyz[None, :, :]) ** 2, axis=-1))

    clashes = int(np.sum(dmat < 2.2))
    min_dist = float(np.min(dmat)) if dmat.size else 99.9
    heavy_pairs = int(np.sum(dmat <= 5.0))

    b_polar = np.array([a["elem"] in ("N", "O") for a in b_atoms], dtype=bool)
    ag_polar = np.array([a["elem"] in ("N", "O") for a in ag_atoms], dtype=bool)
    polar_hbonds = int(np.sum((dmat <= 3.5) & b_polar[:, None] & ag_polar[None, :]))

    b_ca = np.array([a["xyz"] for a in b_atoms if a["aname"] == "CA"], dtype=float)
    ag_ca = np.array([a["xyz"] for a in ag_atoms if a["aname"] == "CA"], dtype=float)
    if len(b_ca) and len(ag_ca):
        ca_d = np.sqrt(np.sum((b_ca[:, None, :] - ag_ca[None, :, :]) ** 2, axis=-1))
        ca_contacts = int(np.sum(ca_d <= 8.0))
    else:
        ca_contacts = 0

    b_elems = [a["elem"] for a in b_atoms]
    ag_elems = [a["elem"] for a in ag_atoms]
    s_b = sasa_atoms(b_xyz, b_elems)
    s_ag = sasa_atoms(ag_xyz, ag_elems)
    s_both = sasa_atoms(np.vstack((b_xyz, ag_xyz)), b_elems + ag_elems)
    bsa = max(0.0, 0.5 * (s_b + s_ag - s_both))

    contact_b_res = {b_atoms[i]["reskey"] for i in np.where(np.any(dmat <= 5.0, axis=1))[0]}
    cdr3_span = find_cdr3("".join(b_res_seq))
    if cdr3_span and contact_b_res:
        c_start, c_len = cdr3_span
        cdr3_res = set(b_res_order[c_start:c_start + c_len])
        cdr3_share = len(contact_b_res & cdr3_res) / len(contact_b_res)
    elif cdr3_span:
        c_len = cdr3_span[1]
        cdr3_share = 0.0
    else:
        c_len = 0
        cdr3_share = 0.0

    return {
        "name": pdb_path.stem,
        "path": str(pdb_path),
        "status": "ok",
        "binder_chain": b_ch,
        "antigen_chain": ag_ch,
        "binder_len": len(b_res_order),
        "antigen_len": len(ag_res_set),
        "clashes_lt_2_2A": clashes,
        "min_dist_A": f"{min_dist:.2f}",
        "heavy_contacts_5A": heavy_pairs,
        "ca_contacts_8A": ca_contacts,
        "polar_hbonds_3_5A": polar_hbonds,
        "bsa_A2": f"{bsa:.1f}",
        "cdr3_len": c_len,
        "cdr3_contact_share": f"{cdr3_share:.3f}",
    }


def resolve_inputs(patterns: list[str]) -> list[Path]:
    """Expand globs safely, skipping unreadable/corrupted filesystem entries."""
    out: list[Path] = []
    seen: set[str] = set()
    for pat in patterns:
        matches = glob.glob(pat)
        for item in sorted(matches):
            p = Path(item)
            key = os.path.normcase(str(p))
            if key in seen or p.is_dir():
                continue
            seen.add(key)
            out.append(p)
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    """Write interface audit table atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".tmp_", suffix=".csv", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp_name, path)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def selftest() -> int:
    """Verify clash detection, H-bond contact count, positive BSA, and single-chain skip."""
    # Construct minimal VHH-like sequence (>=35 aa ending ...C...WGQG...) + 2-atom antigen near CDR3.
    seq_res = ["ALA"] * 18 + ["CYS"] + ["TYR"] * 6 + ["TRP", "GLY", "GLN", "GLY"] + ["SER"] * 7
    lines: list[str] = []
    atom_id = 1
    for idx, res in enumerate(seq_res, start=1):
        # Put CDR3 residue 22 near the origin; rest far away along X.
        x = 0.0 if idx == 22 else idx * 4.0 + 10.0
        lines.append(f"ATOM  {atom_id:5d}  CA  {res} H{idx:4d}    {x:8.3f}   0.000   0.000  1.00 50.00           C\n")
        atom_id += 1
        if idx == 22:
            lines.append(f"ATOM  {atom_id:5d}  OH  {res} H{idx:4d}     0.000   3.000   0.000  1.00 50.00           O\n")
            atom_id += 1
    # Antigen residue making one 3.0 A H-bond (no clash >=2.2 A) plus one clashing version.
    lines_ok = lines + [
        f"ATOM  {atom_id:5d}  N   LYS T   1     0.000   6.000   0.000  1.00 50.00           N\n",
        f"ATOM  {atom_id+1:5d}  CA  LYS T   1     0.000   7.400   0.000  1.00 50.00           C\n",
    ]
    lines_clash = lines + [
        f"ATOM  {atom_id:5d}  N   LYS T   1     0.000   4.500   0.000  1.00 50.00           N\n",
    ]

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        p_ok = root / "ok.pdb"
        p_cl = root / "clash.pdb"
        p_single = root / "single.pdb"
        p_ok.write_text("".join(lines_ok), encoding="utf-8")
        p_cl.write_text("".join(lines_clash), encoding="utf-8")
        p_single.write_text("".join(lines), encoding="utf-8")

        r_ok = analyze_complex(p_ok)
        assert r_ok["status"] == "ok" and r_ok["clashes_lt_2_2A"] == 0, r_ok
        assert r_ok["polar_hbonds_3_5A"] == 1, r_ok
        assert float(r_ok["bsa_A2"]) > 20.0, r_ok
        assert float(r_ok["cdr3_contact_share"]) == 1.0, r_ok

        # Verify column-27 insertion codes (e.g. 22A, 22B) are counted as distinct residues.
        lines_icode = list(lines_ok) + [
            f"ATOM  {atom_id+2:5d}  CA  TYR H  22A   80.000   0.000   0.000  1.00 50.00           C\n",
            f"ATOM  {atom_id+3:5d}  CA  SER H  22B   84.000   0.000   0.000  1.00 50.00           C\n",
        ]
        p_ic = root / "icode.pdb"
        p_ic.write_text("".join(lines_icode), encoding="utf-8")
        r_ic = analyze_complex(p_ic)
        assert r_ic["binder_len"] == r_ok["binder_len"] + 2, r_ic

        r_cl = analyze_complex(p_cl)
        assert r_cl["clashes_lt_2_2A"] >= 1, r_cl

        r_single = analyze_complex(p_single)
        assert r_single["status"] == "single_chain_skipped", r_single

    print("interface_check selftest: pass")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Static heavy-atom interface audit for nanobody complexes.")
    ap.add_argument("--inputs", nargs="+", help="Complex PDB paths or glob patterns")
    ap.add_argument("--out", type=Path, help="Output CSV table path")
    ap.add_argument("--binder-chain", default=None, help="Explicit binder chain ID (default auto: H)")
    ap.add_argument("--antigen-chain", default=None, help="Explicit antigen chain ID (default auto: T/A)")
    ap.add_argument("--selftest", action="store_true", help="Run built-in check and exit")
    args = ap.parse_args()

    if args.selftest:
        return selftest()
    if not args.inputs or args.out is None:
        ap.error("--inputs and --out are required")

    paths = resolve_inputs(args.inputs)
    if not paths:
        print("No complex PDB files matched.", file=sys.stderr)
        return 1

    rows = [analyze_complex(p, args.binder_chain, args.antigen_chain) for p in paths]
    write_csv(args.out, rows)
    ok_rows = [r for r in rows if r["status"] == "ok"]
    print(f"Audited files : {len(rows)} (complexes: {len(ok_rows)})")
    print(f"Out           : {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
