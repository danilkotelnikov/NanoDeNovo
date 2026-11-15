#!/usr/bin/env python3
# version: 0.3.0-a43ea1-4 | date: 15.XI.2026
"""Score complex-prediction PDBs with anchored REMARK 250 parsing and ipTM gating.

Fixes Finding 12 visibility gaps:
- Anchors `Predicted lDDT-Ca score:` vs `Antigen Predicted lDDT-Ca score:` so
  antigen pLDDT never shadows nanobody pLDDT.
- Writes an atomic summary CSV reporting every predicted complex plus explicit
  `iptm_pass` flags at `--min-iptm` (default 0.60, reported beside >=0.7 and >=0.8
  shares so low-ipTM poses are never passed silently to simulation).

Usage:
    py -3 -I fold_ag.py --complex-dir ./tfold_ag_20 \
        --out-csv results\\complexes_20.csv --min-iptm 0.60
    py -3 -I fold_ag.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import statistics as st
import sys
import tempfile
from pathlib import Path

PLDDT_RE = re.compile(r"^REMARK\s+250\s+Predicted lDDT-Ca score:\s*([0-9.]+)", re.M)
IPTM_RE = re.compile(r"^REMARK\s+250\s+Predicted ipTM score:\s*([0-9.]+)", re.M)
PTM_RE = re.compile(r"^REMARK\s+250\s+Predicted pTM score:\s*([0-9.]+)", re.M)
AG_PLDDT_RE = re.compile(r"^REMARK\s+250\s+Antigen Predicted lDDT-Ca score:\s*([0-9.]+)", re.M)

FIELDS = ["name", "plddt", "iptm", "ptm", "antigen_plddt", "iptm_gate", "pass_iptm"]


def parse_complex_remarks(pdb_path: Path) -> dict | None:
    """Parse all four REMARK 250 metrics from one complex PDB header."""
    head: list[str] = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM"):
                break
            head.append(line)
    text = "".join(head)
    plddt = PLDDT_RE.search(text)
    iptm = IPTM_RE.search(text)
    ptm = PTM_RE.search(text)
    ag_plddt = AG_PLDDT_RE.search(text)
    if not (plddt and iptm and ptm):
        return None
    return {
        "plddt": float(plddt.group(1)),
        "iptm": float(iptm.group(1)),
        "ptm": float(ptm.group(1)),
        "antigen_plddt": float(ag_plddt.group(1)) if ag_plddt else None,
    }


def atomic_write(path: Path, rows: list[dict]) -> None:
    """Write complex metrics CSV atomically."""
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


def score_complexes(complex_dir: Path, min_iptm: float) -> tuple[list[dict], list[str]]:
    """Score all `*.pdb` files inside `complex_dir`."""
    rows: list[dict] = []
    unread: list[str] = []
    for pdb in sorted(complex_dir.glob("*.pdb")):
        parsed = parse_complex_remarks(pdb)
        if parsed is None:
            unread.append(pdb.name)
            continue
        passed = parsed["iptm"] >= min_iptm
        rows.append(
            {
                "name": pdb.stem,
                "plddt": f"{parsed['plddt']:.4f}",
                "iptm": f"{parsed['iptm']:.4f}",
                "ptm": f"{parsed['ptm']:.4f}",
                "antigen_plddt": f"{parsed['antigen_plddt']:.4f}" if parsed["antigen_plddt"] is not None else "",
                "iptm_gate": f"{min_iptm:.2f}",
                "pass_iptm": "true" if passed else "false",
            }
        )
    return rows, unread


def summarize(rows: list[dict]) -> dict[str, float | int]:
    """Return cohort ipTM distribution metrics."""
    if not rows:
        return {"n": 0, "mean_iptm": 0.0, "ge_0_6": 0, "ge_0_7": 0, "ge_0_8": 0}
    iptms = [float(r["iptm"]) for r in rows]
    return {
        "n": len(rows),
        "mean_iptm": st.mean(iptms),
        "ge_0_6": sum(v >= 0.6 for v in iptms),
        "ge_0_7": sum(v >= 0.7 for v in iptms),
        "ge_0_8": sum(v >= 0.8 for v in iptms),
    }


def run(args: argparse.Namespace) -> int:
    rows, unread = score_complexes(args.complex_dir, args.min_iptm)
    if not rows:
        print(f"No readable complex PDBs in {args.complex_dir}", file=sys.stderr)
        return 1
    atomic_write(args.out_csv, rows)
    stats = summarize(rows)
    pass_n = sum(1 for r in rows if r["pass_iptm"] == "true")
    print(f"Complexes read : {stats['n']} (unread: {len(unread)})")
    print(f"Mean ipTM      : {stats['mean_iptm']:.4f}")
    print(f"ipTM >= 0.60   : {stats['ge_0_6']} | >= 0.70: {stats['ge_0_7']} | >= 0.80: {stats['ge_0_8']}")
    print(f"Pass @ {args.min_iptm:.2f}    : {pass_n}")
    if args.strict_gate and pass_n == 0:
        print("Strict ipTM gate blocked pipeline: 0 complexes passed.", file=sys.stderr)
        return 2
    return 0


def selftest() -> int:
    """Verify anchored REMARK 250 extraction, low-ipTM flagging, and strict gate behavior."""
    sample = (
        "REMARK 250\n"
        "REMARK 250 Predicted lDDT-Ca score: 0.9058\n"
        "REMARK 250 Structure predicted by tFold\n"
        "REMARK 250 Predicted ipTM score: 0.2355\n"
        "REMARK 250 Predicted pTM score: 0.6667\n"
        "REMARK 250 Antigen Predicted lDDT-Ca score: 0.9109\n"
        "ATOM      1  CA  SER H   1       0.000   0.000   0.000  1.00  0.90           C\n"
    )
    high = (
        "REMARK 250 Antigen Predicted lDDT-Ca score: 0.8800\n"
        "REMARK 250 Predicted pTM score: 0.8400\n"
        "REMARK 250 Predicted ipTM score: 0.7820\n"
        "REMARK 250 Predicted lDDT-Ca score: 0.9120\n"
        "ATOM      1  CA  SER H   1       0.000   0.000   0.000  1.00  0.91           C\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "low.pdb").write_text(sample, encoding="utf-8")
        (root / "high.pdb").write_text(high, encoding="utf-8")
        rows, unread = score_complexes(root, 0.60)
        assert not unread and len(rows) == 2
        by_name = {r["name"]: r for r in rows}
        assert by_name["low"]["plddt"] == "0.9058", by_name["low"]
        assert by_name["low"]["antigen_plddt"] == "0.9109", by_name["low"]
        assert by_name["low"]["pass_iptm"] == "false"
        assert by_name["high"]["pass_iptm"] == "true"
        out_csv = root / "complexes.csv"
        atomic_write(out_csv, rows)
        assert len(out_csv.read_text(encoding="utf-8").splitlines()) == 3

    print("fold_ag selftest: pass")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Score antibody-antigen complex predictions and ipTM gates.")
    ap.add_argument("--complex-dir", type=Path, help="Folder of complex PDB predictions")
    ap.add_argument("--out-csv", type=Path, help="Output CSV path")
    ap.add_argument("--min-iptm", type=float, default=0.60, help="Advisory/gate ipTM cutoff")
    ap.add_argument("--strict-gate", action="store_true", help="Exit 2 when zero complexes pass --min-iptm")
    ap.add_argument("--selftest", action="store_true", help="Run built-in check and exit")
    args = ap.parse_args()

    if args.selftest:
        return selftest()
    if args.complex_dir is None or args.out_csv is None:
        ap.error("--complex-dir and --out-csv are required")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
