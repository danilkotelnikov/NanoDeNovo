#!/usr/bin/env python3
# version: 0.3.0-a43ea1-4 | date: 15.XI.2026
"""Score and manage single-chain nanobody folding outputs without silent row drift.

Fixes Findings 3, 5, 8, and 9 from AUDIT.md:
- Inclusive index range `start..end` (`1..1000`), never stopping at `999`.
- Consistent folder naming (`tfold_{tag}/pdbs`), matching `select_sequences.py`.
- Named `REMARK 250` regex parsing instead of fixed line-number indexing.
- Atomic full CSV rewrite, so re-running never duplicates existing rows.

Usage:
    py -3 -I fold_ab.py --work-dir ./tfold_20 --tag 20 \
        --out-csv w06-nanodenovo\\results\\folded_20.csv --score-only
    py -3 -I fold_ab.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

THREE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
PLDDT_RE = re.compile(r"^REMARK\s+250\s+Predicted lDDT-Ca score:\s*([0-9.]+)", re.M)
PTM_RE = re.compile(r"^REMARK\s+250\s+Predicted pTM score:\s*([0-9.]+)", re.M)
FIELDS = ["Name", "Backbone", "Rank", "Sequence", "pLDDT", "pTM", "Confident"]
NAME_RE = re.compile(r"^(?P<bid>\d+)_top\.(?P<rank>\d+)\.(?P<tag>\w+)$")


def read_ca_sequence(pdb_path: Path) -> str:
    """Extract one-letter amino acid sequence from ATOM CA records."""
    seq: list[str] = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                seq.append(THREE.get(line[17:20].strip(), "X"))
    return "".join(seq)


def parse_remarks(pdb_path: Path) -> tuple[float | None, float | None]:
    """Parse pLDDT and pTM from named REMARK 250 lines regardless of header order."""
    head: list[str] = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM"):
                break
            head.append(line)
    text = "".join(head)
    plddt_m = PLDDT_RE.search(text)
    ptm_m = PTM_RE.search(text)
    return (
        float(plddt_m.group(1)) if plddt_m else None,
        float(ptm_m.group(1)) if ptm_m else None,
    )


def atomic_write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    """Write a complete CSV atomically so re-runs replace rather than duplicate rows."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".tmp_", suffix=".csv", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp_name, path)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def iter_expected_stems(start: int, end: int, tag: str, top: int) -> list[str]:
    """Return inclusive `(start..end) x (1..top)` stem identifiers."""
    stems: list[str] = []
    for idx in range(start, end + 1):
        for rank in range(1, top + 1):
            stems.append(f"{idx:04d}_top.{rank}.{tag}")
    return stems


def score_folder(
    pdbs_dir: Path,
    tag: str,
    start: int,
    end: int,
    top: int,
    plddt_min: float,
    ptm_min: float,
) -> tuple[list[dict], list[str]]:
    """Collect deduplicated rows across PDBs in inclusive range `start..end`."""
    rows: list[dict] = []
    missing: list[str] = []
    seen: set[str] = set()
    for stem in iter_expected_stems(start, end, tag, top):
        pdb = pdbs_dir / f"{stem}.pdb"
        if not pdb.exists():
            missing.append(stem)
            continue
        if stem in seen:
            continue
        seen.add(stem)
        plddt, ptm = parse_remarks(pdb)
        if plddt is None or ptm is None:
            missing.append(f"{stem}:missing-remarks")
            continue
        match = NAME_RE.match(stem)
        bid = match.group("bid") if match else ""
        rank = match.group("rank") if match else ""
        confident = (plddt >= plddt_min) and (ptm >= ptm_min)
        rows.append(
            {
                "Name": stem,
                "Backbone": bid,
                "Rank": rank,
                "Sequence": read_ca_sequence(pdb),
                "pLDDT": f"{plddt:.4f}",
                "pTM": f"{ptm:.4f}",
                "Confident": "true" if confident else "false",
            }
        )
    return rows, missing


def run(args: argparse.Namespace) -> int:
    work_dir: Path = args.work_dir
    pdbs_dir = work_dir / "pdbs"
    pdbs_dir.mkdir(parents=True, exist_ok=True)

    if not args.score_only and args.fold_cmd:
        for stem in iter_expected_stems(args.start, args.end, args.tag, args.top):
            fasta = work_dir / f"{stem}.fa"
            out_pdb = pdbs_dir / f"{stem}.pdb"
            if not fasta.exists() or out_pdb.exists():
                continue
            cmd = args.fold_cmd.format(fasta=str(fasta), out=str(out_pdb), stem=stem)
            proc = subprocess.run(cmd, shell=True, check=False)
            if proc.returncode != 0:
                print(f"Folding failed on {stem}", file=sys.stderr)
                return proc.returncode

    rows, missing = score_folder(
        pdbs_dir, args.tag, args.start, args.end, args.top, args.plddt_min, args.ptm_min
    )
    atomic_write_csv(args.out_csv, FIELDS, rows)
    confident_n = sum(1 for r in rows if r["Confident"] == "true")
    print(f"Scored designs : {len(rows)}")
    print(f"Confident      : {confident_n}")
    print(f"Missing/unread : {len(missing)}")
    print(f"Table          : {args.out_csv}")
    return 0


def selftest() -> int:
    """Verify inclusive range 1..1000, unordered REMARK 250 parsing, and idempotent CSV overwrite."""
    pdb_reordered = (
        "REMARK 250 Structure predicted by open folder\n"
        "REMARK 250 Extra banner line\n"
        "REMARK 250 Predicted pTM score: 0.8712\n"
        "REMARK 250 Predicted lDDT-Ca score: 0.9145\n"
        "ATOM      1  CA  SER H   1       0.000   0.000   0.000  1.00  0.91           C\n"
        "ATOM      2  CA  ALA H   2       3.800   0.000   0.000  1.00  0.92           C\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        pdbs = root / "tfold_20" / "pdbs"
        pdbs.mkdir(parents=True)
        for bid in (1, 1000):
            (pdbs / f"{bid:04d}_top.1.20.pdb").write_text(pdb_reordered, encoding="utf-8")

        stems = iter_expected_stems(1, 1000, "20", 3)
        assert len(stems) == 3000, len(stems)
        assert stems[-1] == "1000_top.3.20", stems[-1]

        out_csv = root / "folded_20.csv"
        rows1, miss1 = score_folder(pdbs, "20", 1, 1000, 3, 0.90, 0.85)
        assert len(rows1) == 2 and len(miss1) == 2998
        assert rows1[0]["pLDDT"] == "0.9145" and rows1[0]["pTM"] == "0.8712", rows1[0]
        assert rows1[1]["Name"] == "1000_top.1.20"
        atomic_write_csv(out_csv, FIELDS, rows1)
        # Second write must keep exact row count (no doubling).
        atomic_write_csv(out_csv, FIELDS, rows1)
        lines = out_csv.read_text(encoding="utf-8").splitlines()
        assert len(lines) == 3, lines

    print("fold_ab selftest: pass")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Inclusive, atomic single-chain nanobody fold scorer.")
    ap.add_argument("--work-dir", type=Path, help="Folder containing `<work-dir>/pdbs` and FASTA inputs")
    ap.add_argument("--out-csv", type=Path, help="Output CSV table path (overwritten atomically)")
    ap.add_argument("--tag", default="20", help="Dataset tag (20 or 30)")
    ap.add_argument("--start", type=int, default=1, help="First backbone index (inclusive)")
    ap.add_argument("--end", type=int, default=1000, help="Last backbone index (inclusive)")
    ap.add_argument("--top", type=int, default=3, help="Ranks per backbone (inclusive 1..top)")
    ap.add_argument("--plddt-min", type=float, default=0.90, help="Advisory pLDDT threshold")
    ap.add_argument("--ptm-min", type=float, default=0.85, help="Advisory pTM threshold")
    ap.add_argument("--fold-cmd", default=None, help="Optional shell template with {fasta}, {out}, {stem}")
    ap.add_argument("--score-only", action="store_true", help="Do not invoke external folding; score existing PDBs")
    ap.add_argument("--selftest", action="store_true", help="Run built-in check and exit")
    args = ap.parse_args()

    if args.selftest:
        return selftest()
    if args.work_dir is None or args.out_csv is None:
        ap.error("--work-dir and --out-csv are required")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
