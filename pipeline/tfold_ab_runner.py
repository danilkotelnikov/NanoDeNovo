#!/usr/bin/env python3
"""Repaired monomer refolding scorer & launcher shared by tfold_ab_20.py and tfold_ab_30.py.

Repairs applied:
- Inclusive range `start..end` (`1..1000`), eliminating the off-by-one drop of backbone `1000`.
- Unified folder naming `tfold_{tag}/pdbs` matching `process_{tag}.py`.
- Anchored `REMARK 250` score parsing regardless of header line positions.
- Atomic CSV rewrite keyed by `Name` so reruns never duplicate rows.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

PLDDT_RE = re.compile(r"^REMARK\s+250\s+Predicted lDDT-Ca score:\s*([0-9.]+)", re.M)
PTM_RE = re.compile(r"^REMARK\s+250\s+Predicted pTM score:\s*([0-9.]+)", re.M)
THREE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
FIELDS = ["Name", "Sequence", "pLDDT", "pTM"]


def read_fasta_sequence(path: Path) -> str:
    if not path.exists():
        return ""
    seq: list[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith(">"):
                seq.append(s)
    return "".join(seq)


def read_pdb_ca_sequence(path: Path) -> str:
    seq: list[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                seq.append(THREE.get(line[17:20].strip(), "X"))
    return "".join(seq)


def parse_remarks(pdb_path: Path) -> tuple[float | None, float | None]:
    head: list[str] = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
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


def atomic_write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".tmp_", suffix=".csv", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp_name, path)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def run_dataset(
    base_dir: Path,
    dataset: str,
    start: int = 1,
    end: int = 1000,
    top: int = 3,
    plddt_min: float = 0.90,
    ptm_min: float = 0.85,
    fold_script: Path = Path("projects/tfold_ab/predict.py"),
    python_bin: str = sys.executable,
    score_only: bool = True,
) -> dict:
    tfold_dir = base_dir / f"tfold_{dataset}"
    output_pdb_dir = tfold_dir / "pdbs"
    confident_pdb_dir = tfold_dir / "pdbs_confident"
    csv_path = base_dir / f"folded_sequences_{dataset}.csv"

    output_pdb_dir.mkdir(parents=True, exist_ok=True)
    confident_pdb_dir.mkdir(parents=True, exist_ok=True)

    rows_by_name: dict[str, dict] = {}
    confident_count = 0

    for idx in range(start, end + 1):
        file_num = f"{idx:04d}"
        for n in range(1, top + 1):
            stem = f"{file_num}_top.{n}.{dataset}"
            fa_path = tfold_dir / f"{stem}.fa"
            pdb_path = output_pdb_dir / f"{stem}.pdb"

            if not score_only and fa_path.exists() and not pdb_path.exists():
                cmd = [python_bin, str(fold_script), "--fasta", str(fa_path), "--output", str(pdb_path)]
                proc = subprocess.run(cmd, check=False, encoding="utf-8", errors="replace")
                if proc.returncode != 0:
                    continue

            if not pdb_path.exists():
                continue

            plddt, ptm = parse_remarks(pdb_path)
            if plddt is None or ptm is None:
                continue

            seq = read_fasta_sequence(fa_path) or read_pdb_ca_sequence(pdb_path)
            rows_by_name[stem] = {
                "Name": stem,
                "Sequence": seq,
                "pLDDT": f"{plddt:.4f}",
                "pTM": f"{ptm:.4f}",
            }
            if plddt > plddt_min and ptm > ptm_min:
                confident_count += 1
                shutil.copy2(pdb_path, confident_pdb_dir / f"{stem}.pdb")

    ordered_rows = [rows_by_name[k] for k in sorted(rows_by_name)]
    atomic_write_csv(csv_path, ordered_rows)
    return {"dataset": dataset, "scored": len(ordered_rows), "confident": confident_count, "csv": str(csv_path)}


def build_parser(default_dataset: str) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=f"Repaired tFold-Ab runner/scorer (dataset {default_dataset}).")
    ap.add_argument("--base-dir", type=Path, default=Path("production"), help="Working root containing tfold_<dataset>")
    ap.add_argument("--dataset", default=default_dataset, help="Dataset tag (20 or 30)")
    ap.add_argument("--start", type=int, default=1, help="First backbone index (inclusive)")
    ap.add_argument("--end", type=int, default=1000, help="Last backbone index (inclusive)")
    ap.add_argument("--top", type=int, default=3, help="Designs per backbone")
    ap.add_argument("--plddt-min", type=float, default=0.90, help="Confidence filter pLDDT minimum")
    ap.add_argument("--ptm-min", type=float, default=0.85, help="Confidence filter pTM minimum")
    ap.add_argument("--fold-script", type=Path, default=Path("projects/tfold_ab/predict.py"), help="tFold-Ab prediction script")
    ap.add_argument("--python-bin", default=sys.executable, help="Python executable for external folding")
    ap.add_argument("--score-only", action="store_true", help="Score existing PDB outputs without invoking folding")
    ap.add_argument("--selftest", action="store_true", help="Run built-in check and exit")
    return ap


def selftest(default_dataset: str = "20") -> int:
    pdb_txt = (
        "REMARK 250 Header noise line\n"
        "REMARK 250 Predicted pTM score: 0.8800\n"
        "REMARK 250 Predicted lDDT-Ca score: 0.9400\n"
        "ATOM      1  CA  ALA H   1       0.000   0.000   0.000  1.00  0.94           C\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        pdbs = root / f"tfold_{default_dataset}" / "pdbs"
        pdbs.mkdir(parents=True)
        for bid in (1, 1000):
            stem = f"{bid:04d}_top.1.{default_dataset}"
            (root / f"tfold_{default_dataset}" / f"{stem}.fa").write_text(">H\nQVQL\n", encoding="utf-8")
            (pdbs / f"{stem}.pdb").write_text(pdb_txt, encoding="utf-8")
        res1 = run_dataset(root, default_dataset, start=1, end=1000, score_only=True)
        res2 = run_dataset(root, default_dataset, start=1, end=1000, score_only=True)
        assert res1["scored"] == 2 and res2["scored"] == 2, (res1, res2)
        lines = (root / f"folded_sequences_{default_dataset}.csv").read_text(encoding="utf-8").splitlines()
        assert len(lines) == 3 and "1000_top.1" in lines[-1], lines
        assert (root / f"tfold_{default_dataset}" / "pdbs_confident" / f"1000_top.1.{default_dataset}.pdb").exists()
    print(f"pipeline tfold_ab_{default_dataset} selftest: pass")
    return 0


def run_cli(default_dataset: str) -> int:
    ap = build_parser(default_dataset)
    ns = ap.parse_args()
    if ns.selftest:
        return selftest(default_dataset)
    res = run_dataset(
        base_dir=ns.base_dir,
        dataset=ns.dataset,
        start=ns.start,
        end=ns.end,
        top=ns.top,
        plddt_min=ns.plddt_min,
        ptm_min=ns.ptm_min,
        fold_script=ns.fold_script,
        python_bin=ns.python_bin,
        score_only=ns.score_only,
    )
    print(f"tfold_ab_{ns.dataset}: scored={res['scored']} confident={res['confident']} -> {res['csv']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run_cli("20"))
