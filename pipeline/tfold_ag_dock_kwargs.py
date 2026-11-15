#!/usr/bin/env python3
"""Repaired complex prediction & REMARK 250 scorer for tFold-Ag.

Repairs applied relative to legacy `tfold_ag_dock_kwargs.py`:
- Fixes critical substring bug where `'Predicted lDDT-Ca score:' in line` shadowed
  `'Antigen Predicted lDDT-Ca score:'` and left `len(stats) == 3 != 4`.
  Uses anchored regexes from `tools/fold_ag.py` so line ordering never breaks parsing.
- Replaces external host paths (`../antigen/...`) with neutral `--base-dir` (`./production`).
- Writes summary CSVs atomically via temp file + `os.replace`.
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
IPTM_RE = re.compile(r"^REMARK\s+250\s+Predicted ipTM score:\s*([0-9.]+)", re.M)
PTM_RE = re.compile(r"^REMARK\s+250\s+Predicted pTM score:\s*([0-9.]+)", re.M)
AG_PLDDT_RE = re.compile(r"^REMARK\s+250\s+Antigen Predicted lDDT-Ca score:\s*([0-9.]+)", re.M)

CSV_FIELDS = [
    "Name",
    "Predicted lDDT-Ca score",
    "Predicted ipTM score",
    "Predicted pTM score",
    "Antigen Predicted lDDT-Ca score",
]


def parse_pdb_remarks(pdb_file: str | Path) -> dict | None:
    """Parse all 4 REMARK 250 complex metrics using anchored line regexes."""
    head: list[str] = []
    with open(pdb_file, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ATOM"):
                break
            head.append(line)
    text = "".join(head)
    plddt = PLDDT_RE.search(text)
    iptm = IPTM_RE.search(text)
    ptm = PTM_RE.search(text)
    ag_plddt = AG_PLDDT_RE.search(text)
    if not (plddt and iptm and ptm and ag_plddt):
        return None
    return {
        "lddt_ca": float(plddt.group(1)),
        "iptm": float(iptm.group(1)),
        "ptm": float(ptm.group(1)),
        "ag_lddt_ca": float(ag_plddt.group(1)),
    }


def atomic_write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".tmp_", suffix=".csv", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
            writer.writeheader()
            for entry in rows:
                writer.writerow({k: entry[k] for k in CSV_FIELDS})
        os.replace(tmp_name, path)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def run_dataset(
    base_dir: Path,
    out_dir: Path,
    dataset: str,
    ag_chain: str,
    iptm_min: float | None = 0.60,
    ptm_min: float | None = None,
    lddt_min: float | None = None,
    ag_lddt_min: float | None = None,
    predict_script: Path = Path("projects/tfold_ag/predict.py"),
    python_bin: str = sys.executable,
    score_only: bool = True,
) -> dict:
    complexes_file = out_dir / f"complexes_v_{dataset}.txt"
    tfold_ag_dir = base_dir / f"tfold_ag_{dataset}"
    tfold_ag_dir.mkdir(parents=True, exist_ok=True)

    entities: list[str] = []
    if complexes_file.exists():
        entities = [ln.strip() for ln in complexes_file.read_text(encoding="utf-8").splitlines() if ln.strip()]
    elif score_only:
        entities = [p.stem + ".fa" for p in sorted(tfold_ag_dir.glob("*.pdb"))]

    all_stats: list[dict] = []
    for entity in entities:
        stem = Path(entity).stem
        fasta_path = base_dir / "filtered_nanobodies" / f"v_{dataset}" / "complexes" / ag_chain / "fastas" / entity
        msa_path = base_dir / "antigens" / "mmseq" / ag_chain / f"{ag_chain}.a3m"
        output_pdb = tfold_ag_dir / f"{stem}.pdb"

        if not score_only and fasta_path.exists() and not output_pdb.exists():
            cmd = [python_bin, str(predict_script), "--fasta", str(fasta_path), "--msa", str(msa_path), "--output", str(output_pdb)]
            proc = subprocess.run(cmd, check=False, encoding="utf-8", errors="replace")
            if proc.returncode != 0:
                continue

        if not output_pdb.exists():
            continue
        stats = parse_pdb_remarks(output_pdb)
        if stats is None:
            continue
        all_stats.append(
            {
                "Name": stem,
                "Predicted lDDT-Ca score": stats["lddt_ca"],
                "Predicted ipTM score": stats["iptm"],
                "Predicted pTM score": stats["ptm"],
                "Antigen Predicted lDDT-Ca score": stats["ag_lddt_ca"],
                "pdb_path": output_pdb,
            }
        )

    all_csv = out_dir / f"all_complexes_stats_{dataset}.csv"
    atomic_write_csv(all_csv, all_stats)

    conf_dir = tfold_ag_dir / "confident_complexes"
    conf_dir.mkdir(parents=True, exist_ok=True)
    confident_stats: list[dict] = []
    for entry in all_stats:
        ok = True
        if lddt_min is not None and entry["Predicted lDDT-Ca score"] < lddt_min:
            ok = False
        if ptm_min is not None and entry["Predicted pTM score"] < ptm_min:
            ok = False
        if iptm_min is not None and entry["Predicted ipTM score"] < iptm_min:
            ok = False
        if ag_lddt_min is not None and entry["Antigen Predicted lDDT-Ca score"] < ag_lddt_min:
            ok = False
        if ok:
            confident_stats.append(entry)
            shutil.copy2(entry["pdb_path"], conf_dir / entry["pdb_path"].name)

    conf_csv = out_dir / f"confident_complexes_stats_{dataset}.csv"
    atomic_write_csv(conf_csv, confident_stats)
    return {"dataset": dataset, "scored": len(all_stats), "confident": len(confident_stats), "csv": str(all_csv)}


def selftest() -> int:
    # Verify anchored regex prevents 'Antigen Predicted lDDT-Ca score:' from overwriting nanobody lDDT-Ca
    pdb_reversed = (
        "REMARK 250 Antigen Predicted lDDT-Ca score: 0.9100\n"
        "REMARK 250 Predicted ipTM score: 0.7200\n"
        "REMARK 250 Predicted pTM score: 0.8100\n"
        "REMARK 250 Predicted lDDT-Ca score: 0.9350\n"
        "ATOM      1  CA  ALA H   1       0.000   0.000   0.000  1.00 0.93           C\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        pdb_dir = root / "tfold_ag_20"
        pdb_dir.mkdir(parents=True)
        (pdb_dir / "0001_top.1.20_8e8z_C_complex.pdb").write_text(pdb_reversed, encoding="utf-8")
        parsed = parse_pdb_remarks(pdb_dir / "0001_top.1.20_8e8z_C_complex.pdb")
        assert parsed is not None
        assert parsed["lddt_ca"] == 0.9350 and parsed["ag_lddt_ca"] == 0.9100, parsed
        assert parsed["iptm"] == 0.7200 and parsed["ptm"] == 0.8100, parsed
        res = run_dataset(root, root, "20", "8e8z_C", iptm_min=0.60, score_only=True)
        assert res["scored"] == 1 and res["confident"] == 1, res
    print("pipeline tfold_ag_dock_kwargs selftest: pass")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Repaired tFold-Ag complex docking scorer.")
    ap.add_argument("--base-dir", type=Path, default=Path("production"), help="Working root with tfold_ag_<dataset>")
    ap.add_argument("--out-dir", type=Path, default=Path("."), help="Output CSV/log directory")
    ap.add_argument("--datasets", nargs="+", default=["20", "30"], help="Datasets to process")
    ap.add_argument("--ag-chain", default="8e8z_C", help="Target antigen chain identifier")
    ap.add_argument("--lddt-ca", type=float, default=None)
    ap.add_argument("--ptm", type=float, default=None)
    ap.add_argument("--iptm", type=float, default=0.60)
    ap.add_argument("--ag-lddt-ca", type=float, default=None)
    ap.add_argument("--score-only", action="store_true", help="Score existing PDBs without invoking folding")
    ap.add_argument("--selftest", action="store_true", help="Run built-in check and exit")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    for ds in args.datasets:
        res = run_dataset(
            base_dir=args.base_dir,
            out_dir=args.out_dir,
            dataset=str(ds),
            ag_chain=args.ag_chain,
            iptm_min=args.iptm,
            ptm_min=args.ptm,
            lddt_min=args.lddt_ca,
            ag_lddt_min=args.ag_lddt_ca,
            score_only=args.score_only,
        )
        print(f"tfold_ag_{ds}: scored={res['scored']} confident={res['confident']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
