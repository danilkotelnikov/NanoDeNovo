#!/usr/bin/env python3
"""Repaired confident-structure CSV filter using O(1) set lookup and atomic writes."""
from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


def atomic_write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".tmp_", suffix=".csv", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp_name, path)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def process_dataset(base_dir: Path, dataset_number: str) -> dict:
    txt_path = base_dir / f"confident_proteins_{dataset_number}.txt"
    csv_path = base_dir / f"folded_sequences_{dataset_number}.csv"
    out_path = base_dir / f"confident_proteins_{dataset_number}.csv"

    if not txt_path.exists():
        return {"dataset": dataset_number, "status": "skipped_missing_txt", "matched": 0}
    if not csv_path.exists():
        return {"dataset": dataset_number, "status": "skipped_missing_csv", "matched": 0}

    with open(txt_path, "r", encoding="utf-8", errors="replace") as fh:
        names = {line.strip().removesuffix(".pdb") for line in fh if line.strip()}

    matched_rows: list[dict] = []
    seen_names: set[str] = set()
    with open(csv_path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            return {"dataset": dataset_number, "status": "empty_header", "matched": 0}
        for row in reader:
            name = (row.get("Name") or "").strip()
            if name in names and name not in seen_names:
                seen_names.add(name)
                matched_rows.append(row)

    atomic_write_csv(out_path, fieldnames, matched_rows)
    return {"dataset": dataset_number, "status": "ok", "matched": len(matched_rows), "out": str(out_path)}


def selftest() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "confident_proteins_20.txt").write_text("0001_top.1.20.pdb\n1000_top.1.20.pdb\n", encoding="utf-8")
        (root / "folded_sequences_20.csv").write_text(
            "Name,Sequence,pLDDT,pTM\n"
            "0001_top.1.20,QVQL,0.95,0.90\n"
            "0001_top.2.20,AAAA,0.80,0.70\n"
            "1000_top.1.20,EVQL,0.96,0.91\n",
            encoding="utf-8",
        )
        res = process_dataset(root, "20")
        assert res["status"] == "ok" and res["matched"] == 2, res
        out_lines = (root / "confident_proteins_20.csv").read_text(encoding="utf-8").splitlines()
        assert len(out_lines) == 3, out_lines
    print("pipeline parse_confident selftest: pass")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Filter folded sequence CSVs against confident PDB lists.")
    ap.add_argument("--base-dir", type=Path, default=Path("."), help="Directory holding CSV and TXT files")
    ap.add_argument("--datasets", nargs="+", default=["20", "30"], help="Dataset tags")
    ap.add_argument("--selftest", action="store_true", help="Run built-in check and exit")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    for ds in args.datasets:
        res = process_dataset(args.base_dir, str(ds))
        print(f"Dataset {ds}: {res['status']} (matched={res['matched']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
