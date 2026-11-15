#!/usr/bin/env python3
# version: w06 | date: 14.XI.2026
"""Align every folded design with its own backbone, kept or rejected.

The original run aligned only the designs that had passed the confidence
filter. That left the filter's misses unmeasurable. This tool aligns all of
them and records both match scores, each read by name from the USalign text.

USalign prints two match scores. One is normalised by the length of structure
1 and one by the length of structure 2. This tool puts the DESIGN first, so
`tm_design` is normalised by the design and `tm_scaffold` by the backbone.

Usage:
    py -3 -I align_all.py --usalign PATH --scaffolds DIR --models DIR \
        --out results/align_all_20.csv [--workers 4]
    py -3 -I align_all.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

NAME_RE = re.compile(r"^(?P<bid>\d+)_top\.(?P<rank>\d+)\.(?P<tag>\w+)\.pdb$")
RMSD_RE = re.compile(r"Aligned length=\s*(\d+),\s*RMSD=\s*([0-9.]+)")
TM_RE = re.compile(r"TM-score=\s*([0-9.]+)\s*\(normalized by length of Structure_([12])")
LEN_RE = re.compile(r"Length of Structure_([12]):\s*(\d+)")
FIELDS = ["name", "scaffold", "len_design", "len_scaffold", "aligned_length",
          "rmsd", "tm_design", "tm_scaffold"]


def parse(text: str) -> dict | None:
    """Read one USalign report. Return None when a field is missing."""
    rmsd = RMSD_RE.search(text)
    tms = {int(num): float(val) for val, num in TM_RE.findall(text)}
    lens = {int(num): int(val) for num, val in LEN_RE.findall(text)}
    if rmsd is None or set(tms) != {1, 2} or set(lens) != {1, 2}:
        return None
    return {
        "len_design": lens[1], "len_scaffold": lens[2],
        "aligned_length": int(rmsd.group(1)), "rmsd": float(rmsd.group(2)),
        "tm_design": tms[1], "tm_scaffold": tms[2],
    }


def align_one(usalign: str, model: Path, scaffold: Path) -> dict | None:
    done = subprocess.run([usalign, str(model), str(scaffold)],
                          capture_output=True, text=True, timeout=120)
    return parse(done.stdout)


def run(usalign: str, scaffolds: Path, models: Path, out: Path, workers: int) -> int:
    jobs, orphans = [], []
    for model in sorted(models.glob("*.pdb")):
        match = NAME_RE.match(model.name)
        if match is None:
            continue
        scaffold = scaffolds / f"gen_{match['bid']}.pdb"
        (jobs if scaffold.exists() else orphans).append((model, scaffold))

    def work(job):
        model, scaffold = job
        row = align_one(usalign, model, scaffold)
        return model, scaffold, row

    if not jobs:
        # An empty result must never look like a finished run.
        print(f"No design with a backbone file in {models}", file=sys.stderr)
        return 1

    rows, failed = [], []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for model, scaffold, row in pool.map(work, jobs):
            if row is None:
                failed.append(model.name)
                continue
            rows.append({"name": model.stem, "scaffold": scaffold.name, **row})

    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(out)  # a rerun replaces the file; it never appends to it

    print(f"Models found     : {len(jobs) + len(orphans)}")
    print(f"Aligned          : {len(rows)}")
    print(f"No backbone file : {len(orphans)}")
    print(f"USalign failed   : {len(failed)}")
    print(f"Written          : {out}")
    return 0 if not failed and not orphans else 2


SAMPLE = """\
Name of Structure_1: a.pdb:H (to be superimposed onto Structure_2)
Name of Structure_2: b.pdb:H
Length of Structure_1: 117 residues
Length of Structure_2: 116 residues

Aligned length= 110, RMSD=   1.25, Seq_ID=n_identical/n_aligned= 0.095
TM-score= 0.91000 (normalized by length of Structure_1: L=117, d0=3.97)
TM-score= 0.92000 (normalized by length of Structure_2: L=116, d0=3.97)
"""


def selftest() -> int:
    got = parse(SAMPLE)
    assert got is not None
    # The value after "Structure_1" is the design's; after "Structure_2" the backbone's.
    assert got["tm_design"] == 0.91 and got["tm_scaffold"] == 0.92, got
    assert got["rmsd"] == 1.25 and got["aligned_length"] == 110, got
    assert (got["len_design"], got["len_scaffold"]) == (117, 116), got
    # A report with a missing field must be refused, not guessed.
    assert parse(SAMPLE.replace("RMSD=", "RMS=")) is None
    assert NAME_RE.match("0866_top.3.20.pdb")["bid"] == "0866"
    print("selftest: pass")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Align every folded design with its backbone.")
    parser.add_argument("--usalign", help="Path to the USalign program")
    parser.add_argument("--scaffolds", type=Path, help="Folder with gen_NNNN.pdb backbones")
    parser.add_argument("--models", type=Path, help="Folder with folded designs, NNNN_top.K.TAG.pdb")
    parser.add_argument("--out", type=Path, help="CSV file to write")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        return selftest()
    if not (args.usalign and args.scaffolds and args.models and args.out):
        parser.error("--usalign, --scaffolds, --models and --out are required")
    return run(args.usalign, args.scaffolds, args.models, args.out, args.workers)


if __name__ == "__main__":
    sys.exit(main())
