#!/usr/bin/env python3
# version: 0.3.0-a43ea1-4 | date: 15.XI.2026
"""USalign structural alignment and scoring of refolded nanobodies against backbones.

Loads refold_bundle.json.gz and harvested refold_ca.npz, pre-writes 1,000 reference
scaffolds into /dev/shm, runs parallel USalign across 6,200 CA traces, computes
normalized TM-scores (by design and scaffold), RMSD, and aligned length, and
writes refold_usalign.csv and cohort summary refold_usalign_summary.json.

Usage:
    py -3 -I score_refold_usalign.py --selftest
    py -3 -I score_refold_usalign.py --bundle PATH --ca PATH --metrics PATH \\
        --out-csv PATH --out-summary PATH [--workers 12]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shlex
import shutil
import statistics as st
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

RMSD_RE = re.compile(r"Aligned length=\s*(\d+),\s*RMSD=\s*([0-9.]+)")
TM_RE = re.compile(r"TM-score=\s*([0-9.]+)\s*\(normalized by length of Structure_([12])")
LEN_RE = re.compile(r"Length of Structure_([12]):\s*(\d+)")

FIELDS = [
    "id",
    "group",
    "name",
    "backbone",
    "rank",
    "nll_score",
    "length",
    "len_scaffold",
    "aligned_length",
    "rmsd",
    "tm_design",
    "tm_scaffold",
    "pass_common",
    "pass_strict",
]


def parse_usalign(text: str) -> dict | None:
    """Read one USalign report. Return None when a required field is missing."""
    rmsd = RMSD_RE.search(text)
    tms = {int(num): float(val) for val, num in TM_RE.findall(text)}
    lens = {int(num): int(val) for num, val in LEN_RE.findall(text)}
    if rmsd is None or set(tms) != {1, 2} or set(lens) != {1, 2}:
        return None
    return {
        "len_design": lens[1],
        "len_scaffold": lens[2],
        "aligned_length": int(rmsd.group(1)),
        "rmsd": float(rmsd.group(2)),
        "tm_design": tms[1],
        "tm_scaffold": tms[2],
    }


def format_ca_pdb(coords: np.ndarray, chain: str = "A") -> str:
    """Format C-alpha coordinates as minimal PDB text with valid ATOM records."""
    lines = []
    for i, row in enumerate(coords, start=1):
        x, y, z = float(row[0]), float(row[1]), float(row[2])
        lines.append(
            f"ATOM  {i:5d}  CA  ALA {chain}{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
        )
    lines.append("TER\nEND\n")
    return "".join(lines)


def summarize_cohorts(rows: list[dict], elapsed_s: float) -> dict:
    """Compute summary statistics across groups for USalign alignments."""
    groups = sorted({r["group"] for r in rows})
    out: dict = {
        "elapsed_s": round(elapsed_s, 2),
        "total": len(rows),
        "groups": {},
    }
    for g in groups:
        sub = [r for r in rows if r["group"] == g]
        rmsds = [float(r["rmsd"]) for r in sub]
        tms_des = [float(r["tm_design"]) for r in sub]
        tms_scaf = [float(r["tm_scaffold"]) for r in sub]
        lens = [int(r["aligned_length"]) for r in sub]
        out["groups"][g] = {
            "n": len(sub),
            "mean_rmsd": round(st.mean(rmsds), 4),
            "median_rmsd": round(st.median(rmsds), 4),
            "mean_tm_design": round(st.mean(tms_des), 4),
            "median_tm_design": round(st.median(tms_des), 4),
            "mean_tm_scaffold": round(st.mean(tms_scaf), 4),
            "median_tm_scaffold": round(st.median(tms_scaf), 4),
            "mean_aligned_length": round(st.mean(lens), 2),
            "median_aligned_length": round(st.median(lens), 2),
            "pass_common": sum(r["pass_common"] == "true" for r in sub),
            "pass_strict": sum(r["pass_strict"] == "true" for r in sub),
        }
    return out


def win_to_wsl(path_str: str) -> str:
    """Convert Windows path C:\\... to WSL /mnt/c/... path."""
    p = Path(path_str).resolve()
    drive = p.drive.rstrip(":").lower()
    posix_tail = p.as_posix().split(":", 1)[-1]
    return f"/mnt/{drive}{posix_tail}"


def resolve_wsl_python() -> str:
    """Resolve WSL Python executable from NANODENOVO_WSL_PYTHON or default PATH python3."""
    return os.environ.get("NANODENOVO_WSL_PYTHON", "python3")


def build_wsl_dispatch_command(args: argparse.Namespace, script_path: str = __file__) -> list[str]:
    """Build a safe WSL invocation using quoted paths so spaces and Cyrillic directories survive."""
    wsl_python = resolve_wsl_python()
    wsl_script = win_to_wsl(script_path)

    py_args = [
        wsl_python,
        wsl_script,
        "--bundle",
        win_to_wsl(args.bundle),
        "--ca",
        win_to_wsl(args.ca),
        "--out-csv",
        win_to_wsl(args.out_csv),
        "--out-summary",
        win_to_wsl(args.out_summary),
        "--usalign",
        args.usalign,
        "--workers",
        str(args.workers),
    ]
    if args.metrics:
        py_args.extend(["--metrics", win_to_wsl(args.metrics)])
    if args.shm:
        py_args.extend(["--shm", args.shm])

    quoted = " ".join(shlex.quote(tok) for tok in py_args)
    return ["wsl", "-e", "bash", "-lc", quoted]


def run_wsl_dispatch(args: argparse.Namespace) -> int:
    """Dispatch execution to WSL python environment with mapped arguments."""
    cmd = build_wsl_dispatch_command(args)
    print(f"Windows host dispatching to WSL: {' '.join(cmd)}", flush=True)
    res = subprocess.run(cmd, encoding="utf-8", errors="replace")
    return res.returncode


def run_alignments(
    bundle_path: Path,
    ca_path: Path,
    metrics_path: Path | None,
    usalign_bin: str,
    shm_root: Path,
    out_csv: Path,
    out_summary: Path,
    workers: int,
) -> int:
    """Execute parallel USalign structural alignments in memory."""
    t0 = time.time()
    if not Path(usalign_bin).exists():
        raise FileNotFoundError(f"USalign executable not found at: {usalign_bin}")

    print(f"Reading reference bundle from {bundle_path}...", flush=True)
    with gzip.open(bundle_path, "rt", encoding="utf-8") as fh:
        bundle = json.load(fh)
    scaffolds_raw = bundle["scaffolds"]
    designs = bundle["designs"]

    print(f"Reading CA coordinates from {ca_path}...", flush=True)
    with np.load(ca_path) as npz:
        ca_dict = dict(npz)

    # Optional metrics lookup for nll_score or cross-validation
    metrics_map: dict[str, dict] = {}
    if metrics_path and metrics_path.exists():
        print(f"Reading metrics from {metrics_path}...", flush=True)
        with open(metrics_path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for r in reader:
                metrics_map[r["id"]] = r

    # Create temporary scratch directory in /dev/shm
    shm_dir = shm_root / f"usalign_refold_{os.getpid()}"
    shm_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Pre-write 1,000 reference scaffold PDBs
        print(f"Pre-writing {len(scaffolds_raw)} reference scaffolds to {shm_dir}...", flush=True)
        t_scaf0 = time.time()
        for bid, sdata in scaffolds_raw.items():
            scaf_ca = np.array(sdata["ca"], dtype=float)
            scaf_pdb = shm_dir / f"scaf_{bid}.pdb"
            scaf_pdb.write_text(format_ca_pdb(scaf_ca), encoding="ascii")
        print(f"Scaffolds ready in {time.time() - t_scaf0:.2f}s", flush=True)

        print(f"Aligning {len(designs)} designs using {workers} worker threads...", flush=True)
        rows: list[dict] = []
        failed: list[str] = []

        def align_worker(d: dict) -> dict | None:
            did = d["id"]
            bid = d["bid"]
            if did not in ca_dict:
                return None
            ca_arr = ca_dict[did]
            tid = threading.get_ident()
            des_pdb = shm_dir / f"des_{tid}.pdb"
            des_pdb.write_text(format_ca_pdb(ca_arr), encoding="ascii")
            scaf_pdb = shm_dir / f"scaf_{bid}.pdb"

            proc = subprocess.run(
                [usalign_bin, str(des_pdb), str(scaf_pdb)],
                capture_output=True,
                text=True,
                timeout=60,
            )
            parsed = parse_usalign(proc.stdout)
            if parsed is None:
                return None

            nll = d.get("nll")
            if nll is None and did in metrics_map:
                nll_str = metrics_map[did].get("nll_score", "")
            else:
                nll_str = f"{nll:.4f}" if nll is not None else ""

            pass_comm = (
                "true" if (parsed["tm_scaffold"] >= 0.90 and parsed["rmsd"] <= 2.0) else "false"
            )
            pass_strc = (
                "true" if (parsed["tm_scaffold"] >= 0.95 and parsed["rmsd"] <= 1.0) else "false"
            )

            return {
                "id": did,
                "group": d["group"],
                "name": d["name"],
                "backbone": bid,
                "rank": d["rank"],
                "nll_score": nll_str,
                "length": len(d.get("seq", "")) or parsed["len_design"],
                "len_scaffold": parsed["len_scaffold"],
                "aligned_length": parsed["aligned_length"],
                "rmsd": f"{parsed['rmsd']:.3f}",
                "tm_design": f"{parsed['tm_design']:.5f}",
                "tm_scaffold": f"{parsed['tm_scaffold']:.5f}",
                "pass_common": pass_comm,
                "pass_strict": pass_strc,
            }

        with ThreadPoolExecutor(max_workers=workers) as pool:
            for i, result in enumerate(pool.map(align_worker, designs), start=1):
                if result is None:
                    failed.append(designs[i - 1]["id"])
                else:
                    rows.append(result)
                if i % 1000 == 0 or i == len(designs):
                    dt = time.time() - t0
                    print(
                        f"Scored {i}/{len(designs)} alignments in {dt:.1f}s ({i/dt:.1f} designs/s)",
                        flush=True,
                    )

        if failed:
            print(f"WARNING: {len(failed)} alignments failed parsing!", file=sys.stderr)

        # Write atomic CSV
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        tmp_csv = out_csv.with_suffix(".tmp")
        with open(tmp_csv, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        tmp_csv.replace(out_csv)
        print(f"Saved {len(rows)} alignments -> {out_csv}", flush=True)

        # Compute summary
        elapsed = time.time() - t0
        summary = summarize_cohorts(rows, elapsed)
        out_summary.parent.mkdir(parents=True, exist_ok=True)
        tmp_sum = out_summary.with_suffix(".tmp")
        tmp_sum.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        tmp_sum.replace(out_summary)
        print(f"Saved summary -> {out_summary}", flush=True)
        print(json.dumps(summary, indent=2), flush=True)

    finally:
        shutil.rmtree(shm_dir, ignore_errors=True)

    return 0 if not failed else 1


SAMPLE_USALIGN = """\
Name of Structure_1: /dev/shm/test_des.pdb:A (to be superimposed onto Structure_2)
Name of Structure_2: /dev/shm/test_scaf.pdb:A
Length of Structure_1: 118 residues
Length of Structure_2: 118 residues

Aligned length= 116, RMSD=   2.09, Seq_ID=n_identical/n_aligned= 1.000
TM-score= 0.80349 (normalized by length of Structure_1: L=118, d0=4.01)
TM-score= 0.80349 (normalized by length of Structure_2: L=118, d0=4.01)
(You should use TM-score normalized by length of the reference structure)

(":" denotes residue pairs of d < 5.0 Angstrom, "." denotes other aligned residues)
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-A-
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: .:::::: :::::: : 
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-AAAAAAA-AAAAAAAAA

#Total CPU time is  0.02 seconds
"""


def selftest() -> int:
    """Run offline unit tests for parser, PDB formatting, and cohort summarizer."""
    parsed = parse_usalign(SAMPLE_USALIGN)
    assert parsed is not None, "Failed to parse valid USalign output"
    assert parsed["len_design"] == 118
    assert parsed["len_scaffold"] == 118
    assert parsed["aligned_length"] == 116
    assert abs(parsed["rmsd"] - 2.09) < 1e-4
    assert abs(parsed["tm_design"] - 0.80349) < 1e-5
    assert abs(parsed["tm_scaffold"] - 0.80349) < 1e-5

    # Test malformed output
    assert parse_usalign("RMSD=invalid") is None
    assert parse_usalign(SAMPLE_USALIGN.replace("TM-score=", "TM=")) is None

    # Test PDB formatting
    coords = np.array([[1.0, 2.0, 3.0], [4.5, 5.5, 6.5]])
    pdb_txt = format_ca_pdb(coords)
    assert "ATOM      1  CA  ALA A   1       1.000   2.000   3.000" in pdb_txt
    assert "ATOM      2  CA  ALA A   2       4.500   5.500   6.500" in pdb_txt
    assert pdb_txt.endswith("TER\nEND\n")

    # Test cohort summarizer
    dummy_rows = [
        {
            "id": "1",
            "group": "corr_20",
            "name": "d1",
            "backbone": "0001",
            "rank": "1",
            "nll_score": "0.5",
            "length": 118,
            "len_scaffold": 118,
            "aligned_length": 116,
            "rmsd": "1.500",
            "tm_design": "0.92000",
            "tm_scaffold": "0.92000",
            "pass_common": "true",
            "pass_strict": "false",
        },
        {
            "id": "2",
            "group": "corr_20",
            "name": "d2",
            "backbone": "0001",
            "rank": "2",
            "nll_score": "0.6",
            "length": 118,
            "len_scaffold": 118,
            "aligned_length": 110,
            "rmsd": "0.800",
            "tm_design": "0.96000",
            "tm_scaffold": "0.96000",
            "pass_common": "true",
            "pass_strict": "true",
        },
        {
            "id": "3",
            "group": "orig_20",
            "name": "d3",
            "backbone": "0002",
            "rank": "1",
            "nll_score": "1.2",
            "length": 120,
            "len_scaffold": 120,
            "aligned_length": 100,
            "rmsd": "4.500",
            "tm_design": "0.60000",
            "tm_scaffold": "0.60000",
            "pass_common": "false",
            "pass_strict": "false",
        },
    ]
    summary = summarize_cohorts(dummy_rows, 1.23)
    assert summary["total"] == 3
    assert summary["groups"]["corr_20"]["n"] == 2
    assert summary["groups"]["corr_20"]["pass_common"] == 2
    assert summary["groups"]["corr_20"]["pass_strict"] == 1
    assert abs(summary["groups"]["corr_20"]["mean_rmsd"] - 1.15) < 1e-4
    assert summary["groups"]["orig_20"]["n"] == 1
    assert summary["groups"]["orig_20"]["pass_common"] == 0

    # Test WSL dispatch quoting with spaces and non-ASCII directory names plus env override
    prev_py = os.environ.get("NANODENOVO_WSL_PYTHON")
    try:
        os.environ["NANODENOVO_WSL_PYTHON"] = "/opt/venv/bin/python3"
        dummy_ns = argparse.Namespace(
            bundle=r"C:\work\Desktop Рабочий\NanoDeNovo\bundle.json.gz",
            ca=r"C:\work\Desktop Рабочий\NanoDeNovo\ca.npz",
            out_csv=r"C:\work\Desktop Рабочий\NanoDeNovo\out.csv",
            out_summary=r"C:\work\Desktop Рабочий\NanoDeNovo\sum.json",
            usalign="USalign",
            workers=4,
            metrics=None,
            shm="/dev/shm",
        )
        disp = build_wsl_dispatch_command(dummy_ns, __file__)
        assert disp[:4] == ["wsl", "-e", "bash", "-lc"]
        assert "/opt/venv/bin/python3" in disp[4]
        assert "'/mnt/c/work/Desktop Рабочий/NanoDeNovo/bundle.json.gz'" in disp[4], disp[4]
    finally:
        if prev_py is None:
            os.environ.pop("NANODENOVO_WSL_PYTHON", None)
        else:
            os.environ["NANODENOVO_WSL_PYTHON"] = prev_py

    print("score_refold_usalign selftest: pass")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="USalign structural alignment and scoring of refolded nanobodies."
    )
    parser.add_argument("--selftest", action="store_true", help="Run offline self-test and exit")
    parser.add_argument("--bundle", type=str, help="Path to refold_bundle.json.gz")
    parser.add_argument("--ca", type=str, help="Path to refold_ca.npz")
    parser.add_argument("--metrics", type=str, default=None, help="Path to refold_metrics.csv")
    parser.add_argument(
        "--usalign",
        type=str,
        default=os.environ.get("NANODENOVO_USALIGN", "USalign"),
        help="Path or command name for USalign executable",
    )
    parser.add_argument(
        "--shm",
        type=str,
        default="/dev/shm",
        help="Directory for shared memory temporary structures",
    )
    parser.add_argument("--out-csv", type=str, help="Path to output CSV")
    parser.add_argument("--out-summary", type=str, help="Path to output JSON summary")
    parser.add_argument(
        "--workers", type=int, default=12, help="Number of parallel alignment workers"
    )

    args = parser.parse_args()
    if args.selftest:
        return selftest()

    if not (args.bundle and args.ca and args.out_csv and args.out_summary):
        parser.error("--bundle, --ca, --out-csv, and --out-summary are required")

    # If running on Windows host, dispatch transparently to WSL
    if sys.platform == "win32":
        return run_wsl_dispatch(args)

    return run_alignments(
        bundle_path=Path(args.bundle),
        ca_path=Path(args.ca),
        metrics_path=Path(args.metrics) if args.metrics else None,
        usalign_bin=args.usalign,
        shm_root=Path(args.shm),
        out_csv=Path(args.out_csv),
        out_summary=Path(args.out_summary),
        workers=args.workers,
    )


if __name__ == "__main__":
    raise SystemExit(main())
