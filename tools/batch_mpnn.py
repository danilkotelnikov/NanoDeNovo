#!/usr/bin/env python3
# version: 0.3.0-a43ea1-4 | date: 15.XI.2026
"""Batch ProteinMPNN sequence design with one model load per dataset.

Fixes Findings 4 and 11 from AUDIT.md:
- Writes and reads the exact same folder (`pmpnn_{tag}/seqs`), eliminating the
  `pmpnn_030` vs `pmpnn_30` rename mismatch.
- Parses and scores a folder of scaffolds in one JSONL batch invocation instead of
  reloading model weights 1,000 times per dataset.

Usage:
    py -3 -I batch_mpnn.py --scaffolds ./scaffolds --out-dir pmpnn_020 \
        --tag 20 --mpnn-dir /path/to/ProteinMPNN [--dry-run]
    py -3 -I batch_mpnn.py --selftest
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ID_RE = re.compile(r"(\d{3,})")


def discover_scaffolds(scaffolds_dir: Path, start: int = 1, end: int | None = None) -> list[tuple[int, Path]]:
    """Return sorted (backbone_index, path) pairs within inclusive [start, end]."""
    found: list[tuple[int, Path]] = []
    for pdb in sorted(scaffolds_dir.glob("*.pdb")):
        match = ID_RE.search(pdb.stem)
        if not match:
            continue
        idx = int(match.group(1))
        if idx < start or (end is not None and idx > end):
            continue
        found.append((idx, pdb))
    return sorted(found, key=lambda item: item[0])


def normalize_fasta_names(seqs_dir: Path, scaffolds: list[tuple[int, Path]], tag: str, temp: float) -> tuple[int, list[str]]:
    """Standardize raw ProteinMPNN FASTA names to `{bid:04d}_{tag}.t.{temp}.fa` inside the same folder."""
    seqs_dir.mkdir(parents=True, exist_ok=True)
    renamed = 0
    missing: list[str] = []
    temp_str = f"{temp:g}"
    for idx, pdb in scaffolds:
        bid = f"{idx:04d}"
        canonical = seqs_dir / f"{bid}_{tag}.t.{temp_str}.fa"
        if canonical.exists():
            renamed += 1
            continue
        candidates = [
            seqs_dir / f"{pdb.stem}.fa",
            seqs_dir / f"{pdb.stem}.fasta",
            seqs_dir / f"gen_{bid}.fa",
        ]
        src = next((c for c in candidates if c.exists()), None)
        if src is None:
            missing.append(bid)
            continue
        os.replace(src, canonical)
        renamed += 1
    return renamed, missing


def build_commands(
    mpnn_dir: Path,
    scaffolds_dir: Path,
    out_dir: Path,
    jsonl_path: Path,
    num_seqs: int,
    temp: float,
    batch_size: int,
    seed: int,
    python_bin: str = sys.executable,
) -> list[list[str]]:
    """Construct the two single-session commands: parse all PDBs, then sample in one model load."""
    parse_cmd = [
        python_bin,
        str(mpnn_dir / "helper_scripts" / "parse_multiple_chains.py"),
        "--input_path", str(scaffolds_dir),
        "--output_path", str(jsonl_path),
    ]
    run_cmd = [
        python_bin,
        str(mpnn_dir / "protein_mpnn_run.py"),
        "--jsonl_path", str(jsonl_path),
        "--out_folder", str(out_dir),
        "--num_seq_per_target", str(num_seqs),
        "--sampling_temp", str(temp),
        "--batch_size", str(batch_size),
        "--seed", str(seed),
    ]
    return [parse_cmd, run_cmd]


def run(args: argparse.Namespace) -> int:
    scaffolds = discover_scaffolds(args.scaffolds, args.start, args.end)
    if not scaffolds:
        print(f"No scaffold PDBs found in {args.scaffolds}", file=sys.stderr)
        return 1

    out_dir: Path = args.out_dir
    seqs_dir = out_dir / "seqs"
    out_dir.mkdir(parents=True, exist_ok=True)
    jsonl_path = out_dir / f"parsed_chains_{args.tag}.jsonl"
    cmds = build_commands(
        args.mpnn_dir, args.scaffolds, out_dir, jsonl_path,
        args.num_seqs, args.temp, args.batch_size, args.seed, args.python_bin,
    )

    plan = {
        "tag": args.tag,
        "scaffold_count": len(scaffolds),
        "first_id": scaffolds[0][0],
        "last_id": scaffolds[-1][0],
        "out_dir": str(out_dir),
        "seqs_dir": str(seqs_dir),
        "model_loads": 1,
        "commands": [" ".join(c) for c in cmds],
    }
    if args.dry_run:
        print(json.dumps(plan, indent=2))
        return 0

    for cmd in cmds:
        proc = subprocess.run(cmd, check=False)
        if proc.returncode != 0:
            print(f"Command failed ({proc.returncode}): {' '.join(cmd)}", file=sys.stderr)
            return proc.returncode

    done, missing = normalize_fasta_names(seqs_dir, scaffolds, args.tag, args.temp)
    print(f"Backbones    : {len(scaffolds)}")
    print(f"Normalized   : {done} in {seqs_dir}")
    if missing:
        print(f"Missing FASTA: {len(missing)} (first: {missing[0]})", file=sys.stderr)
        return 1
    return 0


def selftest() -> int:
    """Verify inclusive indexing, one-load command plan, and consistent tag 20/30 seqs dir naming."""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        scaf = root / "scaffolds"
        scaf.mkdir()
        for idx in (1, 500, 1000):
            (scaf / f"gen_{idx:04d}.pdb").write_text("ATOM      1  CA  ALA H   1       0.0   0.0   0.0  1.00  0.00           C\n", encoding="utf-8")

        found = discover_scaffolds(scaf, 1, 1000)
        assert [i for i, _ in found] == [1, 500, 1000], found

        for tag in ("20", "30"):
            out = root / f"pmpnn_0{tag}"
            seqs = out / "seqs"
            seqs.mkdir(parents=True)
            for idx, pdb in found:
                (seqs / f"{pdb.stem}.fa").write_text(f">sample=1, score=0.90\nAAAA\n", encoding="utf-8")
            cmds = build_commands(root / "ProteinMPNN", scaf, out, out / f"parsed_{tag}.jsonl", 128, 0.15, 32, 42)
            assert len(cmds) == 2
            assert "--jsonl_path" in cmds[1] and "--out_folder" in cmds[1]
            assert str(out) in cmds[1]
            done, missing = normalize_fasta_names(seqs, found, tag, 0.15)
            assert done == 3 and not missing, (done, missing)
            assert (seqs / f"1000_{tag}.t.0.15.fa").exists()

    print("batch_mpnn selftest: pass")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Batch ProteinMPNN design with one model load per dataset.")
    ap.add_argument("--scaffolds", type=Path, help="Folder of scaffold PDB files (gen_0001.pdb..)")
    ap.add_argument("--out-dir", type=Path, help="Output folder (reads and writes `<out-dir>/seqs`)")
    ap.add_argument("--tag", default="20", help="Dataset temperature/run tag (e.g. 20 or 30)")
    ap.add_argument("--mpnn-dir", type=Path, default=Path("ProteinMPNN"), help="Path to ProteinMPNN checkout")
    ap.add_argument("--python-bin", default=sys.executable, help="Python interpreter for ProteinMPNN")
    ap.add_argument("--start", type=int, default=1, help="First scaffold index (inclusive)")
    ap.add_argument("--end", type=int, default=None, help="Last scaffold index (inclusive)")
    ap.add_argument("--num-seqs", type=int, default=128, help="Sequences sampled per scaffold")
    ap.add_argument("--temp", type=float, default=0.15, help="ProteinMPNN sampling temperature")
    ap.add_argument("--batch-size", type=int, default=32, help="ProteinMPNN batch size")
    ap.add_argument("--seed", type=int, default=42, help="Random seed")
    ap.add_argument("--dry-run", action="store_true", help="Print command plan without executing")
    ap.add_argument("--selftest", action="store_true", help="Run built-in verification and exit")
    args = ap.parse_args()

    if args.selftest:
        return selftest()
    if args.scaffolds is None or args.out_dir is None:
        ap.error("--scaffolds and --out-dir are required")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
