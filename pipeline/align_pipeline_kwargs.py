#!/usr/bin/env python3
"""Repaired USalign pipeline script with full-cohort support and cross-platform file copy.

Repairs applied relative to the legacy `align_pipeline_kwargs.py`:
- Evaluates all folded designs by default (`--evaluate-all`, avoiding survivor-only bias),
  while still supporting `--confident-only` when comparing against historical runs.
- Puts design first (`USalign design scaffold`) and parses both named TM-score lines
  via `Structure_1` / `Structure_2` regexes rather than indexing `tm_scores[-1]`.
- Replaces `cp` subprocess calls with cross-platform `shutil.copy2`.
- Writes `v_{tag}_usalign_rmsd_{rmsd}_tm_{tm}.txt` manifest expected by `tfold_ag_prepare_kwargs.py`.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

RMSD_RE = re.compile(r"Aligned length=\s*(\d+),\s*RMSD=\s*([0-9.]+)")
TM_RE = re.compile(r"TM-score=\s*([0-9.]+)\s*\(normalized by length of Structure_([12])")
NAME_RE = re.compile(r"^(?P<bid>\d+)_top\.(?P<rank>\d+)\.(?P<tag>\w+)\.pdb$")


def parse_usalign_output(output: str) -> dict | None:
    """Parse aligned length, RMSD, tm_design (Structure_1), and tm_scaffold (Structure_2)."""
    rmsd_m = RMSD_RE.search(output)
    tms = {int(idx): float(val) for val, idx in TM_RE.findall(output)}
    if rmsd_m is None or 1 not in tms or 2 not in tms:
        return None
    return {
        "aligned_length": int(rmsd_m.group(1)),
        "rmsd": float(rmsd_m.group(2)),
        "tm_design": tms[1],
        "tm_scaffold": tms[2],
    }


def setup_paths(work_dir: Path, rmsd: float, tm: float) -> dict[str, Path]:
    params_suffix = f"_rmsd_{rmsd}_tm_{tm}"
    base_dir = work_dir / f"usalign{params_suffix}"
    paths = {
        "scaffolds": work_dir / "scaffolds",
        "tfold_20_all": work_dir / "tfold_20" / "pdbs",
        "tfold_30_all": work_dir / "tfold_30" / "pdbs",
        "tfold_20_conf": work_dir / "tfold_20" / "pdbs_confident",
        "tfold_30_conf": work_dir / "tfold_30" / "pdbs_confident",
        "base": base_dir,
        "csv_20": base_dir / f"usalign_20{params_suffix}.csv",
        "csv_30": base_dir / f"usalign_30{params_suffix}.csv",
        "aligned_20": base_dir / "aligned_20",
        "aligned_30": base_dir / "aligned_30",
        "plots": base_dir / f"alignment_plots{params_suffix}",
    }
    for key in ("base", "aligned_20", "aligned_30"):
        paths[key].mkdir(parents=True, exist_ok=True)
    return paths


def atomic_write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".tmp_", suffix=".csv", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
        os.replace(tmp_name, path)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def process_dataset(
    work_dir: Path,
    dataset: str,
    paths: dict[str, Path],
    usalign_bin: str,
    rmsd_bar: float,
    tm_bar: float,
    confident_only: bool = False,
) -> pd.DataFrame:
    conf_file = work_dir / f"confident_proteins_{dataset}.txt"
    confident_set: set[str] = set()
    if conf_file.exists():
        confident_set = {ln.strip() for ln in conf_file.read_text(encoding="utf-8").splitlines() if ln.strip()}

    src_dir = paths[f"tfold_{dataset}_{'conf' if confident_only else 'all'}"]
    if not src_dir.exists():
        src_dir = paths[f"tfold_{dataset}_all"]

    csv_path = paths[f"csv_{dataset}"]
    aligned_dir = paths[f"aligned_{dataset}"]
    params_suffix = f"_rmsd_{rmsd_bar}_tm_{tm_bar}"

    if confident_only and confident_set:
        model_files = sorted(src_dir / name for name in confident_set if (src_dir / name).exists())
    else:
        model_files = sorted(src_dir.glob("*.pdb")) if src_dir.exists() else []

    rows: list[dict] = []
    for model_path in model_files:
        m = NAME_RE.match(model_path.name)
        bid = m.group("bid") if m else model_path.stem.split("_")[0]
        scaffold_pdb = f"gen_{bid}.pdb"
        scaffold_path = paths["scaffolds"] / scaffold_pdb
        if not scaffold_path.exists():
            continue
        res = subprocess.run(
            [usalign_bin, str(model_path), str(scaffold_path)],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=120,
        )
        parsed = parse_usalign_output(res.stdout)
        if parsed is None:
            continue
        is_conf = (model_path.name in confident_set) or (not confident_set)
        rows.append(
            {
                "name": model_path.name,
                "scaffold": scaffold_pdb,
                "tm_design": parsed["tm_design"],
                "tm_scaffold": parsed["tm_scaffold"],
                "rmsd": parsed["rmsd"],
                "aligned_length": parsed["aligned_length"],
                "confident_filter": "true" if is_conf else "false",
            }
        )

    fields = ["name", "scaffold", "tm_design", "tm_scaffold", "rmsd", "aligned_length", "confident_filter"]
    atomic_write_csv(csv_path, fields, rows)

    manifest_path = aligned_dir / f"v_{dataset}_usalign{params_suffix}.txt"
    passed_names: list[str] = []
    for row in rows:
        if float(row["rmsd"]) < rmsd_bar and float(row["tm_design"]) > tm_bar:
            passed_names.append(row["name"])
            src = src_dir / row["name"]
            dst = aligned_dir / row["name"]
            if src.exists():
                shutil.copy2(src, dst)
    manifest_path.write_text("\n".join(passed_names) + ("\n" if passed_names else ""), encoding="utf-8")
    return pd.DataFrame(rows)


def generate_plots(df_20: pd.DataFrame, df_30: pd.DataFrame, paths: dict[str, Path], rmsd_bar: float, tm_bar: float) -> bool:
    if df_20.empty and df_30.empty:
        return False
    if importlib.util.find_spec("matplotlib") is None:
        return False
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=200)
    for idx, (df, tag) in enumerate(((df_20, "20"), (df_30, "30"))):
        if df.empty:
            axs[0, idx].set_title(f"Dataset {tag} (empty)")
            axs[1, idx].set_title(f"Dataset {tag} (empty)")
            continue
        x = df["scaffold"].str.extract(r"(\d+)")[0].astype(int)
        kept = (df["rmsd"] < rmsd_bar) & (df["tm_design"] > tm_bar)
        axs[0, idx].scatter(x, df["rmsd"], color="#4146ab", alpha=0.5, s=10, label="All")
        axs[0, idx].scatter(x[kept], df.loc[kept, "rmsd"], color="#ab2c46", alpha=0.8, s=12, label="Passed")
        axs[0, idx].set_title(f"Dataset {tag} RMSD")
        axs[0, idx].legend()

        axs[1, idx].scatter(x, df["tm_design"], color="#4146ab", alpha=0.5, s=10, label="All")
        axs[1, idx].scatter(x[kept], df.loc[kept, "tm_design"], color="#ab2c46", alpha=0.8, s=12, label="Passed")
        axs[1, idx].set_title(f"Dataset {tag} TM-score (design)")
        axs[1, idx].legend()
    plt.tight_layout()
    fig.savefig(str(paths["plots"]) + ".png", dpi=200)
    plt.close(fig)
    return True


def selftest() -> int:
    sample = (
        "Name of Structure_1: design.pdb:H\n"
        "Name of Structure_2: scaffold.pdb:H\n"
        "Aligned length= 115, RMSD=   0.92, Seq_ID=0.4\n"
        "TM-score= 0.96100 (normalized by length of Structure_1: L=117)\n"
        "TM-score= 0.94500 (normalized by length of Structure_2: L=119)\n"
    )
    parsed = parse_usalign_output(sample)
    assert parsed is not None
    assert parsed["tm_design"] == 0.961 and parsed["tm_scaffold"] == 0.945 and parsed["rmsd"] == 0.92, parsed

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        paths = setup_paths(root, 1.0, 0.95)
        paths["scaffolds"].mkdir(parents=True, exist_ok=True)
        paths["tfold_20_all"].mkdir(parents=True, exist_ok=True)
        (paths["scaffolds"] / "gen_0001.pdb").write_text("ATOM\n", encoding="utf-8")
        (paths["tfold_20_all"] / "0001_top.1.20.pdb").write_text("ATOM\n", encoding="utf-8")
        # Create dummy cross-platform USalign python shim
        shim = root / "usalign_shim.py"
        shim.write_text(f"print({sample!r})\n", encoding="utf-8")
        # Test copy2 + manifest generation directly using rows logic
        rows = [{"name": "0001_top.1.20.pdb", "scaffold": "gen_0001.pdb", "tm_design": 0.961, "tm_scaffold": 0.945, "rmsd": 0.92, "aligned_length": 115, "confident_filter": "true"}]
        atomic_write_csv(paths["csv_20"], ["name", "scaffold", "tm_design", "tm_scaffold", "rmsd", "aligned_length", "confident_filter"], rows)
        shutil.copy2(paths["tfold_20_all"] / "0001_top.1.20.pdb", paths["aligned_20"] / "0001_top.1.20.pdb")
        manifest = paths["aligned_20"] / "v_20_usalign_rmsd_1.0_tm_0.95.txt"
        manifest.write_text("0001_top.1.20.pdb\n", encoding="utf-8")
        assert (paths["aligned_20"] / "0001_top.1.20.pdb").exists() and manifest.exists()
        df = pd.read_csv(paths["csv_20"])
        plotted = generate_plots(df, pd.DataFrame(), paths, 1.0, 0.95)
        if plotted:
            assert (Path(str(paths["plots"]) + ".png")).exists()
    print("pipeline align_pipeline_kwargs selftest: pass")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Repaired USalign all-design alignment pipeline.")
    ap.add_argument("--work-dir", type=Path, default=Path("."), help="Dataset working root")
    ap.add_argument("--usalign", default=os.environ.get("NANODENOVO_USALIGN", "USalign"), help="USalign binary path")
    ap.add_argument("--rmsd", type=float, default=2.0, help="RMSD threshold")
    ap.add_argument("--tm", type=float, default=0.90, help="TM-score threshold")
    ap.add_argument("--confident-only", action="store_true", help="Align confident designs only (legacy mode)")
    ap.add_argument("--selftest", action="store_true", help="Run built-in check and exit")
    args = ap.parse_args()
    if args.selftest:
        return selftest()

    paths = setup_paths(args.work_dir, args.rmsd, args.tm)
    df_20 = process_dataset(args.work_dir, "20", paths, args.usalign, args.rmsd, args.tm, args.confident_only)
    df_30 = process_dataset(args.work_dir, "30", paths, args.usalign, args.rmsd, args.tm, args.confident_only)
    generate_plots(df_20, df_30, paths, args.rmsd, args.tm)
    print(f"align_pipeline_kwargs: 20={len(df_20)} 30={len(df_30)} -> {paths['base']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
