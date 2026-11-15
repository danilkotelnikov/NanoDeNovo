#!/usr/bin/env python3
"""Repaired confidence distribution plotter for NanoDeNovo datasets 20 and 30.

Repairs applied relative to legacy `build_plots.py`:
- Wraps execution inside `main()` so importing `build_plots` causes zero side effects.
- Uses headless Matplotlib `Agg` backend and Matplotlib 3.9+ `tick_labels=` compatibility.
- Guards against missing or empty confident CSV tables cleanly.
"""
from __future__ import annotations

import argparse
import importlib.util
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def save_plot(fig, stem_path: Path, plt_mod, formats: tuple[str, ...] = ("png", "svg")) -> None:
    stem_path.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(f"{stem_path}.{fmt}", bbox_inches="tight", dpi=200)
    plt_mod.close(fig)


def generate_summary_plots(work_dir: Path, out_dir: Path, dataset_types: tuple[str, ...] = ("20", "30")) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    has_mpl = importlib.util.find_spec("matplotlib") is not None
    if has_mpl:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    else:
        plt = None

    processed: dict[str, dict] = {}
    stats: dict[str, dict] = {}

    for dtype in dataset_types:
        gen_csv = work_dir / f"folded_sequences_{dtype}.csv"
        conf_csv = work_dir / f"confident_proteins_{dtype}.csv"
        if not gen_csv.exists():
            continue
        general_df = pd.read_csv(gen_csv)
        if general_df.empty:
            continue
        name_to_idx = {name: i for i, name in enumerate(general_df["Name"].tolist())}
        if conf_csv.exists():
            conf_df = pd.read_csv(conf_csv)
            c_idx = conf_df["Name"].map(name_to_idx)
            mask = c_idx.notna()
            conf_df = conf_df[mask]
            conf_indices = c_idx[mask].astype(int)
        else:
            conf_df = general_df.iloc[0:0]
            conf_indices = pd.Series([], dtype=int)

        processed[dtype] = {
            "general_plddt": general_df["pLDDT"],
            "general_ptm": general_df["pTM"],
            "confident_plddt": conf_df["pLDDT"],
            "confident_ptm": conf_df["pTM"],
            "confident_indices": conf_indices,
        }
        stats[dtype] = {
            "General pLDDT": general_df["pLDDT"].describe().to_dict(),
            "General pTM": general_df["pTM"].describe().to_dict(),
        }

        if plt is not None:
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 6), dpi=200)
            fig.subplots_adjust(hspace=0.4)
            ax1.scatter(range(len(general_df)), general_df["pLDDT"], color="#41ab89", s=6, label="All designs")
            if not conf_df.empty:
                ax1.scatter(conf_indices, conf_df["pLDDT"], color="#4146ab", s=8, label="Confident subset")
            ax1.set_title(f"v_48_{dtype}: pLDDT distribution")
            ax1.set_ylabel("pLDDT")
            ax1.set_xlabel("Order")
            ax1.legend(loc="lower right", fontsize=7)

            ax2.scatter(range(len(general_df)), general_df["pTM"], color="#41ab89", s=6, label="All designs")
            if not conf_df.empty:
                ax2.scatter(conf_indices, conf_df["pTM"], color="#4146ab", s=8, label="Confident subset")
            ax2.set_title(f"v_48_{dtype}: pTM distribution")
            ax2.set_ylabel("pTM")
            ax2.set_xlabel("Order")
            ax2.legend(loc="lower right", fontsize=7)
            save_plot(fig, out_dir / f"scatter_plots_{dtype}", plt)

    if processed and plt is not None:
        fig, axes = plt.subplots(len(processed), 1, figsize=(7, 4 * len(processed)), dpi=200, squeeze=False)
        for ax, (dtype, d) in zip(axes[:, 0], processed.items()):
            box_data = [
                d["general_plddt"].to_numpy(),
                d["confident_plddt"].to_numpy() if len(d["confident_plddt"]) else np.array([np.nan]),
                d["general_ptm"].to_numpy(),
                d["confident_ptm"].to_numpy() if len(d["confident_ptm"]) else np.array([np.nan]),
            ]
            labels = ["All pLDDT", "Kept pLDDT", "All pTM", "Kept pTM"]
            positions = [1, 2, 4, 5]
            bplot = ax.boxplot(
                box_data,
                positions=positions,
                patch_artist=True,
                tick_labels=labels,
                widths=0.55,
                showfliers=False,
            )
            for patch in bplot["boxes"]:
                patch.set_facecolor("none")
                patch.set_edgecolor("black")
            ax.set_title(f"Dataset v_48_{dtype}")
            ax.set_ylabel("Score")
        fig.tight_layout()
        save_plot(fig, out_dir / "box_plots", plt)

    return {"datasets_plotted": list(processed.keys()), "matplotlib_active": has_mpl, "out_dir": str(out_dir)}


def selftest() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "folded_sequences_20.csv").write_text(
            "Name,Sequence,pLDDT,pTM\n"
            "0001_top.1.20,QVQL,0.95,0.89\n"
            "0001_top.2.20,AAAA,0.78,0.65\n"
            "0002_top.1.20,EVQL,0.94,0.88\n",
            encoding="utf-8",
        )
        (root / "confident_proteins_20.csv").write_text(
            "Name,Sequence,pLDDT,pTM\n"
            "0001_top.1.20,QVQL,0.95,0.89\n"
            "0002_top.1.20,EVQL,0.94,0.88\n",
            encoding="utf-8",
        )
        res = generate_summary_plots(root, root / "plots", ("20",))
        assert res["datasets_plotted"] == ["20"], res
        if res["matplotlib_active"]:
            assert (root / "plots" / "scatter_plots_20.png").exists()
            assert (root / "plots" / "box_plots.png").exists()
    print("pipeline build_plots selftest: pass")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Repaired pLDDT/pTM scatter and box plotter.")
    ap.add_argument("--work-dir", type=Path, default=Path("."), help="Folder with folded_sequences_*.csv")
    ap.add_argument("--out-dir", type=Path, default=Path("plots"), help="Plot output folder")
    ap.add_argument("--datasets", nargs="+", default=["20", "30"], help="Dataset tags")
    ap.add_argument("--selftest", action="store_true", help="Run built-in check and exit")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    res = generate_summary_plots(args.work_dir, args.out_dir, tuple(args.datasets))
    print(f"build_plots: plotted={res['datasets_plotted']} -> {res['out_dir']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
