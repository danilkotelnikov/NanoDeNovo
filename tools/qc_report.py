#!/usr/bin/env python3
# version: w06 | date: 14.XI.2026
"""Check the quality of a NanoDeNovo result set.

The tool answers three questions about one dataset.

1. Integrity. Are designs missing? Are rows repeated?
2. False positives. How many designs pass the confidence filter but fail an
   independent structural check?
3. Thresholds. Which confidence thresholds remove the most false positives
   for the smallest loss of designs?

Reference standard
------------------
There is no experimental measurement here, so the tool uses structural
self-consistency as the reference. A design is self-consistent when its
refolded model matches the backbone it was designed on: TM-score at or above
`--tm-min` and RMSD at or below `--rmsd-max`.

The screen under test is the confidence filter: pLDDT above `--plddt-min` and
pTM above `--ptm-min`.

A design that passes the confidence filter and fails the structural check is
counted as a failure of the filter. The share of such designs among the
designs that passed is a false discovery rate, not a false positive rate: a
false positive rate needs the rejected designs too, and those were never
measured. The tool therefore reports the self-consistency failure rate among
the designs the filter kept, and says plainly that misses cannot be counted.

The structural bar matters. A nanobody has a long and flexible third loop, so a
strict bar fails designs whose framework is correct and whose loop moved. The
tool reports several bars side by side, so a reader can see how much of the
failure rate comes from the choice of bar.

Usage:
    py -3 -I qc_report.py --root E:\\nanodenovo --dataset 20 --out report_20.md
    py -3 -I qc_report.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import statistics as st
import sys
from pathlib import Path

ALIGN_DIRS = (
    "usalign_rmsd_2.0_tm_0.90",
    "usalign_rmsd_1.5_tm_0.95",
    "usalign_rmsd_1.0_tm_0.95",
)


def read_folded(path: Path) -> list[dict]:
    """Read Name, Sequence, pLDDT, pTM for every folded design."""
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                rows.append(
                    {
                        "name": row["Name"].strip(),
                        "sequence": row.get("Sequence", "").strip(),
                        "plddt": float(row["pLDDT"]),
                        "ptm": float(row["pTM"]),
                    }
                )
            except (KeyError, TypeError, ValueError):
                continue
    return rows


def read_alignment(path: Path, tm_column: str = "tm_design") -> dict[str, dict]:
    """Read TM-score and RMSD per design name.

    Two table layouts are accepted. The original run wrote `Confident Protein`,
    `TM-score` and `RMSD`. `align_all.py` writes `name`, `tm_design`,
    `tm_scaffold` and `rmsd`; `tm_column` picks which of the two scores is `tm`.
    """
    out: dict[str, dict] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            name = (row.get("name") or row.get("Confident Protein") or "").strip()
            if not name:
                continue
            if name.endswith(".pdb"):
                name = name[:-4]
            try:
                if "name" in row:
                    out[name] = {"tm": float(row[tm_column]), "rmsd": float(row["rmsd"])}
                else:
                    out[name] = {"tm": float(row["TM-score"]), "rmsd": float(row["RMSD"])}
            except (KeyError, TypeError, ValueError):
                continue
    return out


BARS = (
    (0.95, 1.0, "strict"),
    (0.95, 2.0, "strict match, relaxed RMSD"),
    (0.90, 2.0, "common design bar"),
    (0.85, 3.0, "loose"),
)


def ratio(num: int, den: int) -> float | None:
    return num / den if den else None


def confusion(folded: list[dict], align: dict[str, dict], plddt_min: float,
              ptm_min: float, tm_bar: float, rmsd_bar: float) -> dict:
    """Count the filter against the reference over every measured design."""
    tp = fp = fn = tn = 0
    for row in folded:
        measure = align.get(row["name"])
        if measure is None:
            continue
        kept = row["plddt"] > plddt_min and row["ptm"] > ptm_min
        good = measure["tm"] >= tm_bar and measure["rmsd"] <= rmsd_bar
        if kept and good:
            tp += 1
        elif kept:
            fp += 1
        elif good:
            fn += 1
        else:
            tn += 1
    return {
        "tp": tp, "fp": fp, "fn": fn, "tn": tn,
        "fdr": ratio(fp, tp + fp),      # failures among the kept
        "miss": ratio(fn, tp + fn),     # good designs the filter threw away
        "fpr": ratio(fp, fp + tn),      # bad designs the filter kept
        "recall": ratio(tp, tp + fn),
        "f1": ratio(2 * tp, 2 * tp + fp + fn),
        # How common correct designs are. A filter that keeps everything scores
        # an F1 of 2p/(1+p) at prevalence p, so read F1 with this beside it.
        "prevalence": ratio(tp + fn, tp + fp + fn + tn),
        # Youden's J: recall minus the false positive rate. It is 0 for a filter
        # that keeps everything or nothing, whatever the prevalence.
        "j": (ratio(tp, tp + fn) or 0.0) - (ratio(fp, fp + tn) or 0.0)
        if (tp + fn) and (fp + tn) else None,
    }


def auc(scores_good: list[float], scores_bad: list[float]) -> float | None:
    """Chance that a good design scores above a bad one (Mann-Whitney)."""
    if not scores_good or not scores_bad:
        return None
    ordered = sorted(scores_bad)
    from bisect import bisect_left, bisect_right
    wins = 0.0
    for value in scores_good:
        wins += bisect_left(ordered, value) + 0.5 * (bisect_right(ordered, value) - bisect_left(ordered, value))
    return wins / (len(scores_good) * len(scores_bad))


def describe(values: list[float]) -> dict:
    if not values:
        return {"n": 0}
    return {
        "n": len(values),
        "mean": st.mean(values),
        "median": st.median(values),
        "std": st.pstdev(values) if len(values) > 1 else 0.0,
        "min": min(values),
        "max": max(values),
    }


def find_gaps(names: list[str], top: int, expect: int | None = None) -> list[str]:
    """Report backbones that do not have the expected number of designs.

    `expect` is the number of backbones the run should cover. Without it the
    check cannot see a run that stopped early, because it has nothing to
    compare the highest number against.
    """
    counts: dict[str, int] = {}
    for name in names:
        counts[name.split("_")[0]] = counts.get(name.split("_")[0], 0) + 1
    if not counts:
        return []
    numbers = sorted(int(k) for k in counts if k.isdigit())
    gaps = []
    if numbers:
        last = expect if expect is not None else max(numbers)
        expected = set(range(1 if expect is not None else min(numbers), last + 1))
        missing = sorted(expected - set(numbers))
        if missing:
            gaps.append(
                f"{len(missing)} backbone(s) absent up to {last:04d}; "
                f"absent: {', '.join(f'{m:04d}' for m in missing[:5])}"
                + (" ..." if len(missing) > 5 else "")
            )
        if expect is not None and len(names) != expect * top:
            gaps.append(
                f"{expect * top} designs expected, {len(names)} present "
                f"({expect * top - len(names)} absent)"
            )
    short = sorted(k for k, v in counts.items() if v != top)
    if short:
        gaps.append(
            f"{len(short)} backbone(s) do not have {top} designs; "
            f"first: {short[0]} has {counts[short[0]]}"
        )
    return gaps


def read_loops(path: Path) -> dict[str, dict]:
    """Read the third-loop table written by cdr3_check.py. Unplaced designs are left out."""
    out: dict[str, dict] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                out[row["name"].strip()] = {
                    "fw": float(row["fw_rmsd"]), "cdr3": float(row["cdr3_rmsd"]),
                    "len": int(row["cdr3_len"]),
                }
            except (KeyError, TypeError, ValueError):
                continue  # an empty field means the loop could not be placed
    return out


LOOP_BARS = (2.0, 3.0, 4.0)


def loop_analysis(folded: list[dict], loops: dict[str, dict], plddt_min: float,
                  ptm_min: float) -> dict:
    """Score the filter against the third loop, the part a whole-chain match hides."""
    kept_all = [r for r in folded if r["plddt"] > plddt_min and r["ptm"] > ptm_min]
    placed = [r for r in folded if r["name"] in loops]
    kept = [r for r in placed if r["plddt"] > plddt_min and r["ptm"] > ptm_min]
    rejected = [r for r in placed if r not in kept]

    def med(rows: list[dict], key: str) -> float | None:
        values = [loops[r["name"]][key] for r in rows]
        return st.median(values) if values else None

    bars = []
    for bar in LOOP_BARS:
        tp = fp = fn = tn = 0
        for r in placed:
            is_kept = r["plddt"] > plddt_min and r["ptm"] > ptm_min
            good = loops[r["name"]]["cdr3"] <= bar
            if is_kept and good:
                tp += 1
            elif is_kept:
                fp += 1
            elif good:
                fn += 1
            else:
                tn += 1
        bars.append({
            "bar": bar, "tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "fdr": ratio(fp, tp + fp), "miss": ratio(fn, tp + fn), "fpr": ratio(fp, fp + tn),
            "prevalence": ratio(tp + fn, len(placed)),
        })
    mid = [r for r in placed]
    good3 = [r for r in mid if loops[r["name"]]["cdr3"] <= 3.0]
    bad3 = [r for r in mid if loops[r["name"]]["cdr3"] > 3.0]
    return {
        "folded": len(folded), "placed": len(placed),
        "kept_total": len(kept_all), "kept_placed": len(kept),
        "rejected_placed": len(rejected),
        "rejected_total": len(folded) - len(kept_all),
        "fw_kept": med(kept, "fw"), "fw_rejected": med(rejected, "fw"),
        "cdr3_kept": med(kept, "cdr3"), "cdr3_rejected": med(rejected, "cdr3"),
        "len_kept": med(kept, "len"), "len_rejected": med(rejected, "len"),
        "bars": bars,
        "plddt_auc": auc([r["plddt"] for r in good3], [r["plddt"] for r in bad3]),
        "ptm_auc": auc([r["ptm"] for r in good3], [r["ptm"] for r in bad3]),
    }


def analyse(folded: list[dict], align: dict[str, dict], plddt_min: float,
            ptm_min: float, tm_min: float, rmsd_max: float, top: int,
            expect: int | None = None) -> dict:
    names = [row["name"] for row in folded]
    duplicates = len(names) - len(set(names))

    confident = [r for r in folded if r["plddt"] > plddt_min and r["ptm"] > ptm_min]
    rejected = [r for r in folded if not (r["plddt"] > plddt_min and r["ptm"] > ptm_min)]

    # Designs that have an independent structural measurement.
    measured = [r for r in confident if r["name"] in align]
    tp, fp = [], []
    for row in measured:
        a = align[row["name"]]
        if a["tm"] >= tm_min and a["rmsd"] <= rmsd_max:
            tp.append(row)
        else:
            fp.append(row)

    rejected_measured = [r for r in rejected if r["name"] in align]

    # Misses can be counted only when every rejected design has a measurement.
    full = bool(rejected) and len(rejected_measured) == len(rejected)
    everything = [r for r in folded if r["name"] in align]
    confusions = [
        {"label": label, "tm_min": tm_bar, "rmsd_max": rmsd_bar,
         **confusion(folded, align, plddt_min, ptm_min, tm_bar, rmsd_bar)}
        for tm_bar, rmsd_bar, label in BARS
    ] if everything else []

    # Does each score separate designs that pass the reference from those that fail?
    good = [r for r in everything
            if align[r["name"]]["tm"] >= tm_min and align[r["name"]]["rmsd"] <= rmsd_max]
    bad_all = [r for r in everything if r not in good]
    separation = {
        "plddt_auc": auc([r["plddt"] for r in good], [r["plddt"] for r in bad_all]),
        "ptm_auc": auc([r["ptm"] for r in good], [r["ptm"] for r in bad_all]),
    }

    # How much of the failure rate comes from the choice of structural bar.
    bars = []
    for tm_bar, rmsd_bar, label in BARS:
        bad = [
            r for r in measured
            if not (align[r["name"]]["tm"] >= tm_bar
                    and align[r["name"]]["rmsd"] <= rmsd_bar)
        ]
        bars.append(
            {
                "label": label,
                "tm_min": tm_bar,
                "rmsd_max": rmsd_bar,
                "failures": len(bad),
                "share": len(bad) / len(measured) if measured else None,
            }
        )

    # Threshold sweep. With every design measured the grid reaches below the
    # old cut and reports misses as well; otherwise it covers only the survivors.
    sweep = []
    p_grid = (0.80, 0.85, 0.88, 0.90, 0.92, 0.94, 0.96) if full else (plddt_min, 0.92, 0.94, 0.95, 0.96)
    t_grid = (0.70, 0.75, 0.80, 0.85, 0.87, 0.89) if full else (ptm_min, 0.87, 0.89, 0.90)
    pool = everything if full else measured
    for p_cut in p_grid:
        for t_cut in t_grid:
            kept = [r for r in pool if r["plddt"] > p_cut and r["ptm"] > t_cut]
            if not kept:
                continue
            bad = [r for r in kept
                   if not (align[r["name"]]["tm"] >= tm_min
                           and align[r["name"]]["rmsd"] <= rmsd_max)]
            row = {
                "plddt_min": p_cut,
                "ptm_min": t_cut,
                "kept": len(kept),
                "false_positives": len(bad),
                "fp_share": len(bad) / len(kept),
            }
            if full:
                c = confusion(folded, align, p_cut, t_cut, tm_min, rmsd_max)
                row.update({"false_negatives": c["fn"], "miss": c["miss"], "j": c["j"]})
            sweep.append(row)

    return {
        "folded": len(folded),
        "duplicates": duplicates,
        "gaps": find_gaps(names, top, expect),
        "confident": len(confident),
        "rejected": len(rejected),
        "measured": len(measured),
        "rejected_measured": len(rejected_measured),
        "full_coverage": full,
        "confusion": confusions,
        "separation": separation,
        "true_positives": len(tp),
        "false_positives": len(fp),
        "fp_share": (len(fp) / len(measured)) if measured else None,
        "bars": bars,
        "fp_examples": sorted(
            ({"name": r["name"], "plddt": r["plddt"], "ptm": r["ptm"],
              **align[r["name"]]} for r in fp),
            key=lambda d: d["tm"],
        )[:10],
        "plddt_all": describe([r["plddt"] for r in folded]),
        "ptm_all": describe([r["ptm"] for r in folded]),
        "plddt_confident": describe([r["plddt"] for r in confident]),
        "ptm_confident": describe([r["ptm"] for r in confident]),
        "tm_measured": describe([align[r["name"]]["tm"] for r in measured]),
        "rmsd_measured": describe([align[r["name"]]["rmsd"] for r in measured]),
        "sweep": sorted(sweep, key=(lambda d: -(d["j"] or 0.0)) if full
                        else (lambda d: (d["fp_share"], -d["kept"]))),
    }


def fmt(stats: dict) -> str:
    if not stats.get("n"):
        return "no data"
    return (
        f"n={stats['n']}, mean={stats['mean']:.4f}, median={stats['median']:.4f}, "
        f"sd={stats['std']:.4f}, min={stats['min']:.4f}, max={stats['max']:.4f}"
    )


def render(dataset: str, res: dict, args) -> str:
    lines = [
        f"# Quality report, dataset {dataset}",
        "",
        "## What this report measures",
        "",
        f"The screen under test keeps a design when pLDDT is above {args.plddt_min} "
        f"and pTM is above {args.ptm_min}.",
        f"The reference standard calls a design correct when its refolded model "
        f"matches its backbone: TM-score at or above {args.tm_min} and RMSD at or "
        f"below {args.rmsd_max} angstrom.",
        "",
        "The reference is a calculation, not a laboratory measurement. So the "
        "share below is a self-consistency failure rate among the designs the "
        "filter kept. It is a false discovery rate. It is not a false positive "
        "rate, because that would need the rejected designs as well.",
        "",
        "## Integrity",
        "",
        f"- Folded designs in the table: {res['folded']}",
        f"- Repeated rows: {res['duplicates']}",
    ]
    for gap in res["gaps"] or ["No absent backbones found."]:
        lines.append(f"- {gap}")

    lines += [
        "",
        "## Confidence filter",
        "",
        f"- Designs that pass: {res['confident']}",
        f"- Designs that fail: {res['rejected']}",
        f"- pLDDT, all designs: {fmt(res['plddt_all'])}",
        f"- pTM, all designs: {fmt(res['ptm_all'])}",
        f"- pLDDT, designs that pass: {fmt(res['plddt_confident'])}",
        f"- pTM, designs that pass: {fmt(res['ptm_confident'])}",
        "",
        "## Designs that pass the filter but fail the structure check",
        "",
        f"- Designs that pass the filter and also have a structural measurement: "
        f"{res['measured']}",
        f"- Correct by the reference standard: {res['true_positives']}",
        f"- Failures: {res['false_positives']}",
    ]
    if res["fp_share"] is not None:
        lines.append(
            f"- Failure share of the measured set: {res['fp_share']:.3f}"
        )
    lines += [
        f"- TM-score of the measured set: {fmt(res['tm_measured'])}",
        f"- RMSD of the measured set: {fmt(res['rmsd_measured'])}",
    ]

    if res["bars"]:
        lines += [
            "",
            "### How much the structural bar decides",
            "",
            "A nanobody has a long and flexible third loop. A strict bar fails a "
            "design whose framework is right and whose loop moved, so the bar is "
            "shown at four settings.",
            "",
            "| Bar | TM-score at or above | RMSD at or below | Failures | Share |",
            "|---|---|---|---|---|",
        ]
        for bar in res["bars"]:
            share = f"{bar['share']:.3f}" if bar["share"] is not None else "n/a"
            lines.append(
                f"| {bar['label']} | {bar['tm_min']} | {bar['rmsd_max']} | "
                f"{bar['failures']} | {share} |"
            )

    if res["fp_examples"]:
        lines += [
            "",
            "The ten weakest failures at the chosen bar:",
            "",
            "| Design | pLDDT | pTM | TM-score | RMSD |",
            "|---|---|---|---|---|",
        ]
        for row in res["fp_examples"]:
            lines.append(
                f"| {row['name']} | {row['plddt']:.4f} | {row['ptm']:.4f} | "
                f"{row['tm']:.4f} | {row['rmsd']:.2f} |"
            )

    lines += [
        "",
        "## False negatives",
        "",
    ]
    if res["full_coverage"]:
        lines += [
            f"Every one of the {res['rejected']} rejected designs has a structural "
            "measurement, so the filter can be scored in both directions. A miss is "
            "a design that the reference calls correct and the filter rejected.",
            "",
            "| Bar | Kept and correct | Kept and wrong | Missed | Rejected and wrong "
            "| Failures among kept | Miss rate | Wrong designs kept | Correct share of all |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        for c in res["confusion"]:
            def pc(v):
                return f"{v:.3f}" if v is not None else "n/a"
            lines.append(
                f"| {c['label']} ({c['tm_min']}/{c['rmsd_max']}) | {c['tp']} | {c['fp']} | "
                f"{c['fn']} | {c['tn']} | {pc(c['fdr'])} | {pc(c['miss'])} | {pc(c['fpr'])} | "
                f"{pc(c['prevalence'])} |"
            )
        lines += [
            "",
            "Failures among kept is the false discovery rate. Miss rate is the share of "
            "correct designs that the filter threw away. Wrong designs kept is the false "
            "positive rate. All three are measured against structural self-consistency, "
            "which is a calculation and not a laboratory result.",
        ]
        sep = res["separation"]
        if sep["plddt_auc"] is not None:
            lines += [
                "",
                f"Separation at the chosen bar: a correct design scores higher pLDDT than a "
                f"wrong one with probability {sep['plddt_auc']:.3f}, and higher pTM with "
                f"probability {sep['ptm_auc']:.3f}. A value of 0.5 means the score tells "
                "nothing; 1.0 means it separates them perfectly. The whole chain is "
                "dominated by the conserved framework, so this value is generous. "
                "Read the third-loop check for the harder test.",
            ]
    elif res["rejected_measured"]:
        lines.append(
            f"Structural measurements exist for {res['rejected_measured']} rejected "
            "design(s). A full false negative count needs a measurement for every "
            "rejected design."
        )
    else:
        lines += [
            f"False negatives cannot be counted. The structural check was run only on "
            f"designs that had already passed the confidence filter, so the "
            f"{res['rejected']} rejected designs have no measurement.",
            "",
            "To close this gap, align the rejected designs against their backbones. "
            "The alignment runs on a processor and needs no graphics card.",
        ]

    if res["sweep"]:
        lines += [
            "",
            "## Threshold choice",
            "",
            "Each row shows how many measured designs a threshold pair keeps, and how "
            "many of those are false positives.",
        ]
        if res["full_coverage"]:
            lines += [
                "Every threshold pair is scored on all designs. Rows are ordered by Youden's "
                "J, the share of correct designs kept minus the share of wrong designs kept. "
                "J is 0 for a filter that keeps everything, however common correct designs are.",
                "",
                "| pLDDT above | pTM above | Kept | False positives | Missed | Failure share | Miss rate | J |",
                "|---|---|---|---|---|---|---|---|",
            ]
            for row in res["sweep"][:12]:
                lines.append(
                    f"| {row['plddt_min']} | {row['ptm_min']} | {row['kept']} | "
                    f"{row['false_positives']} | {row['false_negatives']} | "
                    f"{row['fp_share']:.3f} | {row['miss']:.3f} | {row['j']:.3f} |"
                )
            study = [r for r in res["sweep"]
                     if r["plddt_min"] == args.plddt_min and r["ptm_min"] == args.ptm_min]
            if study:
                s = study[0]
                lines += [
                    "",
                    f"The study's own cut (pLDDT above {args.plddt_min}, pTM above "
                    f"{args.ptm_min}) scores J = {s['j']:.3f}. It keeps {s['kept']} designs, "
                    f"keeps {s['false_positives']} wrong ones and misses {s['false_negatives']}.",
                ]
            best = res["sweep"][0]
            at_edge = (best["plddt_min"] == min(r["plddt_min"] for r in res["sweep"])
                       or best["ptm_min"] == min(r["ptm_min"] for r in res["sweep"]))
            if at_edge:
                lines += [
                    "",
                    "The best row sits on the edge of the grid. A wider grid could do "
                    "better, so this table shows a direction and not an optimum.",
                ]
            lines += [
                "",
                "J ignores what the later stages cost. A looser filter sends more "
                "designs to docking and simulation. Pick a row with that cost in mind.",
            ]
        else:
            lines += [
                "",
                "| pLDDT above | pTM above | Kept | False positives | Share |",
                "|---|---|---|---|---|",
            ]
            for row in res["sweep"][:12]:
                lines.append(
                    f"| {row['plddt_min']} | {row['ptm_min']} | {row['kept']} | "
                    f"{row['false_positives']} | {row['fp_share']:.3f} |"
                )

    lines.append("")
    return "\n".join(lines)


def render_loops(res: dict) -> list[str]:
    """Text for the third-loop check."""
    def n(v, digits=2):
        return f"{v:.{digits}f}" if v is not None else "n/a"

    lines = [
        "",
        "## Third-loop check",
        "",
        "A whole-chain match is dominated by the conserved framework, so most designs "
        "pass it. This check fits the framework alone and measures how far the third "
        "loop (CDR3) moved. The loop is found by a sequence heuristic, so some designs "
        "cannot be placed.",
        "",
        f"- Designs placed: {res['placed']} of {res['folded']}",
        f"- Kept designs placed: {res['kept_placed']} of {res['kept_total']}",
        f"- Rejected designs placed: {res['rejected_placed']} of {res['rejected_total']}",
        "",
        "| Measure | Kept | Rejected |",
        "|---|---|---|",
        f"| Median framework RMSD (angstrom) | {n(res['fw_kept'])} | {n(res['fw_rejected'])} |",
        f"| Median CDR3 RMSD (angstrom) | {n(res['cdr3_kept'])} | {n(res['cdr3_rejected'])} |",
        f"| Median CDR3 length (residues) | {n(res['len_kept'], 0)} | {n(res['len_rejected'], 0)} |",
        "",
        "Loop length matters. If kept designs have shorter loops than rejected ones, "
        "the filter favours short loops, and the long loops that many antigens need "
        "are dropped before docking.",
        "",
        "| CDR3 RMSD at or below (angstrom) | Kept and correct | Kept and wrong | Missed "
        "| Rejected and wrong | Failures among kept | Miss rate | Wrong designs kept | Correct share |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for b in res["bars"]:
        lines.append(
            f"| {b['bar']} | {b['tp']} | {b['fp']} | {b['fn']} | {b['tn']} | {n(b['fdr'], 3)} "
            f"| {n(b['miss'], 3)} | {n(b['fpr'], 3)} | {n(b['prevalence'], 3)} |"
        )
    lines += [
        "",
        f"Separation at 3 angstrom: a design with a good loop has a higher pLDDT than one "
        f"with a bad loop with probability {n(res['plddt_auc'], 3)}, and a higher pTM with "
        f"probability {n(res['ptm_auc'], 3)}.",
        "",
        "The placed designs are mostly the kept ones. Designs that cannot be placed are "
        "mostly rejected, so the miss counts here are a floor for the rejected side.",
    ]
    return lines


def selftest() -> int:
    """One runnable check on a small synthetic dataset."""
    import tempfile

    folded = [
        # name, plddt, ptm
        ("0001_top.1.20", 0.95, 0.90),  # passes, self-consistent -> true positive
        ("0001_top.2.20", 0.96, 0.91),  # passes, poor TM -> false positive
        ("0001_top.3.20", 0.80, 0.70),  # fails the filter
        ("0002_top.1.20", 0.97, 0.92),  # passes, poor RMSD -> false positive
    ]
    align = {
        "0001_top.1.20": (0.97, 0.80),
        "0001_top.2.20": (0.70, 0.90),
        "0002_top.1.20": (0.98, 3.50),
    }
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        f_path = root / "folded.csv"
        with open(f_path, "w", newline="", encoding="utf-8") as handle:
            w = csv.writer(handle)
            w.writerow(["Name", "Sequence", "pLDDT", "pTM"])
            for name, plddt, ptm in folded:
                w.writerow([name, "AAAA", plddt, ptm])

        a_path = root / "align.csv"
        with open(a_path, "w", newline="", encoding="utf-8") as handle:
            w = csv.writer(handle)
            w.writerow(["Confident Protein", "Scaffold ID", "TM-score", "RMSD"])
            for name, (tm, rmsd) in align.items():
                w.writerow([name + ".pdb", "gen.pdb", tm, rmsd])

        rows = read_folded(f_path)
        assert len(rows) == 4, rows
        got = read_alignment(a_path)
        assert set(got) == set(align), got

        res = analyse(rows, got, 0.90, 0.85, 0.95, 1.0, top=3)
        assert res["confident"] == 3, res["confident"]
        assert res["rejected"] == 1, res["rejected"]
        assert res["measured"] == 3, res["measured"]
        assert res["true_positives"] == 1, res["true_positives"]
        assert res["false_positives"] == 2, res["false_positives"]
        assert abs(res["fp_share"] - 2 / 3) < 1e-9
        # The rejected design has no measurement, so no false negative count.
        assert res["rejected_measured"] == 0
        # 0002 has only one design, so the gap check must notice.
        assert any("do not have 3 designs" in g for g in res["gaps"]), res["gaps"]

        # With the expected count given, a run that stopped early is noticed.
        res3 = analyse(rows, got, 0.90, 0.85, 0.95, 1.0, top=3, expect=3)
        assert any("0003" in g for g in res3["gaps"]), res3["gaps"]
        assert any("9 designs expected" in g for g in res3["gaps"]), res3["gaps"]
        # A stricter threshold must never keep more designs.
        kept = [r["kept"] for r in res["sweep"]]
        assert max(kept) <= 3, kept

        # Full coverage: the new table layout, with every design measured.
        full_path = root / "align_all.csv"
        with open(full_path, "w", newline="", encoding="utf-8") as handle:
            w = csv.writer(handle)
            w.writerow(["name", "scaffold", "len_design", "len_scaffold",
                        "aligned_length", "rmsd", "tm_design", "tm_scaffold"])
            w.writerow(["0001_top.1.20", "g", 117, 116, 110, 0.80, 0.97, 0.50])
            w.writerow(["0001_top.2.20", "g", 117, 116, 110, 0.90, 0.70, 0.99])
            # Rejected by the filter, yet correct by the reference: a miss.
            w.writerow(["0001_top.3.20", "g", 117, 116, 110, 0.60, 0.98, 0.98])
            w.writerow(["0002_top.1.20", "g", 117, 116, 110, 3.50, 0.98, 0.98])
        # The column choice must change which score is read.
        by_design = read_alignment(full_path, "tm_design")
        by_scaffold = read_alignment(full_path, "tm_scaffold")
        assert by_design["0001_top.1.20"]["tm"] == 0.97, by_design
        assert by_scaffold["0001_top.1.20"]["tm"] == 0.50, by_scaffold

        res4 = analyse(rows, by_design, 0.90, 0.85, 0.95, 1.0, top=3)
        assert res4["full_coverage"] is True
        c = confusion(rows, by_design, 0.90, 0.85, 0.95, 1.0)
        # 0001.1 kept+good, 0001.2 kept+bad, 0001.3 rejected+good, 0002.1 kept+bad.
        assert (c["tp"], c["fp"], c["fn"], c["tn"]) == (1, 2, 1, 0), c
        assert c["miss"] == 0.5 and abs(c["fdr"] - 2 / 3) < 1e-9, c
        # With every design measured, a miss must show in the sweep too.
        assert any(r["false_negatives"] >= 1 for r in res4["sweep"]), res4["sweep"]

        # AUC: perfect separation is 1.0, no separation is 0.5, empty is undefined.
        assert auc([0.9, 0.8], [0.1, 0.2]) == 1.0
        assert auc([0.5, 0.5], [0.5, 0.5]) == 0.5
        assert auc([], [0.5]) is None

        # Youden's J is 0 for a filter that keeps everything, even when 3 in 4 are correct.
        keep_all = [{"name": f"n{i}", "plddt": 0.99, "ptm": 0.99} for i in range(4)]
        table = {f"n{i}": {"tm": 0.99 if i < 3 else 0.5, "rmsd": 1.0} for i in range(4)}
        c_all = confusion(keep_all, table, 0.90, 0.85, 0.90, 2.0)
        assert c_all["prevalence"] == 0.75 and c_all["j"] == 0.0, c_all
        # F1 of that useless filter is 0.857, which is why J orders the sweep.
        assert abs(c_all["f1"] - 6 / 7) < 1e-9, c_all

        # Third loop: three designs, two placed, one rejected with a bad loop.
        folded3 = [
            {"name": "a", "plddt": 0.95, "ptm": 0.90},  # kept, good loop
            {"name": "b", "plddt": 0.95, "ptm": 0.90},  # kept, bad loop
            {"name": "c", "plddt": 0.70, "ptm": 0.60},  # rejected, good loop -> a miss
            {"name": "d", "plddt": 0.70, "ptm": 0.60},  # rejected, not placed
        ]
        loops3 = {"a": {"fw": 0.5, "cdr3": 1.0, "len": 9},
                  "b": {"fw": 0.6, "cdr3": 5.0, "len": 12},
                  "c": {"fw": 0.7, "cdr3": 2.0, "len": 10}}
        la = loop_analysis(folded3, loops3, 0.90, 0.85)
        assert la["placed"] == 3 and la["rejected_total"] == 2 and la["rejected_placed"] == 1, la
        three = [b for b in la["bars"] if b["bar"] == 3.0][0]
        assert (three["tp"], three["fp"], three["fn"], three["tn"]) == (1, 1, 1, 0), three
        assert la["cdr3_kept"] == 3.0 and la["len_kept"] == 10.5, la
        assert render_loops(la), "the loop section must render"

    print("selftest: pass")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check the quality of a NanoDeNovo result set.")
    parser.add_argument("--root", type=Path, help="Folder with the result tables")
    parser.add_argument("--dataset", default="20", help="Dataset tag, for example 20")
    parser.add_argument("--align-csv", type=Path, default=None,
                        help="Alignment table; found under --root when not given")
    parser.add_argument("--cdr3-csv", type=Path, default=None,
                        help="Third-loop table from cdr3_check.py; adds the loop section")
    parser.add_argument("--tm-column", default="tm_design",
                        choices=("tm_design", "tm_scaffold"),
                        help="Which match score to use from an align_all.py table")
    parser.add_argument("--plddt-min", type=float, default=0.90)
    parser.add_argument("--ptm-min", type=float, default=0.85)
    parser.add_argument("--tm-min", type=float, default=0.90,
                        help="Reference bar: the common design bar, not the strict one")
    parser.add_argument("--rmsd-max", type=float, default=2.0)
    parser.add_argument("--top", type=int, default=3,
                        help="Designs expected per backbone")
    parser.add_argument("--expect-backbones", type=int, default=None,
                        help="How many backbones the run should cover, "
                             "so a run that stopped early is noticed")
    parser.add_argument("--out", type=Path, default=None, help="Write the report here")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()

    if args.selftest:
        return selftest()
    if args.root is None:
        parser.error("--root is required")

    folded_path = args.root / f"folded_sequences_{args.dataset}.csv"
    if not folded_path.exists():
        print(f"Not found: {folded_path}", file=sys.stderr)
        return 1

    align_path = args.align_csv
    if align_path is None:
        for name in ALIGN_DIRS:
            candidate = args.root / name / f"usalign_{args.dataset}.csv"
            if candidate.exists():
                align_path = candidate
                break
    if align_path is None or not align_path.exists():
        print("No alignment table found; the structural check cannot run.",
              file=sys.stderr)
        return 1

    folded = read_folded(folded_path)
    align = read_alignment(align_path, args.tm_column)
    res = analyse(folded, align, args.plddt_min, args.ptm_min, args.tm_min,
                  args.rmsd_max, args.top, args.expect_backbones)
    report = render(args.dataset, res, args)
    if args.cdr3_csv is not None:
        loops = read_loops(args.cdr3_csv)
        report += "\n".join(render_loops(loop_analysis(
            folded, loops, args.plddt_min, args.ptm_min))) + "\n"

    if args.out:
        args.out.write_text(report, encoding="utf-8")
        print(f"Report written: {args.out}")
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
