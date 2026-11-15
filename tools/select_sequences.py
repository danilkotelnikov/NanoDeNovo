#!/usr/bin/env python3
# version: w06 | date: 14.XI.2026
"""Select the best ProteinMPNN sequences for each backbone.

ProteinMPNN writes a score in every FASTA header. The score is the mean
negative log-likelihood of the sequence. A LOW score means a likely sequence.
The original pipeline sorted the score downwards and kept the highest values,
so it kept the least likely sequences of each batch. This tool keeps the
lowest values.

Usage:
    py -3 -I select_sequences.py --seqs-dir pmpnn_020/seqs --out-dir tfold_20 \
        --tag 20 --top 3 --log selections_20.csv
    py -3 -I select_sequences.py --selftest
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# A sampled record has "sample=" in its header. The first record of a
# ProteinMPNN file is the input sequence and has no "sample=".
SAMPLE_RE = re.compile(r"\bsample=")
# "\bscore=" does not match inside "global_score=", because "_" and "s" are
# both word characters, so no word boundary sits between them.
# The number pattern accepts scientific notation, for example 1e-3.
SCORE_RE = re.compile(r"\bscore=([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)")
ID_RE = re.compile(r"(\d{3,})")


def read_samples(path: Path) -> list[tuple[float, str]]:
    """Return (score, sequence) for every sampled record in one FASTA file."""
    records: list[tuple[float, str]] = []
    header: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        if header is None:
            return
        if not SAMPLE_RE.search(header):
            return
        match = SCORE_RE.search(header)
        if match is None:
            return
        sequence = "".join(chunks)
        if sequence:
            records.append((float(match.group(1)), sequence))

    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                flush()
                header = line
                chunks = []
            elif line:
                chunks.append(line)
    flush()
    return records


def select(records: list[tuple[float, str]], top: int) -> list[tuple[int, float, str]]:
    """Return the `top` most likely records, best first, as (rank, score, seq).

    Sampling with replacement can return the same sequence more than once.
    Folding the same sequence three times costs three runs and gives one
    answer, so only the first copy of a sequence is kept.
    """
    # Sort upwards: a low negative log-likelihood is a likely sequence.
    ordered = sorted(records, key=lambda item: item[0])
    chosen: list[tuple[int, float, str]] = []
    seen: set[str] = set()
    for score, seq in ordered:
        if seq in seen:
            continue
        seen.add(seq)
        chosen.append((len(chosen) + 1, score, seq))
        if len(chosen) == top:
            break
    return chosen


def backbone_id(path: Path) -> str:
    match = ID_RE.search(path.stem)
    return match.group(1) if match else path.stem


def run(seqs_dir: Path, out_dir: Path, tag: str, top: int, log_path: Path | None) -> int:
    files = sorted(seqs_dir.glob("*.fa")) + sorted(seqs_dir.glob("*.fasta"))
    if not files:
        print(f"No FASTA files in {seqs_dir}", file=sys.stderr)
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    written = 0
    skipped: list[str] = []

    for path in files:
        records = read_samples(path)
        if len(records) < top:
            skipped.append(f"{path.name}: {len(records)} samples")
            continue
        scores = [score for score, _ in records]
        best, worst = min(scores), max(scores)
        bid = backbone_id(path)

        for rank, score, sequence in select(records, top):
            name = f"{bid}_top.{rank}.{tag}"
            with open(out_dir / f"{name}.fa", "w", encoding="utf-8") as out:
                out.write(f">H\n{sequence}\n")
            written += 1
            rows.append(
                {
                    "name": name,
                    "backbone": bid,
                    "rank": rank,
                    "score": score,
                    "batch_samples": len(records),
                    "batch_best": best,
                    "batch_worst": worst,
                    # Share of the batch this sequence beats. The best sample of
                    # a batch of n scores (n-1)/n, so a value near 1 is good.
                    "batch_percentile": round(
                        sum(1 for s in scores if s > score) / len(scores), 4
                    ),
                    "length": len(sequence),
                }
            )

    if log_path is not None and rows:
        with open(log_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    print(f"Backbones read : {len(files)}")
    print(f"Sequences kept : {written} (top {top} of each batch)")
    if rows:
        mean_pct = sum(float(r["batch_percentile"]) for r in rows) / len(rows)
        mean_score = sum(float(r["score"]) for r in rows) / len(rows)
        print(f"Mean score     : {mean_score:.4f} (lower is better)")
        print(f"Mean rank share: {mean_pct:.4f} (near 1 = best of its batch)")
    if skipped:
        print(f"Skipped {len(skipped)} file(s) with too few samples; first: {skipped[0]}")
    if log_path is not None:
        print(f"Log            : {log_path}")
    return 0


def selftest() -> int:
    """One runnable check: the tool must keep the low scores, not the high ones."""
    import tempfile

    fasta = (
        ">gen_0001, score=3.0537, global_score=3.0537, designed_chains=['H']\n"
        "VVVV\n"
        ">T=0.15, sample=1, score=0.9000, seq_recovery=0.08\n"
        "AAAA\n"
        ">T=0.15, sample=2, score=1.5000, seq_recovery=0.08\n"
        "BBBB\n"
        ">T=0.15, sample=3, score=0.5000, seq_recovery=0.08\n"
        "CCCC\n"
        ">T=0.15, sample=4, score=1.2000, seq_recovery=0.08\n"
        "DDDD\n"
        # The same sequence sampled twice, with a worse score the second time.
        ">T=0.15, sample=5, score=0.6000, seq_recovery=0.08\n"
        "CCCC\n"
        # A score in scientific notation must still be read.
        ">T=0.15, sample=6, score=2.5e+00, seq_recovery=0.08\n"
        "EEEE\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        seqs = root / "seqs"
        seqs.mkdir()
        (seqs / "0001_20.t.0.15.fa").write_text(fasta, encoding="utf-8")

        records = read_samples(seqs / "0001_20.t.0.15.fa")
        # The input record (score 3.0537, sequence VVVV) must not be a candidate.
        assert len(records) == 6, records
        assert all(seq != "VVVV" for _, seq in records)
        # A score in scientific notation is read, and "global_score" is ignored.
        assert (2.5, "EEEE") in records, records

        chosen = select(records, 3)
        assert [seq for _, _, seq in chosen] == ["CCCC", "AAAA", "DDDD"], chosen
        assert [score for _, score, _ in chosen] == [0.5, 0.9, 1.2], chosen
        # The worst sample must never be kept when better ones exist.
        assert "BBBB" not in [seq for _, _, seq in chosen]
        # The repeated sequence must be kept once, not twice.
        assert [seq for _, _, seq in chosen].count("CCCC") == 1, chosen

        out = root / "out"
        assert run(seqs, out, "20", 3, root / "log.csv") == 0
        assert (out / "0001_top.1.20.fa").read_text(encoding="utf-8") == ">H\nCCCC\n"

        with open(root / "log.csv", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        # The best of six samples beats five of them.
        assert rows[0]["batch_percentile"] == "0.8333", rows[0]
        assert rows[0]["batch_samples"] == "6", rows[0]

    print("selftest: pass")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Select the best ProteinMPNN sequences for each backbone.")
    parser.add_argument("--seqs-dir", type=Path, help="Folder with ProteinMPNN FASTA files")
    parser.add_argument("--out-dir", type=Path, help="Folder for the selected FASTA files")
    parser.add_argument("--tag", default="20", help="Dataset tag used in the output names")
    parser.add_argument("--top", type=int, default=3, help="How many sequences to keep per backbone")
    parser.add_argument("--log", type=Path, default=None, help="CSV file recording every choice")
    parser.add_argument("--selftest", action="store_true", help="Run the built-in check and exit")
    args = parser.parse_args()

    if args.selftest:
        return selftest()
    if args.seqs_dir is None or args.out_dir is None:
        parser.error("--seqs-dir and --out-dir are required")
    return run(args.seqs_dir, args.out_dir, args.tag, args.top, args.log)


if __name__ == "__main__":
    raise SystemExit(main())
