#!/usr/bin/env python3
"""Repaired ProteinMPNN sequence selection shared by process_020.py and process_030.py.

Repairs applied:
- Ascending negative log-likelihood (NLL) sort with deterministic `(score, sample_idx, seq)` tie-break.
- Safe multiline FASTA parser that excludes the native seed record (`sample=` missing).
- Exact amino-acid sequence deduplication before taking top-K candidates.
- Consistent `tfold_{tag}/` directory naming matching `tfold_ab_{tag}.py`.
- Atomic log and CSV writes via temporary file + `os.replace`.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import tempfile
from pathlib import Path

SAMPLE_RE = re.compile(r"\bsample=(\d+)\b")
SCORE_RE = re.compile(r"\bscore=([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)")
ID_RE = re.compile(r"(\d{3,})")


def read_sampled_records(path: Path) -> list[tuple[float, int, str]]:
    """Parse multiline ProteinMPNN FASTA records containing `sample=` and `score=`."""
    records: list[tuple[float, int, str]] = []
    header: str | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        if header is None:
            return
        sm = SAMPLE_RE.search(header)
        if sm is None:
            return
        scm = SCORE_RE.search(header)
        if scm is None:
            return
        seq = "".join(seq_parts).strip()
        if seq:
            records.append((float(scm.group(1)), int(sm.group(1)), seq))

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith(">"):
                flush()
                header = line
                seq_parts = []
            elif line:
                seq_parts.append(line)
    flush()
    return records


def select_top_unique(records: list[tuple[float, int, str]], top: int) -> list[tuple[int, float, int, str]]:
    """Return lowest-NLL unique sequences with deterministic `(score, sample_id, seq)` tie-break."""
    ordered = sorted(records, key=lambda item: (item[0], item[1], item[2]))
    chosen: list[tuple[int, float, int, str]] = []
    seen: set[str] = set()
    for score, sample_id, seq in ordered:
        if seq in seen:
            continue
        seen.add(seq)
        chosen.append((len(chosen) + 1, score, sample_id, seq))
        if len(chosen) == top:
            break
    return chosen


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".tmp_", suffix=path.suffix or ".txt", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        os.replace(tmp_name, path)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def atomic_write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
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


def process_tag(
    tag: str,
    seqs_dir: Path,
    out_dir: Path,
    log_txt: Path,
    log_csv: Path | None = None,
    top: int = 3,
    start: int = 1,
    end: int = 1000,
) -> dict:
    """Process one dataset tag (`20` or `30`) across inclusive `[start, end]`."""
    out_dir.mkdir(parents=True, exist_ok=True)
    txt_lines: list[str] = []
    csv_rows: list[dict] = []
    backbones = 0
    written = 0

    for idx in range(start, end + 1):
        bid = f"{idx:04d}"
        candidates = [
            seqs_dir / f"{bid}_{tag}.t.0.15.fa",
            seqs_dir / f"gen_{bid}.fa",
            seqs_dir / f"{bid}_{tag}.fa",
        ]
        in_path = next((p for p in candidates if p.exists()), None)
        if in_path is None:
            continue
        backbones += 1
        records = read_sampled_records(in_path)
        chosen = select_top_unique(records, top)
        for rank, score, sample_id, seq in chosen:
            out_fa = out_dir / f"{bid}_top.{rank}.{tag}.fa"
            atomic_write_text(out_fa, f">H\n{seq}\n")
            written += 1
            txt_lines.append(f">{bid}_top_{rank}_sample_{sample_id}_score_{score:.6f}\n{seq}\n")
            csv_rows.append(
                {
                    "name": f"{bid}_top.{rank}.{tag}",
                    "backbone": bid,
                    "rank": rank,
                    "sample": sample_id,
                    "score": f"{score:.6f}",
                    "length": len(seq),
                }
            )

    atomic_write_text(log_txt, "".join(txt_lines))
    if log_csv is not None:
        atomic_write_csv(log_csv, csv_rows, ["name", "backbone", "rank", "sample", "score", "length"])

    return {"tag": tag, "backbones_found": backbones, "sequences_written": written, "out_dir": str(out_dir)}


def build_parser(default_tag: str) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=f"Repaired ProteinMPNN top-K selector (dataset {default_tag}).")
    ap.add_argument("--tag", default=default_tag, help="Dataset tag (20 or 30)")
    ap.add_argument("--seqs-dir", type=Path, default=None, help="Directory of ProteinMPNN FASTA files")
    ap.add_argument("--out-dir", type=Path, default=None, help="Output tfold_<tag> directory")
    ap.add_argument("--log-txt", type=Path, default=None, help="Parsed sequence log text path")
    ap.add_argument("--log-csv", type=Path, default=None, help="Optional CSV log path")
    ap.add_argument("--top", type=int, default=3, help="Sequences to keep per backbone")
    ap.add_argument("--start", type=int, default=1, help="First backbone index (inclusive)")
    ap.add_argument("--end", type=int, default=1000, help="Last backbone index (inclusive)")
    ap.add_argument("--selftest", action="store_true", help="Run built-in check and exit")
    return ap


def selftest(default_tag: str = "20") -> int:
    """Verify ascending NLL sort, tie-breaking, multiline FASTA parsing, and deduplication."""
    sample_fa = (
        ">gen_0001, score=3.1000, global_score=3.1000\n"
        "MKTAY\nNATIVE\n"
        ">T=0.15, sample=1, score=0.9500\n"
        "AAAA\nAAAA\n"
        ">T=0.15, sample=2, score=0.8200\n"
        "BBBB\nCCCC\n"
        ">T=0.15, sample=3, score=0.8200\n"
        "AAAA\nDDDD\n"
        ">T=0.15, sample=4, score=0.8200\n"
        "BBBB\nCCCC\n"
        ">T=0.15, sample=5, score=1.4000\n"
        "ZZZZ\nZZZZ\n"
    )
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        seqs = root / f"pmpnn_0{default_tag}" / "seqs"
        seqs.mkdir(parents=True)
        (seqs / f"0001_{default_tag}.t.0.15.fa").write_text(sample_fa, encoding="utf-8")
        out_dir = root / f"tfold_{default_tag}"
        res = process_tag(
            tag=default_tag,
            seqs_dir=seqs,
            out_dir=out_dir,
            log_txt=root / f"parsed_sequences_{default_tag}.txt",
            log_csv=root / f"selections_{default_tag}.csv",
            top=3,
            start=1,
            end=5,
        )
        assert res["sequences_written"] == 3, res
        top1 = (out_dir / f"0001_top.1.{default_tag}.fa").read_text(encoding="utf-8")
        top2 = (out_dir / f"0001_top.2.{default_tag}.fa").read_text(encoding="utf-8")
        top3 = (out_dir / f"0001_top.3.{default_tag}.fa").read_text(encoding="utf-8")
        assert top1 == ">H\nBBBBCCCC\n", top1
        assert top2 == ">H\nAAAADDDD\n", top2
        assert top3 == ">H\nAAAAAAAA\n", top3
    print(f"pipeline process_{default_tag.zfill(3)} selftest: pass")
    return 0


def run_cli(default_tag: str) -> int:
    ap = build_parser(default_tag)
    ns = ap.parse_args()
    if ns.selftest:
        return selftest(default_tag)
    tag = ns.tag
    seqs_dir = ns.seqs_dir or Path(f"pmpnn_0{tag}") / "seqs"
    out_dir = ns.out_dir or Path(f"tfold_{tag}")
    log_txt = ns.log_txt or Path(f"parsed_sequences_{tag}.txt")
    res = process_tag(tag, seqs_dir, out_dir, log_txt, ns.log_csv, ns.top, ns.start, ns.end)
    print(f"Dataset {tag}: backbones={res['backbones_found']} kept={res['sequences_written']} -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run_cli("20"))
