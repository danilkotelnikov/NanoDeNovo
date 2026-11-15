#!/usr/bin/env python3
"""Repaired complex FASTA generator for tFold-Ag docking.

Repairs applied relative to legacy `tfold_ag_prepare_kwargs.py`:
- Preserves the target antigen chain identifier (`>{chain}`, e.g. `>C` for `8e8z_C`)
  instead of overwriting it to `>A`.
- Adds missing `continue` guard when the antigen PDB file does not exist (instead of crashing
  later opening an uncreated `raw_fa_path`).
- Includes a pure-Python CA-to-FASTA extractor fallback when external `./pdb2fasta` is absent.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

THREE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def extract_chains_from_pdb(pdb_path: Path) -> dict[str, str]:
    """Pure-Python PDB chain sequence reader from CA records preserving insertion order."""
    chains: dict[str, list[str]] = {}
    seen: set[tuple[str, str]] = set()
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            ch = line[21].strip() or "A"
            res_id = line[22:27].strip()
            key = (ch, res_id)
            if key in seen:
                continue
            seen.add(key)
            chains.setdefault(ch, []).append(THREE.get(line[17:20].strip(), "X"))
    return {ch: "".join(seq) for ch, seq in chains.items()}


def process_datasets(
    work_dir: Path,
    datasets: list[str],
    antigens: list[str],
    chains: list[str],
    rmsd: float = 1.0,
    tm: float = 0.95,
    pdb2fasta_bin: Path | None = None,
) -> dict:
    created = 0
    skipped_missing_ag = 0
    suffix = f"_rmsd_{rmsd}_tm_{tm}"

    for dataset in datasets:
        aligned_dir = work_dir / f"usalign{suffix}" / f"aligned_{dataset}"
        txt_path = aligned_dir / f"v_{dataset}_usalign{suffix}.txt"
        if not txt_path.exists():
            continue
        entities = [ln.strip() for ln in txt_path.read_text(encoding="utf-8").splitlines() if ln.strip()]

        tfold_dir = work_dir / f"tfold_{dataset}"
        dest_fastas_dir = work_dir / "filtered_nanobodies" / f"v_{dataset}" / "fastas"
        dest_fastas_dir.mkdir(parents=True, exist_ok=True)

        for entity_pdb in entities:
            base = Path(entity_pdb).stem
            src_fa = tfold_dir / f"{base}.fa"
            if src_fa.exists():
                shutil.copy2(src_fa, dest_fastas_dir / f"{base}.fa")

        for antigen in antigens:
            antigen_pdb_path = work_dir / "antigens" / "pdbs" / f"{antigen}.pdb"
            raw_fa_path = work_dir / "antigens" / "fastas" / f"{antigen}.fa"
            edited_fa_path = work_dir / "antigens" / "fastas" / f"{antigen}_edited.fa"

            antigen_sequences: dict[str, str] = {}
            if not edited_fa_path.exists():
                if not raw_fa_path.exists():
                    if not antigen_pdb_path.exists():
                        skipped_missing_ag += 1
                        continue
                    raw_fa_path.parent.mkdir(parents=True, exist_ok=True)
                    if pdb2fasta_bin and pdb2fasta_bin.exists():
                        with open(raw_fa_path, "w", encoding="utf-8") as f_out:
                            subprocess.run([str(pdb2fasta_bin), str(antigen_pdb_path)], stdout=f_out, check=True)
                    else:
                        extracted = extract_chains_from_pdb(antigen_pdb_path)
                        raw_lines = [f">{antigen}:{ch}\n{seq}\n" for ch, seq in extracted.items()]
                        raw_fa_path.write_text("".join(raw_lines), encoding="utf-8")

                with open(raw_fa_path, "r", encoding="utf-8") as f_in, open(edited_fa_path, "w", encoding="utf-8") as f_out:
                    for line in f_in:
                        if line.startswith(">"):
                            parts = line.strip().split(":")
                            chain_part = parts[-1].split()[0] if len(parts) >= 2 else line[1:].strip().split()[0]
                            f_out.write(f">{chain_part or 'A'}\n")
                        else:
                            f_out.write(line)

            cur_h: str | None = None
            cur_s: list[str] = []
            for raw in edited_fa_path.read_text(encoding="utf-8").splitlines():
                line = raw.strip()
                if line.startswith(">"):
                    if cur_h is not None:
                        antigen_sequences[cur_h] = "".join(cur_s)
                    cur_h = line[1:].strip()
                    cur_s = []
                elif line:
                    cur_s.append(line)
            if cur_h is not None:
                antigen_sequences[cur_h] = "".join(cur_s)

            for chain in chains:
                if chain not in antigen_sequences:
                    continue
                ag_seq = antigen_sequences[chain]
                complex_dir = work_dir / "filtered_nanobodies" / f"v_{dataset}" / "complexes" / f"{antigen}_{chain}" / "fastas"
                complex_dir.mkdir(parents=True, exist_ok=True)

                for entity_pdb in entities:
                    base = Path(entity_pdb).stem
                    nb_fa_path = dest_fastas_dir / f"{base}.fa"
                    if not nb_fa_path.exists():
                        continue
                    h_seq_parts: list[str] = []
                    read_h = False
                    for raw in nb_fa_path.read_text(encoding="utf-8").splitlines():
                        line = raw.strip()
                        if line.startswith(">"):
                            read_h = (line[1:].strip() == "H")
                        elif read_h and line:
                            h_seq_parts.append(line)
                    h_seq = "".join(h_seq_parts)
                    if not h_seq:
                        continue
                    c_path = complex_dir / f"{base}_{antigen}_{chain}_complex.fa"
                    c_path.write_text(f">H\n{h_seq}\n>{chain}\n{ag_seq}\n", encoding="utf-8")
                    created += 1

    return {"created": created, "skipped_missing_antigen": skipped_missing_ag}


def selftest() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        aligned = root / "usalign_rmsd_1.0_tm_0.95" / "aligned_20"
        aligned.mkdir(parents=True)
        (aligned / "v_20_usalign_rmsd_1.0_tm_0.95.txt").write_text("0001_top.1.20.pdb\n", encoding="utf-8")
        tfold = root / "tfold_20"
        tfold.mkdir(parents=True)
        (tfold / "0001_top.1.20.fa").write_text(">H\nQVQLVES\n", encoding="utf-8")
        pdbs = root / "antigens" / "pdbs"
        pdbs.mkdir(parents=True)
        (pdbs / "8e8z.pdb").write_text(
            "ATOM      1  CA  ALA C   1       0.000   0.000   0.000  1.00 50.00           C\n"
            "ATOM      2  CA  GLY C   2       3.800   0.000   0.000  1.00 50.00           C\n",
            encoding="utf-8",
        )
        res = process_datasets(root, ["20"], ["8e8z", "missing_ag"], ["C"])
        assert res["created"] == 1 and res["skipped_missing_antigen"] == 1, res
        out_fa = root / "filtered_nanobodies" / "v_20" / "complexes" / "8e8z_C" / "fastas" / "0001_top.1.20_8e8z_C_complex.fa"
        txt = out_fa.read_text(encoding="utf-8")
        assert ">H\nQVQLVES\n>C\nAG\n" == txt, txt
    print("pipeline tfold_ag_prepare_kwargs selftest: pass")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Repaired complex FASTA generator preserving chain IDs.")
    ap.add_argument("--work-dir", type=Path, default=Path("."), help="Working directory")
    ap.add_argument("--datasets", nargs="+", default=["20", "30"], help="Dataset tags")
    ap.add_argument("--antigen", nargs="+", default=["8e8z"], help="Antigen PDB stems")
    ap.add_argument("--chains", nargs="+", default=["C"], help="Target chain IDs (e.g. C)")
    ap.add_argument("--rmsd", type=float, default=1.0, help="Manifest RMSD suffix")
    ap.add_argument("--tm", type=float, default=0.95, help="Manifest TM suffix")
    ap.add_argument("--selftest", action="store_true", help="Run built-in check and exit")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    res = process_datasets(args.work_dir, args.datasets, args.antigen, args.chains, args.rmsd, args.tm)
    print(f"tfold_ag_prepare_kwargs: created={res['created']} skipped_missing_ag={res['skipped_missing_antigen']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
