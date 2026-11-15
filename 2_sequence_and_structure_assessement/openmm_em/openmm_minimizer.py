#!/usr/bin/env python3
"""Clean non-interactive OpenMM energy minimizer for nanobody-antigen complexes.

Repairs the original interactive script in-place:
- Replaces top-level `input()` blocking calls with CLI flags (retaining TTY fallback
  only when `--interactive` is given or stdin is an interactive terminal without args).
- Uses `tempfile.TemporaryDirectory()` instead of hardcoded `/tmp` for Windows/WSL safety.
- Defaults to implicit solvent (`amber14/protein.ff14SB.xml` + `implicit/gbn2.xml`)
  and harmonic C-alpha positional restraints (`CustomExternalForce`), with explicit
  `--solvent vacuum` legacy override.
- Preserves original chain IDs and residue numbering (`keepIds=True`), disables
  cross-gap loop hallucination (`fixer.missingResidues = {}` by default), and records
  before/after heavy-atom clashes (`< 2.2` A) and energies atomically.
- Returns non-zero exit code when any complex fails or no PDBs match.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import tempfile
from pathlib import Path
from typing import Any


def parse_heavy_atoms(pdb_path: Path) -> list[tuple[str, int, str, str, float, float, float]]:
    """Return `(chain, resseq, icode, elem, x, y, z)` for heavy ATOM records."""
    atoms: list[tuple[str, int, str, str, float, float, float]] = []
    seen: set[tuple[str, int, str, str]] = set()
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            altloc = line[16]
            if altloc not in (" ", "A", "1"):
                continue
            atom_name = line[12:16].strip()
            elem = (line[76:78].strip() or atom_name[0]).upper()
            if elem == "H" or atom_name.startswith("H"):
                continue
            chain = line[21].strip() or "_"
            resseq = int(line[22:26])
            icode = line[26].strip()
            key = (chain, resseq, icode, atom_name)
            if key in seen:
                continue
            seen.add(key)
            atoms.append(
                (
                    chain,
                    resseq,
                    icode,
                    elem,
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
            )
    return atoms


def audit_heavy_clashes(pdb_path: Path, cutoff: float = 2.2) -> dict[str, Any]:
    """Count non-adjacent heavy-atom clashes `< cutoff` A and report min distance."""
    atoms = parse_heavy_atoms(pdb_path)
    if len(atoms) < 2:
        return {"heavy_atoms": len(atoms), "clashes": 0, "min_distance_a": None}

    cutoff_sq = cutoff * cutoff
    clashes = 0
    inter_clashes = 0
    min_dist_sq = float("inf")
    n = len(atoms)
    for i in range(n):
        ci, ri, ii, _, xi, yi, zi = atoms[i]
        for j in range(i + 1, n):
            cj, rj, ij, _, xj, yj, zj = atoms[j]
            same_chain = ci == cj
            if same_chain and ii == ij and abs(ri - rj) <= 1:
                continue
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            d2 = dx * dx + dy * dy + dz * dz
            if d2 < min_dist_sq:
                min_dist_sq = d2
            if d2 < cutoff_sq:
                clashes += 1
                if not same_chain:
                    inter_clashes += 1

    return {
        "heavy_atoms": n,
        "clashes": clashes,
        "inter_chain_clashes": inter_clashes,
        "min_distance_a": round(math.sqrt(min_dist_sq), 4) if math.isfinite(min_dist_sq) else None,
    }


def atomic_write_text(path: Path, text: str) -> None:
    """Atomically write text file via temporary sibling and `os.replace`."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".tmp_", suffix=path.suffix, dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
        os.replace(tmp_name, path)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def sanitize_input_pdb(raw_text: str) -> str:
    """Drop extra `MODEL`/`ENDMDL` delimiters and keep clean ATOM/HETATM/TER records."""
    out: list[str] = []
    for line in raw_text.splitlines():
        if line.startswith(("ATOM  ", "TER", "REMARK", "CRYST1")):
            out.append(line)
        elif line.startswith("END") and not line.startswith("ENDMDL"):
            continue
    out.append("END")
    return "\n".join(out) + "\n"


def _add_ca_restraints(system: Any, topology: Any, positions: Any, restraint_k_kcal: float) -> int:
    """Attach harmonic positional restraints (kcal/mol/A^2) on all C-alpha atoms."""
    if restraint_k_kcal <= 0.0:
        return 0
    from openmm import CustomExternalForce, unit

    k_kj = restraint_k_kcal * 4.184 / (0.1 * 0.1)  # kcal/(mol*A^2) -> kJ/(mol*nm^2)
    force = CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    force.addGlobalParameter("k", k_kj)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    restrained = 0
    for atom in topology.atoms():
        if atom.name == "CA":
            pos = positions[atom.index].value_in_unit(unit.nanometers)
            force.addParticle(atom.index, [pos[0], pos[1], pos[2]])
            restrained += 1
    if restrained:
        system.addForce(force)
    return restrained


def minimize_pdb(
    input_pdb: Path,
    output_pdb: Path,
    log_path: Path,
    solvent: str = "gbn2",
    restraint_k: float = 10.0,
    tolerance_kj: float = 10.0,
    max_iterations: int = 1000,
    keep_missing_residues: bool = False,
) -> dict[str, Any]:
    """Run PDBFixer + OpenMM minimization on one complex and return audit metrics."""
    try:
        from openmm import CustomExternalForce, LocalEnergyMinimizer, VerletIntegrator, app, unit  # noqa: F401
        from pdbfixer.pdbfixer import PDBFixer
    except ImportError as exc:
        raise RuntimeError(
            "OpenMM and pdbfixer are required for live minimization. "
            "Run inside the configured WSL/conda environment (`immunebuilder` or `openmm`)."
        ) from exc

    before_clashes = audit_heavy_clashes(input_pdb)
    log_lines = [
        f"Minimization log for {input_pdb.stem}",
        "=" * 40,
        f"Input PDB: {input_pdb}",
        f"Solvent: {solvent} | CA restraint k: {restraint_k} kcal/(mol*A^2) | Tolerance: {tolerance_kj} kJ/(mol*nm)",
        f"Pre-relax clashes (<2.2 A): {before_clashes['clashes']} (min_dist={before_clashes['min_distance_a']} A)",
    ]

    with tempfile.TemporaryDirectory(prefix="openmm_em_") as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        cleaned_pdb = tmp_dir / "cleaned.pdb"
        cleaned_pdb.write_text(
            sanitize_input_pdb(input_pdb.read_text(encoding="utf-8", errors="replace")),
            encoding="utf-8",
        )

        fixer = PDBFixer(filename=str(cleaned_pdb))
        fixer.findMissingResidues()
        if not keep_missing_residues:
            fixer.missingResidues = {}
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        fixed_pdb = tmp_dir / "fixed.pdb"
        with open(fixed_pdb, "w", encoding="utf-8") as handle:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

        pdb = app.PDBFile(str(fixed_pdb))
        ff_files = ["amber14/protein.ff14SB.xml"]
        if solvent == "gbn2":
            ff_files.append("implicit/gbn2.xml")
        forcefield = app.ForceField(*ff_files)

        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(forcefield, pH=7.4)
        modeller.addExtraParticles(forcefield)

        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
        )
        n_ca = _add_ca_restraints(system, modeller.topology, modeller.positions, restraint_k)

        integrator = VerletIntegrator(0.002 * unit.picoseconds)
        simulation = app.Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

        e_init = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        tol = tolerance_kj * unit.kilojoule_per_mole / unit.nanometer
        LocalEnergyMinimizer.minimize(simulation.context, tolerance=tol, maxIterations=max_iterations)

        state = simulation.context.getState(getEnergy=True, getPositions=True)
        e_final = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

        out_tmp = tmp_dir / "minimized.pdb"
        with open(out_tmp, "w", encoding="utf-8") as handle:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), handle, keepIds=True)

        output_pdb.parent.mkdir(parents=True, exist_ok=True)
        atomic_write_text(output_pdb, out_tmp.read_text(encoding="utf-8"))

    after_clashes = audit_heavy_clashes(output_pdb)
    log_lines.extend(
        [
            f"Restrained CA atoms: {n_ca}",
            f"Initial potential energy: {e_init:.2f} kJ/mol",
            f"Final potential energy: {e_final:.2f} kJ/mol",
            f"Post-relax clashes (<2.2 A): {after_clashes['clashes']} (min_dist={after_clashes['min_distance_a']} A)",
        ]
    )
    atomic_write_text(log_path, "\n".join(log_lines) + "\n")

    return {
        "name": input_pdb.stem,
        "status": "ok",
        "output_pdb": str(output_pdb),
        "log_path": str(log_path),
        "solvent": solvent,
        "restraint_k": restraint_k,
        "restrained_ca": n_ca,
        "initial_energy_kj": round(e_init, 3),
        "final_energy_kj": round(e_final, 3),
        "clashes_before": before_clashes["clashes"],
        "clashes_after": after_clashes["clashes"],
        "min_distance_before_a": before_clashes["min_distance_a"],
        "min_distance_after_a": after_clashes["min_distance_a"],
    }


def resolve_io_paths(args: argparse.Namespace) -> list[tuple[Path, Path, Path]]:
    """Resolve input/output PDB tuples from CLI flags or optional TTY prompts."""
    if args.input:
        in_pdb = Path(args.input)
        out_pdb = Path(args.output) if args.output else in_pdb.with_name(f"{in_pdb.stem}_em.pdb")
        log_path = out_pdb.with_suffix(".log")
        return [(in_pdb, out_pdb, log_path)]

    in_dir_str = args.input_dir
    out_dir_str = args.output_dir
    if (not in_dir_str or not out_dir_str) and (args.interactive or sys.stdin.isatty()):
        if not in_dir_str:
            in_dir_str = input("Enter the path to the input directory containing PDB files: ").strip()
        if not out_dir_str:
            out_dir_str = input("Enter the path to the output directory: ").strip()

    if not in_dir_str or not out_dir_str:
        raise ValueError("Provide --input-dir and --output-dir (or --input for a single PDB).")

    in_dir = Path(in_dir_str)
    out_dir = Path(out_dir_str)
    if not in_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {in_dir}")

    pdbs = sorted(in_dir.glob("*.pdb"))
    if not pdbs:
        raise FileNotFoundError(f"No .pdb files found in {in_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)
    return [(p, out_dir / f"{p.stem}_em.pdb", out_dir / f"{p.stem}_em.log") for p in pdbs]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Clean OpenMM implicit-solvent minimizer with heavy-atom clash audit.")
    p.add_argument("--input-dir", help="Input directory containing complex .pdb files")
    p.add_argument("--output-dir", help="Output directory for minimized PDBs and logs")
    p.add_argument("--input", help="Single input .pdb file (alternative to --input-dir)")
    p.add_argument("--output", help="Single output .pdb path (used with --input)")
    p.add_argument("--solvent", choices=["gbn2", "vacuum"], default="gbn2", help="Solvent model (default: gbn2)")
    p.add_argument("--restraint-k", type=float, default=10.0, help="CA harmonic restraint in kcal/(mol*A^2)")
    p.add_argument("--tolerance", type=float, default=10.0, help="Force convergence tolerance in kJ/(mol*nm)")
    p.add_argument("--max-iterations", type=int, default=1000, help="Maximum OpenMM minimization iterations")
    p.add_argument("--keep-missing-residues", action="store_true", help="Allow PDBFixer loop insertion across gaps")
    p.add_argument("--interactive", action="store_true", help="Prompt interactively when IO flags are omitted")
    p.add_argument("--dry-run", action="store_true", help="Run CPU clash audit & IO validation without OpenMM")
    p.add_argument("--selftest", action="store_true", help="Run CPU self-tests (no OpenMM required) and exit")
    return p


def selftest() -> int:
    """Verify clash detection, sanitization, atomic writes, and non-TTY argument validation."""
    with tempfile.TemporaryDirectory(prefix="openmm_em_test_") as tmp_str:
        tmp = Path(tmp_str)
        pdb_dir = tmp / "in"
        out_dir = tmp / "out"
        pdb_dir.mkdir()
        clash_pdb = pdb_dir / "complex_a.pdb"
        clash_pdb.write_text(
            "MODEL        1\n"
            "ATOM      1  CA  ALA H   1       0.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      2  CA  GLY C  50       1.200   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  CA  VAL C  51       4.900   0.000   0.000  1.00  0.00           C\n"
            "ENDMDL\n"
            "END\n",
            encoding="utf-8",
        )
        cleaned = sanitize_input_pdb(clash_pdb.read_text(encoding="utf-8"))
        assert "ENDMDL" not in cleaned
        assert cleaned.strip().endswith("END")

        audit = audit_heavy_clashes(clash_pdb, cutoff=2.2)
        assert audit["clashes"] == 1
        assert audit["inter_chain_clashes"] == 1
        assert audit["min_distance_a"] == 1.2

        ns = build_parser().parse_args(["--input-dir", str(pdb_dir), "--output-dir", str(out_dir), "--dry-run"])
        jobs = resolve_io_paths(ns)
        assert len(jobs) == 1
        assert jobs[0][1].name == "complex_a_em.pdb"

    print("openmm_minimizer selftest: pass")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.selftest:
        return selftest()

    try:
        jobs = resolve_io_paths(args)
    except (ValueError, FileNotFoundError) as exc:
        print(f"Argument error: {exc}", file=sys.stderr)
        return 2

    records: list[dict[str, Any]] = []
    failed = 0
    for in_pdb, out_pdb, log_path in jobs:
        if args.dry_run:
            pre = audit_heavy_clashes(in_pdb)
            records.append({"name": in_pdb.stem, "status": "dry-run", "input": str(in_pdb), "pre_audit": pre})
            continue
        try:
            rec = minimize_pdb(
                input_pdb=in_pdb,
                output_pdb=out_pdb,
                log_path=log_path,
                solvent=args.solvent,
                restraint_k=args.restraint_k,
                tolerance_kj=args.tolerance,
                max_iterations=args.max_iterations,
                keep_missing_residues=args.keep_missing_residues,
            )
            records.append(rec)
            print(f"[ok] {in_pdb.name} -> {out_pdb.name} (clashes {rec['clashes_before']} -> {rec['clashes_after']})")
        except Exception as exc:
            failed += 1
            atomic_write_text(log_path, f"Minimization failed for {in_pdb.name}: {exc}\n")
            records.append({"name": in_pdb.stem, "status": "failed", "error": str(exc)})
            print(f"[fail] {in_pdb.name}: {exc}", file=sys.stderr)

    summary_dir = Path(args.output_dir) if args.output_dir else jobs[0][1].parent
    atomic_write_text(summary_dir / "minimization_summary.json", json.dumps({"jobs": records, "failed": failed}, indent=2))
    if failed:
        print(f"Completed with {failed} failed complex(es).", file=sys.stderr)
        return 1
    print(f"Processed {len(records)} complex(es) successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
