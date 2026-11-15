#!/usr/bin/env python3
# version: 0.3.0-a43ea1-4 | date: 15.XI.2026
"""Clean RFantibody design wrapper with pre-run target hotspot & SASA checks.

Fixes Finding 7 from AUDIT.md:
- Replaces conflicting duplicate `--config` options with `--load-config` and
  `--save-config`, preserving working CLI defaults without `required=not config`.
- Streams subprocess stdout line-by-line rather than buffering until completion.
- Validates that every hotspot residue (`<chain><resseq>`, e.g. `T120`) exists in
  the target antigen PDB and clears `--min-hotspot-sasa` before launching GPU work.
- Adds `--resume` to skip stages whose output PDB/quiver already exists.

Usage:
    py -3 -I design_rfantibody.py --target ./rosetta_dock\\scFv-0389-459-5-Ag.pdb \
        --hotspots T42,T45 --validate-only
    py -3 -I design_rfantibody.py --selftest
"""
from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "SE": 1.90}
HOTSPOT_RE = re.compile(r"^([A-Za-z_])(-?\d+)([A-Za-z]?)$")


def fibonacci_sphere(n: int = 92) -> np.ndarray:
    """Return unit sphere sample vectors for Shrake-Rupley SASA estimates."""
    idx = np.arange(0, n, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * idx / n)
    theta = math.pi * (1.0 + math.sqrt(5.0)) * idx
    return np.column_stack((np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)))


def read_heavy_atoms(pdb_path: Path) -> tuple[np.ndarray, np.ndarray, list[tuple[str, int, str, str]]]:
    """Read heavy atom coords, expanded VDW radii, and `(chain, resseq, icode, resname)` tags."""
    coords, radii, tags = [], [], []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            altloc = line[16].strip() if len(line) > 16 else ""
            if altloc not in ("", "A", "1"):
                continue
            atom_name = line[12:16].strip()
            elem = (line[76:78].strip() or atom_name[0]).upper()
            if elem == "H" or atom_name.startswith("H"):
                continue
            ch = line[21].strip() or "_"
            resseq = int(line[22:26])
            icode = line[26].strip() if len(line) > 26 else ""
            resname = line[17:20].strip()
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            radii.append(VDW.get(elem, 1.70) + 1.40)
            tags.append((ch, resseq, icode, resname))
    if not coords:
        return np.zeros((0, 3), dtype=float), np.zeros((0,), dtype=float), []
    return np.array(coords, dtype=float), np.array(radii, dtype=float), tags


def parse_hotspots(spec: str) -> list[tuple[str, int, str]]:
    """Parse `T42,T100A` hotspot tokens into `(chain, resseq, icode)` tuples."""
    out: list[tuple[str, int, str]] = []
    for raw in [tok.strip() for tok in spec.split(",") if tok.strip()]:
        m = HOTSPOT_RE.match(raw)
        if not m:
            raise ValueError(f"Invalid hotspot token `{raw}`; expected e.g. T42 or A108A")
        out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


def validate_hotspots(pdb_path: Path, spec: str, min_sasa: float) -> tuple[bool, list[dict]]:
    """Check that all hotspots exist in `pdb_path` and meet solvent exposure `min_sasa`."""
    coords, radii, tags = read_heavy_atoms(pdb_path)
    if len(coords) == 0:
        raise ValueError(f"No heavy atoms found in {pdb_path}")

    sphere = fibonacci_sphere(92)
    by_res: dict[tuple[str, int, str], list[int]] = {}
    for idx, (ch, resseq, icode, _) in enumerate(tags):
        by_res.setdefault((ch, resseq, icode), []).append(idx)

    diff = coords[:, None, :] - coords[None, :, :]
    dist_mat = np.sqrt(np.sum(diff * diff, axis=-1))

    report: list[dict] = []
    ok = True
    for ch, resseq, icode in parse_hotspots(spec):
        indices = by_res.get((ch, resseq, icode))
        label = f"{ch}{resseq}{icode}"
        if not indices:
            ok = False
            report.append({"hotspot": label, "present": False, "resname": "", "sasa_a2": 0.0, "exposed": False})
            continue
        resname = tags[indices[0]][3]
        sasa = 0.0
        for ai in indices:
            r = radii[ai]
            # Neighbor atoms excluding self
            nb = np.where((dist_mat[ai] > 1e-6) & (dist_mat[ai] < r + radii + 1e-5))[0]
            pts = coords[ai] + r * sphere
            if len(nb):
                d_nb = np.sqrt(np.sum((pts[:, None, :] - coords[nb][None, :, :]) ** 2, axis=-1))
                free = np.all(d_nb >= radii[nb][None, :] - 1e-6, axis=1)
                frac = float(np.mean(free))
            else:
                frac = 1.0
            sasa += 4.0 * math.pi * r * r * frac
        exposed = bool(sasa >= min_sasa)
        if not exposed:
            ok = False
        report.append({"hotspot": label, "present": True, "resname": resname, "sasa_a2": round(float(sasa), 2), "exposed": exposed})
    return bool(ok), report


def stream_command(cmd: list[str], log_path: Path) -> int:
    """Run command with true line-by-line streaming to stdout and log file."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a", encoding="utf-8") as logf:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            logf.write(line)
        return proc.wait()


def build_stage_commands(cfg: dict) -> list[tuple[str, Path, list[str]]]:
    """Return `(stage_name, expected_output, cmd)` tuples for RFdiffusion, MPNN, and RF2."""
    out_dir = Path(cfg["output_dir"])
    rfd_out = out_dir / "rfdiffusion" / "design.pdb"
    mpnn_out = out_dir / "proteinmpnn" / "design_mpnn.pdb"
    rf2_out = out_dir / "rf2" / "design_rf2.pdb"
    rf_dir = Path(cfg["rfantibody_dir"])

    return [
        (
            "rfdiffusion",
            rfd_out,
            [
                cfg["python_bin"],
                str(rf_dir / "scripts" / "rfdiffusion_inference.py"),
                "--config-name", "antibody",
                f"antibody.target_pdb={cfg['target']}",
                f"antibody.framework_pdb={cfg['framework']}",
                f"inference.num_designs={cfg['num_designs']}",
                f"ppi.hotspot_res=[{cfg['hotspots']}]",
                f"inference.output_prefix={rfd_out.with_suffix('')}",
            ],
        ),
        (
            "proteinmpnn",
            mpnn_out,
            [
                cfg["python_bin"],
                str(rf_dir / "scripts" / "proteinmpnn_interface_design.py"),
                "-pdbdir", str(rfd_out.parent),
                "-outpdbdir", str(mpnn_out.parent),
                "-seqs_per_struct", str(cfg["seqs_per_struct"]),
                "-temperature", str(cfg["temperature"]),
            ],
        ),
        (
            "rf2",
            rf2_out,
            [
                cfg["python_bin"],
                str(rf_dir / "scripts" / "rf2_predict.py"),
                f"input.pdb_dir={mpnn_out.parent}",
                f"output.pdb_dir={rf2_out.parent}",
            ],
        ),
    ]


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Clean RFantibody pipeline wrapper with hotspot validation.")
    ap.add_argument("--load-config", type=Path, help="Load JSON defaults before CLI overrides")
    ap.add_argument("--save-config", type=Path, help="Write resolved configuration JSON and continue/exit")
    ap.add_argument("--target", type=Path, help="Target antigen PDB (e.g. 8e8z_C)")
    ap.add_argument("--framework", type=Path, default=Path("framework.pdb"), help="VHH framework PDB")
    ap.add_argument("--rfantibody-dir", type=Path, default=Path("RFantibody"), help="MIT RFantibody directory")
    ap.add_argument("--output-dir", type=Path, default=Path("rfantibody_out"), help="Output working directory")
    ap.add_argument("--hotspots", default="", help="Comma-separated hotspots, e.g. T42,T45")
    ap.add_argument("--min-hotspot-sasa", type=float, default=15.0, help="Minimum SASA (A^2) per hotspot residue")
    ap.add_argument("--num-designs", type=int, default=100, help="RFdiffusion backbone count")
    ap.add_argument("--seqs-per-struct", type=int, default=4, help="ProteinMPNN sequences per backbone")
    ap.add_argument("--temperature", type=float, default=0.15, help="MPNN temperature")
    ap.add_argument("--python-bin", default=sys.executable, help="Python executable for RFantibody scripts")
    ap.add_argument("--resume", action="store_true", help="Skip stages whose output path already exists")
    ap.add_argument("--validate-only", action="store_true", help="Run target hotspot SASA check only")
    ap.add_argument("--dry-run", action="store_true", help="Print resolved stages without executing")
    ap.add_argument("--selftest", action="store_true", help="Run built-in verification and exit")
    return ap


def resolve_config(ns: argparse.Namespace) -> dict:
    """Merge JSON file config (if given) with explicit CLI flags."""
    cfg = {}
    if ns.load_config and ns.load_config.exists():
        cfg.update(json.loads(ns.load_config.read_text(encoding="utf-8")))
    for key in (
        "target", "framework", "rfantibody_dir", "output_dir", "hotspots",
        "min_hotspot_sasa", "num_designs", "seqs_per_struct", "temperature", "python_bin",
    ):
        val = getattr(ns, key)
        if val is not None:
            cfg[key] = str(val) if isinstance(val, Path) else val
    return cfg


def selftest() -> int:
    """Verify non-conflicting config options, hotspot existence/SASA check, and resume plan."""
    ap = build_parser()
    ns = ap.parse_args(["--target", "t.pdb", "--hotspots", "T1", "--save-config", "c.json", "--resume", "--dry-run"])
    assert ns.resume and ns.dry_run

    # Two residues + insertion code residues: exposed T1, buried T2, T100A and T100B.
    lines = [
        "ATOM      1  CA  ARG T   1      20.000   0.000   0.000  1.00 50.00           C\n",
        "ATOM      2  CA  ALA T   2       0.000   0.000   0.000  1.00 50.00           C\n",
        "ATOM     10  CA  TYR T 100A     40.000   0.000   0.000  1.00 50.00           C\n",
        "ATOM     11  CA  SER T 100B     50.000   0.000   0.000  1.00 50.00           C\n",
    ]
    for idx, (dx, dy, dz) in enumerate(((2,0,0),(-2,0,0),(0,2,0),(0,-2,0),(0,0,2),(0,0,-2)), start=3):
        lines.append(f"ATOM    {idx:3d}  C   GLY T   9     {dx:7.3f} {dy:7.3f} {dz:7.3f}  1.00 50.00           C\n")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        pdb = root / "ag.pdb"
        pdb.write_text("".join(lines), encoding="utf-8")

        ok1, rep1 = validate_hotspots(pdb, "T1", 15.0)
        assert ok1 and rep1[0]["sasa_a2"] > 100.0, rep1
        assert json.dumps({"ok": ok1, "hotspots": rep1})

        ok_ic, rep_ic = validate_hotspots(pdb, "T100A,T100B", 15.0)
        assert ok_ic and [r["resname"] for r in rep_ic] == ["TYR", "SER"], rep_ic

        ok2, rep2 = validate_hotspots(pdb, "T2", 15.0)
        assert not ok2 and rep2[0]["sasa_a2"] == 0.0, rep2

        ok3, rep3 = validate_hotspots(pdb, "T999", 15.0)
        assert not ok3 and not rep3[0]["present"], rep3

        cfg = {
            "target": str(pdb), "framework": "fw.pdb", "rfantibody_dir": str(root),
            "output_dir": str(root / "out"), "hotspots": "T1", "min_hotspot_sasa": 15.0,
            "num_designs": 10, "seqs_per_struct": 2, "temperature": 0.15, "python_bin": sys.executable,
        }
        stages = build_stage_commands(cfg)
        assert [s[0] for s in stages] == ["rfdiffusion", "proteinmpnn", "rf2"]

    print("design_rfantibody selftest: pass")
    return 0


def main() -> int:
    ap = build_parser()
    ns = ap.parse_args()
    if ns.selftest:
        return selftest()

    cfg = resolve_config(ns)
    if ns.save_config:
        ns.save_config.parent.mkdir(parents=True, exist_ok=True)
        ns.save_config.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")

    if not cfg.get("target"):
        ap.error("--target is required")
    target = Path(cfg["target"])
    if not target.exists():
        print(f"Missing target PDB: {target}", file=sys.stderr)
        return 1

    hotspot_report: list[dict] = []
    if cfg.get("hotspots"):
        ok, hotspot_report = validate_hotspots(target, cfg["hotspots"], float(cfg["min_hotspot_sasa"]))
        print(json.dumps({"target": str(target), "hotspots_ok": ok, "hotspots": hotspot_report}, indent=2))
        if not ok:
            print("Hotspot validation failed before GPU stage.", file=sys.stderr)
            return 2
    if ns.validate_only:
        return 0

    stages = build_stage_commands(cfg)
    if ns.dry_run:
        print(json.dumps([{"stage": s, "out": str(p), "cmd": " ".join(c)} for s, p, c in stages], indent=2))
        return 0

    log_path = Path(cfg["output_dir"]) / "rfantibody.log"
    for stage_name, expected_out, cmd in stages:
        if ns.resume and expected_out.exists():
            print(f"[resume] skipping {stage_name}: {expected_out} exists")
            continue
        rc = stream_command(cmd, log_path)
        if rc != 0:
            print(f"Stage {stage_name} exited with {rc}", file=sys.stderr)
            return rc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
