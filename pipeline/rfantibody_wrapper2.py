#!/usr/bin/env python3
"""Repaired RFantibody wrapper for NanoDeNovo Stage 1.

Repairs applied relative to the original script:
- Separates `--load-config` and `--save-config` (removes conflicting duplicate `--config` options).
- Removes contradictory `required=not config` flags while providing TTY prompts only when interactive.
- Streams subprocess output line by line to stdout and log file with UTF-8 encoding.
- Validates hotspot presence and solvent-accessible surface area (SASA) prior to GPU launch.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TOOLS_DIR = ROOT / "tools"
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

from design_rfantibody import validate_hotspots  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Repaired RFantibody pipeline wrapper with hotspot SASA check.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--load-config", type=Path, help="Load JSON defaults before CLI overrides")
    ap.add_argument("--save-config", type=Path, help="Optional JSON path to save resolved configuration")
    ap.add_argument("--target_pdb", "--target", dest="target_pdb", help="Path to antigen PDB file (e.g. 8e8z_C)")
    ap.add_argument("--framework_pdb", "--framework", dest="framework_pdb", help="Path to VHH framework PDB file")
    ap.add_argument("--ckpt_override_path", default="", help="Path to RFdiffusion checkpoint")
    ap.add_argument("--hotspot_res", "--hotspots", dest="hotspot_res", default="", help="Comma-separated hotspot residues (e.g. C42,C45)")
    ap.add_argument("--design_loops", default="H3:7-15", help="Comma-separated design loop ranges (e.g. H3:7-15)")
    ap.add_argument("--num_designs", type=int, default=100, help="Number of backbones to generate")
    ap.add_argument("--final_step", type=int, default=48, help="Final diffusion step")
    ap.add_argument("--loop_string", default="H3", help="Loop string passed to ProteinMPNN interface design")
    ap.add_argument("--seqs_per_struct", type=int, default=4, help="Sequences sampled per structure")
    ap.add_argument("--num_connections", type=int, default=48, help="ProteinMPNN k-NN connections")
    ap.add_argument("--checkpoint_name", default="", help="Optional ProteinMPNN checkpoint override")
    ap.add_argument("--rfantibody_dir", type=Path, default=Path("RFantibody"), help="MIT RFantibody checkout path")
    ap.add_argument("--output_dir", type=Path, default=None, help="Output directory (defaults beside framework PDB)")
    ap.add_argument("--python_bin", default=sys.executable, help="Python executable for RFantibody scripts")
    ap.add_argument("--min_hotspot_sasa", type=float, default=15.0, help="Minimum hotspot SASA (A^2)")
    ap.add_argument("--resume", action="store_true", help="Skip stages whose output PDB already exists")
    ap.add_argument("--validate_only", "--validate-only", dest="validate_only", action="store_true", help="Run hotspot validation only")
    ap.add_argument("--dry_run", "--dry-run", dest="dry_run", action="store_true", help="Print resolved stage commands without running")
    ap.add_argument("--selftest", action="store_true", help="Run built-in check and exit")
    return ap


def resolve_config(ns: argparse.Namespace) -> dict:
    cfg: dict = {}
    if ns.load_config and ns.load_config.exists():
        cfg.update(json.loads(ns.load_config.read_text(encoding="utf-8")))
    keys = (
        "target_pdb", "framework_pdb", "ckpt_override_path", "hotspot_res",
        "design_loops", "num_designs", "final_step", "loop_string",
        "seqs_per_struct", "num_connections", "checkpoint_name",
        "rfantibody_dir", "output_dir", "python_bin", "min_hotspot_sasa",
    )
    for k in keys:
        val = getattr(ns, k, None)
        if val is not None and val != "":
            cfg[k] = str(val) if isinstance(val, Path) else val
        elif k not in cfg:
            cfg[k] = "" if isinstance(val, str) else val
    return cfg


def stream_command(cmd: list[str], log_file: Path) -> int:
    """Execute command with real-time UTF-8 streaming to console and log file."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    with open(log_file, "a", encoding="utf-8", errors="replace") as logf:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            logf.write(line)
            logf.flush()
        return proc.wait()


def build_commands(cfg: dict) -> list[tuple[str, Path, list[str]]]:
    fw_path = Path(cfg["framework_pdb"])
    fw_name = fw_path.stem.replace(".", "_")
    out_dir = Path(cfg["output_dir"]) if cfg.get("output_dir") else (fw_path.parent / fw_name)
    rf_dir = Path(cfg.get("rfantibody_dir") or "RFantibody")
    py_bin = str(cfg.get("python_bin") or sys.executable)

    rfdiff_dir = out_dir / "rfdiffusion"
    pmpnn_dir = out_dir / "pmpnn"
    rf2_dir = out_dir / "rf2"
    rfdiff_prefix = rfdiff_dir / f"{fw_name}_rfdiff"

    rfdiff_cmd = [
        py_bin,
        str(rf_dir / "scripts" / "rfdiffusion_inference.py"),
        "--config-name", "antibody",
        f"antibody.target_pdb={cfg['target_pdb']}",
        f"antibody.framework_pdb={cfg['framework_pdb']}",
        f"ppi.hotspot_res=[{cfg.get('hotspot_res', '')}]",
        f"antibody.design_loops=[{cfg.get('design_loops', 'H3:7-15')}]",
        f"inference.num_designs={cfg.get('num_designs', 100)}",
        f"inference.final_step={cfg.get('final_step', 48)}",
        "inference.deterministic=True",
        f"inference.output_prefix={rfdiff_prefix}",
    ]
    if cfg.get("ckpt_override_path"):
        rfdiff_cmd.append(f"inference.ckpt_override_path={cfg['ckpt_override_path']}")

    pmpnn_cmd = [
        py_bin,
        str(rf_dir / "scripts" / "proteinmpnn_interface_design.py"),
        "-pdbdir", str(rfdiff_dir),
        "-outpdbdir", str(pmpnn_dir),
        "-loop_string", str(cfg.get("loop_string", "H3")),
        "-seqs_per_struct", str(cfg.get("seqs_per_struct", 4)),
        "-num_connections", str(cfg.get("num_connections", 48)),
    ]
    if cfg.get("checkpoint_name"):
        pmpnn_cmd.extend(["-checkpoint_name", str(cfg["checkpoint_name"])])

    rf2_cmd = [
        py_bin,
        str(rf_dir / "scripts" / "rf2_predict.py"),
        f"input.pdb_dir={pmpnn_dir}",
        f"output.pdb_dir={rf2_dir}",
    ]

    return [
        ("rfdiffusion", rfdiff_dir / f"{fw_name}_rfdiff_0.pdb", rfdiff_cmd),
        ("proteinmpnn", pmpnn_dir / f"{fw_name}_mpnn.pdb", pmpnn_cmd),
        ("rf2", rf2_dir / f"{fw_name}_rf2.pdb", rf2_cmd),
    ]


def selftest() -> int:
    ap = build_parser()
    ns = ap.parse_args(["--target_pdb", "t.pdb", "--framework_pdb", "f.pdb", "--save-config", "out.json", "--dry-run"])
    assert ns.dry_run and ns.save_config == Path("out.json")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        pdb = root / "ag.pdb"
        pdb.write_text(
            "ATOM      1  CA  ARG C   1      20.000   0.000   0.000  1.00 50.00           C\n",
            encoding="utf-8",
        )
        ok, rep = validate_hotspots(pdb, "C1", 15.0)
        assert ok and rep[0]["present"], rep
        cfg = {
            "target_pdb": str(pdb),
            "framework_pdb": str(root / "fw.pdb"),
            "output_dir": str(root / "out"),
            "hotspot_res": "C1",
            "design_loops": "H3:7-15",
            "num_designs": 10,
        }
        cmds = build_commands(cfg)
        assert [c[0] for c in cmds] == ["rfdiffusion", "proteinmpnn", "rf2"]
    print("pipeline rfantibody_wrapper2 selftest: pass")
    return 0


def main() -> int:
    ap = build_parser()
    ns = ap.parse_args()
    if ns.selftest:
        return selftest()

    cfg = resolve_config(ns)
    if not cfg.get("target_pdb") and sys.stdin.isatty():
        cfg["target_pdb"] = input("Path to target antigen PDB: ").strip()
    if not cfg.get("framework_pdb") and sys.stdin.isatty():
        cfg["framework_pdb"] = input("Path to VHH framework PDB: ").strip()

    if not cfg.get("target_pdb") or not cfg.get("framework_pdb"):
        print("Error: --target_pdb and --framework_pdb are required.", file=sys.stderr)
        return 2

    target = Path(cfg["target_pdb"])
    if not target.exists():
        print(f"Error: target PDB not found: {target}", file=sys.stderr)
        return 1

    fw_path = Path(cfg["framework_pdb"])
    out_dir = Path(cfg["output_dir"]) if cfg.get("output_dir") else (fw_path.parent / fw_path.stem.replace(".", "_"))
    out_dir.mkdir(parents=True, exist_ok=True)
    save_path = Path(ns.save_config) if ns.save_config else (out_dir / "config.json")
    save_path.parent.mkdir(parents=True, exist_ok=True)
    save_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")

    if cfg.get("hotspot_res"):
        ok, rep = validate_hotspots(target, str(cfg["hotspot_res"]), float(cfg.get("min_hotspot_sasa") or 15.0))
        print(json.dumps({"target": str(target), "hotspots_ok": ok, "hotspots": rep}, indent=2))
        if not ok:
            print("Hotspot validation failed prior to diffusion.", file=sys.stderr)
            return 2
    if ns.validate_only:
        return 0

    stages = build_commands(cfg)
    if ns.dry_run:
        print(json.dumps([{"stage": s, "out": str(p), "cmd": " ".join(c)} for s, p, c in stages], indent=2))
        return 0

    for name, out_pdb, cmd in stages:
        out_pdb.parent.mkdir(parents=True, exist_ok=True)
        if ns.resume and out_pdb.exists():
            print(f"[resume] skipping {name}: {out_pdb} exists")
            continue
        rc = stream_command(cmd, out_pdb.parent / f"{name}.log")
        if rc != 0:
            print(f"Stage {name} exited with code {rc}", file=sys.stderr)
            return rc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
