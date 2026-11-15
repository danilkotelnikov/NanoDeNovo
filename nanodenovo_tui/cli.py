"""CLI entrypoint and non-destructive self-test suite for NanoDeNovo TUI."""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path
from typing import Sequence

from .deps import render_doctor_rich, run_doctor
from .instructions import list_sections, load_instructions, render_instructions_rich
from .rich_tui import render_dashboard_text, run_rich_tui
from .textual_tui import is_textual_available, run_textual_tui
from .workflow import WorkflowEngine, configure_utf8


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="NanoDeNovo automated 8-stage workflow CLI and TUI.")
    parser.add_argument("--repo", type=Path, default=None, help="Repository root directory (defaults to parent of nanodenovo.py).")
    parser.add_argument("--ui", choices=("auto", "rich", "textual"), default="auto", help="Interactive UI backend.")
    parser.add_argument("--headless", action="store_true", help="Run non-interactively without starting the TUI.")
    parser.add_argument("--mode", choices=("replay", "dry-run", "live"), default="replay", help="Stage execution mode.")
    parser.add_argument("--stage", choices=tuple(f"S{i}" for i in range(1, 9)), default=None, help="Optional single stage ID in headless mode.")
    parser.add_argument("--resume", action="store_true", help="Skip stages whose recorded SHA-256 ledger hashes still match.")
    parser.add_argument("--doctor", action="store_true", help="Inspect Python packages, binaries, WSL2, and repository paths.")
    parser.add_argument("--instructions", action="store_true", help="Display built-in operator instructions from docs/INSTRUCTIONS.md.")
    parser.add_argument("--verify-ledger", action="store_true", help="Verify results/ledger.jsonl against current artifacts fail-closed.")
    parser.add_argument("--json", action="store_true", help="Print machine-readable JSON output.")
    parser.add_argument("--selftest", action="store_true", help="Run internal non-destructive checks and exit.")
    return parser


def _default_repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _run_selftest(repo_root: Path) -> int:
    engine = WorkflowEngine(repo_root)
    doc = run_doctor(repo_root)
    if not doc.get("ok"):
        sys.stderr.write(f"selftest failed: doctor reported not ok: {doc}\n")
        return 1

    # Dry-run across all 8 stages without modifying results/ledger.jsonl
    for spec in engine.specs:
        rec = engine.run_stage(spec.stage_id, mode="dry-run")
        if rec["status"] != "done":
            sys.stderr.write(f"selftest failed on dry-run {spec.stage_id}: {rec}\n")
            return 1

    # Verify non-TTY dashboard frame renders all 8 stages cleanly
    frame = render_dashboard_text(engine, mode="replay", ascii_only=True)
    for spec in engine.specs:
        if spec.stage_id not in frame:
            sys.stderr.write(f"selftest failed: {spec.stage_id} missing from dashboard render\n")
            return 1

    # Verify fail-closed replay check inside an isolated temporary directory
    with tempfile.TemporaryDirectory() as tmp:
        tmp_engine = WorkflowEngine(Path(tmp))
        bad = tmp_engine.run_stage("S1", mode="replay")
        if bad["status"] != "failed":
            sys.stderr.write(f"selftest failed: missing input did not fail-closed: {bad}\n")
            return 1

    print(f"nanodenovo_tui selftest: pass (stages={len(engine.specs)}, python={doc['python']['version']})")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    ascii_only = configure_utf8()
    parser = build_parser()
    args = parser.parse_args(argv)
    repo_root = (args.repo or _default_repo_root()).resolve()

    if args.selftest:
        return _run_selftest(repo_root)

    if args.doctor:
        report = run_doctor(repo_root)
        if args.json:
            print(json.dumps(report, indent=2, sort_keys=True))
        else:
            from rich.console import Console

            render_doctor_rich(Console(), report, ascii_only=ascii_only)
        return 0 if report.get("ok") else 1

    if args.instructions:
        text = load_instructions(repo_root)
        if args.json:
            sections = [{"title": t, "body": b} for t, b in list_sections(text)]
            print(json.dumps({"path": str(repo_root / "docs" / "INSTRUCTIONS.md"), "sections": sections}, indent=2))
        else:
            from rich.console import Console

            render_instructions_rich(Console(), repo_root, ascii_only=ascii_only)
        return 0

    engine = WorkflowEngine(repo_root)

    if args.verify_ledger:
        res = engine.verify_ledger()
        if args.json:
            print(json.dumps(res, indent=2, sort_keys=True))
        else:
            if res.get("ok"):
                print(f"Ledger OK: {res['stages_checked']} stages verified against {engine.ledger_path}")
            else:
                print(f"Ledger FAILED: {res.get('error') or res.get('mismatches')}")
        return 0 if res.get("ok") else 1

    if args.headless:
        def _printer(line: str) -> None:
            if not args.json:
                print(line, flush=True)

        if args.stage:
            rec = engine.run_stage(args.stage, mode=args.mode, resume=args.resume, on_line=_printer)
            out = {"ok": rec["status"] == "done", "mode": args.mode, "records": [rec], "ledger": str(engine.ledger_path)}
        else:
            out = engine.run_all(mode=args.mode, resume=args.resume, on_line=_printer)

        if args.json:
            print(json.dumps(out, indent=2, sort_keys=True))
        elif not out["ok"]:
            sys.stderr.write("Workflow failed in headless mode.\n")
        return 0 if out["ok"] else 1

    ui_choice = args.ui
    if ui_choice == "auto":
        ui_choice = "textual" if is_textual_available() else "rich"
    if ui_choice == "textual":
        return run_textual_tui(repo_root, initial_mode=args.mode, resume=args.resume)
    return run_rich_tui(repo_root, initial_mode=args.mode, resume=args.resume)
