"""Preflight dependency doctor reading pinned contracts in `config/dependencies.json`."""
from __future__ import annotations

import importlib.metadata
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


def parse_version_tuple(ver: str) -> tuple[int, ...]:
    parts: list[int] = []
    for tok in ver.replace("-", ".").split("."):
        digits = "".join(ch for ch in tok if ch.isdigit())
        if digits:
            parts.append(int(digits))
        else:
            break
    return tuple(parts) if parts else (0,)


def probe_python_package(import_name: str, dist_name: str, min_version: str) -> dict[str, Any]:
    spec_ok = importlib.util.find_spec(import_name) is not None
    version_str = ""
    if spec_ok:
        try:
            version_str = importlib.metadata.version(dist_name)
        except importlib.metadata.PackageNotFoundError:
            try:
                mod = __import__(import_name)
                version_str = str(getattr(mod, "__version__", "unknown"))
            except Exception:
                version_str = "unknown"
    meets_min = False
    if spec_ok and version_str and version_str != "unknown":
        meets_min = parse_version_tuple(version_str) >= parse_version_tuple(min_version)
    elif spec_ok:
        meets_min = True
    return {
        "name": dist_name,
        "import_name": import_name,
        "present": spec_ok,
        "version": version_str,
        "min_version": min_version,
        "meets_min": meets_min,
    }


def probe_wsl_status() -> dict[str, Any]:
    if sys.platform != "win32":
        return {"applicable": False, "available": False, "detail": "Non-Windows host"}
    wsl_bin = shutil.which("wsl.exe") or shutil.which("wsl")
    if not wsl_bin:
        return {"applicable": True, "available": False, "detail": "wsl.exe not found on PATH"}
    try:
        proc = subprocess.run(
            [wsl_bin, "--status"],
            capture_output=True,
            timeout=4,
            check=False,
        )
        raw = proc.stdout or proc.stderr or b""
        # Windows wsl.exe often emits UTF-16LE even in modern terminals
        if b"\x00" in raw:
            text = raw.decode("utf-16-le", errors="replace").replace("\x00", "")
        else:
            text = raw.decode("utf-8", errors="replace")
        summary = " | ".join(line.strip() for line in text.splitlines() if line.strip())[:160]
        return {"applicable": True, "available": proc.returncode == 0, "detail": summary or f"exit {proc.returncode}"}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"applicable": True, "available": False, "detail": f"{type(exc).__name__}"}


def run_doctor(repo_root: Path) -> dict[str, Any]:
    """Inspect Python runtime, packages, binaries, WSL2, and repository artifacts."""
    manifest_path = repo_root / "config" / "dependencies.json"
    if not manifest_path.exists():
        return {"ok": False, "error": f"Missing dependency manifest: {manifest_path}"}

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    py_ok = sys.version_info >= (3, 11)

    host_results: list[dict[str, Any]] = []
    for entry in manifest.get("host_packages", []):
        res = probe_python_package(entry["import_name"], entry["name"], entry["min_version"])
        res["required_for"] = entry.get("required_for", [])
        res["locked_version"] = entry.get("locked_version", "")
        host_results.append(res)

    optional_results: list[dict[str, Any]] = []
    for entry in manifest.get("optional_packages", []):
        res = probe_python_package(entry["import_name"], entry["name"], entry["min_version"])
        res["purpose"] = entry.get("purpose", "")
        optional_results.append(res)

    bin_results: list[dict[str, Any]] = []
    for entry in manifest.get("external_binaries", []):
        if entry.get("platform") and entry["platform"] != sys.platform:
            continue
        env_key = entry.get("env_override", "")
        env_val = os.environ.get(env_key, "") if env_key else ""
        resolved = env_val or shutil.which(entry["name"]) or ""
        bin_results.append(
            {
                "name": entry["name"],
                "env_override": env_key,
                "present": bool(resolved),
                "resolved": resolved,
                "required_for": entry.get("required_for", []),
            }
        )

    expected_dirs = [
        "1_virtual_screening",
        "2_sequence_and_structure_assessement",
        "pipeline",
        "tools",
        "results",
        "docs/INSTRUCTIONS.md",
        "AUDIT.md",
    ]
    path_checks = {rel: (repo_root / rel).exists() for rel in expected_dirs}

    # Replay + Dry-run require numpy, scipy, pandas, rich, and all repo directories
    replay_core = {"numpy", "scipy", "pandas", "rich"}
    core_ok = all(r["present"] and r["meets_min"] for r in host_results if r["name"] in replay_core)
    paths_ok = all(path_checks.values())
    live_host_ok = all(r["present"] and r["meets_min"] for r in host_results if "live" in r["required_for"])

    return {
        "ok": bool(py_ok and core_ok and paths_ok),
        "live_host_ready": bool(py_ok and live_host_ok and paths_ok),
        "isolated_mode": bool(getattr(sys.flags, "isolated", 0)),
        "python": {
            "version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
            "executable": sys.executable,
            "ok": py_ok,
        },
        "host_packages": host_results,
        "optional_packages": optional_results,
        "binaries": bin_results,
        "wsl": probe_wsl_status(),
        "repo_paths": path_checks,
        "manifest_version": manifest.get("schema_version", "unknown"),
    }


def render_doctor_rich(console: Any, report: dict[str, Any], ascii_only: bool = False) -> None:
    """Render structured doctor report using Rich tables."""
    from rich import box
    from rich.table import Table

    box_style = box.ASCII if ascii_only else box.ROUNDED
    status_tag = "READY" if report.get("ok") else "BLOCKED"
    console.print(f"[bold]NanoDeNovo Preflight Doctor[/bold] | Status: {status_tag} | Python {report.get('python', {}).get('version', '?')}")
    if report.get("isolated_mode"):
        console.print("[dim]Note: Python -I (isolated mode) skips user-site packages such as matplotlib if installed per-user.[/dim]")

    tbl = Table(title="Host & Optional Python Packages", box=box_style)
    tbl.add_column("Package")
    tbl.add_column("Version")
    tbl.add_column("Min Required")
    tbl.add_column("Status")
    tbl.add_column("Role")

    for pkg in report.get("host_packages", []):
        st = "OK" if (pkg["present"] and pkg["meets_min"]) else "MISSING"
        tbl.add_row(pkg["name"], pkg["version"] or "-", pkg["min_version"], st, ",".join(pkg.get("required_for", [])))
    for pkg in report.get("optional_packages", []):
        st = "OK" if pkg["present"] else "OPTIONAL-ABSENT"
        tbl.add_row(pkg["name"], pkg["version"] or "-", pkg["min_version"], st, pkg.get("purpose", "optional")[:42])
    console.print(tbl)

    bin_tbl = Table(title="External Tools & WSL2", box=box_style)
    bin_tbl.add_column("Component")
    bin_tbl.add_column("Status")
    bin_tbl.add_column("Resolved / Detail")
    for b in report.get("binaries", []):
        bin_tbl.add_row(b["name"], "OK" if b["present"] else "NOT-FOUND", b["resolved"] or b.get("env_override", ""))
    wsl = report.get("wsl", {})
    if wsl.get("applicable"):
        bin_tbl.add_row("WSL2", "OK" if wsl.get("available") else "UNAVAILABLE", wsl.get("detail", ""))
    console.print(bin_tbl)
