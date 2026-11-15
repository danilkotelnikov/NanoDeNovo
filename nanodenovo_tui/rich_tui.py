"""Cross-platform Rich TUI dashboard and raw single-key controller for NanoDeNovo."""
from __future__ import annotations

import ctypes
import io
import sys
import threading
from pathlib import Path
from typing import Any

from .deps import render_doctor_rich, run_doctor
from .instructions import render_instructions_rich
from .workflow import WorkflowEngine, configure_utf8


def is_interactive_tty() -> bool:
    """Return True only when stdin is a genuine interactive console (Win32 NUL-safe)."""
    try:
        if not sys.stdin.isatty():
            return False
        if sys.platform == "win32":
            import msvcrt

            handle = msvcrt.get_osfhandle(sys.stdin.fileno())
            mode = ctypes.c_uint32()
            return ctypes.windll.kernel32.GetConsoleMode(handle, ctypes.byref(mode)) != 0
        return True
    except (OSError, ValueError, io.UnsupportedOperation):
        return False


def read_single_key() -> str:
    """Read one keypress without requiring Enter; restores terminal attributes in finally."""
    if not is_interactive_tty():
        ch = sys.stdin.read(1)
        return ch if ch else "q"

    if sys.platform == "win32":
        import msvcrt

        ch = msvcrt.getwch()
        if ch in ("\x00", "\xe0"):
            _ = msvcrt.getwch()
            return ""
        if ch == "\x03":
            raise KeyboardInterrupt
        return ch.lower()

    import termios
    import tty

    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setcbreak(fd)
        ch = sys.stdin.read(1)
        if ch == "\x03":
            raise KeyboardInterrupt
        return ch.lower()
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)


def _status_badge(status: str, ascii_only: bool) -> str:
    labels = {
        "pending": "IDLE" if ascii_only else "[dim]IDLE[/dim]",
        "running": "RUN" if ascii_only else "[bold yellow]RUN[/bold yellow]",
        "done": "DONE" if ascii_only else "[bold green]DONE[/bold green]",
        "failed": "FAIL" if ascii_only else "[bold red]FAIL[/bold red]",
        "cancelled": "CANCEL" if ascii_only else "[magenta]CANCEL[/magenta]",
        "timeout": "TIMEOUT" if ascii_only else "[red]TIMEOUT[/red]",
    }
    return labels.get(status, status.upper())


def build_dashboard_renderable(
    engine: WorkflowEngine,
    mode: str,
    resume: bool,
    stage_states: dict[str, dict[str, Any]],
    log_lines: list[str],
    ledger_status: str,
    ascii_only: bool = False,
) -> Any:
    """Build a Rich renderable representing the full 8-stage NanoDeNovo dashboard."""
    from rich import box
    from rich.console import Group
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text

    box_style = box.ASCII if ascii_only else box.ROUNDED
    header_text = Text.assemble(
        ("NanoDeNovo Workflow Console", "bold cyan"),
        f"  |  Mode: {mode}  |  Resume: {'on' if resume else 'off'}  |  Ledger: {ledger_status}\n",
        f"Workspace: {engine.repo_root}",
    )
    header = Panel(header_text, box=box_style)

    table = Table(title="Pipeline Stages (S1-S8)", box=box_style, expand=True)
    table.add_column("Key", width=4, justify="center")
    table.add_column("Stage", width=5)
    table.add_column("Name", min_width=24)
    table.add_column("Status", width=9, justify="center")
    table.add_column("Elapsed", width=9, justify="right")
    table.add_column("Detail / Artifact Check", ratio=1)

    for idx, spec in enumerate(engine.specs, start=1):
        st = stage_states.get(spec.stage_id, {})
        raw_status = str(st.get("status", "pending"))
        elapsed = st.get("elapsed_s")
        elapsed_str = f"{elapsed:.2f}s" if isinstance(elapsed, (int, float)) else "-"
        detail = str(st.get("detail", spec.short_desc))
        if st.get("resumed"):
            detail = f"[resumed] {detail}"
        table.add_row(
            str(idx),
            spec.stage_id,
            spec.title,
            _status_badge(raw_status, ascii_only=ascii_only),
            elapsed_str,
            detail[:110],
        )

    tail = log_lines[-8:] if log_lines else ["Ready. Press 1-8 to run one stage, or 'a' to run all stages."]
    log_panel = Panel("\n".join(tail), title="Event Log", box=box_style)
    keys_panel = Panel(
        "[1-8] Run Stage   [a] Run All   [m] Cycle Mode   [r] Toggle Resume   "
        "[d] Doctor   [i] Instructions   [v] Verify Ledger   [q] Quit",
        title="Controls",
        box=box_style,
    )
    return Group(header, table, log_panel, keys_panel)


def render_dashboard_text(
    engine: WorkflowEngine,
    mode: str = "replay",
    resume: bool = False,
    stage_states: dict[str, dict[str, Any]] | None = None,
    log_lines: list[str] | None = None,
    width: int = 120,
    ascii_only: bool = True,
) -> str:
    """Render one static dashboard frame as a string for non-TTY output and unit tests."""
    from io import StringIO
    from rich.console import Console

    buf = StringIO()
    console = Console(file=buf, width=width, force_terminal=False, no_color=True)
    ledger_check = engine.verify_ledger()
    ledger_status = "VERIFIED" if ledger_check.get("ok") else "UNVERIFIED"
    renderable = build_dashboard_renderable(
        engine=engine,
        mode=mode,
        resume=resume,
        stage_states=stage_states or {},
        log_lines=log_lines or [],
        ledger_status=ledger_status,
        ascii_only=ascii_only,
    )
    console.print(renderable)
    return buf.getvalue()


def run_rich_tui(repo_root: Path, initial_mode: str = "replay", resume: bool = False) -> int:
    """Launch the Rich interactive workflow console (or print a single frame if stdin is not a TTY)."""
    from rich.console import Console

    ascii_only = configure_utf8()
    engine = WorkflowEngine(repo_root)
    console = Console()

    stage_states: dict[str, dict[str, Any]] = engine._load_existing_ledger()
    log_lines: list[str] = []
    mode_order = ("replay", "dry-run", "live")
    mode = initial_mode if initial_mode in mode_order else "replay"

    if not is_interactive_tty():
        console.print(
            render_dashboard_text(
                engine=engine,
                mode=mode,
                resume=resume,
                stage_states=stage_states,
                log_lines=["Non-interactive stdin detected; use --headless for automated execution."],
                ascii_only=ascii_only,
            )
        )
        return 0

    def _draw() -> None:
        ledger_check = engine.verify_ledger()
        ledger_status = "VERIFIED" if ledger_check.get("ok") else "UNVERIFIED"
        console.clear()
        console.print(
            build_dashboard_renderable(
                engine=engine,
                mode=mode,
                resume=resume,
                stage_states=stage_states,
                log_lines=log_lines,
                ledger_status=ledger_status,
                ascii_only=ascii_only,
            )
        )

    def _wait_any_key() -> None:
        console.print("\nPress any key to return to the dashboard...")
        try:
            read_single_key()
        except KeyboardInterrupt:
            pass

    while True:
        _draw()
        try:
            key = read_single_key()
        except KeyboardInterrupt:
            return 0

        if key == "q":
            return 0
        if key == "m":
            mode = mode_order[(mode_order.index(mode) + 1) % len(mode_order)]
            log_lines.append(f"Switched execution mode to '{mode}'.")
            continue
        if key == "r":
            resume = not resume
            log_lines.append(f"Toggled --resume={'on' if resume else 'off'}.")
            continue
        if key == "d":
            console.clear()
            render_doctor_rich(console, run_doctor(repo_root), ascii_only=ascii_only)
            _wait_any_key()
            continue
        if key == "i":
            console.clear()
            render_instructions_rich(console, repo_root, ascii_only=ascii_only)
            _wait_any_key()
            continue
        if key == "v":
            res = engine.verify_ledger()
            if res.get("ok"):
                log_lines.append(f"Ledger verified cleanly ({res.get('stages_checked')} stages).")
            else:
                log_lines.append(f"Ledger verification FAILED: {res.get('error') or res.get('mismatches')}")
            continue

        if key in {"1", "2", "3", "4", "5", "6", "7", "8"}:
            stage_id = f"S{key}"
            stage_states[stage_id] = {"status": "running", "detail": f"Executing {stage_id} ({mode})..."}
            _draw()
            cancel_ev = threading.Event()
            try:
                rec = engine.run_stage(
                    stage_id,
                    mode=mode,
                    resume=resume,
                    on_line=log_lines.append,
                    cancel_event=cancel_ev,
                )
            except KeyboardInterrupt:
                cancel_ev.set()
                rec = {"stage": stage_id, "status": "cancelled", "detail": "Interrupted by operator."}
            stage_states[stage_id] = rec
            log_lines.append(f"{stage_id} -> {rec['status']}: {rec.get('detail', '')}")
            continue

        if key == "a":
            cancel_ev = threading.Event()
            try:
                out = engine.run_all(
                    mode=mode,
                    resume=resume,
                    on_line=log_lines.append,
                    cancel_event=cancel_ev,
                )
                for rec in out["records"]:
                    stage_states[rec["stage"]] = rec
                log_lines.append(f"Run-all ({mode}) completed: ok={out['ok']}")
            except KeyboardInterrupt:
                cancel_ev.set()
                log_lines.append("Run-all cancelled by operator.")
            continue
