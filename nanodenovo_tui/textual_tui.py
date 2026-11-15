"""Optional Textual full-screen TUI wrapper for NanoDeNovo (lazy-imported)."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

from .deps import run_doctor
from .instructions import load_instructions
from .workflow import WorkflowEngine, configure_utf8


def is_textual_available() -> bool:
    """Return True when the optional `textual` UI dependency is importable."""
    return importlib.util.find_spec("textual") is not None


def run_textual_tui(repo_root: Path, initial_mode: str = "replay", resume: bool = False) -> int:
    """Run the optional Textual dashboard, or explain clearly when `textual` is absent."""
    configure_utf8()
    if not is_textual_available():
        sys.stderr.write(
            "Optional UI dependency 'textual' is not installed in this Python environment.\n"
            "Use the built-in Rich dashboard instead: py -3 nanodenovo.py --ui rich\n"
            "Or install optional Textual: py -3 -m pip install 'textual>=0.85.0'\n"
        )
        return 2

    from textual.app import App, ComposeResult
    from textual.binding import Binding
    from textual.widgets import DataTable, Footer, Header, Log, Static

    engine = WorkflowEngine(repo_root)
    modes = ("replay", "dry-run", "live")

    class NanoDeNovoTextualApp(App[int]):
        TITLE = "NanoDeNovo Automated Pipeline"
        BINDINGS = [
            Binding("a", "run_all", "Run All"),
            Binding("m", "cycle_mode", "Cycle Mode"),
            Binding("d", "show_doctor", "Doctor"),
            Binding("i", "show_instructions", "Instructions"),
            Binding("v", "verify_ledger", "Verify Ledger"),
            Binding("q", "quit", "Quit"),
        ]

        def __init__(self) -> None:
            super().__init__()
            self.mode = initial_mode if initial_mode in modes else "replay"
            self.resume = resume
            self.states = engine._load_existing_ledger()

        def compose(self) -> ComposeResult:
            yield Header()
            yield Static(id="status_bar")
            yield DataTable(id="stages_table")
            yield Log(id="event_log")
            yield Footer()

        def on_mount(self) -> None:
            tbl = self.query_one("#stages_table", DataTable)
            tbl.add_columns("Stage", "Name", "Status", "Elapsed", "Detail")
            self._refresh_table()

        def _refresh_table(self) -> None:
            bar = self.query_one("#status_bar", Static)
            bar.update(f"Mode: {self.mode} | Resume: {self.resume} | Repo: {repo_root}")
            tbl = self.query_one("#stages_table", DataTable)
            tbl.clear()
            for spec in engine.specs:
                st = self.states.get(spec.stage_id, {})
                el = st.get("elapsed_s")
                el_s = f"{el:.2f}s" if isinstance(el, (int, float)) else "-"
                tbl.add_row(spec.stage_id, spec.title, str(st.get("status", "pending")), el_s, str(st.get("detail", spec.short_desc)))

        def action_cycle_mode(self) -> None:
            self.mode = modes[(modes.index(self.mode) + 1) % len(modes)]
            self.query_one("#event_log", Log).write_line(f"Mode set to {self.mode}")
            self._refresh_table()

        def action_run_all(self) -> None:
            log = self.query_one("#event_log", Log)
            res = engine.run_all(mode=self.mode, resume=self.resume, on_line=log.write_line)
            for rec in res["records"]:
                self.states[rec["stage"]] = rec
            self._refresh_table()

        def action_show_doctor(self) -> None:
            rep = run_doctor(repo_root)
            self.query_one("#event_log", Log).write_line(f"Doctor ok={rep.get('ok')} live_ready={rep.get('live_host_ready')}")

        def action_show_instructions(self) -> None:
            txt = load_instructions(repo_root)
            self.query_one("#event_log", Log).write_line(txt[:800])

        def action_verify_ledger(self) -> None:
            res = engine.verify_ledger()
            self.query_one("#event_log", Log).write_line(f"Ledger verify: {res}")

    NanoDeNovoTextualApp().run()
    return 0
