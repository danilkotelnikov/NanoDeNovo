"""Operator instructions reader and Rich formatter."""
from __future__ import annotations

from pathlib import Path
from typing import Any


def load_instructions(repo_root: Path) -> str:
    doc = repo_root / "docs" / "INSTRUCTIONS.md"
    if not doc.exists():
        return "# NanoDeNovo Operator Guide\n\nInstructions file docs/INSTRUCTIONS.md was not found."
    return doc.read_text(encoding="utf-8", errors="replace")


def list_sections(markdown_text: str) -> list[tuple[str, str]]:
    sections: list[tuple[str, str]] = []
    cur_title = "Overview"
    buf: list[str] = []
    for line in markdown_text.splitlines():
        if line.startswith("## "):
            if buf:
                sections.append((cur_title, "\n".join(buf).strip()))
            cur_title = line[3:].strip()
            buf = []
        else:
            buf.append(line)
    if buf:
        sections.append((cur_title, "\n".join(buf).strip()))
    return sections


def render_instructions_rich(console: Any, repo_root: Path, ascii_only: bool = False) -> None:
    from rich import box
    from rich.panel import Panel

    box_style = box.ASCII if ascii_only else box.ROUNDED
    text = load_instructions(repo_root)
    console.print(Panel(text, title="NanoDeNovo Built-In Operator Instructions (docs/INSTRUCTIONS.md)", box=box_style))
