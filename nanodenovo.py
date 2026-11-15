#!/usr/bin/env python3
"""Root launcher for the automated 8-stage NanoDeNovo workflow CLI and TUI."""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("PYTHONUTF8", "1")
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from nanodenovo_tui.cli import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
