#!/usr/bin/env python3
"""Repaired entrypoint for dataset 20 sequence selection."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from process_sequences import run_cli  # noqa: E402


def main() -> int:
    return run_cli("20")


if __name__ == "__main__":
    raise SystemExit(main())
