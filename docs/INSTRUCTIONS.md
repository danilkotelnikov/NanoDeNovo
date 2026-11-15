# NanoDeNovo Operator Guide

This guide explains how to check predetermined dependencies and run the automated 8-stage NanoDeNovo workflow on Windows and WSL2.
The date of this specification is 15.XI.2026.

## Workflow Architecture

The workflow runs eight stages against target structure `PDB 8e8z_C` (`DOI 10.3390/ijms26199262`).

| Stage | Purpose | Core Script |
|---|---|---|
| S1 | Check backbone hotspots and SASA exposure | `tools/design_rfantibody.py` |
| S2 | Batch inverse folding with one model load | `tools/batch_mpnn.py` |
| S3 | Select lowest-NLL unique sequences per backbone | `tools/select_sequences.py` |
| S4 | Score inclusive monomer refolds across 1..1000 | `tools/fold_ab.py` |
| S5 | Align every folded design and measure CDR3 | `tools/align_all.py` and `tools/cdr3_check.py` |
| S6 | Parse complex REMARK metrics and interface geometry | `tools/fold_ag.py` and `tools/interface_check.py` |
| S7 | Run OpenMM minimization and clash audit | `2_sequence_and_structure_assessement/openmm_em/openmm_minimizer.py` |
| S8 | Score Youden thresholds and write SHA-256 ledger | `tools/qc_report.py` |

## Predetermined Dependencies

Host Python requires version 3.11 or newer.
The dependency manifest sits in `config/dependencies.json` and `requirements-lock.txt`.
GPU checkpoints and OpenMM live inside WSL2 (`environment.wsl.yml`).

Run preflight inspection without downloading packages.

```powershell
py -3 -I nanodenovo.py --doctor
```

Run dependency installation only when an operator approves downloads.

```powershell
powershell -ExecutionPolicy Bypass -File scripts\setup_deps.ps1 -InstallHost
```

## Execution Modes

Choose one of three execution modes before starting a run.

- `replay`: Verify repository artifacts on CPU and build a portable SHA-256 ledger (`results/ledger.jsonl`). Missing or corrupt files fail closed.
- `dry-run`: Validate paths, contracts, and quoted WSL commands without writing output structures.
- `live`: Run stages locally or across WSL2 with `--resume` input/output hash verification.

## Interactive TUI Controls

Start the interactive terminal interface from the repository root.

```powershell
py -3 nanodenovo.py
```

Pass `--ui rich` for the built-in Rich dashboard or `--ui textual` when the optional Textual package is installed.

- Keys `1` to `8`: Select one stage (`S1` through `S8`).
- Key `r`: Run the selected stage under the active mode.
- Key `a`: Run all eight stages in dependency order.
- Key `m`: Cycle active mode (`replay`, `dry-run`, `live`).
- Key `d`: Open the Preflight Doctor panel.
- Key `i`: Open the Operator Instructions panel.
- Key `v`: Verify every SHA-256 hash inside `results/ledger.jsonl`.
- Key `q`: Restore terminal settings and exit cleanly.

## Headless Automation

Use `--headless` for non-interactive shells, scheduled jobs, and continuous integration tests.

```powershell
py -3 -I nanodenovo.py --headless --mode replay
py -3 -I nanodenovo.py --verify-ledger
py -3 -I nanodenovo.py --selftest
```

## Failure Recovery

Press `Ctrl+C` during a long run to terminate child process groups cleanly.
Cancelled runs record `status: cancelled` in the stage summary.
Inspect the missing file or mismatched checksum reported by `--verify-ledger`.
Run the built-in test suite to verify every script after changes.

```powershell
py -3 -m pytest tests
```
