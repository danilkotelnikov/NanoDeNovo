#!/usr/bin/env bash
set -euo pipefail

# Single-session batch ProteinMPNN runner (eliminates 1,000x checkpoint reloads
# and fixes the pmpnn_030 vs pmpnn_30 folder mismatch).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

TAG="${1:-20}"
SCAFFOLDS_DIR="${2:-${REPO_ROOT}/scaffolds}"
OUT_DIR="${3:-${REPO_ROOT}/pmpnn_0${TAG}}"
MPNN_DIR="${MPNN_DIR:-${REPO_ROOT}/third_party/ProteinMPNN}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
NUM_SEQS="${NUM_SEQS:-128}"
TEMP="${TEMP:-0.15}"
BATCH_SIZE="${BATCH_SIZE:-32}"
SEED="${SEED:-47832}"
EXTRA_ARGS=("${@:4}")

exec "${PYTHON_BIN}" "${REPO_ROOT}/tools/batch_mpnn.py" \
  --scaffolds "${SCAFFOLDS_DIR}" \
  --out-dir "${OUT_DIR}" \
  --tag "${TAG}" \
  --mpnn-dir "${MPNN_DIR}" \
  --python-bin "${PYTHON_BIN}" \
  --num-seqs "${NUM_SEQS}" \
  --temp "${TEMP}" \
  --batch-size "${BATCH_SIZE}" \
  --seed "${SEED}" \
  "${EXTRA_ARGS[@]}"
