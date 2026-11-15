#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
MODE="check"

for arg in "$@"; do
  case "$arg" in
    --check-only) MODE="check" ;;
    --install) MODE="install" ;;
    --install-conda) MODE="conda" ;;
    *)
      echo "Unknown option: $arg" >&2
      echo "Usage: $0 [--check-only | --install | --install-conda]" >&2
      exit 2
      ;;
  esac
done

cd "${REPO_ROOT}"

if [ "${MODE}" = "check" ]; then
  exec python3 -I "${REPO_ROOT}/nanodenovo.py" --doctor
elif [ "${MODE}" = "install" ]; then
  python3 -m pip install -r "${REPO_ROOT}/requirements-lock.txt"
  exec python3 -I "${REPO_ROOT}/nanodenovo.py" --doctor
elif [ "${MODE}" = "conda" ]; then
  conda env update -f "${REPO_ROOT}/environment.wsl.yml" --prune || conda env create -f "${REPO_ROOT}/environment.wsl.yml"
  exec conda run -n nanodenovo-wsl python -I "${REPO_ROOT}/nanodenovo.py" --doctor
fi
