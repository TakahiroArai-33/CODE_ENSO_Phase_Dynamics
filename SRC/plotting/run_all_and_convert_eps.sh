#!/usr/bin/env bash
# Generate only paper panels. PNG and EPS are saved directly by Matplotlib.
if [ -z "${BASH_VERSION:-}" ]; then
    printf 'Run this script with bash.\n' >&2
    return 2 2>/dev/null || exit 2
fi
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
    printf 'Run this script with bash, not source.\n' >&2
    return 2
fi
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENV="${CONDA_ENV-virtual_env}"
if [[ -n "${CONDA_ENV}" ]]; then
    python_cmd=(conda run --no-capture-output -n "${CONDA_ENV}" python)
else
    python_cmd=(python)
fi
exec "${python_cmd[@]}" "${SCRIPT_DIR}/run_paper_figures.py" "$@"
