#!/usr/bin/env bash
# Recompute the paper estimates, then publish the manuscript's PNG/EPS figures.
if [ -z "${BASH_VERSION:-}" ]; then
  printf 'Run with bash.\n' >&2
  return 2 2>/dev/null || exit 2
fi
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  printf 'Run with bash, not source.\n' >&2
  return 2
fi
set -euo pipefail

DATASET_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${DATASET_DIR}/../.." && pwd)"
FIGURE_RUNNER="${ROOT_DIR}/SRC/plotting/run_all_and_convert_eps.sh"
figures_only=false
show_help=false
figure_args=()
for arg in "$@"; do
  case "${arg}" in
    --figures-only) figures_only=true ;;
    --collect-from|--collect-from=*|--dry-run|--preflight-estimation)
      figures_only=true
      figure_args+=("${arg}")
      ;;
    -h|--help)
      show_help=true
      figure_args+=("${arg}")
      ;;
    *) figure_args+=("${arg}") ;;
  esac
done
# Bash 3.2 (macOS) treats an empty array expansion as unset under set -u.
if [[ "${#figure_args[@]}" -gt 0 ]]; then
  set -- "${figure_args[@]}"
else
  set --
fi

if [[ "${show_help}" == true ]]; then
  cat <<'HELP'
Usage: bash SRC/enso_ERSSTv5_2026/run_all.sh [--figures-only] [figure options]

Default: recompute the 15 paper estimates from myFolder/*.mat, then generate
the manuscript figures in SRC/figure/. MATLAB data are not changed.
Download NOAA SST separately first: bash SRC/download_sst_data.sh

  --figures-only  Plot using the existing intermediate NPZ files.
  --output-dir   Use a dedicated output directory (relative to current cwd).
  --collect-from Collect existing plots without recomputing or plotting.

CONDA_ENV defaults to virtual_env; set CONDA_ENV="" to use the current python.
SDE_TRIAL_COUNT defaults to 5000. Figure-runner options follow:
HELP
  exec bash "${FIGURE_RUNNER}" --help
fi

if [[ "${figures_only}" == false ]]; then
  # Check plotting inputs, dependencies, and options before changing estimates.
  bash "${FIGURE_RUNNER}" --preflight-estimation "$@"
  bash "${DATASET_DIR}/SRC_py/run_all_fixed.sh"
fi
exec bash "${FIGURE_RUNNER}" "$@"
