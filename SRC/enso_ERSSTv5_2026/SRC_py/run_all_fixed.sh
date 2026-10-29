#!/usr/bin/env bash
# Internal estimator-only entry point; ../run_all.sh also creates paper figures.
if [ -z "${BASH_VERSION:-}" ]; then
  printf 'Run with bash.\n' >&2
  return 2 2>/dev/null || exit 2
fi
if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  printf 'Run with bash, not source.\n' >&2
  return 2
fi
set -euo pipefail
check_inputs=false
case "${1:-}" in
  -h|--help)
    printf 'Usage: bash SRC_py/run_all_fixed.sh [--check-inputs]\n'
    printf 'Compute only the 15 paper NPZ estimates. For figures too, use ../run_all.sh.\n'
    exit 0
    ;;
  --check-inputs) check_inputs=true; shift ;;
esac
if [[ "$#" -ne 0 ]]; then
  printf 'Unknown estimator option: %s\nUse --help for usage.\n' "$1" >&2
  exit 2
fi
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
CONDA_ENV="${CONDA_ENV-virtual_env}"
if [[ -n "${CONDA_ENV}" ]]; then
  python_cmd=(conda run --no-capture-output -n "${CONDA_ENV}" python)
else
  python_cmd=(python)
fi
missing=false
for input in Z.mat Z_left.mat mu.mat Nino34_data.mat Lambda.mat; do
  if [[ ! -s "../myFolder/${input}" ]]; then
    printf 'Missing MATLAB data: myFolder/%s\n' "${input}" >&2
    missing=true
  fi
done
for script in Estimate_Q_2D_fixed.py Averaging.py Estimate_DirectAveraging_fixed.py \
  common_fixed.py Module_Bayes_Estimation.py Module_Bayes_Routine.py \
  DirectAvaraging/Module_DirectAveraging_Core.py DirectAvaraging/Module_Estimator.py \
  DirectAvaraging/Module_Langevin.py DirectAvaraging/Module_routine.py; do
  if [[ ! -s "${script}" ]]; then
    printf 'Missing estimator source: SRC_py/%s\n' "${script}" >&2
    missing=true
  fi
done
[[ "${missing}" == false ]] || exit 1
# Import every estimator and its runtime dependencies before the first write.
"${python_cmd[@]}" -B -c '
import numpy, scipy, netCDF4
import common_fixed, Module_Bayes_Routine, Estimate_Q_2D_fixed, Averaging
from Estimate_DirectAveraging_fixed import load_direct_module
load_direct_module()
'
if [[ "${check_inputs}" == true ]]; then
  printf 'Estimator inputs and dependencies are ready.\n'
  exit 0
fi

printf '[1/4] Q estimation (M=3; amplitude filters 10-50%%)\n'
"${python_cmd[@]}" Estimate_Q_2D_fixed.py --M-list 3 --amp-ranks 10 20 30 40 50
printf '[2/4] Q averaging (3:1, amplitude filter 30%%)\n'
"${python_cmd[@]}" Averaging.py --orders 3 --amps 30 --ratio 1:3
printf '[3/4] Direct averaging (five ratios, amplitude filter 30%%)\n'
"${python_cmd[@]}" Estimate_DirectAveraging_fixed.py \
  --orders 3 --amp-ranks 30 --ratios 1:2 1:3 1:4 1:5 2:5 --log10-lambda-max 8
printf '[4/4] Direct averaging (3:1, remaining amplitude filters)\n'
"${python_cmd[@]}" Estimate_DirectAveraging_fixed.py \
  --orders 3 --amp-ranks 10 20 40 50 --ratios 1:3 --log10-lambda-max 8
printf 'Estimation complete; MATLAB data were kept unchanged.\n'
