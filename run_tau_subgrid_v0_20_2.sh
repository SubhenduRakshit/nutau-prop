#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
printf '==== v0.20.2 tau subgrid diagnostic started %s ====\n' "$(date)"
python -m py_compile scripts/diagnose_tau_subgrid_v0_20_2.py
python scripts/diagnose_tau_subgrid_v0_20_2.py
printf '==== v0.20.2 tau subgrid diagnostic finished %s ====\n' "$(date)"
