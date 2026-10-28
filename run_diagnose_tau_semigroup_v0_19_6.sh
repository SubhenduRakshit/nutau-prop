#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
python -m py_compile scripts/diagnose_tau_semigroup_v0_19_6.py
printf '==== v0.19.6 tau semigroup diagnostic started %s ====\n' "$(date)"
python scripts/diagnose_tau_semigroup_v0_19_6.py
printf '==== v0.19.6 tau semigroup diagnostic finished %s ====\n' "$(date)"
