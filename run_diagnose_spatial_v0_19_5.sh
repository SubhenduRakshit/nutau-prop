#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export NUTAU_NPROC="${NUTAU_NPROC:-2}"
python -m py_compile scripts/diagnose_spatial_convergence_v0_19_5.py
printf '==== v0.19.5 spatial diagnostic started %s ====\n' "$(date)"
printf 'workers=%s; internal numeric threads capped at 1\n' "$NUTAU_NPROC"
python scripts/diagnose_spatial_convergence_v0_19_5.py
printf '==== v0.19.5 spatial diagnostic finished %s ====\n' "$(date)"
