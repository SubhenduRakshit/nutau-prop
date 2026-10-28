#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
printf '==== v0.20.4 normalization-corrected smoke started %s ====\n' "$(date)"
python scripts/smoke_mixed_expm_v0_20_4.py
printf '==== v0.20.4 smoke finished %s ====\n' "$(date)"
