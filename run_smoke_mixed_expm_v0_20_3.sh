#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
printf '==== v0.20.3 mixed-grid smoke started %s ====\n' "$(date)"
python scripts/smoke_mixed_expm_v0_20_3.py
printf '==== v0.20.3 mixed-grid smoke finished %s ====\n' "$(date)"
