#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUTAU_NPROC="${NUTAU_NPROC:-6}"

echo "==== v0.20.9 reuse-first coupled gate started $(date) ===="
echo "workers=$NUTAU_NPROC; internal numeric threads=1"
python scripts/validate_coupled_v0_20_9.py
echo "==== v0.20.9 reuse-first coupled gate finished $(date) ===="
