#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "==== v0.20.7 bounded coupled transport gate started $(date) ===="
python scripts/validate_coupled_v0_20_7.py
echo "==== v0.20.7 bounded coupled transport gate finished $(date) ===="
