#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "==== v0.20.5 decay-z quadrature diagnostic started $(date) ===="
python scripts/diagnose_decay_quadrature_v0_20_5.py
echo "==== v0.20.5 decay-z quadrature diagnostic finished $(date) ===="
