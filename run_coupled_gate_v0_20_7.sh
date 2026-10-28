#!/bin/bash
set -euo pipefail

# Independent configurations run in separate processes.  Keep each process
# single-threaded internally to avoid BLAS/OpenMP oversubscription.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUTAU_NPROC="${NUTAU_NPROC:-6}"

echo "==== v0.20.7 multicore coupled transport gate started $(date) ===="
echo "workers=$NUTAU_NPROC; internal numeric threads=1"
python scripts/validate_coupled_v0_20_7.py
echo "==== v0.20.7 multicore coupled transport gate finished $(date) ===="
