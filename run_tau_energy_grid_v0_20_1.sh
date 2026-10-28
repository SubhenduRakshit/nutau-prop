#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUTAU_NPROC="${NUTAU_NPROC:-3}"
echo "==== v0.20.1 tau energy-grid diagnostic started $(date) ===="
echo "workers=$NUTAU_NPROC; internal numeric threads capped at 1"
python -m py_compile scripts/diagnose_tau_energy_grid_v0_20_1.py
python scripts/diagnose_tau_energy_grid_v0_20_1.py
echo "==== v0.20.1 tau energy-grid diagnostic finished $(date) ===="
