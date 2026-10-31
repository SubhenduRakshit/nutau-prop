#!/usr/bin/env bash
set -euo pipefail
PKR_ROOT="$(cd "$(dirname "$0")" && pwd)"
PKR_PYTHON="${PKR_PYTHON:-/Users/subhendurakshit/miniforge3/bin/python}"
export PYTHONPATH="$PKR_ROOT/src"
export NUTAU_NPROC=14
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$PKR_ROOT"
mkdir -p results/v0.87_constraints
"$PKR_PYTHON" scripts/pkr_exact_constraints_v0_87.py > results/v0.87_constraints/run.txt 2>&1
"$PKR_PYTHON" scripts/check_pkr_transport_identity_v0_87.py > results/v0.87_constraints/transport_identity_run.txt
printf '%s\n' 'v0.87 exact-input collider envelope and IceCube input audit saved; no full IceCube likelihood claimed.'
