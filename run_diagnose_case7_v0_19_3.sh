#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUTAU_NPROC=${NUTAU_NPROC:-5}
python -m py_compile scripts/diagnose_case7_old_vs_expm_v0_19_3.py
python scripts/diagnose_case7_old_vs_expm_v0_19_3.py
