#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "==== v0.20 combined-generator smoke started $(date) ===="
python -m py_compile src/nutau_prop/chirality_transport_full_expm.py scripts/smoke_full_expm_v0_20.py
python scripts/smoke_full_expm_v0_20.py
echo "==== v0.20 combined-generator smoke finished $(date) ===="
