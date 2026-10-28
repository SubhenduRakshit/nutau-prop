#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "==== v0.19 expm candidate validation started $(date) ===="
"$PYTHON" -m py_compile src/nutau_prop/chirality_transport_expm.py scripts/validate_expm_transport_v0_19.py
"$PYTHON" scripts/validate_expm_transport_v0_19.py --resume
echo "==== v0.19 expm candidate validation PASSED $(date) ===="
