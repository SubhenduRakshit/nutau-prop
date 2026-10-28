#!/usr/bin/env bash
set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUTAU_NPROC=${NUTAU_NPROC:-12}

echo "==== v0.19.2 expm candidate validation started $(date) ===="
echo "workers=${NUTAU_NPROC}; internal numeric threads capped at 1"
python -m py_compile src/nutau_prop/chirality_transport_expm.py scripts/validate_expm_transport_v0_19_2.py
python scripts/validate_expm_transport_v0_19_2.py --nproc "${NUTAU_NPROC}"
