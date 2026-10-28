#!/bin/bash
set -euo pipefail
NPROC=${NUTAU_NPROC:-5}
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "==== v0.19.4 transport diagnostics started $(date) ===="
echo "workers=$NPROC; internal numeric threads capped at 1"
python -m py_compile \
  scripts/diagnose_transport_v0_19_4.py \
  scripts/build_case7_refined_grid_v0_19_4.py \
  scripts/diagnose_energy_grid_v0_19_4.py

python scripts/diagnose_transport_v0_19_4.py --nproc "$NPROC"

if [ "${NUTAU_RUN_ENERGY_GRID:-0}" = "1" ]; then
  GRID_NPROC=${NUTAU_GRID_NPROC:-12}
  echo "==== optional 181-point energy-grid stage; workers=$GRID_NPROC ===="
  python scripts/build_case7_refined_grid_v0_19_4.py --workers "$GRID_NPROC"
  python scripts/diagnose_energy_grid_v0_19_4.py
else
  echo "Energy-grid stage skipped. Set NUTAU_RUN_ENERGY_GRID=1 only if Stage 1-2 leave an unresolved plateau."
fi

echo "==== v0.19.4 diagnostics finished $(date) ===="
