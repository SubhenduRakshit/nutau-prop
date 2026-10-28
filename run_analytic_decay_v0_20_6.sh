#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
printf '==== v0.20.6 analytic-decay diagnostic started %s ====\n' "$(date)"
python -m py_compile src/nutau_prop/analytic_decay_projection.py scripts/diagnose_analytic_decay_v0_20_6.py
python scripts/diagnose_analytic_decay_v0_20_6.py
printf '==== v0.20.6 analytic-decay diagnostic finished %s ====\n' "$(date)"
