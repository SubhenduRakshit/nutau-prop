#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"
NPROC="${NUTAU_NPROC:-12}"

# Bound BLAS/OpenMP threading inside each multiprocessing worker.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "==== v0.18.1 certified transport started $(date) ===="
echo "workers=$NPROC; internal numeric threads capped at 1"

# Reassess the superseded v0.17 refinement-trend heuristic from saved results.
"$PYTHON" scripts/reassess_tau_transport_v0_17_1.py

# Validate the exact production solver classes with Ehlert and pointwise UHE checks.
"$PYTHON" scripts/validate_chirality_transport_v0_18_1.py --workers "$NPROC"

# Composite gate: historical + v0.17 substantive + v0.17.1 + chirality solver.
"$PYTHON" scripts/certify_transport_v0_18_1.py

# Only now allow production transport. No detector folding is done here.
"$PYTHON" scripts/run_certified_transport_v0_18_1.py --workers "$NPROC"

echo "==== v0.18.1 certified transport completed $(date) ===="
