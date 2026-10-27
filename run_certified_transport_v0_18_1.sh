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

echo "[0/5] syntax + editable-package metadata preflight"
"$PYTHON" -m py_compile \
  scripts/reassess_tau_transport_v0_17_1.py \
  scripts/validate_chirality_transport_v0_18_1.py \
  scripts/certify_transport_v0_18_1.py \
  scripts/run_certified_transport_v0_18_1.py
"$PYTHON" -m pip install -e . --no-deps -q
"$PYTHON" -c 'import nutau_prop; print("nutau_prop version", nutau_prop.__version__)'

# Reassess the superseded v0.17 refinement-trend heuristic from saved results.
echo "[1/5] v0.17.1 refinement reassessment"
"$PYTHON" scripts/reassess_tau_transport_v0_17_1.py

# Validate the exact production solver classes with Ehlert and pointwise UHE checks.
echo "[2/5] chirality-aware production-solver validation"
"$PYTHON" scripts/validate_chirality_transport_v0_18_1.py --workers "$NPROC"

# Composite gate: historical + v0.17 substantive + v0.17.1 + chirality solver.
echo "[3/5] composite numerical certification"
"$PYTHON" scripts/certify_transport_v0_18_1.py

# Only now allow production transport. No detector folding is done here.
echo "[4/5] certified contact + resolved-mediator transport"
"$PYTHON" scripts/run_certified_transport_v0_18_1.py --workers "$NPROC"

echo "[5/5] complete; detector weighting intentionally not run"
echo "==== v0.18.1 certified transport completed $(date) ===="
