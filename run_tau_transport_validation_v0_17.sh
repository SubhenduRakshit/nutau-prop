#!/usr/bin/env bash
set -euo pipefail
PYTHON="${PYTHON:-python3}"
NPROC="${NUTAU_NPROC:-12}"

echo "==== v0.17 finite-beta_tau validation started $(date) ===="
echo "workers=$NPROC"

# Historical beta_tau=0 regression. Reuse cached historical cross-section tables
# when present; only rebuild them if the cache is missing.
if [[ ! -f cache/historical_nu.npz ]]; then
  "$PYTHON" scripts/precompute_xsec.py --config config/historical.yaml \
    --workers "$NPROC" --out cache/historical_nu.npz --kind sm
fi
if [[ ! -f cache/historical_antinu.npz ]]; then
  "$PYTHON" scripts/precompute_xsec.py --config config/historical.yaml \
    --workers "$NPROC" --out cache/historical_antinu.npz --kind sm --antinu
fi

echo "[1/4] Rakshit-Reya beta_tau=0 regression"
"$PYTHON" scripts/validate_rakshit_reya.py --workers "$NPROC"
echo "[2/4] Dutta-Reno-Sarcevic beta_tau=0 regression"
"$PYTHON" scripts/validate_dutta_reno_sarcevic.py --workers "$NPROC"
echo "[3/4] combined historical gate"
"$PYTHON" scripts/run_historical_gate.py

echo "[4/4] independent finite-beta_tau analytic + convergence gate"
"$PYTHON" scripts/validate_tau_transport_v0_17.py --workers "$NPROC"

echo "==== v0.17 finite-beta_tau validation PASSED $(date) ===="
