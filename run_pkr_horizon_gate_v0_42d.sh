#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${NUTAU_NPROC:=6}"
python scripts/pkr_horizon_gate_mac_v0_42d.py | tee results/v0.42d/run.log

git add results/v0.42d
git commit -m "results: v0.42d Mac Studio PKR absorption-horizon gate" || true
git push
