#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p results/v0.45
python scripts/pkr_width_br_dijet_v0_45.py | tee results/v0.45/run.log
git add results/v0.45
git commit -m "results: v0.45 PKR width BR dijet prescreen" || true
git push
echo ""
echo "v0.45 complete; shell remains open."
