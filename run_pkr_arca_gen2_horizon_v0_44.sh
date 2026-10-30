#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p results/v0.44

python scripts/pkr_arca_gen2_horizon_v0_44.py | tee results/v0.44/run.log

git add results/v0.44
git commit -m "results: v0.44 ARCA Gen2 absorption-horizon forecast" || true
git push

echo ""
echo "v0.44 complete; shell remains open."
