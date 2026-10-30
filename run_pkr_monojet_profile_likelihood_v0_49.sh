#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p results/v0.49

python scripts/pkr_monojet_profile_likelihood_v0_49.py | tee results/v0.49/run.log

git add results/v0.49
git commit -m "results: v0.49 independent CMS monojet profile likelihood" || true
git push

echo ""
echo "v0.49 complete; shell remains open."
