#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p results/v0.49c

python scripts/pkr_monojet_profile_likelihood_v0_49c.py | tee results/v0.49c/run.log

git add results/v0.49c
git commit -m "results: v0.49c PKR CMS monojet profile likelihood" || true
git push

echo ""
echo "v0.49c complete; shell remains open."
