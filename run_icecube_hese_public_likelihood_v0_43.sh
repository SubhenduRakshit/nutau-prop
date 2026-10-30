#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p cache results/v0.43

if [ ! -d cache/HESE-7-year-data-release/.git ]; then
  git clone --depth 1 https://github.com/icecube/HESE-7-year-data-release.git cache/HESE-7-year-data-release
else
  git -C cache/HESE-7-year-data-release pull --ff-only
fi

python scripts/icecube_hese_public_likelihood_v0_43b.py | tee results/v0.43/run.log

git add results/v0.43
git commit -m "results: v0.43b IceCube HESE public-data likelihood" || true
git push
