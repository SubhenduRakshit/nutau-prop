#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p results/v0.47

python scripts/pkr_monojet_hepdata_v0_47.py | tee results/v0.47/run.log

git add results/v0.47
git commit -m "results: v0.47 exact CMS HEPData 50 GeV monojet gate" || true
git push

echo ""
echo "v0.47 complete; shell remains open."
