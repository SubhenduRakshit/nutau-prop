#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p cache results/v0.50

MIRROR="cache/LLP-CDFO"
if [ ! -d "$MIRROR/.git" ]; then
  git clone --depth 1 --filter=blob:none --sparse https://github.com/andlessa/LLP-CDFO.git "$MIRROR"
fi

git -C "$MIRROR" pull --ff-only
git -C "$MIRROR" sparse-checkout set \
  CMS-EXO-20-004 \
  statisticalTools \
  results_dataFrames

python scripts/pkr_monojet_official_chain_validation_v0_50.py | tee results/v0.50/run.log

git add results/v0.50
git commit -m "results: v0.50 official CMS monojet recast-chain validation" || true
git push

echo ""
echo "v0.50 complete; shell remains open."
