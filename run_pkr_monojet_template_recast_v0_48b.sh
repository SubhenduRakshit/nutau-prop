#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p cache results/v0.48b

MIRROR="cache/LLP-CDFO"
if [ ! -d "$MIRROR/.git" ]; then
  git clone --depth 1 --filter=blob:none --sparse https://github.com/andlessa/LLP-CDFO.git "$MIRROR"
  git -C "$MIRROR" sparse-checkout set CMS-EXO-20-004/CMS_data/HEPData-ins1894408-v2-csv statisticalTools
else
  git -C "$MIRROR" pull --ff-only
  git -C "$MIRROR" sparse-checkout set CMS-EXO-20-004/CMS_data/HEPData-ins1894408-v2-csv statisticalTools
fi

python scripts/pkr_monojet_template_recast_v0_48b.py | tee results/v0.48b/run.log

git add results/v0.48b
git commit -m "results: v0.48b PKR CMS monojet numerical-stability recast" || true
git push

echo ""
echo "v0.48b complete; shell remains open."
