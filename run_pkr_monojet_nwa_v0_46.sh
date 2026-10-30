#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p results/v0.46
python scripts/pkr_monojet_nwa_v0_46.py | tee results/v0.46/run.log

git add results/v0.46
git commit -m "results: v0.46 PKR monojet NWA public-limit screen" || true
git push

echo ""
echo "v0.46 complete; shell remains open."
