#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export BROWSER=true
mkdir -p results/v0.64_mg5_norm
python scripts/pkr_recertify_mg5_normalization_v0_64.py | tee results/v0.64_mg5_norm/run.log
RC=${PIPESTATUS[0]}
git add results/v0.64_mg5_norm
git commit -m "results: v0.64 MG5 normalization re-certification" || true
git push || true
exit "$RC"
