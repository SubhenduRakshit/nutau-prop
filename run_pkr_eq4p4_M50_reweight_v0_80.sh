#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/v0.80_eq4p4_M50_reweight cache/v0.80_eq4p4_M50_reweight

eval "$(conda shell.bash hook)"
set +u
conda activate "${PKR_CM_ENV:-nutau-checkmate-v051}"
set -u
source cache/v0.51/env.sh
export BROWSER=true

echo "v0.80: same-sample M50 reweight; no MG5/Pythia/new events"
python scripts/pkr_eq4p4_M50_reweight_checkmate_v0_80.py | tee results/v0.80_eq4p4_M50_reweight/run.txt

git add results/v0.80_eq4p4_M50_reweight
git commit -m "results: v0.80 Eq4.4 M50 same-sample reweight" || true
git push || true
