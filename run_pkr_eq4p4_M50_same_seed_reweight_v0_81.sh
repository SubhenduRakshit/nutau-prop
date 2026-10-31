#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/v0.81_eq4p4_M50_same_seed_reweight cache/v0.81_eq4p4_M50_same_seed_reweight

eval "$(conda shell.bash hook)"
set +u
conda activate "${PKR_CM_ENV:-nutau-checkmate-v051}"
set -u
source cache/v0.51/env.sh
export BROWSER=true

echo "v0.81: deterministic same-seed M50 reweight; no MG5/Pythia/new events"
python scripts/pkr_eq4p4_M50_same_seed_reweight_v0_81.py | tee results/v0.81_eq4p4_M50_same_seed_reweight/run.txt

git add results/v0.81_eq4p4_M50_same_seed_reweight docs/PKR_EQ4P4_V080_SEED_AUDIT_v0.81.md
git commit -m "results: v0.81 deterministic Eq4.4 M50 same-seed reweight" || true
git push || true
