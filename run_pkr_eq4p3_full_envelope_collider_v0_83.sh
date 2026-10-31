#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/v0.83_eq4p3_full_envelope_collider

echo "v0.83: Eq4.3 full-a UHE + collider screen; no MG5/Pythia/Delphes/CheckMATE"
python scripts/pkr_eq4p3_full_envelope_collider_v0_83.py | tee results/v0.83_eq4p3_full_envelope_collider/run.txt

git add results/v0.83_eq4p3_full_envelope_collider
git commit -m "results: v0.83 Eq4.3 full-a no-new-events collider screen" || true
git push || true
