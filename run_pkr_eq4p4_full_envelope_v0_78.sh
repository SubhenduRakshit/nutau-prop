#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p results/v0.78_eq4p4_full_envelope

echo "v0.78: full-a Eq4.4 envelope audit; no MG5/Pythia/CheckMATE"
python scripts/pkr_eq4p4_full_envelope_v0_78.py | tee results/v0.78_eq4p4_full_envelope/run.txt

git add results/v0.78_eq4p4_full_envelope
git commit -m "results: v0.78 corrected Eq4.4 full-a envelope" || true
git push || true
