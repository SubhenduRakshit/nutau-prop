#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python scripts/pkr_eq4p4_width_resolution_v0_73.py
git add results/v0.73_eq4p4_width_resolution
git commit -m "results: v0.73 Eq4.4 width normalization resolution" || true
git push || true
echo "v0.73 complete: no MG5/Pythia/CheckMATE rerun"
