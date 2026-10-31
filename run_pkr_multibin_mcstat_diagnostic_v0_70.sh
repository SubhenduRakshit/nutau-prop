#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/v0.70_multibin_mcstat
python scripts/pkr_multibin_mcstat_diagnostic_v0_70.py | tee results/v0.70_multibin_mcstat/run.txt
git add results/v0.70_multibin_mcstat
git commit -m "results: v0.70 multibin MC-stat diagnostic" || true
git push || true
