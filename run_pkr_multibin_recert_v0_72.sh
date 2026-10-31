#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/v0.72_multibin_recert
python scripts/pkr_multibin_recert_v0_72.py | tee results/v0.72_multibin_recert/run.txt
git add results/v0.72_multibin_recert
git commit -m "results: v0.72 multibin recertification" || true
git push || true
