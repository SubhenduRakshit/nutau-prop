#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/v0.71_multibin_output_audit
python scripts/pkr_multibin_output_audit_v0_71.py | tee results/v0.71_multibin_output_audit/run.txt
git add results/v0.71_multibin_output_audit
git commit -m "results: v0.71 multibin output audit" || true
git push || true
