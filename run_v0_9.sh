#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python scripts/run_public_ehe_audit_v0_9.py | tee results/v0.9/public_ehe_audit_v0_9_summary.txt
