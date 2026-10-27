#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
python - <<'PY'
try:
 import lhapdf
except Exception as e:
 raise SystemExit("LHAPDF is required. On macOS: conda install -c conda-forge lhapdf, then: lhapdf install CT10nlo")
PY
python scripts/run_case6_pdf_cut_v0_8.py | tee results/v0.8/case6_pdf_cut_v0_8_summary.txt
python scripts/run_uhe_events_v0_8.py | tee results/v0.8/uhe_event_counts_v0_8_summary.txt
