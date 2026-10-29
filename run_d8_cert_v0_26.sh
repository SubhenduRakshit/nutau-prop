#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/v0.26

echo '=== v0.26 analytic normalization / positivity / unitarity ==='
python scripts/check_d8_analytic_v0_25.py | tee results/v0.26/analytic_gate.txt

echo '=== v0.26 exact dynamic-Q CT18NLO fold, strict C=-8pi/3 ==='
python scripts/d8_uhe_cert_v0_25.py --pdf CT18NLO --lambda-tev 4 --unitarity-mode to_lambda \
  --xis 1,0.5,0.25 --energies-pev 1,10,100,300 --nx 300 --ny 48 --outdir results/v0.26 \
  | tee results/v0.26/uhe_gate.txt

echo '=== v0.26 exact HighPT Run-II di-tau recast ==='
if command -v wolframscript >/dev/null 2>&1; then
  wolframscript -file scripts/highpt_d8_tautau_v0_26.wls | tee results/v0.26/highpt_gate.txt
else
  echo 'SKIP HighPT: wolframscript not found' | tee results/v0.26/highpt_gate.txt
fi

echo '=== v0.26 provenance ==='
python - <<'PY' > results/v0.26/provenance.txt
import sys
print('python',sys.version.replace('\n',' '))
try:
 import lhapdf
 print('lhapdf paths',lhapdf.paths())
 p=lhapdf.mkPDF('CT18NLO',0)
 print('CT18NLO q2min/q2max',p.q2Min,p.q2Max)
except Exception as e:
 print('lhapdf/CT18NLO ERROR',repr(e))
PY

echo '=== auto-commit results ==='
git add results/v0.26
if ! git diff --cached --quiet; then
  git commit -m 'results: v0.26 d8 certification gates'
  git push
else
  echo 'No new results to commit.'
fi

echo '=== v0.26 finished; Earth propagation NOT run ==='
