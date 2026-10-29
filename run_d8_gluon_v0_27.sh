#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results/v0.27

echo '=== v0.27 Q_l2G2D exact CT18NLO UHE screen ==='
python scripts/d8_gluon_uhe_v0_27.py \
  --pdf CT18NLO --lambda-tev 4 --unitarity-mode to_lambda \
  --xis 1,0.5,0.25 --energies-pev 1,10,100,300 \
  --nx 300 --ny 48 --outdir results/v0.27 \
  | tee results/v0.27/gluon_uhe_gate.txt

echo '=== provenance ==='
python - <<'PY' > results/v0.27/provenance.txt
import sys
print('python',sys.version.replace('\n',' '))
import lhapdf
print('lhapdf paths',lhapdf.paths())
p=lhapdf.mkPDF('CT18NLO',0)
print('CT18NLO q2min/q2max',p.q2Min,p.q2Max)
PY

echo '=== auto-commit results ==='
git add results/v0.27
if ! git diff --cached --quiet; then
  git commit -m 'results: v0.27 exact gluon d8 UHE screen'
  git push
else
  echo 'No new results to commit.'
fi

echo '=== v0.27 finished; Earth propagation NOT run ==='
