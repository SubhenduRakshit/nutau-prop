#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
OUT=results/v0.29; mkdir -p "$OUT" build/v0.29
N=${GLUON_D8_NEVENTS:-50000}
for x in python lhapdf rivet contur pythia8-config yodastack; do command -v "$x" >/dev/null || { echo "MISSING: $x"; exit 2; }; done
rivet --list-analyses | grep -q 'ATLAS_2025_I2905252' || { echo 'MISSING Rivet ATLAS_2025_I2905252'; exit 3; }
python - <<'PY'
import lhapdf
try: lhapdf.mkPDF('NNPDF31_lo_as_0130',0)
except Exception: raise SystemExit('Install NNPDF31_lo_as_0130 with: lhapdf install NNPDF31_lo_as_0130')
PY
python scripts/generate_gluon_d8_lhe_v0_29.py --events "$N" --out "$OUT/gluon_d8_C8pi.lhe" | tee "$OUT/generation.log"
c++ -O2 $(pythia8-config --cxxflags) scripts/shower_gluon_d8_pythia8_v0_29.cc -o build/v0.29/shower_d8 $(pythia8-config --libs --hepmc3)
build/v0.29/shower_d8 "$OUT/gluon_d8_C8pi.lhe" "$OUT/gluon_d8_C8pi.hepmc" | tee "$OUT/pythia.log"
rivet -a ATLAS_2025_I2905252 -o "$OUT/gluon_d8_C8pi.yoda" "$OUT/gluon_d8_C8pi.hepmc" | tee "$OUT/rivet.log"
contur "$OUT/gluon_d8_C8pi.yoda" -o "$OUT/contur_C8pi" --nomultip 2>&1 | tee "$OUT/contur_C8pi.log"
bash scripts/scan_gluon_d8_contur_v0_29.sh "$OUT/gluon_d8_C8pi.yoda" "$OUT/contur_scan"
echo '=== v0.29 collider recast complete; Earth propagation NOT run ==='
