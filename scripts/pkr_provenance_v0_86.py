#!/usr/bin/env python3
"""Record source/input identities for the v0.86 reproducibility checkpoint."""
import hashlib
import json
import platform
import subprocess
from pathlib import Path
import lhapdf
import numpy
import scipy
import sympy

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/v0.86_normalization'
files=[ROOT/p for p in (
    'src/nutau_prop/pkr_resolved.py','src/nutau_prop/case6.py','src/nutau_prop/pdfs.py',
    'src/nutau_prop/constants.py','scripts/check_pkr_traces_v0_86.py',
    'scripts/pkr_normalization_fold_v0_86.py','scripts/pkr_refine_screens_v0_86.py',
    'results/v0.76_eq4p4_mass_prescreen/base_mass_envelopes.csv',
    'results/v0.79_eq4p4_mass_acceptance_pilot/mass_ranking.csv',
    'results/v0.81_eq4p4_M50_same_seed_reweight/M50_same_seed_reweight.json')]
pdf_files=[]
for directory in lhapdf.paths():
    for name in ('CT18NLO.info','CT18NLO_0000.dat'):
        p=Path(directory)/'CT18NLO'/name
        if p.exists(): pdf_files.append(p)
    if pdf_files: break
records=[]
for p in files+pdf_files:
    records.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
result=dict(version='v0.86',pre_commit_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
    host=platform.node(),python=platform.python_version(),numpy=numpy.__version__,scipy=scipy.__version__,
    sympy=sympy.__version__,lhapdf=lhapdf.version(),lhapdf_module=lhapdf.__file__,files=records,
    analytic_status='PASS: specified distinct-state massless chiral vertices',
    numerical_status='PASS: targeted quadrature and interpolation gates',
    collider_status='conditional NWA screens, no new experimental exclusion',
    unresolved=['original PKR UFO and field/multiplicity mapping','exact mass-specific dijet data',
                'operator-specific collider validity and off-shell effects','tau-specific IceCube likelihood',
                'UV completion and physical cutoff','updated coupled transport/detector prediction'])
(OUT/'provenance.json').write_text(json.dumps(result,indent=2)+'\n')
print('Recorded',len(records),'source/PDF identities in',OUT/'provenance.json')
