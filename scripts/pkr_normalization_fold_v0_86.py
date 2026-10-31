#!/usr/bin/env python3
"""Fresh dynamic-Q CT18 folds and analytic coupling screens, v0.86.

No detector generation or Git side effects. All results are labelled screens.
"""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import json
import math
import platform
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.polynomial.legendre import leggauss
from nutau_prop.case6 import Case6Config, Case6DIS
from nutau_prop.constants import M_N, GEV2_TO_CM2
from nutau_prop.pkr_resolved import (nucleon_d2_from_xpdf, invisible_width,
    operator_factor, quark_width_coefficient, allowed_gq_squared)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/v0.86_normalization'
CALCS = {}
SMS = {}

def calc_for(n):
    if n not in CALCS:
        import lhapdf
        lhapdf.setVerbosity(0)
        CALCS[n] = Case6DIS(Case6Config(pdf_set='CT18NLO', nx=n))
    return CALCS[n]

def yrule(ymin,n):
    nodes, weights = leggauss(n)
    # Both integrals have smooth endpoints; logarithmic y resolves low Q^2.
    edges = (ymin, 0.08, 1.0)
    for lo,hi,log in [(edges[0],edges[1],True),(edges[1],edges[2],False)]:
        if lo >= hi:
            continue
        if log:
            aa,bb = math.log(lo),math.log(hi)
            yy = np.exp((nodes+1)*(bb-aa)/2+aa)
            ww = weights*(bb-aa)/2*yy
        else:
            yy = (nodes+1)*(hi-lo)/2+lo
            ww = weights*(hi-lo)/2
        yield from zip(yy,ww)

def integrate_np(E,M,a,cutoff,n=48,legacy_domain=False):
    calc=calc_for(n); s=2*M_N*E; qmin2=calc.cfg.qmin_gev**2
    xmax=min(1-1e-10,cutoff**2/s)
    ymin=max(qmin2/(s*xmax),1e-5 if legacy_domain else 0)
    if ymin>=1: return 0.0
    value=0.0
    for y,w in yrule(ymin,n):
        xmin=max(qmin2/(s*y),1e-9)
        if xmin>=xmax: continue
        def integrand(x):
            # Universal sum is isospin invariant: proton sum = isoscalar sum.
            q2=x*y*s
            xf=sum(calc.pdf.xfx_proton(pid,x,q2) for pid in range(-5,6) if pid)
            return nucleon_d2_from_xpdf(s,x,y,xf,M,a,1.0)*GEV2_TO_CM2
        value+=w*calc._integrate_logx(integrand,xmin,xmax)
    return value

def integrate_sm(E,n=48,legacy_domain=False):
    key=(E,n,legacy_domain)
    if key in SMS: return SMS[key]
    calc=calc_for(n); s=2*M_N*E
    ymin=max(calc.cfg.qmin_gev**2/s,1e-5 if legacy_domain else 0)
    value=np.zeros(2)
    for y,w in yrule(ymin,n):
        for anti in (False,True):
            value[int(anti)]+=w*sum(calc.dsdy_sm(E,y,antinu=anti))
    SMS[key]=value
    return value

def fold(task):
    E,M,a,cut,n=task
    sm=integrate_sm(E,n)
    np44=integrate_np(E,M,a,cut,n)
    return dict(E_GeV=E,M_GeV=M,a_GeV_inv=a,cutoff_GeV=cut,nquad=n,
                sigma44_gq1_cm2=np44,SM_nu_cm2=sm[0],SM_antinu_cm2=sm[1],
                ratio44_gq1=float(0.5*np44*(1/sm[0]+1/sm[1])))

def collider_refs():
    rank=pd.read_csv(ROOT/'results/v0.79_eq4p4_mass_acceptance_pilot/mass_ranking.csv')
    j=json.loads((ROOT/'results/v0.81_eq4p4_M50_same_seed_reweight/M50_same_seed_reweight.json').read_text())
    refs={}
    for M in (50.,100.,220.):
        r=rank[np.isclose(rank.M_GeV,M)].iloc[0]
        P=float(r.coupling_factor); rr=float(r.pilot_multibin_r)
        if M==50:
            P=j['representative']['coupling_factor_gq2BR']; rr=j['checkmate']['multibin_r']
        refs[M]=dict(K=rr/P,P_ref=P,r_ref=rr,source='v0.81 100k' if M==50 else 'v0.79 20k pilot')
    return refs

def screens(base):
    refs=collider_refs(); out=[]
    # None: relaxed necessary envelope (ignore dijet); .03: sufficient floor screen.
    # Approximate mass-specific values are sensitivity diagnostics, NOT measured limits.
    for operator in ('4.4','4.3'):
        fac=operator_factor(operator)
        for M in (50.,100.,220.):
            q=base[np.isclose(base.M_GeV,M)]
            for label,dijet in [('no_dijet_relaxed',None),('global_floor_sufficient',.03),
                                ('approx_mass_value',{50.:.038,100.:.083,220.:.116}[M])]:
                for rmax in (0.8,1.0,1.2):
                    candidates=[]
                    for row in q.itertuples():
                        a=row.a_GeV_inv
                        z=allowed_gq_squared(M,a,operator,rmax/refs[M]['K'],dijet)
                        I=invisible_width(M,a,operator); B=quark_width_coefficient(M)
                        P=z*I/(B*z+I) if z else 0.0
                        candidates.append(dict(operator=operator,M_GeV=M,dijet_screen=label,
                            r_cap=rmax,max_ratio100=z*fac*row.ratio44_gq1,
                            gq=math.sqrt(z),a_GeV_inv=a,cutoff_GeV=1/a,
                            Gamma_over_M=(B*z+I)/M,monojet_r=refs[M]['K']*P,
                            gq_equiv_dijet=math.sqrt(max(0,z-P)),
                            reference=refs[M]['source']))
                    out.append(max(candidates,key=lambda v:v['max_ratio100']))
    return pd.DataFrame(out)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--stage',choices=['initial','verify'],default='initial')
    args=parser.parse_args(); OUT.mkdir(exist_ok=True,parents=True)
    start=time.time(); workers=min(14,max(1,int(os.environ.get('NUTAU_NPROC','14'))))
    if args.stage=='initial':
        # Recompute all masses, rather than assuming old quadrature was adequate.
        tasks=[(1e8,M,float(a),1/float(a),48) for M in (50.,100.,220.)
               for a in np.geomspace(2e-5,.999/M,96)]
        output='base_folds.csv'
    else:
        envelope=pd.read_csv(OUT/'screen_maxima.csv')
        tasks=[]
        for row in envelope[np.isclose(envelope.r_cap,1.0)].itertuples():
            for n in (72,112):
                tasks.append((1e8,row.M_GeV,row.a_GeV_inv,1/row.a_GeV_inv,n))
        # Legacy named 5% points plus the most recent representatives.
        benchmarks=[('v0.73_44',50.,.03,.0035630660806766665,'4.4'),
                    ('v0.42_43',50.,.03,.0019725821667,'4.3'),
                    ('v0.81_44',50.,.038,.002379776,'4.4'),
                    ('v0.83_43',50.,.017,.006743708483124297,'4.3')]
        for name,M,g,a,op in benchmarks:
            for E in (1e8,2.2e8,3e8):
                for cut in (1/a,.5/a):
                    for n in (72,112): tasks.append((E,M,a,cut,n))
        # Exact equal-amplitude/common-cutoff equivalence needs no second PDF integral.
        tasks=list(dict.fromkeys(tasks)); output='verification_folds.csv'
    print(f'{args.stage}: {len(tasks)} targeted CT18 folds, {workers} workers; BLAS threads=1',flush=True)
    rows=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i,row in enumerate(pool.map(fold,tasks),1):
            rows.append(row)
            if i%24==0 or i==len(tasks): print(f'{i}/{len(tasks)} folds complete',flush=True)
    frame=pd.DataFrame(rows); frame.to_csv(OUT/output,index=False)
    if args.stage=='initial':
        summary=screens(frame); summary.to_csv(OUT/'screen_maxima.csv',index=False)
        print(summary[(summary.r_cap==1)&(summary.dijet_screen=='no_dijet_relaxed')].to_string(index=False))
    else:
        reports=[]
        for name,M,g,a,op in benchmarks:
            q=frame[np.isclose(frame.a_GeV_inv,a,rtol=1e-10,atol=0)]
            for r in q[q.nquad==112].itertuples():
                low=q[(q.E_GeV==r.E_GeV)&np.isclose(q.cutoff_GeV,r.cutoff_GeV)&(q.nquad==72)].iloc[0]
                reports.append(dict(benchmark=name,operator=op,E_GeV=r.E_GeV,gq=g,a_GeV_inv=a,
                    cutoff_GeV=r.cutoff_GeV,cutoff_factor=r.cutoff_GeV*a,
                    sigma_NP_cm2=g*g*operator_factor(op)*r.sigma44_gq1_cm2,
                    ratio_to_SM=g*g*operator_factor(op)*r.ratio44_gq1,
                    quadrature_relative_change=abs(r.ratio44_gq1/low.ratio44_gq1-1)))
        pd.DataFrame(reports).to_csv(OUT/'benchmarks.csv',index=False)
        print(pd.DataFrame(reports).to_string(index=False))
    metadata=dict(version='v0.86',stage=args.stage,workers=workers,thread_count_per_worker=1,
                  host=platform.node(),platform=platform.platform(),seconds=time.time()-start,
                  pdf='CT18NLO member 0',hard_process='LO',y_min='Qmin^2/min(s,Lambda^2)',
                  ratio_definition='mean of NP/SM_nu and NP/SM_antinu; SM=CC+NC',
                  collider='calibrated NWA screen, not an all-constraints certificate')
    (OUT/f'{args.stage}_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'Done in {time.time()-start:.1f}s',flush=True)

if __name__=='__main__': main()
