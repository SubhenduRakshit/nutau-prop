#!/usr/bin/env python3
"""v0.84 — high-resolution PKR Eq.4.3 M=50 threshold-edge certification.

Purpose
-------
v0.83 used dgq=0.001 and found its best 50-GeV point at gq=0.017 with
r_MB=0.8234.  Because the 5%-UHE capacity threshold is near gq=0.0165, that
coarse grid can miss the true minimum collider visibility.

This script performs NO MG5/Pythia/Delphes/CheckMATE generation.  It directly
re-evaluates the CT18 resolved Eq.4.4 base fold near its M=50 maximum; Eq.4.3 is
exactly twice that cross section at equal (M,gq,a).  At fixed a, the smallest
collider signal compatible with sigma_NP/sigma_SM >= 5% occurs at the equality
boundary

    gq_min(a)^2 = 0.05 / [2 f43base(a)],

because P=gq^2 BR_inv increases monotonically with gq^2 at fixed a.  We then
compute the physical O+h.c. dipole width Gamma_inv=a^2 M^3/(6 pi), the exact
quark width, the global-minimum dijet sufficient screen, and rescale the
canonical deterministic v0.81 M50 multibin result linearly in P=gq^2 BR_inv.

The principal question is whether there is a *finite* local region with the
predeclared MC-screen upper edge r_hi < 0.8.  A threshold-edge point within
sub-percent numerical distance of 0.8 is reported as FRAGILE, not promoted to
a robust physical benchmark.
"""
from __future__ import annotations

import importlib.util, json, math, os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd

from nutau_prop.case6 import Case6Config, Case6DIS
from nutau_prop.constants import GEV2_TO_CM2, M_N

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'v0.84_eq4p3_M50_threshold_refine'
OUT.mkdir(parents=True,exist_ok=True)
V81=ROOT/'results'/'v0.81_eq4p4_M50_same_seed_reweight'/'M50_same_seed_reweight.json'
R79=ROOT/'results'/'v0.79_eq4p4_mass_acceptance_pilot'/'mass_ranking.csv'

PDF='CT18NLO'; M=50.0; TARGET=0.05
# Local interval brackets the v0.74 base-envelope maximum at a~9.106e-3.
A_MIN=0.0065; A_MAX=0.0125; NA=int(os.environ.get('PKR_V084_NA','320'))
GLOBAL_DIJET_FLOOR=0.03
QUARK_MASS={'u':0.00216,'d':0.00467,'s':0.0934,'c':1.27,'b':4.18,'t':172.69}

class Eq44Base(Case6DIS):
    def _d2(self,E,x,y,a):
        s=2*M_N*E; q2=x*y*s
        if q2<self.cfg.qmin_gev**2: return 0.0
        p=self._pdf_components(x,q2)
        c4=a/(M*M)  # gq=1
        val=(4/math.pi)*c4*c4*x*s*s*y*(1-y)*sum(p.values())
        val*= (M*M/(M*M+q2))**2
        return val*GEV2_TO_CM2
    def dsnp(self,E,y,a):
        s=2*M_N*E
        xmin=max(self.cfg.qmin_gev**2/(s*y),1e-9)
        xmax=min(1-1e-10,(1/a)**2/s)
        if xmin>=xmax: return 0.0
        return self._integrate_logx(lambda x:self._d2(E,x,y,a),xmin,xmax)

def fold_chunk(avals):
    calc=Eq44Base(Case6Config(pdf_set=PDF,nx=64,y_log_bins=28,y_linear_bins=28))
    _,ys,dy=calc.y_grid(); E=1e8
    smn=sma=0.0
    for y,d in zip(ys,dy):
        cc,nc=calc.dsdy_sm(E,y,antinu=False); smn+=(cc+nc)*d
        cc,nc=calc.dsdy_sm(E,y,antinu=True);  sma+=(cc+nc)*d
    out=[]
    for a in avals:
        snp=sum(calc.dsnp(E,float(y),float(a))*d for y,d in zip(ys,dy))
        f44=0.5*(snp/max(smn,1e-300)+snp/max(sma,1e-300))
        out.append((float(a),float(f44)))
    return out

def gamma_q(gq):
    s=0.0
    for mq in QUARK_MASS.values():
        if M<=2*mq: continue
        r=(mq/M)**2
        s += 3*gq*gq*M/(12*math.pi)*(1+2*r)*math.sqrt(1-4*r)
    return s

def metrics(a,f44,rref,Pref,mc):
    # Eq4.3 ratio = 2*gq^2*f44; choose equality boundary for minimum P at fixed a.
    if f44<=0: return None
    g2=TARGET/(2*f44); gq=math.sqrt(g2)
    ginv=a*a*M**3/(6*math.pi)
    gqwid=gamma_q(gq); gtot=ginv+gqwid
    br=ginv/gtot; brqq=gqwid/gtot
    P=g2*br
    r=rref*P/Pref; rhi=r*(1+mc); rlo=r*(1-mc)
    gqeq=gq*math.sqrt(brqq)
    return dict(a_GeV_inv=a,f44_base=f44,eq43_base=2*f44,gq_threshold=gq,
                Gamma_inv_GeV=ginv,Gamma_q_GeV=gqwid,Gamma_over_M=gtot/M,
                BR_inv=br,BR_qq=brqq,coupling_factor=P,
                monojet_r=r,monojet_r_mc_lo=rlo,monojet_r_mc_hi=rhi,
                gq_equiv_dijet=gqeq,dijet_global_floor_pass=bool(gqeq<GLOBAL_DIJET_FLOOR),
                vertex_scale_GeV=1/a,vertex_pass=bool(1/a>M),narrow_width_pass=bool(gtot/M<0.10))

def main():
    if not V81.exists() or not R79.exists(): raise RuntimeError('missing v0.81 or v0.79 input')
    j=json.loads(V81.read_text()); rank=pd.read_csv(R79)
    rref=float(j['checkmate']['multibin_r']); Pref=float(j['representative']['coupling_factor_gq2BR'])
    mc=float(rank[np.isclose(rank.M_GeV,50.0)].iloc[0].pilot_mc_scale)

    agrid=np.geomspace(A_MIN,A_MAX,NA)
    nproc=min(max(1,int(os.environ.get('NUTAU_NPROC','14'))),NA)
    chunks=[x for x in np.array_split(agrid,nproc) if len(x)]
    print(f'v0.84: direct CT18 local refinement, {NA} a-points on {nproc} processes',flush=True)
    vals=[]
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        for z in ex.map(fold_chunk,chunks): vals.extend(z)
    vals.sort()
    rows=[]
    for a,f in vals:
        m=metrics(a,f,rref,Pref,mc)
        if m: rows.append(m)
    D=pd.DataFrame(rows)
    D.to_csv(OUT/'local_refinement.csv',index=False)
    valid=D[D.vertex_pass & D.narrow_width_pass & D.dijet_global_floor_pass]
    if valid.empty: raise RuntimeError('no locally valid threshold points')
    b=valid.loc[valid.monojet_r.idxmin()]
    robust=valid[valid.monojet_r_mc_hi<0.8]

    # Quantify how finite the sub-0.8 region is rather than relying on a single grid point.
    finite_width=False; alo=ahi=np.nan
    if len(robust)>=2:
        alo=float(robust.a_GeV_inv.min()); ahi=float(robust.a_GeV_inv.max())
        finite_width=(ahi/alo-1)>1e-3

    margin=0.8-float(b.monojet_r_mc_hi)
    # Sub-0.8 by <1% in r_hi is deliberately labelled fragile.
    if margin<=0:
        status='NO_ROBUST_SURVIVAL_AT_CONTINUOUS_THRESHOLD'
    elif margin<0.008 or not finite_width:
        status='FRAGILE_THRESHOLD_EDGE_SUB08_NOT_CERTIFIED_ROBUST'
    else:
        status='FINITE_ROBUST_SURVIVAL_REGION_FOUND'

    summary={
      'version':'v0.84_eq4p3_M50_threshold_refine',
      'method':'direct CT18 local fold; Eq4.3=2xEq4.4; analytic gq threshold at each a; no new collider events',
      'reference_r_v081':rref,'reference_P_v081':Pref,'mc_scale':mc,
      'best':{k:(bool(v) if isinstance(v,(bool,np.bool_)) else float(v)) for k,v in b.to_dict().items()},
      'n_points_rhi_lt_0p8':int(len(robust)),
      'robust_a_interval':[alo,ahi] if len(robust) else None,
      'robust_interval_fractional_width':float(ahi/alo-1) if len(robust)>=2 else 0.0,
      'margin_0p8_minus_best_rhi':margin,
      'status':status,
      'interpretation':'A sub-0.8 threshold-edge point is not promoted if the margin is <0.008 (1% of the threshold) or the surviving a-region is not finite at >0.1% width.'
    }
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('=== v0.84 Eq4.3 M50 threshold-edge refinement ===')
    print(f"a_best={b.a_GeV_inv:.10g}  gq_crit={b.gq_threshold:.10g}")
    print(f"BRinv={b.BR_inv:.9f}  Gamma/M={b.Gamma_over_M:.6g}  gqeq_dijet={b.gq_equiv_dijet:.6g}")
    print(f"P={b.coupling_factor:.10g}  r={b.monojet_r:.6f}  r_mc_hi={b.monojet_r_mc_hi:.6f}")
    print(f"points with r_mc_hi<0.8: {len(robust)}; a interval={None if len(robust)==0 else (alo,ahi)}")
    print(f"margin 0.8-r_mc_hi(best)={margin:.6f}")
    print(f"status: {status}")
    print('new collider events: NONE')
    return 0

if __name__=='__main__': raise SystemExit(main())
