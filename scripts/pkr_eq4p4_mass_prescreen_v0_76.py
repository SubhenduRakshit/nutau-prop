#!/usr/bin/env python3
"""v0.76 — Eq.4.4 resolved-mediator mass prescreen.

Scans M = 50, 100, 220 GeV. For each mass and gq, solve the smallest a=c/Lambda
that gives sigma_NP/sigma_SM = 5% at 100 PeV using the same resolved CT18NLO
kernel and sqrt(shat)<1/a validity cut as v0.74. Then apply the physical
O+h.c. invisible width, Gamma_inv=a^2 M^3/(24 pi), narrow-width and
Lambda_vertex=1/a > M requirements.

The dijet input is a digitized SCREEN from the official CMS EXO-24-007
Figure 2a observed curve (BR(Z'->qq)=1), not an exact HEPData likelihood.
Values are intentionally rounded and accompanied by a +/-0.005 coupling band:
  M=50:  gq95 ~ 0.038
  M=100: gq95 ~ 0.083
  M=220: gq95 ~ 0.116
A point is called robust-dijet-pass only if gq*sqrt(BRqq) < gq95-0.005.

No cross-mass monojet exclusion is inferred. We only rank surviving points by
P_monojet = gq^2 BRinv, since production x invisible BR scales this way in NWA;
acceptance is mass-dependent and requires one fresh CheckMATE validation.
"""
from __future__ import annotations

import json, math, os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from nutau_prop.case6 import Case6Config, Case6DIS
from nutau_prop.constants import GEV2_TO_CM2, M_N

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.76_eq4p4_mass_prescreen"
OUT.mkdir(parents=True, exist_ok=True)

PDF = "CT18NLO"
TARGET_E = 1.0e8
TARGET_RATIO = 0.05
MASSES = (50.0, 100.0, 220.0)
CMS_GQ95 = {50.0: 0.038, 100.0: 0.083, 220.0: 0.116}
CMS_DIGITIZATION_ERR = 0.005
GQ_VALUES = np.arange(0.005, 0.16001, 0.001)
NA = 144
QUARK_MASS = {"u":0.00216,"d":0.00467,"s":0.0934,"c":1.27,"b":4.18,"t":172.69}


def gamma_qq_one(M, gq, mq):
    if M <= 2*mq: return 0.0
    r=(mq/M)**2
    return 3*gq*gq*M/(12*math.pi)*(1+2*r)*math.sqrt(1-4*r)


class Eq4p4DIS(Case6DIS):
    def __init__(self, M):
        super().__init__(Case6Config(pdf_set=PDF, nx=52, y_log_bins=22, y_linear_bins=22))
        self.M=float(M)

    def _d2_base(self,E,x,y,a):
        s=2*M_N*E; q2=x*y*s
        if q2 < self.cfg.qmin_gev**2: return 0.0
        p=self._pdf_components(x,q2)
        c4=a/(self.M*self.M)  # gq=1, kappa=a
        val=(4/math.pi)*c4*c4*x*s*s*y*(1-y)*sum(p.values())
        val*= (self.M*self.M/(self.M*self.M+q2))**2
        return val*GEV2_TO_CM2

    def dsnp_base(self,E,y,a):
        s=2*M_N*E
        xmin=max(self.cfg.qmin_gev**2/(s*y),1e-9)
        xmax=min(1-1e-10,(1/a)**2/s)
        if xmin>=xmax: return 0.0
        return self._integrate_logx(lambda x:self._d2_base(E,x,y,a),xmin,xmax)


def envelope_for_mass(M):
    calc=Eq4p4DIS(M)
    agrid=np.geomspace(2e-5,0.999/M,NA)
    _,ys,dy=calc.y_grid()
    sm_n=sm_a=0.0
    for y,d in zip(ys,dy):
        cc,nc=calc.dsdy_sm(TARGET_E,y,antinu=False); sm_n+=(cc+nc)*d
        cc,nc=calc.dsdy_sm(TARGET_E,y,antinu=True);  sm_a+=(cc+nc)*d
    vals=[]
    for a in agrid:
        snp=sum(calc.dsnp_base(TARGET_E,y,float(a))*d for y,d in zip(ys,dy))
        vals.append(0.5*(snp/max(sm_n,1e-300)+snp/max(sm_a,1e-300)))
    return M,agrid,np.asarray(vals)


def first_crossing(a,f,target):
    idx=np.where((f[:-1]<target)&(f[1:]>=target))[0]
    if not len(idx): return None
    i=int(idx[0]); x0,x1=np.log(a[i]),np.log(a[i+1]); y0,y1=np.log(f[i]),np.log(f[i+1])
    t=(math.log(target)-y0)/(y1-y0)
    return float(math.exp(x0+t*(x1-x0)))


def main():
    nproc=min(max(1,int(os.environ.get("NUTAU_NPROC","14"))),len(MASSES))
    print(f"[1/2] CT18 mass envelopes: {len(MASSES)} masses on {nproc} processes",flush=True)
    env={}
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        for M,a,f in ex.map(envelope_for_mass,MASSES): env[M]=(a,f)

    base=[]; rows=[]
    for M in MASSES:
        agrid,f=env[M]
        for a,r in zip(agrid,f): base.append(dict(M_GeV=M,a_GeV_inv=a,ratio100_at_gq1=r))
        fmax=float(f.max()); gqcrit=math.sqrt(TARGET_RATIO/fmax)
        for gq in GQ_VALUES:
            a=first_crossing(agrid,f,TARGET_RATIO/(gq*gq))
            if a is None: continue
            gq_width=sum(gamma_qq_one(M,gq,mq) for mq in QUARK_MASS.values())
            ginv=a*a*M**3/(24*math.pi); gtot=gq_width+ginv
            brinv=ginv/gtot; brqq=gq_width/gtot; gqeq=gq*math.sqrt(brqq)
            lim=CMS_GQ95[M]
            rows.append(dict(M_GeV=M,gq=gq,a_GeV_inv=a,kappa_GeV_inv=gq*a,
                Lambda_vertex_GeV=1/a,Gamma_over_M=gtot/M,BR_inv=brinv,BR_qq=brqq,
                gq_equiv_dijet=gqeq,cms_gq95_digitized=lim,
                dijet_margin=lim-gqeq,
                dijet_central_pass=bool(gqeq<lim),
                dijet_robust_pass=bool(gqeq<lim-CMS_DIGITIZATION_ERR),
                narrow_width_pass=bool(gtot/M<0.10),vertex_scale_pass=bool(1/a>M),
                monojet_nwa_proxy=gq*gq*brinv,gq_5pct_capacity_threshold=gqcrit))

    pd.DataFrame(base).to_csv(OUT/"base_mass_envelopes.csv",index=False)
    df=pd.DataFrame(rows); df.to_csv(OUT/"mass_prescreen.csv",index=False)
    valid=df[df.narrow_width_pass & df.vertex_scale_pass & df.dijet_robust_pass].copy()
    summary=[]
    for M in MASSES:
        q=valid[valid.M_GeV==M]
        allm=df[df.M_GeV==M]
        if len(q):
            b=q.loc[q.monojet_nwa_proxy.idxmin()]
            summary.append(dict(M_GeV=M,status="ROBUST_DIJET_SCREEN_WINDOW",n_valid=len(q),
                gq=float(b.gq),a_GeV_inv=float(b.a_GeV_inv),Lambda_vertex_GeV=float(b.Lambda_vertex_GeV),
                BR_inv=float(b.BR_inv),Gamma_over_M=float(b.Gamma_over_M),
                gq_equiv_dijet=float(b.gq_equiv_dijet),cms_gq95_digitized=float(b.cms_gq95_digitized),
                monojet_nwa_proxy=float(b.monojet_nwa_proxy)))
        elif len(allm):
            b=allm.loc[allm.dijet_margin.idxmax()]
            summary.append(dict(M_GeV=M,status="NO_ROBUST_DIJET_SCREEN_WINDOW",n_valid=0,
                best_dijet_margin=float(b.dijet_margin)))
        else:
            summary.append(dict(M_GeV=M,status="NO_5PCT_CAPACITY",n_valid=0))
    S=pd.DataFrame(summary); S.to_csv(OUT/"mass_summary.csv",index=False)

    eligible=S[S.status=="ROBUST_DIJET_SCREEN_WINDOW"].copy()
    candidate=None
    if len(eligible):
        candidate=eligible.loc[eligible.monojet_nwa_proxy.idxmin()].to_dict()
    meta={"version":"v0.76_eq4p4_mass_prescreen","target":"5% at 100 PeV",
          "cms_source":"CMS EXO-24-007 Figure 2a observed gq curve; digitized/rounded screen only",
          "cms_digitization_error_assumed":CMS_DIGITIZATION_ERR,
          "next_checkmate_candidate":candidate,
          "rule":"Run one fresh CheckMATE point only for the mass with smallest gq^2 BRinv among robust-dijet screen survivors; do not infer cross-mass monojet exclusion from NWA proxy."}
    (OUT/"summary.json").write_text(json.dumps(meta,indent=2)+"\n")

    print("=== v0.76 Eq4.4 mass prescreen ===")
    print(S.to_string(index=False))
    print("next CheckMATE candidate:")
    print(candidate if candidate is not None else "NONE")
    return 0

if __name__=="__main__": raise SystemExit(main())
