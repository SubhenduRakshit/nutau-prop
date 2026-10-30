#!/usr/bin/env python3
"""v0.41 constrained Earth propagation for final PKR Eq.(4.4)/(4.3) options.

Physics choices
---------------
* resolved derivative Eq.(4.4) and transition-dipole Eq.(4.3)
* sterile/Dirac final state => active-flux sink (no BSM regeneration)
* separate quark coupling gq and dimension-5 vertex coefficient
* validity domain sqrt(shat) < Lambda/c = gq/kappa, imposed inside every PDF integral
* exact resolved propagator 1/(M^2+Q^2)^2
* Eq.(4.3) hard kernel = 2 x Eq.(4.4) for equal (M,kappa)
* CT18NLO PDFs
* full PREM nu_tau<->tau transport with tau decay, energy loss, regeneration,
  and polarized secondary nu_mu/antinu_mu using the validated mixed-grid solver.

The current CMS low-mass dijet publication quotes observed gq limits spanning
0.03--0.13 for 50--250 GeV. We therefore use gq=0.03 as a deliberately
conservative common benchmark and solve kappa self-consistently for a 5% NP/SM
cross-section ratio at 100 PeV. If no 5% solution exists, that model/mass is
reported closed for this benchmark and is not propagated.

v0.41a: workflow-trigger checkpoint; physics unchanged.
"""
from __future__ import annotations
import json, math, os
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

from nutau_prop.case6 import Case6Config, Case6DIS, probability_from_dsdy
from nutau_prop.constants import M_N, GEV2_TO_CM2
from nutau_prop.chirality_transport_mixed_expm import TauHelicityFlipTransportMixedExpm

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'v0.41'
OUT.mkdir(parents=True,exist_ok=True)

PDF='CT18NLO'
GQ=0.03
MASSES=(50.0,100.0,220.0)
MODELS={'eq4p4_derivative':1.0,'eq4p3_dipole':2.0}
TARGET_E=1.0e8
TARGET_RATIO=0.05
ENERGIES=np.geomspace(1e4,2e10,101)
ZENITHS=(0.1,0.3,0.5,0.7,1.0)
GAMMAS=(2.0,2.5)
PROBE_PEV=(1.,3.,10.,30.,100.,300.)

class ResolvedPKRDIS(Case6DIS):
    def _d2_resolved(self,E,x,y,M,kappa,factor):
        s=2.0*M_N*E; q2=x*y*s
        if q2 < self.cfg.qmin_gev**2: return 0.0
        p=self._pdf_components(x,q2); sum_xf=sum(p.values())
        c4=kappa/(M*M)
        d2=(4.0/math.pi)*c4*c4*x*s*s*y*(1.0-y)*sum_xf
        d2 *= (M*M/(M*M+q2))**2
        return factor*d2*GEV2_TO_CM2
    def dsdy_resolved(self,E,y,M,kappa,gq,factor):
        lam_vert=gq/kappa
        s=2.0*M_N*E
        xmin=max(self.cfg.qmin_gev**2/(s*y),1e-9)
        xmax=min(1.0-1e-10,lam_vert*lam_vert/s)
        if xmin>=xmax: return 0.0
        return self._integrate_logx(lambda x:self._d2_resolved(E,x,y,M,kappa,factor),xmin,xmax)

CALC=ResolvedPKRDIS(Case6Config(pdf_set=PDF,nx=72,y_log_bins=28,y_linear_bins=28))

def sigma_parts(E,M,kappa,gq,factor,antinu=False):
    _,ys,dy=CALC.y_grid(); snp=0.; ssm=0.
    for y,d in zip(ys,dy):
        cc,nc=CALC.dsdy_sm(E,y,antinu=antinu)
        npv=CALC.dsdy_resolved(E,y,M,kappa,gq,factor)
        ssm+=(cc+nc)*d; snp+=npv*d
    return snp,ssm

def ratio(E,M,kappa,gq,factor):
    a,b=sigma_parts(E,M,kappa,gq,factor,False)
    return a/max(b,1e-300)

def solve_kappa5(M,gq,factor):
    ks=np.geomspace(1e-7,5e-2,180)
    rs=np.array([ratio(TARGET_E,M,k,gq,factor) for k in ks])
    idx=np.where((rs[:-1]<TARGET_RATIO)&(rs[1:]>=TARGET_RATIO))[0]
    if len(idx)==0:
        j=int(np.argmax(rs)); return None,float(ks[j]),float(rs[j])
    lo,hi=float(ks[idx[0]]),float(ks[idx[0]+1])
    for _ in range(28):
        mid=math.sqrt(lo*hi)
        if ratio(TARGET_E,M,mid,gq,factor)>=TARGET_RATIO: hi=mid
        else: lo=mid
    k=hi
    return float(k),float(k),float(ratio(TARGET_E,M,k,gq,factor))

def build_table(M,kappa,gq,factor,antinu=False):
    _,ys,dy=CALC.y_grid(); rec=[]
    for E in ENERGIES:
        cc=np.zeros_like(ys); nc=np.zeros_like(ys); npv=np.zeros_like(ys)
        for i,y in enumerate(ys):
            cc[i],nc[i]=CALC.dsdy_sm(E,y,antinu=antinu)
            npv[i]=CALC.dsdy_resolved(E,y,M,kappa,gq,factor)
        scc,pcc=probability_from_dsdy(cc,dy); snc,pnc=probability_from_dsdy(nc,dy); snp,pnp=probability_from_dsdy(npv,dy)
        rec.append((scc,snc,snp,pcc,pnc,pnp))
    return dict(energies=ENERGIES,y=ys,dy=dy,
                sigma_cc=np.asarray([r[0] for r in rec]),sigma_nc_sm=np.asarray([r[1] for r in rec]),sigma_np=np.asarray([r[2] for r in rec]),
                p_cc=np.asarray([r[3] for r in rec]),p_nc_sm=np.asarray([r[4] for r in rec]),p_np=np.asarray([r[5] for r in rec]))

def interp_log(E,y,e):
    if e<E[0] or e>E[-1]: return float('nan')
    return float(np.interp(np.log(e),np.log(E),y))

def incident(gamma):
    return (ENERGIES/1e6)**(-gamma)

def run_one(model,M,kappa,factor):
    tn=build_table(M,kappa,GQ,factor,False); ta=build_table(M,kappa,GQ,factor,True)
    tr=TauHelicityFlipTransportMixedExpm(tn,ta,tau_refine=8)
    rows=[]
    for gamma in GAMMAS:
        inc=incident(gamma)
        for mu in ZENITHS:
            sm=tr.propagate(mu,inc,inc,mode='dirac_sink',use_np=False,max_step_km=48.0,track_secondary_numu=True)
            bs=tr.propagate(mu,inc,inc,mode='dirac_sink',use_np=True,max_step_km=48.0,track_secondary_numu=True)
            st_sm=sm['nu_tau']+sm['antinu_tau']; st_bs=bs['nu_tau']+bs['antinu_tau']
            smu_sm=sm['nu_mu_secondary']+sm['antinu_mu_secondary']; smu_bs=bs['nu_mu_secondary']+bs['antinu_mu_secondary']
            for ep in PROBE_PEV:
                e=ep*1e6
                a=interp_log(ENERGIES,st_sm,e); b=interp_log(ENERGIES,st_bs,e)
                ma=interp_log(ENERGIES,smu_sm,e); mb=interp_log(ENERGIES,smu_bs,e)
                rows.append(dict(model=model,M_GeV=M,gq=GQ,kappa_GeV_m1=kappa,vertex_scale_GeV=GQ/kappa,gamma=gamma,mu=mu,E_PeV=ep,
                                 tau_flux_ratio=b/max(a,1e-300),tau_frac=(b-a)/max(a,1e-300),
                                 secondary_numu_ratio=mb/max(ma,1e-300) if ma>0 else float('nan'),
                                 secondary_numu_frac=(mb-ma)/max(ma,1e-300) if ma>0 else float('nan')))
    return rows

def main():
    summary=[]; allrows=[]
    for model,factor in MODELS.items():
        for M in MASSES:
            k5,kpeak,rpeak=solve_kappa5(M,GQ,factor)
            rec=dict(model=model,M_GeV=M,gq=GQ,kappa5_GeV_m1=k5,peak_kappa_GeV_m1=kpeak,peak_ratio_100PeV=rpeak,status='OPEN_5PCT' if k5 else 'CLOSED_AT_GQ_0p03')
            summary.append(rec); print(rec,flush=True)
            if k5 is not None:
                allrows.extend(run_one(model,M,k5,factor))
    import pandas as pd
    pd.DataFrame(summary).to_csv(OUT/'benchmark_gate_v0_41.csv',index=False)
    pd.DataFrame(allrows).to_csv(OUT/'earth_propagation_v0_41.csv',index=False)
    if allrows:
        df=pd.DataFrame(allrows); out=[]
        for (model,M,gam),g in df.groupby(['model','M_GeV','gamma']):
            j=g['tau_frac'].abs().idxmax(); r=g.loc[j]
            jm=g['secondary_numu_frac'].abs().idxmax(); rm=g.loc[jm]
            out.append(dict(model=model,M_GeV=M,gamma=gam,max_abs_tau_frac=float(abs(r.tau_frac)),tau_at_mu=float(r.mu),tau_at_E_PeV=float(r.E_PeV),
                            max_abs_secondary_numu_frac=float(abs(rm.secondary_numu_frac)),numu_at_mu=float(rm.mu),numu_at_E_PeV=float(rm.E_PeV)))
        pd.DataFrame(out).to_csv(OUT/'earth_propagation_summary_v0_41.csv',index=False)
    (OUT/'metadata_v0_41.json').write_text(json.dumps(dict(pdf=PDF,gq=GQ,masses=MASSES,models=MODELS,target_ratio=TARGET_RATIO,target_E_GeV=TARGET_E,
        validity='sqrt(shat)<gq/kappa inside PDF fold',transport_mode='dirac_sink',tau_refine=8,max_step_km=48.0,zeniths=ZENITHS,gammas=GAMMAS),indent=2))
    print('WROTE',OUT,flush=True)
if __name__=='__main__': main()
