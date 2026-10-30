#!/usr/bin/env python3
"""v0.42: coupled PKR absorption-horizon transport + public-response detector weighting.

This is the first post-v0.41 calculation centered on the useful region
1 <= tau_SM <= 3 rather than the vertical chord.

Benchmarks
----------
M_Z' = 50 GeV, g_q = 0.03, CT18NLO.
For each PKR branch kappa is solved self-consistently so sigma_NP/sigma_SM=5%
at 100 PeV with sqrt(shat)<g_q/kappa imposed inside every PDF integral.

Transport topologies
--------------------
1. dirac_sink: explicit sterile-N disappearance completion (maximal sink envelope).
2. majorana_flip: active nu_tau <-> antinu_tau redistribution proxy.  This is
   NOT claimed to be the unique literal-PKR flavour completion; it quantifies
   the difference between disappearance and active redistribution in the same
   hard kernel.

Detector weighting
------------------
Uses the calibrated public ARCA Fig.7 tau-CC upgoing response as an energy
weight.  Zenith information comes from the coupled PREM propagation grid.
A controlled log-normal reconstructed-energy migration test with sigma_log10E
=0.30 is also written.  This is a public-response projection, not an official
ARCA likelihood.
"""
from __future__ import annotations

import json, math, os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from nutau_prop.case6 import Case6Config, Case6DIS, probability_from_dsdy
from nutau_prop.constants import M_N, GEV2_TO_CM2, NA
from nutau_prop.prem import chord_profile
from nutau_prop.chirality_transport_mixed_expm import TauHelicityFlipTransportMixedExpm
from nutau_prop.arca import fig7_response, calibrate_fig7_response
from nutau_prop.config import load_config

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'v0.42'
OUT.mkdir(parents=True, exist_ok=True)
TABLEDIR = OUT / 'tables'
TABLEDIR.mkdir(parents=True, exist_ok=True)

PDF = 'CT18NLO'
MZ = 50.0
GQ = 0.03
TARGET_E = 1.0e8
TARGET_RATIO = 0.05
MODELS = {'eq4p4_derivative': 1.0, 'eq4p3_dipole': 2.0}
TOPOLOGIES = {'sterile_sink': 'dirac_sink', 'active_redistribution_proxy': 'majorana_flip'}
# horizon-focused zenith grid; vertical is retained only as a diagnostic endpoint
MUS = np.array([0.08,0.12,0.18,0.25,0.35,0.45,0.55,0.70,1.00])
EGRID = np.geomspace(1e4, 3e8, 61)
GAMMA = 2.5
SIGMA_LOG10E = 0.30
TAU_REFINE = 8
MAX_STEP_KM = 48.0


class ResolvedPKRDIS(Case6DIS):
    def _d2(self, E, x, y, kappa, factor):
        s = 2.0*M_N*E
        q2 = x*y*s
        if q2 < self.cfg.qmin_gev**2:
            return 0.0
        p = self._pdf_components(x, q2)
        c4 = kappa/(MZ*MZ)
        val = (4.0/math.pi)*c4*c4*x*s*s*y*(1.0-y)*sum(p.values())
        val *= (MZ*MZ/(MZ*MZ+q2))**2
        return factor*val*GEV2_TO_CM2

    def dsnp(self, E, y, kappa, factor):
        s = 2.0*M_N*E
        xmin = max(self.cfg.qmin_gev**2/(s*y), 1e-9)
        xmax = min(1.0-1e-10, (GQ/kappa)**2/s)
        if xmin >= xmax:
            return 0.0
        return self._integrate_logx(lambda x: self._d2(E,x,y,kappa,factor), xmin, xmax)


CALC = ResolvedPKRDIS(Case6Config(pdf_set=PDF, nx=52, y_log_bins=22, y_linear_bins=22))


def integrated_sigmas(E, kappa, factor, antinu=False):
    _, ys, dy = CALC.y_grid()
    snp = ssm = 0.0
    for y,d in zip(ys,dy):
        cc,nc = CALC.dsdy_sm(E,y,antinu=antinu)
        snp += CALC.dsnp(E,y,kappa,factor)*d
        ssm += (cc+nc)*d
    return snp, ssm


def ratio100(kappa, factor):
    vals=[]
    for anti in (False,True):
        a,b=integrated_sigmas(TARGET_E,kappa,factor,anti)
        vals.append(a/max(b,1e-300))
    return 0.5*(vals[0]+vals[1])


def solve_kappa(factor):
    ks=np.geomspace(2e-6,2e-3,90)
    rs=np.array([ratio100(k,factor) for k in ks])
    idx=np.where((rs[:-1]<TARGET_RATIO)&(rs[1:]>=TARGET_RATIO))[0]
    if not len(idx):
        j=int(np.argmax(rs)); return None,float(ks[j]),float(rs[j])
    lo,hi=float(ks[idx[0]]),float(ks[idx[0]+1])
    for _ in range(28):
        mid=math.sqrt(lo*hi)
        if ratio100(mid,factor)>=TARGET_RATIO: hi=mid
        else: lo=mid
    return hi,hi,ratio100(hi,factor)


def build_table(kappa, factor, antinu=False):
    _,ys,dy=CALC.y_grid(); rows=[]
    for E in EGRID:
        cc=np.zeros_like(ys); nc=np.zeros_like(ys); npv=np.zeros_like(ys)
        for i,y in enumerate(ys):
            cc[i],nc[i]=CALC.dsdy_sm(E,y,antinu=antinu)
            npv[i]=CALC.dsnp(E,y,kappa,factor)
        scc,pcc=probability_from_dsdy(cc,dy)
        snc,pnc=probability_from_dsdy(nc,dy)
        snp,pnp=probability_from_dsdy(npv,dy)
        rows.append((scc,snc,snp,pcc,pnc,pnp))
    return dict(energies=EGRID,y=ys,dy=dy,
                sigma_cc=np.asarray([r[0] for r in rows]),
                sigma_nc_sm=np.asarray([r[1] for r in rows]),
                sigma_np=np.asarray([r[2] for r in rows]),
                p_cc=np.asarray([r[3] for r in rows]),
                p_nc_sm=np.asarray([r[4] for r in rows]),
                p_np=np.asarray([r[5] for r in rows]))


def save_table(path, table):
    np.savez_compressed(path, **table)


def load_table(path):
    z=np.load(path); return {k:z[k] for k in z.files}


def incident(E):
    return (E/1e6)**(-GAMMA)


def column_depth(mu):
    p=chord_profile(float(mu),max_step_km=24.0)
    return float(np.sum(np.asarray(p['rho'])*np.asarray(p['ds_km'])*1e5))


def sm_tau(E, mu, tnu, tanu):
    # inclusive nu/antinu-averaged SM optical depth guide
    s=0.25*(tnu['sigma_cc']+tnu['sigma_nc_sm']+tanu['sigma_cc']+tanu['sigma_nc_sm'])*2.0
    X=column_depth(mu)
    return NA*X*s


def interp_log(grid, vals, e):
    grid=np.asarray(grid,float); vals=np.asarray(vals,float)
    if e<grid[0] or e>grid[-1]: return np.nan
    if np.all(vals>0):
        return float(np.exp(np.interp(np.log(e),np.log(grid),np.log(vals))))
    return float(np.interp(np.log(e),np.log(grid),vals))


def transport_task(args):
    model, mu = args
    tnu=load_table(TABLEDIR/f'{model}_nu.npz')
    tanu=load_table(TABLEDIR/f'{model}_antinu.npz')
    tr=TauHelicityFlipTransportMixedExpm(tnu,tanu,tau_refine=TAU_REFINE)
    inc=incident(EGRID)
    sm=tr.propagate(float(mu),inc,inc,mode='dirac_sink',use_np=False,
                    max_step_km=MAX_STEP_KM,track_secondary_numu=True)
    out=[]
    for label,mode in TOPOLOGIES.items():
        bs=tr.propagate(float(mu),inc,inc,mode=mode,use_np=True,
                        max_step_km=MAX_STEP_KM,track_secondary_numu=True)
        tau_sm=sm['nu_tau']+sm['antinu_tau']
        tau_bs=bs['nu_tau']+bs['antinu_tau']
        mu_sm=sm['nu_mu_secondary']+sm['antinu_mu_secondary']
        mu_bs=bs['nu_mu_secondary']+bs['antinu_mu_secondary']
        tguide=sm_tau(EGRID,mu,tnu,tanu)
        for i,E in enumerate(EGRID):
            out.append(dict(model=model,topology=label,mu=float(mu),E_GeV=float(E),
                            tau_SM=float(tguide[i]),
                            tau_flux_SM=float(tau_sm[i]),tau_flux_BSM=float(tau_bs[i]),
                            tau_frac=float((tau_bs[i]-tau_sm[i])/max(tau_sm[i],1e-300)),
                            sec_numu_SM=float(mu_sm[i]),sec_numu_BSM=float(mu_bs[i]),
                            sec_numu_frac=float((mu_bs[i]-mu_sm[i])/max(mu_sm[i],1e-300)) if mu_sm[i]>0 else np.nan))
    return out


def lognormal_kernel(Etrue, Erec, sigma=SIGMA_LOG10E):
    lt=np.log10(Etrue)[:,None]; lr=np.log10(Erec)[None,:]
    K=np.exp(-0.5*((lr-lt)/sigma)**2)
    K/=np.maximum(np.sum(K,axis=1,keepdims=True),1e-300)
    return K


def detector_weight(df):
    cfg=load_config(ROOT/'config'/'arca.yaml')
    sf_sh,_=calibrate_fig7_response(cfg)
    R=fig7_response(EGRID)
    # Public ARCA tau-CC upgoing energy response; no fine zenith response is public here.
    A_tau=np.asarray(R['shower']['up_tcc'])*sf_sh
    phi=incident(EGRID)
    # logarithmic integration weights
    dE=np.gradient(EGRID)
    Erec=np.geomspace(1e5,1e8,31)
    K=lognormal_kernel(EGRID,Erec)
    rows=[]
    for (model,topology,mu),q in df.groupby(['model','topology','mu']):
        q=q.sort_values('E_GeV')
        sm=q.tau_flux_SM.to_numpy(float); bs=q.tau_flux_BSM.to_numpy(float)
        base=phi*A_tau*dE
        nsm_true=base*sm; nbs_true=base*bs
        nsm_rec=nsm_true@K; nbs_rec=nbs_true@K
        # tau_SM assigned to reconstructed bin by SM-event-weighted true-energy average
        tau=q.tau_SM.to_numpy(float)
        tau_rec=np.zeros_like(Erec)
        for j in range(len(Erec)):
            w=nsm_true*K[:,j]
            tau_rec[j]=np.sum(w*tau)/max(np.sum(w),1e-300)
        for j,er in enumerate(Erec):
            frac=(nbs_rec[j]-nsm_rec[j])/max(nsm_rec[j],1e-300)
            rows.append(dict(model=model,topology=topology,mu=float(mu),Erec_GeV=float(er),
                             tau_SM_eff=float(tau_rec[j]),N_SM_proxy=float(nsm_rec[j]),
                             N_BSM_proxy=float(nbs_rec[j]),frac_distortion=float(frac),
                             in_horizon_band=bool(1.0<=tau_rec[j]<=3.0)))
    return pd.DataFrame(rows)


def main():
    gates=[]
    for model,factor in MODELS.items():
        k,kp,rp=solve_kappa(factor)
        if k is None:
            gates.append(dict(model=model,status='NO_5PCT',peak_kappa=kp,peak_ratio=rp)); continue
        gates.append(dict(model=model,status='OPEN_CONDITIONAL',M_GeV=MZ,gq=GQ,kappa=k,
                          vertex_scale_GeV=GQ/k,ratio100=rp))
        save_table(TABLEDIR/f'{model}_nu.npz',build_table(k,factor,False))
        save_table(TABLEDIR/f'{model}_antinu.npz',build_table(k,factor,True))
    pd.DataFrame(gates).to_csv(OUT/'benchmark_gate.csv',index=False)

    tasks=[(m,float(mu)) for m in MODELS for mu in MUS if any(g.get('model')==m and g.get('status')=='OPEN_CONDITIONAL' for g in gates)]
    workers=min(max(1,int(os.environ.get('NUTAU_NPROC','4'))),len(tasks))
    allrows=[]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs={ex.submit(transport_task,t):t for t in tasks}
        for fut in as_completed(futs):
            t=futs[fut]; r=fut.result(); allrows.extend(r); print('DONE',t,flush=True)
    df=pd.DataFrame(allrows)
    df.to_csv(OUT/'coupled_true_energy_map.csv',index=False)

    det=detector_weight(df)
    det.to_csv(OUT/'arca_public_response_map.csv',index=False)

    summary=[]
    for (model,topology),q in det.groupby(['model','topology']):
        h=q[q.in_horizon_band].copy()
        if len(h):
            # Poisson-only proxy: |delta N|/sqrt(N), arbitrary common normalization cancels for ranking
            h['Sproxy']=np.abs(h.N_BSM_proxy-h.N_SM_proxy)/np.sqrt(np.maximum(h.N_SM_proxy,1e-300))
            j=h.Sproxy.idxmax(); b=h.loc[j]
            w=np.abs(h.N_BSM_proxy-h.N_SM_proxy)
            weighted=np.sum(w*np.abs(h.frac_distortion))/max(np.sum(w),1e-300)
            summary.append(dict(model=model,topology=topology,
                                max_abs_horizon_distortion=float(h.frac_distortion.abs().max()),
                                weighted_abs_horizon_distortion=float(weighted),
                                best_mu=float(b.mu),best_Erec_GeV=float(b.Erec_GeV),
                                best_tau_SM_eff=float(b.tau_SM_eff),Sproxy=float(b.Sproxy)))
    pd.DataFrame(summary).to_csv(OUT/'horizon_summary.csv',index=False)

    meta=dict(version='v0.42',pdf=PDF,M_GeV=MZ,gq=GQ,target_ratio=TARGET_RATIO,
              validity='sqrt(shat)<gq/kappa inside PDF fold',gamma=GAMMA,mus=MUS.tolist(),
              tau_refine=TAU_REFINE,max_step_km=MAX_STEP_KM,
              topologies=TOPOLOGIES,sigma_log10E=SIGMA_LOG10E,
              detector='ARCA calibrated Fig7 tau-CC upgoing response; public-response projection, not official likelihood')
    (OUT/'metadata.json').write_text(json.dumps(meta,indent=2))
    print(pd.DataFrame(summary).to_string(index=False),flush=True)

if __name__=='__main__':
    main()
