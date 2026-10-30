#!/usr/bin/env python3
"""v0.42b correction: preserve v0.42, fix detector weighting.

Corrections relative to v0.42:
1) propagated spectra already contain the incident E^-gamma spectrum, so the
   detector fold must NOT multiply by the source spectrum a second time;
2) each sampled zenith point is assigned its midpoint |cos(theta)| bin width,
   so the Poisson proxy carries the correct solid-angle measure.
All transport physics is unchanged.
"""
from __future__ import annotations
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
CORE=HERE/'pkr_absorption_horizon_v0_42.py'
spec=importlib.util.spec_from_file_location('pkr_v042_core',CORE)
mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

mod.OUT=mod.ROOT/'results'/'v0.42b'; mod.OUT.mkdir(parents=True,exist_ok=True)
mod.TABLEDIR=mod.OUT/'tables'; mod.TABLEDIR.mkdir(parents=True,exist_ok=True)


def mu_widths(mus):
    m=np.asarray(sorted(mus),float)
    edges=np.empty(len(m)+1,float); edges[0]=0.0; edges[-1]=1.0
    edges[1:-1]=0.5*(m[:-1]+m[1:])
    return {float(x):float(edges[i+1]-edges[i]) for i,x in enumerate(m)}


def detector_weight_fixed(df):
    cfg=mod.load_config(mod.ROOT/'config'/'arca.yaml')
    sf_sh,_=mod.calibrate_fig7_response(cfg)
    R=mod.fig7_response(mod.EGRID)
    A_tau=np.asarray(R['shower']['up_tcc'])*sf_sh
    dE=np.gradient(mod.EGRID)
    Erec=np.geomspace(1e5,1e8,31)
    K=mod.lognormal_kernel(mod.EGRID,Erec)
    widths=mu_widths(df.mu.unique())
    rows=[]
    for (model,topology,mu),q in df.groupby(['model','topology','mu']):
        q=q.sort_values('E_GeV')
        sm=q.tau_flux_SM.to_numpy(float); bs=q.tau_flux_BSM.to_numpy(float)
        # IMPORTANT: sm and bs already include the incident source spectrum.
        domega=2.0*np.pi*widths[float(mu)]
        base=A_tau*dE*domega
        nsm_true=base*sm; nbs_true=base*bs
        nsm_rec=nsm_true@K; nbs_rec=nbs_true@K
        tau=q.tau_SM.to_numpy(float)
        tau_rec=np.zeros_like(Erec)
        for j in range(len(Erec)):
            w=nsm_true*K[:,j]
            tau_rec[j]=np.sum(w*tau)/max(np.sum(w),1e-300)
        for j,er in enumerate(Erec):
            rows.append(dict(model=model,topology=topology,mu=float(mu),dmu=widths[float(mu)],
                Erec_GeV=float(er),tau_SM_eff=float(tau_rec[j]),N_SM_proxy=float(nsm_rec[j]),
                N_BSM_proxy=float(nbs_rec[j]),frac_distortion=float((nbs_rec[j]-nsm_rec[j])/max(nsm_rec[j],1e-300)),
                in_horizon_band=bool(1.0<=tau_rec[j]<=3.0)))
    return pd.DataFrame(rows)

mod.detector_weight=detector_weight_fixed

if __name__=='__main__':
    mod.main()
