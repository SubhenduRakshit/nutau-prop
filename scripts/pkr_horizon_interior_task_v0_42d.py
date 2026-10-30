#!/usr/bin/env python3
"""Boundary-safe interior absorption-horizon convergence task.

Certification point: |cos(theta)|=0.35, reconstructed energy ~=5.01 PeV.
With sigma_log10E=0.30, the two-sigma lower migration scale is ~=1.26 PeV,
so the 1 PeV transport cutoff should not control the reconstructed result.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import pkr_absorption_horizon_v0_42 as mod
from nutau_prop.arca import fig7_response, calibrate_fig7_response
from nutau_prop.config import load_config

ROOT=Path(__file__).resolve().parents[1]
MU=0.35
TARGET_RECO=10**6.7
KAPPA={'eq4p4_derivative':1.0579695766626163e-4,
       'eq4p3_dipole':5.882689479800913e-5}
SETTINGS={'medium':(25,2,96.0),'fine':(31,4,48.0)}

def interp_logx(x,y,x0):
    return float(np.interp(np.log(x0),np.log(np.asarray(x,float)),np.asarray(y,float)))

def fold(q):
    q=q.sort_values('E_GeV'); E=q.E_GeV.to_numpy(float)
    cfg=load_config(ROOT/'config'/'arca.yaml'); sf,_=calibrate_fig7_response(cfg)
    A=np.asarray(fig7_response(E)['shower']['up_tcc'])*sf
    dE=np.gradient(E); Erec=np.geomspace(1e5,1e8,31)
    K=mod.lognormal_kernel(E,Erec,.30)
    sm=q.tau_flux_SM.to_numpy(float); bs=q.tau_flux_BSM.to_numpy(float)
    nsm=(A*dE*sm)@K; nbs=(A*dE*bs)@K
    tau=q.tau_SM.to_numpy(float); te=np.zeros_like(Erec)
    for j in range(len(Erec)):
        w=A*dE*sm*K[:,j]; te[j]=np.sum(w*tau)/max(np.sum(w),1e-300)
    return Erec,te,nsm,nbs,(nbs-nsm)/np.maximum(nsm,1e-300)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--model',required=True,choices=list(KAPPA)); ap.add_argument('--setting',required=True,choices=list(SETTINGS)); ap.add_argument('--outdir',required=True); a=ap.parse_args()
    nE,tr,step=SETTINGS[a.setting]
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    mod.OUT=out; mod.TABLEDIR=out/'tables'; mod.TABLEDIR.mkdir(parents=True,exist_ok=True)
    mod.EGRID=np.geomspace(1e6,3e8,nE); mod.TAU_REFINE=tr; mod.MAX_STEP_KM=step; mod.MUS=np.array([MU])
    model=a.model; factor=mod.MODELS[model]; k=KAPPA[model]
    mod.save_table(mod.TABLEDIR/f'{model}_nu.npz',mod.build_table(k,factor,False))
    mod.save_table(mod.TABLEDIR/f'{model}_antinu.npz',mod.build_table(k,factor,True))
    df=pd.DataFrame(mod.transport_task((model,MU)))
    rows=[]
    for topology,q in df.groupby('topology'):
        q=q.sort_values('E_GeV'); E=q.E_GeV.to_numpy(float)
        er,te,nsm,nbs,frac=fold(q); j=int(np.argmin(np.abs(np.log(er/TARGET_RECO))))
        hor=(te>=1)&(te<=3)
        rows.append(dict(setting=a.setting,model=model,topology=topology,mu=MU,nE=nE,tau_refine=tr,max_step_km=step,
            target_Erec_GeV=float(er[j]),target_tau_SM_eff=float(te[j]),target_reco_frac=float(frac[j]),
            true_tau_frac_3PeV=interp_logx(E,q.tau_frac.to_numpy(float),3e6),
            true_tau_frac_5PeV=interp_logx(E,q.tau_frac.to_numpy(float),5e6),
            sec_numu_frac_3PeV=interp_logx(E,q.sec_numu_frac.to_numpy(float),3e6),
            sec_numu_frac_5PeV=interp_logx(E,q.sec_numu_frac.to_numpy(float),5e6),
            max_abs_horizon_reco=float(np.max(np.abs(frac[hor]))),
            weighted_abs_horizon_reco=float(np.sum(nsm[hor]*np.abs(frac[hor]))/max(np.sum(nsm[hor]),1e-300))))
    pd.DataFrame(rows).to_csv(out/f'summary_{model}_{a.setting}.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False),flush=True)
if __name__=='__main__': main()
