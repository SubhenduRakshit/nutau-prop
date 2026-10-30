#!/usr/bin/env python3
"""Focused convergence check for the v0.42 absorption-horizon result.

Only the dominant horizon chord |cos(theta)|=0.70 is evolved.  The lower
energy boundary stays at 1 PeV because the tested observable is the UHE
1--300 PeV band; the boundary itself is checked separately by comparing the
1.0 and 1.5 PeV true-energy cells and reconstructed-energy folds.

This script reuses the certified v0.41d conditional 5% kappas and compares
three numerical resolutions.  It is a numerical-convergence check, not an
allowed-model likelihood.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import pandas as pd
import pkr_absorption_horizon_v0_42 as mod
from nutau_prop.arca import fig7_response, calibrate_fig7_response
from nutau_prop.config import load_config

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'v0.42c'; OUT.mkdir(parents=True,exist_ok=True)
MU=0.70
KAPPA={'eq4p4_derivative':1.0579695766626163e-4,
       'eq4p3_dipole':5.882689479800913e-5}
SETTINGS={
 'coarse':dict(nE=15,tau_refine=1,max_step_km=192.0),
 'medium':dict(nE=25,tau_refine=2,max_step_km=96.0),
 'fine':dict(nE=31,tau_refine=4,max_step_km=48.0),
}

def interp_logx(x,y,x0):
 x=np.asarray(x,float); y=np.asarray(y,float)
 return float(np.interp(np.log(x0),np.log(x),y))

def fold_group(q):
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
 frac=(nbs-nsm)/np.maximum(nsm,1e-300)
 return Erec,te,nsm,nbs,frac

def run_setting(label,s):
 mod.EGRID=np.geomspace(1e6,3e8,int(s['nE']))
 mod.TAU_REFINE=int(s['tau_refine']); mod.MAX_STEP_KM=float(s['max_step_km'])
 mod.MUS=np.array([MU])
 mod.OUT=OUT/label; mod.OUT.mkdir(parents=True,exist_ok=True)
 mod.TABLEDIR=mod.OUT/'tables'; mod.TABLEDIR.mkdir(parents=True,exist_ok=True)
 rows=[]
 for model,factor in mod.MODELS.items():
  k=KAPPA[model]
  mod.save_table(mod.TABLEDIR/f'{model}_nu.npz',mod.build_table(k,factor,False))
  mod.save_table(mod.TABLEDIR/f'{model}_antinu.npz',mod.build_table(k,factor,True))
  rows.extend(mod.transport_task((model,MU)))
 df=pd.DataFrame(rows); df.to_csv(mod.OUT/'true_energy.csv',index=False)
 out=[]
 for (model,topology),q in df.groupby(['model','topology']):
  q=q.sort_values('E_GeV')
  E=q.E_GeV.to_numpy(float); tf=q.tau_frac.to_numpy(float); sf=q.sec_numu_frac.to_numpy(float)
  Erec,tau_eff,nsm,nbs,frac=fold_group(q)
  horizon=(tau_eff>=1)&(tau_eff<=3)
  weighted=np.sum(nsm[horizon]*np.abs(frac[horizon]))/max(np.sum(nsm[horizon]),1e-300)
  target_rec=10**6.1  # 1.258925 PeV, discovery-grid optimum
  j=int(np.argmin(np.abs(np.log(Erec/target_rec))))
  out.append(dict(setting=label,model=model,topology=topology,nE=int(s['nE']),
    tau_refine=int(s['tau_refine']),max_step_km=float(s['max_step_km']),
    tau_frac_1PeV=interp_logx(E,tf,1e6),tau_frac_1p5PeV=interp_logx(E,tf,1.5e6),
    sec_numu_frac_1PeV=interp_logx(E,sf,1e6),sec_numu_frac_1p5PeV=interp_logx(E,sf,1.5e6),
    reco_frac_1p259PeV=float(frac[j]),reco_tau_eff_1p259PeV=float(tau_eff[j]),
    max_abs_horizon_reco=float(np.max(np.abs(frac[horizon]))),
    SM_event_weighted_abs_horizon_reco=float(weighted)))
 return pd.DataFrame(out)

def main():
 allres=[]
 for label,s in SETTINGS.items():
  print('RUN',label,flush=True); allres.append(run_setting(label,s))
 r=pd.concat(allres,ignore_index=True); r.to_csv(OUT/'convergence.csv',index=False)
 fine=r[r.setting=='fine'].set_index(['model','topology'])
 med=r[r.setting=='medium'].set_index(['model','topology'])
 comps=[]
 for idx in fine.index:
  f=fine.loc[idx]; m=med.loc[idx]
  comps.append(dict(model=idx[0],topology=idx[1],
    rel_change_reco_1p259=abs(f.reco_frac_1p259PeV-m.reco_frac_1p259PeV)/max(abs(f.reco_frac_1p259PeV),1e-12),
    rel_change_weighted_horizon=abs(f.SM_event_weighted_abs_horizon_reco-m.SM_event_weighted_abs_horizon_reco)/max(abs(f.SM_event_weighted_abs_horizon_reco),1e-12),
    rel_change_true_1PeV=abs(f.tau_frac_1PeV-m.tau_frac_1PeV)/max(abs(f.tau_frac_1PeV),1e-12),
    rel_change_true_1p5PeV=abs(f.tau_frac_1p5PeV-m.tau_frac_1p5PeV)/max(abs(f.tau_frac_1p5PeV),1e-12)))
 c=pd.DataFrame(comps); c.to_csv(OUT/'medium_vs_fine.csv',index=False)
 (OUT/'metadata.json').write_text(json.dumps(dict(version='v0.42c',mu=MU,
   energy_min_GeV=1e6,energy_max_GeV=3e8,settings=SETTINGS,
   convergence_target='medium-to-fine relative change <=5% preferred',
   benchmark='v0.41d conditional 5% kappas'),indent=2))
 print('\nCONVERGENCE RESULTS'); print(r.to_string(index=False)); print('\nMEDIUM vs FINE'); print(c.to_string(index=False))
if __name__=='__main__': main()
