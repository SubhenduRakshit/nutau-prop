#!/usr/bin/env python3
"""Reuse-first serial v0.42r discovery pass.

Reuses the certified v0.41d conditional 5% benchmarks and avoids multiprocessing
so the dynamically loaded core module is never pickled. This is a discovery
screen only; selected horizon cells require a finer convergence rerun.
"""
from __future__ import annotations
import importlib.util, json
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('pkr_v042_core',HERE/'pkr_absorption_horizon_v0_42.py')
mod=importlib.util.module_from_spec(spec); assert spec.loader is not None; spec.loader.exec_module(mod)

mod.OUT=mod.ROOT/'results'/'v0.42r'; mod.OUT.mkdir(parents=True,exist_ok=True)
mod.TABLEDIR=mod.OUT/'tables'; mod.TABLEDIR.mkdir(parents=True,exist_ok=True)
mod.MUS=np.array([0.12,0.25,0.45,0.70])
mod.EGRID=np.geomspace(1e4,3e8,25)
mod.TAU_REFINE=2
mod.MAX_STEP_KM=192.0

KAPPA={
 'eq4p4_derivative':1.0579695766626163e-4,
 'eq4p3_dipole':5.882689479800913e-5,
}

def main():
 gates=[]; allrows=[]
 for model,factor in mod.MODELS.items():
  k=KAPPA[model]
  gates.append(dict(model=model,status='OPEN_CONDITIONAL',M_GeV=mod.MZ,gq=mod.GQ,
                    kappa=k,vertex_scale_GeV=mod.GQ/k,ratio100=0.05,
                    source='reused_v0.41d_certified_benchmark'))
  mod.save_table(mod.TABLEDIR/f'{model}_nu.npz',mod.build_table(k,factor,False))
  mod.save_table(mod.TABLEDIR/f'{model}_antinu.npz',mod.build_table(k,factor,True))
 pd.DataFrame(gates).to_csv(mod.OUT/'benchmark_gate.csv',index=False)

 for model in mod.MODELS:
  for mu in mod.MUS:
   rows=mod.transport_task((model,float(mu)))
   allrows.extend(rows)
   print('DONE',model,float(mu),flush=True)

 df=pd.DataFrame(allrows)
 df.to_csv(mod.OUT/'coupled_true_energy_map.csv',index=False)
 det=mod.detector_weight(df)
 det.to_csv(mod.OUT/'arca_public_response_map.csv',index=False)

 summary=[]
 for (model,topology),q in det.groupby(['model','topology']):
  h=q[q.in_horizon_band].copy()
  if not len(h): continue
  h['Sproxy']=np.abs(h.N_BSM_proxy-h.N_SM_proxy)/np.sqrt(np.maximum(h.N_SM_proxy,1e-300))
  j=h.Sproxy.idxmax(); b=h.loc[j]
  w=h.N_SM_proxy.to_numpy(float)
  weighted=np.sum(w*np.abs(h.frac_distortion))/max(np.sum(w),1e-300)
  summary.append(dict(model=model,topology=topology,
                      max_abs_horizon_distortion=float(h.frac_distortion.abs().max()),
                      SM_event_weighted_abs_horizon_distortion=float(weighted),
                      best_mu=float(b.mu),best_Erec_GeV=float(b.Erec_GeV),
                      best_tau_SM_eff=float(b.tau_SM_eff),Sproxy=float(b.Sproxy)))
 pd.DataFrame(summary).to_csv(mod.OUT/'horizon_summary.csv',index=False)
 (mod.OUT/'metadata.json').write_text(json.dumps(dict(
   version='v0.42r',role='coarse discovery only',pdf=mod.PDF,M_GeV=mod.MZ,gq=mod.GQ,
   kappas=KAPPA,validity='sqrt(shat)<gq/kappa inside PDF fold',gamma=mod.GAMMA,
   mus=mod.MUS.tolist(),tau_refine=mod.TAU_REFINE,max_step_km=mod.MAX_STEP_KM,
   sigma_log10E=mod.SIGMA_LOG10E,
   topologies=mod.TOPOLOGIES,
   detector='calibrated public ARCA Fig7 tau-CC upgoing response; not official likelihood'),indent=2))
 print(pd.DataFrame(summary).to_string(index=False),flush=True)

if __name__=='__main__': main()
