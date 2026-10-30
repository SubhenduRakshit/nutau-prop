#!/usr/bin/env python3
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import pkr_absorption_horizon_v0_42 as mod

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--indir',required=True);ap.add_argument('--outdir',required=True);a=ap.parse_args();src=Path(a.indir);out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
 fs=sorted(src.rglob('task_*.csv'))
 if not fs:raise RuntimeError('no task csv')
 dfs=[pd.read_csv(f) for f in fs]
 df=pd.concat(dfs,ignore_index=True);df.to_csv(out/'coupled_true_energy_map.csv',index=False)
 # Every matrix task used the same 15-point np.geomspace grid, but CSV round-trip
 # can give ulp-level differences across jobs.  Use one task as the canonical
 # energy grid rather than unioning float values from all tasks.
 canonical=np.sort(dfs[0].E_GeV.to_numpy(float))
 if any(len(x.E_GeV)!=len(canonical) for x in dfs):raise RuntimeError('inconsistent task energy-grid lengths')
 mod.MUS=np.array(sorted(df.mu.unique()),float);mod.EGRID=canonical;mod.SIGMA_LOG10E=.30
 det=mod.detector_weight(df);det.to_csv(out/'arca_public_response_map.csv',index=False)
 rows=[]
 for (model,topology),q in det.groupby(['model','topology']):
  h=q[q.in_horizon_band].copy()
  if not len(h):continue
  h['Sproxy']=np.abs(h.N_BSM_proxy-h.N_SM_proxy)/np.sqrt(np.maximum(h.N_SM_proxy,1e-300));b=h.loc[h.Sproxy.idxmax()];w=h.N_SM_proxy.to_numpy(float)
  rows.append(dict(model=model,topology=topology,max_abs_horizon_distortion=float(h.frac_distortion.abs().max()),SM_event_weighted_abs_horizon_distortion=float(np.sum(w*np.abs(h.frac_distortion))/max(np.sum(w),1e-300)),best_mu=float(b.mu),best_Erec_GeV=float(b.Erec_GeV),best_tau_SM_eff=float(b.tau_SM_eff),Sproxy=float(b.Sproxy),n_horizon_cells=int(len(h))))
 sm=pd.DataFrame(rows);sm.to_csv(out/'horizon_summary.csv',index=False)
 gs=sorted(src.rglob('gate_*.csv'))
 if gs:pd.concat([pd.read_csv(f) for f in gs],ignore_index=True).drop_duplicates('model').to_csv(out/'benchmark_gate.csv',index=False)
 (out/'metadata.json').write_text(json.dumps(dict(version='v0.42h',role='UHE-only coupled horizon discovery; E_min=1 PeV',energy_grid_points=len(mod.EGRID),energy_min_GeV=float(mod.EGRID.min()),energy_max_GeV=float(mod.EGRID.max()),mus=mod.MUS.tolist(),tau_refine=1,max_step_km=192,sigma_log10E=.30,convergence_required=True,lower_boundary_check_required=True),indent=2))
 print(sm.to_string(index=False),flush=True)
if __name__=='__main__':main()
