#!/usr/bin/env python3
"""Aggregate independent v0.42m coupled tasks and apply public ARCA weighting."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import pkr_absorption_horizon_v0_42 as mod

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--indir',required=True); ap.add_argument('--outdir',required=True); a=ap.parse_args()
    src=Path(a.indir); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    files=sorted(src.rglob('task_*.csv'))
    if not files: raise RuntimeError('no task csv files found')
    df=pd.concat([pd.read_csv(f) for f in files],ignore_index=True)
    df.to_csv(out/'coupled_true_energy_map.csv',index=False)

    mod.MUS=np.array(sorted(df.mu.unique()),float)
    mod.EGRID=np.array(sorted(df.E_GeV.unique()),float)
    mod.SIGMA_LOG10E=.30
    det=mod.detector_weight(df)
    det.to_csv(out/'arca_public_response_map.csv',index=False)

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
            best_tau_SM_eff=float(b.tau_SM_eff),Sproxy=float(b.Sproxy),
            n_horizon_cells=int(len(h))))
    sm=pd.DataFrame(summary); sm.to_csv(out/'horizon_summary.csv',index=False)

    gates=sorted(src.rglob('gate_*.csv'))
    if gates:
        g=pd.concat([pd.read_csv(f) for f in gates],ignore_index=True)
        g.drop_duplicates(subset=['model']).to_csv(out/'benchmark_gate.csv',index=False)
    (out/'metadata.json').write_text(json.dumps(dict(version='v0.42m',role='coarse parallel discovery only',
        energy_grid_points=int(len(mod.EGRID)),mus=mod.MUS.tolist(),tau_refine=2,max_step_km=192,
        sigma_log10E=.30,detector='calibrated public ARCA Fig7 tau-CC upgoing energy response; not official likelihood',
        convergence_required=True),indent=2))
    print(sm.to_string(index=False),flush=True)

if __name__=='__main__': main()
