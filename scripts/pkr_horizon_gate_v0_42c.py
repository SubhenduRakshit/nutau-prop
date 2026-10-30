#!/usr/bin/env python3
"""v0.42c fast absorption-horizon gate.

Purpose: answer the immediate physics question with coupled transport only in
zenith chords likely to overlap 1<=tau_SM<=3. This is a SCREEN, not the final
convergence run. Uses the corrected detector bookkeeping from v0.42b.
"""
from __future__ import annotations
import importlib.util, os
from pathlib import Path
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE=Path(__file__).resolve().parent
CORE=HERE/'pkr_absorption_horizon_v0_42.py'
spec=importlib.util.spec_from_file_location('core42',CORE)
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

m.OUT=m.ROOT/'results'/'v0.42c'; m.OUT.mkdir(parents=True,exist_ok=True)
m.TABLEDIR=m.OUT/'tables'; m.TABLEDIR.mkdir(parents=True,exist_ok=True)
m.MUS=np.array([0.25,0.35,0.45,0.55])
m.EGRID=np.geomspace(3e4,3e8,41)
m.TAU_REFINE=4
m.MAX_STEP_KM=96.0

# corrected detector fold: propagated spectra already include source spectrum.
def mu_widths(mus):
    x=np.asarray(sorted(mus),float); e=np.empty(len(x)+1)
    e[0]=max(0.0,x[0]-(x[1]-x[0])/2); e[-1]=min(1.0,x[-1]+(x[-1]-x[-2])/2)
    e[1:-1]=0.5*(x[:-1]+x[1:])
    return {float(v):float(e[i+1]-e[i]) for i,v in enumerate(x)}

def detector_fixed(df):
    cfg=m.load_config(m.ROOT/'config'/'arca.yaml'); sf,_=m.calibrate_fig7_response(cfg)
    R=m.fig7_response(m.EGRID); A=np.asarray(R['shower']['up_tcc'])*sf
    dE=np.gradient(m.EGRID); Er=np.geomspace(1e5,1e8,25); K=m.lognormal_kernel(m.EGRID,Er)
    dw=mu_widths(df.mu.unique()); rows=[]
    for (model,top,mu),q in df.groupby(['model','topology','mu']):
        q=q.sort_values('E_GeV'); sm=q.tau_flux_SM.to_numpy(float); bs=q.tau_flux_BSM.to_numpy(float)
        base=A*dE*(2*np.pi*dw[float(mu)]); ns=(base*sm)@K; nb=(base*bs)@K
        tau=q.tau_SM.to_numpy(float); te=[]
        for j in range(len(Er)):
            w=base*sm*K[:,j]; te.append(np.sum(w*tau)/max(np.sum(w),1e-300))
        for j,e in enumerate(Er):
            rows.append(dict(model=model,topology=top,mu=float(mu),Erec_GeV=float(e),tau_SM_eff=float(te[j]),
                N_SM_proxy=float(ns[j]),N_BSM_proxy=float(nb[j]),frac_distortion=float((nb[j]-ns[j])/max(ns[j],1e-300)),
                in_horizon_band=bool(1<=te[j]<=3)))
    return pd.DataFrame(rows)
m.detector_weight=detector_fixed

def main():
    gates=[]
    for model,f in m.MODELS.items():
        k,kp,rp=m.solve_kappa(f)
        gates.append(dict(model=model,kappa=k,peak_kappa=kp,ratio100=rp,status='OPEN_CONDITIONAL' if k else 'NO_5PCT'))
        if k:
            m.save_table(m.TABLEDIR/f'{model}_nu.npz',m.build_table(k,f,False)); m.save_table(m.TABLEDIR/f'{model}_antinu.npz',m.build_table(k,f,True))
    pd.DataFrame(gates).to_csv(m.OUT/'benchmark_gate.csv',index=False)
    tasks=[(model,float(mu)) for model in m.MODELS for mu in m.MUS]
    rows=[]; workers=min(4,len(tasks))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        fs={ex.submit(m.transport_task,t):t for t in tasks}
        for f in as_completed(fs): rows.extend(f.result()); print('DONE',fs[f],flush=True)
    df=pd.DataFrame(rows); df.to_csv(m.OUT/'coupled_true_energy_map.csv',index=False)
    det=detector_fixed(df); det.to_csv(m.OUT/'arca_public_response_map.csv',index=False)
    out=[]
    for (model,top),q in det.groupby(['model','topology']):
        h=q[q.in_horizon_band].copy()
        if not len(h): continue
        h['Sproxy']=np.abs(h.N_BSM_proxy-h.N_SM_proxy)/np.sqrt(np.maximum(h.N_SM_proxy,1e-300))
        j=h.Sproxy.idxmax(); b=h.loc[j]
        # event-weighted distortion over the horizon strip
        w=h.N_SM_proxy.to_numpy(float)
        out.append(dict(model=model,topology=top,n_horizon_cells=len(h),
            max_abs_distortion=float(h.frac_distortion.abs().max()),
            event_weighted_abs_distortion=float(np.sum(w*np.abs(h.frac_distortion))/max(np.sum(w),1e-300)),
            best_mu=float(b.mu),best_Erec_GeV=float(b.Erec_GeV),best_tau_SM_eff=float(b.tau_SM_eff),Sproxy=float(b.Sproxy)))
    pd.DataFrame(out).to_csv(m.OUT/'horizon_summary.csv',index=False)
    print(pd.DataFrame(out).to_string(index=False),flush=True)

if __name__=='__main__': main()
