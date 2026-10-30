#!/usr/bin/env python3
"""One independent coupled-horizon task for the v0.42m matrix run.

Uses the certified v0.41d conditional 5% kappas, a coarse discovery grid, and
runs exactly one (operator, zenith) coupled transport task. No multiprocessing
inside the task, so each GitHub matrix job is isolated and deterministic.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import pkr_absorption_horizon_v0_42 as mod

KAPPA={
    'eq4p4_derivative':1.0579695766626163e-4,
    'eq4p3_dipole':5.882689479800913e-5,
}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--model',required=True,choices=list(KAPPA))
    ap.add_argument('--mu',required=True,type=float)
    ap.add_argument('--outdir',required=True)
    a=ap.parse_args()

    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    mod.OUT=out
    mod.TABLEDIR=out/'tables'; mod.TABLEDIR.mkdir(parents=True,exist_ok=True)
    mod.MUS=np.array([0.12,0.25,0.45,0.70])
    mod.EGRID=np.geomspace(1e4,3e8,25)
    mod.TAU_REFINE=2
    mod.MAX_STEP_KM=192.0

    model=a.model; factor=mod.MODELS[model]; k=KAPPA[model]
    mod.save_table(mod.TABLEDIR/f'{model}_nu.npz',mod.build_table(k,factor,False))
    mod.save_table(mod.TABLEDIR/f'{model}_antinu.npz',mod.build_table(k,factor,True))
    rows=mod.transport_task((model,float(a.mu)))
    tag=f"{model}_mu{int(round(a.mu*100)):02d}"
    pd.DataFrame(rows).to_csv(out/f'task_{tag}.csv',index=False)
    pd.DataFrame([dict(model=model,M_GeV=mod.MZ,gq=mod.GQ,kappa=k,
                       vertex_scale_GeV=mod.GQ/k,mu=float(a.mu),
                       role='coarse_discovery')]).to_csv(out/f'gate_{tag}.csv',index=False)
    print('DONE',model,a.mu,flush=True)

if __name__=='__main__': main()
