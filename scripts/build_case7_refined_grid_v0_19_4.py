#!/usr/bin/env python3
"""Optional v0.19.4 refined-energy-grid builder for Case 7, Lambda=3 TeV.

Builds only the failing diagnostic benchmark on a 181-point energy grid.
This is intentionally separate from the default v0.19.4 run because it
requires fresh PDF integrations.
"""
from __future__ import annotations
import argparse, json, multiprocessing as mp, os
from pathlib import Path
import numpy as np
from nutau_prop.case6 import Case6Config, probability_from_dsdy
from nutau_prop.momentum_operators import MomentumOperatorConfig, MomentumOperatorDIS

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'v0.19.4'/'tables'
_CALC=None


def init_worker(cfg):
    global _CALC
    _CALC=MomentumOperatorDIS(MomentumOperatorConfig(base=Case6Config(**cfg)))


def worker(task):
    E,anti=task
    edges,y,dy=_CALC.y_grid()
    cc=np.zeros_like(y); nc=np.zeros_like(y); np7=np.zeros_like(y)
    for i,yy in enumerate(y):
        cc[i],nc[i]=_CALC.dsdy_sm(E,yy,antinu=anti)
        np7[i]=_CALC.dsdy_np(7,E,yy,3.0,cut='shat')
    return dict(E=float(E),anti=bool(anti),y=y,dy=dy,cc=cc,nc=nc,np7=np7)


def build(records,anti):
    q=sorted((r for r in records if r['anti']==anti),key=lambda r:r['E'])
    E=np.asarray([r['E'] for r in q],float); y=np.asarray(q[0]['y'],float); dy=np.asarray(q[0]['dy'],float)
    scc=[]; snc=[]; snp=[]; pcc=[]; pnc=[]; pnp=[]
    for r in q:
        a,b=probability_from_dsdy(np.asarray(r['cc']),dy)
        c,d=probability_from_dsdy(np.asarray(r['nc']),dy)
        e,f=probability_from_dsdy(np.asarray(r['np7']),dy)
        scc.append(a); pcc.append(b); snc.append(c); pnc.append(d); snp.append(e); pnp.append(f)
    return dict(energies=E,y=y,dy=dy,sigma_cc=np.asarray(scc),sigma_nc_sm=np.asarray(snc),
                sigma_np=np.asarray(snp),p_cc=np.asarray(pcc),p_nc_sm=np.asarray(pnc),p_np=np.asarray(pnp))


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--workers',type=int,default=int(os.environ.get('NUTAU_GRID_NPROC','12')))
    a=ap.parse_args(); OUT.mkdir(parents=True,exist_ok=True)
    targets=[OUT/'case7_L3_refined181_nu.npz',OUT/'case7_L3_refined181_antinu.npz']
    if all(p.exists() for p in targets):
        print('REFINED GRID already present; reusing',flush=True); return
    E=np.geomspace(1e4,2e10,181)
    cfg=Case6Config(nx=96,y_log_bins=36,y_linear_bins=36)
    tasks=[(float(e),anti) for anti in (False,True) for e in E]
    ctx=mp.get_context('spawn')
    with ctx.Pool(processes=max(1,a.workers),initializer=init_worker,initargs=(dict(cfg.__dict__),)) as pool:
        records=list(pool.imap_unordered(worker,tasks,chunksize=1))
    for anti,name in ((False,'nu'),(True,'antinu')):
        t=build(records,anti); p=OUT/f'case7_L3_refined181_{name}.npz'
        np.savez_compressed(p,**t)
        p.with_suffix('.json').write_text(json.dumps(dict(version='v0.19.4',operator=7,Lambda_UV_TeV=3.0,
            cut='shat',pdf=cfg.pdf_set,nenergy=181,emin_GeV=1e4,emax_GeV=2e10,purpose='energy-grid diagnostic only'),indent=2))
        print('WROTE',p,flush=True)

if __name__=='__main__': main()
