#!/usr/bin/env python3
"""Build transport-ready exact-PDF tables for Cases 5 and 7.

This is the v0.15 follow-up to the v0.14 prescreen.  It keeps only the
numerically relevant benchmark points: 3 TeV and the nominal c/(4pi)=1 endpoint
for each operator.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path

import numpy as np

from nutau_prop.case6 import Case6Config, mean_y, probability_from_dsdy
from nutau_prop.momentum_operators import (
    MomentumOperatorConfig,
    MomentumOperatorDIS,
    C5_BENCHMARK,
    C7_BENCHMARK,
    lambda_4pi_tev,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.15"
TABLES = OUT / "tables"
_CALC = None
_CASE5_LAMS = None
_CASE7_LAMS = None
_CUT = None


def _init_worker(cfg_dict, l5, l7, cut):
    global _CALC, _CASE5_LAMS, _CASE7_LAMS, _CUT
    base = Case6Config(**cfg_dict)
    _CALC = MomentumOperatorDIS(MomentumOperatorConfig(base=base))
    _CASE5_LAMS = tuple(float(x) for x in l5)
    _CASE7_LAMS = tuple(float(x) for x in l7)
    _CUT = cut


def _worker(task):
    E, antinu = task
    edges, y, dy = _CALC.y_grid()
    sm_cc = np.zeros_like(y); sm_nc = np.zeros_like(y)
    np5 = {lam: np.zeros_like(y) for lam in _CASE5_LAMS}
    np7 = {lam: np.zeros_like(y) for lam in _CASE7_LAMS}
    for i, yy in enumerate(y):
        sm_cc[i], sm_nc[i] = _CALC.dsdy_sm(E, yy, antinu=antinu)
        for lam in _CASE5_LAMS:
            np5[lam][i] = _CALC.dsdy_np(5, E, yy, lam, cut=_CUT)
        for lam in _CASE7_LAMS:
            np7[lam][i] = _CALC.dsdy_np(7, E, yy, lam, cut=_CUT)
    return dict(E=float(E), antinu=bool(antinu), y=y, dy=dy,
                sm_cc=sm_cc, sm_nc=sm_nc, np5=np5, np7=np7)


def tag(x):
    return (f"{float(x):.6g}").replace(".", "p")


def build_table(records, op, lam, antinu):
    q = sorted((r for r in records if r["antinu"] == antinu), key=lambda r: r["E"])
    E = np.asarray([r["E"] for r in q], float)
    y = np.asarray(q[0]["y"], float); dy = np.asarray(q[0]["dy"], float)
    scc=[]; snc=[]; snp=[]; pcc=[]; pnc=[]; pnp=[]; my=[]
    key = "np5" if int(op)==5 else "np7"
    for r in q:
        cc=np.asarray(r["sm_cc"],float); nc=np.asarray(r["sm_nc"],float)
        npds=np.asarray(r[key][float(lam)],float)
        a,b=probability_from_dsdy(cc,dy); c,d=probability_from_dsdy(nc,dy); e,f=probability_from_dsdy(npds,dy)
        scc.append(a); pcc.append(b); snc.append(c); pnc.append(d); snp.append(e); pnp.append(f)
        my.append(mean_y(npds,y,dy))
    return dict(energies=E,y=y,dy=dy,
                sigma_cc=np.asarray(scc),sigma_nc_sm=np.asarray(snc),sigma_np=np.asarray(snp),
                p_cc=np.asarray(pcc),p_nc_sm=np.asarray(pnc),p_np=np.asarray(pnp),
                mean_y_np=np.asarray(my))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--workers",type=int,default=int(os.environ.get("NUTAU_NPROC","14")))
    ap.add_argument("--nenergy",type=int,default=91)
    ap.add_argument("--emin-gev",type=float,default=1e4)
    ap.add_argument("--emax-gev",type=float,default=2e10)
    ap.add_argument("--cut",choices=["shat","Q"],default="shat")
    a=ap.parse_args()
    l5=[3.0,lambda_4pi_tev(C5_BENCHMARK)]
    l7=[3.0,lambda_4pi_tev(C7_BENCHMARK)]
    OUT.mkdir(parents=True,exist_ok=True); TABLES.mkdir(parents=True,exist_ok=True)
    energies=np.geomspace(a.emin_gev,a.emax_gev,a.nenergy)
    cfg=Case6Config(nx=96,y_log_bins=36,y_linear_bins=36)
    tasks=[(float(E),anti) for anti in (False,True) for E in energies]
    ctx=mp.get_context("spawn")
    with ctx.Pool(processes=max(1,a.workers),initializer=_init_worker,
                  initargs=(dict(cfg.__dict__),l5,l7,a.cut)) as pool:
        records=list(pool.imap_unordered(_worker,tasks,chunksize=1))

    rows=[]
    for op,lams,C in ((5,l5,C5_BENCHMARK),(7,l7,C7_BENCHMARK)):
        for lam in lams:
            for anti in (False,True):
                t=build_table(records,op,lam,anti)
                species="antinu" if anti else "nu"
                p=TABLES/f"case{op}_L{tag(lam)}_{species}.npz"
                np.savez_compressed(p,**t)
                p.with_suffix(".json").write_text(json.dumps({
                    "version":"v0.15","operator":op,"C_GeV^-3":C,
                    "Lambda_UV_TeV":lam,"cut":a.cut,"pdf":cfg.pdf_set,
                    "external_state":"transported separately in v0.15 chirality brackets"
                },indent=2))
                r=t["sigma_np"]/np.maximum(t["sigma_cc"]+t["sigma_nc_sm"],1e-300)
                for epev in (1.,10.,100.,220.,300.):
                    rr=float(10**np.interp(np.log10(epev*1e6),np.log10(t["energies"]),np.log10(np.maximum(r,1e-300))))
                    yy=float(np.interp(np.log10(epev*1e6),np.log10(t["energies"]),t["mean_y_np"]))
                    rows.append(dict(operator=op,Lambda_UV_TeV=lam,species=species,E_PeV=epev,r=rr,mean_y=yy))
    import pandas as pd
    D=pd.DataFrame(rows); D.to_csv(OUT/"momentum_tables_audit_v0_15.csv",index=False)
    summary={"version":"v0.15","case5_lambdas_TeV":l5,"case7_lambdas_TeV":l7,
             "cut":a.cut,"energy_grid_GeV":[a.emin_gev,a.emax_gev,a.nenergy]}
    (OUT/"momentum_table_build_v0_15.json").write_text(json.dumps(summary,indent=2))
    print(D[D.species=="nu"].to_string(index=False,float_format=lambda x:f"{x:.5g}"))
    print("WROTE",TABLES)

if __name__=="__main__": main()
