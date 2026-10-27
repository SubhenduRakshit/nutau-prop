#!/usr/bin/env python3
"""Build exact-PDF transport tables for the resolved derivative-Z' Case-6 parent.

Primary v0.16 question: can the physical propagator soften UHE growth while
retaining an interesting 10--300 PeV interaction?  We therefore match the
low-Q coefficient to the Case-6 benchmark and scan mediator mass.  No hard
sqrt(shat) contact cutoff is used in the primary scan.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path

import numpy as np
import pandas as pd

from nutau_prop.case6 import Case6Config, mean_y, probability_from_dsdy
from nutau_prop.resolved_mediator import ResolvedMediatorConfig, ResolvedMediatorDIS

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.16"
TABLES = OUT / "tables"
_CALC = None
_MASSES = None
_SCALE = None
_VERTEX_CUT = None


def tag(x):
    return (f"{float(x):.6g}").replace(".", "p")


def _init_worker(cfg_dict, masses, scale, vertex_cut):
    global _CALC, _MASSES, _SCALE, _VERTEX_CUT
    base = Case6Config(**cfg_dict)
    _CALC = ResolvedMediatorDIS(ResolvedMediatorConfig(base=base))
    _MASSES = tuple(float(x) for x in masses)
    _SCALE = float(scale)
    _VERTEX_CUT = vertex_cut


def _worker(task):
    E, antinu = task
    edges, y, dy = _CALC.y_grid()
    sm_cc = np.zeros_like(y); sm_nc = np.zeros_like(y)
    np_by_m = {m: np.zeros_like(y) for m in _MASSES}
    for i, yy in enumerate(y):
        sm_cc[i], sm_nc[i] = _CALC.dsdy_sm(E, yy, antinu=antinu)
        for m in _MASSES:
            np_by_m[m][i] = _CALC.dsdy_np(
                E, yy, m, lowq_scale=_SCALE, vertex_cut_tev=_VERTEX_CUT
            )
    return dict(E=float(E), antinu=bool(antinu), y=y, dy=dy,
                sm_cc=sm_cc, sm_nc=sm_nc, np=np_by_m)


def build_table(records, mass, antinu):
    q = sorted((r for r in records if r["antinu"] == antinu), key=lambda r: r["E"])
    E = np.asarray([r["E"] for r in q], float)
    y = np.asarray(q[0]["y"], float); dy = np.asarray(q[0]["dy"], float)
    scc=[]; snc=[]; snp=[]; pcc=[]; pnc=[]; pnp=[]; my=[]
    for r in q:
        cc=np.asarray(r["sm_cc"],float); nc=np.asarray(r["sm_nc"],float)
        npds=np.asarray(r["np"][float(mass)],float)
        a,b=probability_from_dsdy(cc,dy)
        c,d=probability_from_dsdy(nc,dy)
        e,f=probability_from_dsdy(npds,dy)
        scc.append(a); pcc.append(b); snc.append(c); pnc.append(d); snp.append(e); pnp.append(f)
        my.append(mean_y(npds,y,dy))
    return dict(energies=E,y=y,dy=dy,
                sigma_cc=np.asarray(scc),sigma_nc_sm=np.asarray(snc),sigma_np=np.asarray(snp),
                p_cc=np.asarray(pcc),p_nc_sm=np.asarray(pnc),p_np=np.asarray(pnp),
                mean_y_np=np.asarray(my))


def interp_log(target, xp, fp):
    return float(10**np.interp(np.log10(target), np.log10(xp),
                               np.log10(np.maximum(fp,1e-300))))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--workers",type=int,default=int(os.environ.get("NUTAU_NPROC","14")))
    ap.add_argument("--nenergy",type=int,default=101)
    ap.add_argument("--emin-gev",type=float,default=1e4)
    ap.add_argument("--emax-gev",type=float,default=2e10)
    ap.add_argument("--masses-tev",nargs="+",type=float,
                    default=[0.1,0.2,0.5,1.0,2.0,3.0])
    ap.add_argument("--lowq-scale",type=float,default=1.0,
                    help="multiply the Case-6 low-Q C4 benchmark")
    ap.add_argument("--vertex-cut-tev",type=float,default=None,
                    help="optional diagnostic sqrt(shat) cut for dim-5 nu-Z' vertex")
    a=ap.parse_args()

    OUT.mkdir(parents=True,exist_ok=True); TABLES.mkdir(parents=True,exist_ok=True)
    energies=np.geomspace(a.emin_gev,a.emax_gev,a.nenergy)
    cfg=Case6Config(nx=96,y_log_bins=36,y_linear_bins=36)
    tasks=[(float(E),anti) for anti in (False,True) for E in energies]
    ctx=mp.get_context("spawn")
    with ctx.Pool(processes=max(1,a.workers),initializer=_init_worker,
                  initargs=(dict(cfg.__dict__),a.masses_tev,a.lowq_scale,a.vertex_cut_tev)) as pool:
        records=list(pool.imap_unordered(_worker,tasks,chunksize=1))

    calc=ResolvedMediatorDIS(ResolvedMediatorConfig(base=cfg))
    rows=[]
    probe_pev=(1.,10.,100.,220.,300.,1000.,10000.)
    for mass in a.masses_tev:
        kappa=calc.kappa_from_lowq_c4(calc.cfg.c4_lowq_gev_m3*a.lowq_scale,mass)
        for anti in (False,True):
            table=build_table(records,mass,anti)
            species="antinu" if anti else "nu"
            p=TABLES/f"resolved_M{tag(mass)}_S{tag(a.lowq_scale)}_{species}.npz"
            np.savez_compressed(p,**table)
            p.with_suffix(".json").write_text(json.dumps({
                "version":"v0.16",
                "model":"derivative-Zprime resolved parent of Case 6",
                "M_Zprime_TeV":mass,
                "lowQ_C4_GeV^-3":calc.cfg.c4_lowq_gev_m3*a.lowq_scale,
                "kappa_c2gq_over_Lambda_GeV^-1":kappa,
                "vertex_cut_TeV":a.vertex_cut_tev,
                "pdf":cfg.pdf_set,
                "note":"Primary scan isolates physical propagator softening; no hard contact cutoff unless vertex_cut_TeV is explicitly set."
            },indent=2))
            r=table["sigma_np"]/np.maximum(table["sigma_cc"]+table["sigma_nc_sm"],1e-300)
            for epev in probe_pev:
                e=epev*1e6
                rr=interp_log(e,table["energies"],r)
                yy=float(np.interp(np.log10(e),np.log10(table["energies"]),table["mean_y_np"]))
                rows.append(dict(M_Zprime_TeV=mass,lowq_scale=a.lowq_scale,species=species,
                                 kappa_GeV_m1=kappa,E_PeV=epev,r_NP_over_SM=rr,mean_y_NP=yy))

    D=pd.DataFrame(rows)
    D.to_csv(OUT/"resolved_mediator_crosssection_v0_16.csv",index=False)
    nu=D[D.species=="nu"].copy()
    piv=nu.pivot(index=["M_Zprime_TeV","kappa_GeV_m1"],columns="E_PeV",values="r_NP_over_SM").reset_index()
    if 220.0 in piv and 10000.0 in piv:
        piv["growth_10EeV_over_220PeV"]=piv[10000.0]/np.maximum(piv[220.0],1e-300)
    piv.to_csv(OUT/"resolved_mediator_shape_gate_v0_16.csv",index=False)
    print("\nResolved mediator exact-PDF shape gate (nu):")
    print(nu.to_string(index=False,float_format=lambda x:f"{x:.5g}"))
    print("\nCompact turnover diagnostic:")
    print(piv.to_string(index=False,float_format=lambda x:f"{x:.5g}"))
    status={
        "version":"v0.16",
        "masses_TeV":a.masses_tev,
        "lowq_scale":a.lowq_scale,
        "vertex_cut_TeV":a.vertex_cut_tev,
        "primary_physics":"Case-6 low-Q matching with exact M^4/(M^2+Q^2)^2 propagator factor",
        "caution":"kappa required for low-Q matching is reported separately; this scan is not yet a current collider-limit recast."
    }
    (OUT/"resolved_mediator_build_status_v0_16.json").write_text(json.dumps(status,indent=2))
    print("WROTE",OUT)

if __name__=="__main__": main()
