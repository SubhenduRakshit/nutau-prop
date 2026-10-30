#!/usr/bin/env python3
"""v0.32 exact dynamic-Q UHE screen for active->sterile d=6 S/T and resolved derivative mediator.

Uses the repository's validated Case6DIS SM/PDF machinery.  The d=6 coefficients are the
SMEFT coefficients G=C/Lambda^2 multiplying the chiral operators O_lNqd^S,T.  The direct
partonic kernels below are in that SMEFT convention; do NOT insert low-energy GNI epsilon
coefficients directly.  In particular the Han-Liao-Liu-Marfatia epsilon_{S,T} convention
carries the overall 2 sqrt(2) G_F normalization, which would otherwise create a factor-16
cross-section mismatch when combined with HighPT Wilson coefficients.

Main d=6 collider anchors used here:
  tau23 scalar: |G_S| <= 0.41 TeV^-2
  tau23 tensor: |G_T| <= 0.25 TeV^-2
These are treated as current collider benchmark bounds for the specified flavor direction,
not universal bounds on all sterile operators.

Resolved derivative mediator anchor:
  kappa <= 2.5e-4 GeV^-1 at M=5 GeV (historical PKR IceCube benchmark).
The fixed-kappa mass scan is a shape diagnostic only away from that mass.
"""
from __future__ import annotations

import argparse, csv, json, math
from pathlib import Path
import numpy as np

from nutau_prop.case6 import Case6Config, Case6DIS
from nutau_prop.constants import M_N, GEV2_TO_CM2

PEV_TO_GEV = 1.0e6

D6_BENCHMARKS = [
    # O_lNqd[tau,2,3] contains the NC component nu_tau s -> N b plus conjugate.
    # Charge-averaged nucleon weight = 1/2(s+sbar+b+bbar).
    dict(name="S_tau23_current", lorentz="S", G_TeV_m2=0.41, pids=(3,-3,5,-5), weight_factor=0.5,
         provenance="current tau+MET/HighPT benchmark for lNqd tau23"),
    dict(name="T_tau23_current", lorentz="T", G_TeV_m2=0.25, pids=(3,-3,5,-5), weight_factor=0.5,
         provenance="current tau+MET/HighPT benchmark for lNqd tensor tau23"),
    # Older first-generation tau+MET sensitivity anchors retained only as comparisons.
    dict(name="S_tau11_older", lorentz="S", G_TeV_m2=0.35, pids=(1,-1), weight_factor=1.0,
         provenance="older tau+MET first-generation comparison; not promoted as current exact bound"),
    dict(name="S_u_tau_older", lorentz="S", G_TeV_m2=0.19, pids=(2,-2), weight_factor=1.0,
         provenance="older tau+MET up-quark comparison; not promoted as current exact bound"),
]


def integrate_y(arr, dy):
    return float(np.sum(np.asarray(arr) * np.asarray(dy)))


def mean_y(arr, y, dy):
    den = integrate_y(arr, dy)
    if den <= 0: return float("nan")
    return float(np.sum(np.asarray(arr)*np.asarray(y)*np.asarray(dy))/den)


def sm_total(calc: Case6DIS, E):
    edges,y,dy=calc.y_grid()
    vals=[]
    for yy in y:
        cc1,nc1=calc.dsdy_sm(E,yy,antinu=False)
        cc2,nc2=calc.dsdy_sm(E,yy,antinu=True)
        vals.append(0.5*(cc1+nc1+cc2+nc2))
    return integrate_y(vals,dy)


def d6_dsdy(calc: Case6DIS, E, yval, bench, lambda_cut_tev=None):
    bounds=calc._x_bounds(E,yval,lambda_tev=lambda_cut_tev,cut="shat" if lambda_cut_tev else "none")
    if bounds is None: return 0.0
    xmin,xmax=bounds
    s=2.0*M_N*E
    G=bench["G_TeV_m2"]*1e-6 # GeV^-2, direct SMEFT chiral-operator coefficient C/Lambda^2
    pids=bench["pids"]
    fac=bench["weight_factor"]
    lor=bench["lorentz"]
    def f(x):
        q2=x*yval*s
        if q2 < calc.cfg.qmin_gev**2: return 0.0
        xf=fac*sum(calc.pdf.xfx_isoscalar(pid,x,q2) for pid in pids)
        # Direct chiral SMEFT normalization.  Relative to the low-energy GNI epsilon
        # kernels these are 1/16 once G=(2 sqrt(2) G_F) C is used consistently.
        if lor=="S":
            d2=(G*G*s/(64.0*math.pi))*yval*yval*xf
        elif lor=="T":
            d2=(G*G*s/(2.0*math.pi))*(2.0-yval)**2*xf
        else:
            raise ValueError(lor)
        return d2*GEV2_TO_CM2
    return calc._integrate_logx(f,xmin,xmax)


def d6_sigma(calc,E,bench,lambda_cut_tev=None):
    edges,y,dy=calc.y_grid()
    vals=[d6_dsdy(calc,E,yy,bench,lambda_cut_tev) for yy in y]
    return integrate_y(vals,dy), mean_y(vals,y,dy)


def resolved_dsdy(calc,E,yval,M_TeV,kappa_GeV_m1,vertex_cut_tev=None):
    bounds=calc._x_bounds(E,yval,lambda_tev=vertex_cut_tev,cut="shat" if vertex_cut_tev else "none")
    if bounds is None:return 0.0
    xmin,xmax=bounds
    s=2.0*M_N*E
    m2=(1e3*M_TeV)**2
    kap=kappa_GeV_m1
    def f(x):
        q2=x*yval*s
        if q2<calc.cfg.qmin_gev**2:return 0.0
        p=calc._pdf_components(x,q2)
        sumxf=sum(p.values())
        d2=(4.0/math.pi)*kap*kap*x*s*s*yval*(1-yval)*sumxf/(m2+q2)**2
        return d2*GEV2_TO_CM2
    return calc._integrate_logx(f,xmin,xmax)


def resolved_sigma(calc,E,M_TeV,kappa,vertex_cut_tev=None):
    edges,y,dy=calc.y_grid()
    vals=[resolved_dsdy(calc,E,yy,M_TeV,kappa,vertex_cut_tev) for yy in y]
    return integrate_y(vals,dy),mean_y(vals,y,dy)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--pdfs",default="CT18NLO,NNPDF31_nnlo_as_0118")
    ap.add_argument("--energies-pev",default="1,10,100,300")
    ap.add_argument("--nx",type=int,default=144)
    ap.add_argument("--outdir",default="results/v0.32")
    a=ap.parse_args()
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    pdfs=[x.strip() for x in a.pdfs.split(',') if x.strip()]
    energies=[float(x) for x in a.energies_pev.split(',')]
    rows=[]
    meta={"pdfs":pdfs,"energies_pev":energies,"d6_benchmarks":D6_BENCHMARKS,
          "d6_normalization":"direct chiral SMEFT C/Lambda^2; S and T kernels are 1/16 of the erroneous first v0.32 GNI-epsilon insertion",
          "resolved_anchor":{"M_TeV":0.005,"kappa_GeV_m1":2.5e-4,
          "note":"historical PKR IceCube benchmark; fixed-kappa mass scan is shape-only away from 5 GeV"}}
    for pdf in pdfs:
        calc=Case6DIS(Case6Config(pdf_set=pdf,nx=a.nx))
        for ep in energies:
            E=ep*PEV_TO_GEV
            sm=sm_total(calc,E)
            # d6 exact folds at full contact, c=1 EFT domain and c=4pi edge.
            for b in D6_BENCHMARKS:
                G=b["G_TeV_m2"]
                lam_c1=1.0/math.sqrt(G)
                lam_4pi=math.sqrt(4.0*math.pi/G)
                for cutname,lam in [("none",None),("c1_shat",lam_c1),("4pi_shat",lam_4pi)]:
                    sig,my=d6_sigma(calc,E,b,lam)
                    rows.append(dict(class_name="d6",pdf=pdf,E_PeV=ep,benchmark=b["name"],lorentz=b["lorentz"],
                        coefficient=b["G_TeV_m2"],coefficient_unit="TeV^-2",cut=cutname,cut_TeV=lam or 0.0,
                        sigma_SM_cm2=sm,sigma_NP_cm2=sig,r_NP_over_SM=sig/sm,mean_y_NP=my,
                        note=b["provenance"]))
            # resolved derivative mediator: historical constrained anchor + fixed-kappa shape scan.
            kappa=2.5e-4
            weak_vertex=1.0/kappa/1e3 # c=1 => Lambda_v=1/kappa =4 TeV
            for mass in [0.005,0.04,0.1,0.22,0.5,1.0,2.0]:
                for cutname,vc in [("none",None),("c1_vertex_shat",weak_vertex)]:
                    sig,my=resolved_sigma(calc,E,mass,kappa,vc)
                    rows.append(dict(class_name="resolved_D",pdf=pdf,E_PeV=ep,
                        benchmark=f"M{mass:g}TeV_k2.5e-4",lorentz="derivative-vector",
                        coefficient=kappa,coefficient_unit="GeV^-1",cut=cutname,cut_TeV=vc or 0.0,
                        sigma_SM_cm2=sm,sigma_NP_cm2=sig,r_NP_over_SM=sig/sm,mean_y_NP=my,
                        note=("M=0.005 TeV is historical PKR IceCube-bound anchor; other masses at fixed kappa are shape diagnostics, not certified allowed points")))
    fields=list(rows[0].keys())
    with open(out/"sterile_d6_resolved_dynamicQ_v0_32.csv","w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    with open(out/"sterile_d6_resolved_dynamicQ_v0_32.json","w") as f:
        json.dump({"meta":meta,"rows":rows},f,indent=2)
    lines=["# v0.32 dynamic-Q numerical summary","", "Ratios use the same-PDF SM baseline from this run.",
           "d=6 S/T numbers use the corrected direct-SMEFT chiral normalization.",""]
    for pdf in pdfs:
        lines += [f"## {pdf}",""]
        for ep in [e for e in energies if e in (100.0,300.0)]:
            lines.append(f"### {ep:g} PeV")
            subset=[r for r in rows if r["pdf"]==pdf and r["E_PeV"]==ep]
            for r in subset:
                if r["class_name"]=="d6" and r["cut"]=="c1_shat":
                    lines.append(f"- {r['benchmark']} c=1 cut: {100*r['r_NP_over_SM']:.4g}% (mean y={r['mean_y_NP']:.3g})")
            anchor=[r for r in subset if r["class_name"]=="resolved_D" and r["benchmark"].startswith("M0.005TeV") and r["cut"]=="c1_vertex_shat"]
            if anchor:
                r=anchor[0];lines.append(f"- resolved D, M=5 GeV, kappa=2.5e-4, c=1 vertex cut: {100*r['r_NP_over_SM']:.4g}% [historical bound, needs current tau-specific recast]")
            rr=[r for r in subset if r["class_name"]=="resolved_D" and r["cut"]=="c1_vertex_shat"]
            if rr:
                r=max(rr,key=lambda z:z["r_NP_over_SM"]);lines.append(f"- resolved fixed-kappa shape maximum: {r['benchmark']}: {100*r['r_NP_over_SM']:.4g}% [not mass-bound-certified]")
            lines.append("")
    (out/"sterile_d6_resolved_summary_v0_32.md").write_text("\n".join(lines)+"\n")
    print("wrote",out)
    print("rows",len(rows))

if __name__=="__main__":
    main()
