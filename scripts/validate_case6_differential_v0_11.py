#!/usr/bin/env python3
"""Validation gate for the v0.11 Case-6 dσ/dy tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
TABLES=ROOT/"results"/"v0.11"/"tables"
OUT=ROOT/"results"/"v0.11"
LAMBDAS=(1.5,2.0,2.44,3.0)

def tag(x): return str(x).replace(".","p")

def load(lam,species):
    z=np.load(TABLES/f"case6_L{tag(lam)}_{species}.npz")
    return {k:z[k] for k in z.files}

def loginterp(x,xp,fp):
    return float(10**np.interp(np.log10(x),np.log10(xp),np.log10(np.maximum(fp,1e-300))))

def main():
    old=pd.read_csv(ROOT/"results"/"v0.8"/"case6_pdf_cut_v0_8.csv")
    rows=[]; pmax=0.0; ccmax=0.0
    ref_cc=None
    for lam in LAMBDAS:
        n=load(lam,"nu"); a=load(lam,"antinu")
        for t in (n,a):
            pmax=max(pmax,float(np.max(np.abs(t["p_cc"].sum(axis=1)-1.0))))
            pmax=max(pmax,float(np.max(np.abs(t["p_nc_sm"].sum(axis=1)-1.0))))
            active=t["sigma_np"]>0
            if np.any(active): pmax=max(pmax,float(np.max(np.abs(t["p_np"][active].sum(axis=1)-1.0))))
            pmax=max(pmax,float(np.max(np.abs(t["p_nc"].sum(axis=1)-1.0))))
        if ref_cc is None: ref_cc=(n["sigma_cc"].copy(),a["sigma_cc"].copy())
        else:
            ccmax=max(ccmax,float(np.max(np.abs(n["sigma_cc"]/ref_cc[0]-1.0))))
            ccmax=max(ccmax,float(np.max(np.abs(a["sigma_cc"]/ref_cc[1]-1.0))))
        E=n["energies"]; smn=n["sigma_cc"]+n["sigma_nc_sm"]; sma=a["sigma_cc"]+a["sigma_nc_sm"]
        r=(n["sigma_np"]+a["sigma_np"])/np.maximum(smn+sma,1e-300)
        q=old[np.isclose(old.Lambda_TeV,lam)].sort_values("E_PeV")
        for epev in (0.1,1.0,10.0,100.0,220.0,300.0):
            if epev*1e6<E.min() or epev*1e6>E.max(): continue
            rr=loginterp(epev*1e6,E,r)
            ro=loginterp(epev,q.E_PeV.to_numpy(float),q.r_shatcut.to_numpy(float))
            rows.append(dict(Lambda_TeV=lam,E_PeV=epev,r_v011=rr,r_v08=ro,
                             fractional_difference=rr/ro-1.0))
    R=pd.DataFrame(rows); R.to_csv(OUT/"case6_v011_vs_v08_regression.csv",index=False)
    maxdiff=float(np.max(np.abs(R.fractional_difference)))
    result={"PASS":bool(pmax<1e-10 and ccmax<1e-12 and maxdiff<0.10),
            "max_probability_normalization_error":pmax,
            "max_CC_lambda_dependence":ccmax,
            "max_abs_ratio_difference_vs_v08":maxdiff,
            "ratio_tolerance":0.10,
            "note":"v0.8 used the same CT10 benchmark but coarser y quadrature; comparison uses nu+antinu averaged SM denominator."}
    (OUT/"case6_differential_gate_v0_11.json").write_text(json.dumps(result,indent=2))
    print(R.to_string(index=False,float_format=lambda x:f"{x:.5g}"))
    print(json.dumps(result,indent=2))
    if not result["PASS"]: raise SystemExit("v0.11 differential Case-6 gate FAILED")

if __name__=="__main__": main()
