#!/usr/bin/env python3
"""v0.14 physics-first prescreen of momentum-growing Cases 5 and 7.

Computes exact CT10nlo PDF folds with the same Qmin and sqrt(shat)<Lambda_UV
prescription as Case 6.  This is deliberately a cross-section/y-kernel gate
before any new Earth-transport implementation.
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
from pathlib import Path

import numpy as np
import pandas as pd

from nutau_prop.case6 import Case6DIS, probability_from_dsdy, mean_y
from nutau_prop.momentum_operators import (
    MomentumOperatorDIS, C5_BENCHMARK, C7_BENCHMARK, lambda_4pi_tev,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.14"
ENERGIES_PEV = (0.1, 1.0, 10.0, 100.0, 220.0, 300.0)
COMMON_LAMBDAS = (1.5, 2.0, 2.44, 3.0, 4.0)


def _integral(dsdy, dy):
    return float(np.sum(np.asarray(dsdy) * np.asarray(dy)))


def worker(epev):
    E = float(epev) * 1e6
    mom = MomentumOperatorDIS()
    c6 = mom.base
    _, y, dy = mom.y_grid()
    smnu = np.zeros_like(y); sma = np.zeros_like(y)
    for i, yy in enumerate(y):
        cc, nc = c6.dsdy_sm(E, yy, antinu=False); smnu[i] = cc + nc
        cc, nc = c6.dsdy_sm(E, yy, antinu=True);  sma[i] = cc + nc
    sigma_sm = 0.5 * (_integral(smnu,dy) + _integral(sma,dy))
    rows=[]
    l5max=lambda_4pi_tev(C5_BENCHMARK); l7max=lambda_4pi_tev(C7_BENCHMARK)
    for op,C,lmax in ((5,C5_BENCHMARK,l5max),(7,C7_BENCHMARK,l7max)):
        lambdas=list(COMMON_LAMBDAS)+[lmax]
        for lam in sorted(set(round(x,6) for x in lambdas)):
            spec=mom.spectra(op,E,lam,cut="shat")
            sig=_integral(spec["np"],spec["dy"])
            my=mean_y(spec["np"],spec["y"],spec["dy"])
            # same-cut Case-6 reference; useful for an exact Case-5 normalization gate
            c6np=np.array([c6.dsdy_np(E,yy,lam,cut="shat") for yy in y])
            sig6=_integral(c6np,dy)
            rows.append(dict(
                operator=op,E_PeV=epev,Lambda_UV_TeV=lam,
                C_GeV_m3=C,c_over_4pi=C*(1e3*lam)**3/(4*math.pi),
                Lambda_4pi_TeV=lmax,sigma_SM_cm2=sigma_sm,
                sigma_NP_cm2=sig,r_NP_over_SM=sig/max(sigma_sm,1e-300),
                mean_y_NP=my,
                ratio_to_Case6_same_cut=sig/max(sig6,1e-300),
            ))
    return rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--workers",type=int,default=min(6,int(os.environ.get("NUTAU_NPROC","6"))))
    a=ap.parse_args(); OUT.mkdir(parents=True,exist_ok=True)
    ctx=mp.get_context("spawn")
    with ctx.Pool(processes=max(1,a.workers)) as pool:
        chunks=list(pool.map(worker,ENERGIES_PEV))
    D=pd.DataFrame([r for c in chunks for r in c]).sort_values(["operator","Lambda_UV_TeV","E_PeV"])
    D.to_csv(OUT/"momentum_operator_prescreen_v0_14.csv",index=False)

    expected=0.5*(C5_BENCHMARK/8.6e-10)**2
    q=D[D.operator==5]
    rel=np.max(np.abs(q.ratio_to_Case6_same_cut/expected-1.0))
    gate={
      "version":"v0.14",
      "case5_expected_ratio_to_case6_same_cut":expected,
      "case5_kernel_regression_max_fractional_error":float(rel),
      "case5_kernel_gate_PASS":bool(rel<1e-8),
      "case5_external_state_caution":"Do not assume an active numu daughter. Fix Dirac/Majorana and charge-conjugation interpretation before transport.",
      "case7_normalization":"L=C7(nu bilinear)GG; alpha_s/(12pi) is absorbed into C7, giving |M|^2=16 C7^2 shat^3 y^3.",
      "cutoff_caution":"sqrt(shat)<Lambda_UV is an EFT-validity boundary, not a physical turnover.",
    }
    (OUT/"momentum_operator_prescreen_gate_v0_14.json").write_text(json.dumps(gate,indent=2))

    key=D[(D.E_PeV.isin([1.0,10.0,100.0,220.0,300.0])) &
          ((np.isclose(D.Lambda_UV_TeV,2.0)) | (np.isclose(D.Lambda_UV_TeV,3.0)) |
           (np.abs(D.Lambda_UV_TeV-D.Lambda_4pi_TeV)<2e-5))]
    print("\nCases 5/7 exact-PDF prescreen")
    print(key[["operator","Lambda_UV_TeV","c_over_4pi","E_PeV","r_NP_over_SM","mean_y_NP","ratio_to_Case6_same_cut"]].to_string(index=False,float_format=lambda x:f"{x:.5g}"))
    print("\nGATE",json.dumps(gate,indent=2))
    print("WROTE",OUT)

if __name__=="__main__": main()
