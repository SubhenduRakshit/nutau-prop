#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"v0.13"


def need(name):
    p=OUT/name
    if not p.exists(): raise SystemExit(f"missing {p}")
    return p


def main():
    A=pd.read_csv(need("case6_arca_joint_v0_13.csv"))
    I=pd.read_csv(need("case6_icecube_joint_v0_13.csv"))
    K=pd.read_csv(need("case6_km3_joint_v0_13.csv"))
    status=json.loads(need("case6_joint_status_v0_13.json").read_text())
    checks={}
    checks["grid_gate"] = bool(status.get("grid_gate_PASS"))
    checks["arca_finite"] = bool(np.isfinite(A.select_dtypes(include=[float,int]).to_numpy()).all())
    checks["icecube_finite"] = bool(np.isfinite(I.select_dtypes(include=[float,int]).to_numpy()).all())
    checks["km3_finite"] = bool(np.isfinite(K.select_dtypes(include=[float,int]).to_numpy()).all())
    checks["nonnegative_counts"] = bool(
        (A[["N_shower_SM","N_shower_BSM","N_track_SM","N_track_BSM"]]>=0).all().all()
        and (I[[c for c in I.columns if c.startswith("IceCube_") and "scale" not in c]]>=0).all().all()
        and (K[[c for c in K.columns if c.startswith("KM3_") and "fractional" not in c and "ratio" not in c]]>=0).all().all()
    )
    # Same published SM baseline must be recovered in both ARCA response modes.
    piv=A.pivot(index="Lambda_TeV",columns="mode",values=["N_shower_SM","N_track_SM"])
    ds=np.max(np.abs(piv[("N_shower_SM","joint")]/piv[("N_shower_SM","propagation_only")]-1.0))
    dt=np.max(np.abs(piv[("N_track_SM","joint")]/piv[("N_track_SM","propagation_only")]-1.0))
    checks["arca_same_sm_baseline"] = bool(ds<1e-10 and dt<1e-10)
    # Public-response envelope should order low <= central/high only weakly; enforce low<=high.
    checks["icecube_envelope_order"] = bool((I.IceCube_joint_low <= I.IceCube_joint_high + 1e-12).all())
    scales=[c for c in I.columns if c.startswith("public_window_flux_scale_")]
    checks["positive_public_scales"] = bool((I[scales]>0).all().all())
    passed=all(checks.values())
    report={"version":"v0.13","checks":checks,"PASS":passed,
            "interpretation":"PASS validates numerical bookkeeping of the public-response joint fold; it is not an IceCube collaboration likelihood."}
    (OUT/"case6_joint_gate_v0_13.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    if not passed: raise SystemExit("v0.13 joint Case-6 QA gate FAILED")

if __name__=="__main__": main()
