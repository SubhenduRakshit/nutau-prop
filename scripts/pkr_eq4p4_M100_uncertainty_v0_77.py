#!/usr/bin/env python3
"""v0.77 companion uncertainty audit for the M=100 GeV Eq4.4 benchmark.

No generation or likelihood rerun. Reads the existing CheckMATE total_results
and reports selected-bin ds/s, an inferred MC-only component after removing the
configured 10% cross-section uncertainty in quadrature, and a diagnostic MC
scale band on the multibin r-value.  This is deliberately coarse; 0.8<r<1.2
remains boundary-unresolved.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"v0.77_eq4p4_M100_checkmate"
RES=OUT/"M100_checkmate.json"
CM=ROOT/"cache"/"v0.77_eq4p4_M100_checkmate"/"checkmate"
ANALYSIS="cms_2107_13201"
XSEC_REL=0.10


def find_total(root):
    hits=list(root.rglob("total_results.txt"))
    if not hits: raise RuntimeError(f"no total_results.txt under {root}")
    return sorted(hits,key=lambda p:p.stat().st_mtime)[-1]


def parse(path):
    out=[]
    for line in path.read_text(errors="replace").splitlines()[1:]:
        p=line.split()
        if len(p)<7 or p[0]!=ANALYSIS: continue
        try: out.append((p[1],float(p[5]),float(p[6])))
        except ValueError: pass
    if not out: raise RuntimeError(f"no {ANALYSIS} rows in {path}")
    return out


def wq(v,w,q):
    v=np.asarray(v,float); w=np.asarray(w,float); o=np.argsort(v); v=v[o]; w=w[o]
    c=np.cumsum(w)/max(np.sum(w),1e-300)
    return float(np.interp(q,c,v))


def main():
    rj=json.loads(RES.read_text()); r=rj["checkmate"]["multibin_r"]
    if r is None: raise RuntimeError("multibin_r unavailable")
    r=float(r); p=find_total(CM); rows=parse(p)
    drel=[]; mrel=[]; weights=[]; mcabs=[]; sig=[]
    for _,s,ds in rows:
        if s<=0: continue
        drel.append(ds/s); weights.append(s)
        m=math.sqrt(max(ds*ds-(XSEC_REL*s)**2,0.0))
        mrel.append(m/s); mcabs.append(m); sig.append(s)
    med=wq(drel,weights,0.5); p84=wq(drel,weights,0.84)
    mmed=wq(mrel,weights,0.5); mp84=wq(mrel,weights,0.84)
    meff=math.sqrt(sum(x*x for x in mcabs))/max(sum(sig),1e-300)
    band=[r*(1-meff),r*(1+meff)]
    if r>=1.2: status="ROBUST_EXCLUSION_THRESHOLD"
    elif r<=0.8: status="ROBUST_SURVIVAL_THRESHOLD"
    else: status="MARGINAL_BOUNDARY_UNRESOLVED_VALIDATION_ONLY"
    result={
      "version":"v0.77_M100_uncertainty_companion","source_total_results":str(p),
      "multibin_r_nominal":r,"decision_status":status,
      "configured_xsec_relative_uncertainty":XSEC_REL,
      "selected_bin_total_ds_over_s_signal_weighted_median":med,
      "selected_bin_total_ds_over_s_signal_weighted_p84":p84,
      "inferred_mc_only_ds_over_s_signal_weighted_median":mmed,
      "inferred_mc_only_ds_over_s_signal_weighted_p84":mp84,
      "inferred_independent_bin_mc_scale_on_summed_signal":meff,
      "approx_r_mc_scale_band":band,
      "caveats":[
        "MC-only subtraction assumes CheckMATE combines the configured 10% XS error and MC error in quadrature.",
        "The r band is an overall-scale diagnostic, not a nuisance-profile likelihood interval.",
        "Shower/matching/detector and surrogate-operator modelling systematics are not reduced by more events.",
        "No cross-mass shape comparison is used because changing mediator mass physically changes acceptance."
      ]}
    (OUT/"uncertainty_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print("=== v0.77 M100 uncertainty audit ===")
    print(f"multibin r: {r:.4f}  status: {status}")
    print(f"selected-bin ds/s weighted median: {med:.2%}  p84: {p84:.2%}")
    print(f"inferred MC-only ds/s weighted median: {mmed:.2%}  p84: {mp84:.2%}")
    print(f"effective independent-bin MC scale: {meff:.2%}")
    print(f"approx r MC-scale band: [{band[0]:.3f}, {band[1]:.3f}]")
    print("interpretation: 0.8<r<1.2 remains boundary-unresolved; do not chase sub-percent shifts")
    return 0

if __name__=="__main__": raise SystemExit(main())
