#!/usr/bin/env python3
"""v0.78 — Corrected full-a Eq.4.4 envelope audit.

This supersedes v0.76 only for envelope interpretation.  v0.76 is preserved.

Corrections relative to v0.76
-----------------------------
1. Do NOT keep only the first a crossing of sigma_NP/sigma_SM = 5%.
   For every (M,gq), find every connected interval in a for which
       sigma_NP/sigma_SM >= 5%
   using the already-computed CT18 base envelopes from v0.76.  This naturally
   includes a second crossing when the validity cutoff shrinks the phase space
   at large a.
2. Apply width, vertex-scale and mass-dependent dijet screens pointwise across
   the complete >=5% interval, rather than only at the lower boundary.
3. Do NOT rank different masses by gq^2 BRinv.  Cross-mass monojet ranking is
   explicitly deferred until production x acceptance information is supplied.

No MG5, Pythia, Delphes or CheckMATE run occurs here.
"""
from __future__ import annotations

import json, math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "results" / "v0.76_eq4p4_mass_prescreen" / "base_mass_envelopes.csv"
OUT = ROOT / "results" / "v0.78_eq4p4_full_envelope"
OUT.mkdir(parents=True, exist_ok=True)

TARGET = 0.05
MASSES = (50.0, 100.0, 220.0)
CMS_GQ95 = {50.0: 0.038, 100.0: 0.083, 220.0: 0.116}
CMS_DIGITIZATION_ERR = 0.005
GQ_VALUES = np.arange(0.005, 0.16001, 0.001)
QUARK_MASS = {"u":0.00216,"d":0.00467,"s":0.0934,"c":1.27,"b":4.18,"t":172.69}


def gamma_qq_one(M, gq, mq):
    if M <= 2*mq:
        return 0.0
    r=(mq/M)**2
    return 3*gq*gq*M/(12*math.pi)*(1+2*r)*math.sqrt(1-4*r)


def log_interp_cross(a0,a1,f0,f1,target):
    # Envelope values are positive. Log interpolation is consistent with v0.76.
    x0,x1=math.log(a0),math.log(a1)
    y0,y1=math.log(max(f0,1e-300)),math.log(max(f1,1e-300))
    yt=math.log(target)
    if y1==y0:
        return math.sqrt(a0*a1)
    t=(yt-y0)/(y1-y0)
    return math.exp(x0+t*(x1-x0))


def intervals_above(a, f, target):
    """Return all connected intervals in a with f>=target.

    Endpoints are log-interpolated crossings.  The input grid endpoints are used
    when an interval reaches the trusted scan boundary.
    """
    a=np.asarray(a,float); f=np.asarray(f,float)
    mask=f>=target
    out=[]
    i=0
    while i<len(a):
        if not mask[i]:
            i+=1; continue
        j=i
        while j+1<len(a) and mask[j+1]:
            j+=1
        if i==0:
            lo=float(a[0])
        else:
            lo=log_interp_cross(a[i-1],a[i],f[i-1],f[i],target)
        if j==len(a)-1:
            hi=float(a[-1])
        else:
            hi=log_interp_cross(a[j],a[j+1],f[j],f[j+1],target)
        out.append((lo,hi,i,j))
        i=j+1
    return out


def point_metrics(M,gq,a):
    gq_width=sum(gamma_qq_one(M,gq,mq) for mq in QUARK_MASS.values())
    ginv=a*a*M**3/(24*math.pi)
    gtot=gq_width+ginv
    brinv=ginv/gtot
    brqq=gq_width/gtot
    gqeq=gq*math.sqrt(brqq)
    lim=CMS_GQ95[M]
    return dict(
        Gamma_over_M=gtot/M, BR_inv=brinv, BR_qq=brqq,
        gq_equiv_dijet=gqeq, cms_gq95_digitized=lim,
        dijet_margin=lim-gqeq,
        dijet_central_pass=bool(gqeq<lim),
        dijet_robust_pass=bool(gqeq<lim-CMS_DIGITIZATION_ERR),
        narrow_width_pass=bool(gtot/M<0.10),
        vertex_scale_pass=bool(1.0/a>M),
        monojet_coupling_factor=gq*gq*brinv,
    )


def sample_interval(lo,hi,n=120):
    if hi<=lo:
        return np.asarray([lo])
    return np.geomspace(lo,hi,n)


def main():
    if not IN.exists():
        raise RuntimeError(f"missing v0.76 base envelope: {IN}")
    base=pd.read_csv(IN)
    interval_rows=[]
    valid_rows=[]
    mass_summary=[]

    for M in MASSES:
        q=base[base.M_GeV==M].sort_values("a_GeV_inv")
        if q.empty:
            raise RuntimeError(f"missing mass {M} in {IN}")
        a=q.a_GeV_inv.to_numpy(float)
        f=q.ratio100_at_gq1.to_numpy(float)
        fmax=float(np.max(f))
        gqcrit=math.sqrt(TARGET/fmax)
        mass_interval_count=0
        mass_valid_count=0
        valid_gqs=set()
        best_margin=-math.inf
        best_valid=None

        for gq in GQ_VALUES:
            target_base=TARGET/(gq*gq)
            ints=intervals_above(a,f,target_base)
            for k,(lo,hi,i0,i1) in enumerate(ints):
                mass_interval_count+=1
                pts=sample_interval(lo,hi)
                valid_pts=[]
                for av in pts:
                    m=point_metrics(M,float(gq),float(av))
                    ok=m["narrow_width_pass"] and m["vertex_scale_pass"] and m["dijet_robust_pass"]
                    if ok:
                        valid_pts.append((float(av),m))
                    best_margin=max(best_margin,float(m["dijet_margin"]))

                row=dict(
                    M_GeV=M,gq=float(gq),interval_index=k,
                    a_lo_GeV_inv=lo,a_hi_GeV_inv=hi,
                    Lambda_vertex_lo_GeV=1/hi,Lambda_vertex_hi_GeV=1/lo,
                    grid_i_lo=i0,grid_i_hi=i1,
                    target_ratio=TARGET,
                    n_probe=len(pts),n_valid_probe=len(valid_pts),
                    has_valid_portion=bool(valid_pts),
                    gq_5pct_capacity_threshold=gqcrit,
                )
                if valid_pts:
                    mass_valid_count+=1
                    valid_gqs.add(round(float(gq),12))
                    # Find the valid point with the smallest monojet coupling factor
                    av,m=min(valid_pts,key=lambda z:z[1]["monojet_coupling_factor"])
                    row.update(dict(
                        valid_a_min=min(x[0] for x in valid_pts),
                        valid_a_max=max(x[0] for x in valid_pts),
                        representative_a=av,
                        representative_BR_inv=m["BR_inv"],
                        representative_Gamma_over_M=m["Gamma_over_M"],
                        representative_gq_equiv_dijet=m["gq_equiv_dijet"],
                        representative_dijet_margin=m["dijet_margin"],
                        representative_monojet_coupling_factor=m["monojet_coupling_factor"],
                    ))
                    cand={"M_GeV":M,"gq":float(gq),"a_GeV_inv":av,**m}
                    if best_valid is None or cand["monojet_coupling_factor"]<best_valid["monojet_coupling_factor"]:
                        best_valid=cand
                    for av2,m2 in valid_pts:
                        valid_rows.append(dict(M_GeV=M,gq=float(gq),a_GeV_inv=av2,**m2))
                interval_rows.append(row)

        if best_valid is None:
            status="NO_ROBUST_DIJET_SCREEN_PORTION"
        else:
            status="FULL_A_ROBUST_DIJET_SCREEN_PORTION_EXISTS"
        mass_summary.append(dict(
            M_GeV=M,status=status,
            gq_capacity_threshold=gqcrit,
            n_ge5pct_intervals=mass_interval_count,
            n_intervals_with_valid_portion=mass_valid_count,
            n_gq_with_valid_portion=len(valid_gqs),
            best_dijet_margin_seen=best_margin,
            representative_gq=(best_valid or {}).get("gq",np.nan),
            representative_a=(best_valid or {}).get("a_GeV_inv",np.nan),
            representative_BR_inv=(best_valid or {}).get("BR_inv",np.nan),
            representative_gq_equiv_dijet=(best_valid or {}).get("gq_equiv_dijet",np.nan),
            representative_monojet_coupling_factor=(best_valid or {}).get("monojet_coupling_factor",np.nan),
        ))

    I=pd.DataFrame(interval_rows)
    V=pd.DataFrame(valid_rows)
    S=pd.DataFrame(mass_summary)
    I.to_csv(OUT/"full_intervals.csv",index=False)
    V.to_csv(OUT/"valid_probe_points.csv",index=False)
    S.to_csv(OUT/"mass_summary.csv",index=False)

    meta={
      "version":"v0.78_eq4p4_full_envelope",
      "source":"v0.76 base CT18 envelopes; no new UHE integration",
      "target":"sigma_NP/sigma_SM >= 5% at 100 PeV",
      "correction":"all connected a intervals retained, including second crossings",
      "dijet":"CMS EXO-24-007 digitized screening values from v0.76 with +/-0.005 robustness band",
      "cross_mass_monojet_ranking":"NOT PERFORMED: requires sigma_prod(M)*A_b(M), not gq^2 BRinv alone",
      "next_step":"obtain mass-dependent production x acceptance proxy before choosing any new CheckMATE benchmark",
    }
    (OUT/"summary.json").write_text(json.dumps(meta,indent=2)+"\n")

    print("=== v0.78 Eq4.4 full-a envelope audit ===")
    print(S.to_string(index=False))
    print("cross-mass CheckMATE candidate: NONE (production x acceptance not yet supplied)")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
