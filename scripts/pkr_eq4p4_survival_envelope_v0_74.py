#!/usr/bin/env python3
"""v0.74: Eq.4.4 M=50 GeV survival-envelope screen after v0.73.

Purpose
-------
At fixed M_Z'=50 GeV, vary g_q and the dimension-5 vertex a=c/Lambda.
For each g_q find the first (smallest-a) solution giving
sigma_NP/sigma_SM = 5% at 100 PeV, with the same resolved propagator and
partonic validity cut sqrt(shat)<1/a used in v0.42d.

Then evaluate the now-fixed O+h.c. invisible width and predict the CMS
monojet multibin r-value by calibrated NWA scaling
    r ~ g_q^2 BR_inv
from v0.73.  This is a cheap survival screen at the SAME mediator mass;
no new detector simulation is performed.

If a robustly surviving point appears, it should be checked with one fresh
CheckMATE sample before any IceCube likelihood or Earth propagation.
"""
from __future__ import annotations

import json
import math
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from nutau_prop.case6 import Case6Config, Case6DIS
from nutau_prop.constants import GEV2_TO_CM2, M_N

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.74_eq4p4_survival_envelope"
OUT.mkdir(parents=True, exist_ok=True)

PDF = "CT18NLO"
MZ = 50.0
TARGET_E = 1.0e8       # 100 PeV
TARGET_RATIO = 0.05

# v0.73 physical O+h.c. calibration at the same mediator mass.
GQ_REF = 0.03
A_REF = 0.0035630660806766665
KAPPA_REF = 1.068919824203e-4
BRINV_REF = 0.5403511562674538
R_MONOJET_REF = 1.238632

# Scan only the vertex-scale-controlled region Lambda_vertex=1/a > M_Z'.
A_MIN = 2.0e-5
A_MAX = 0.999 / MZ
NA = 168
GQ_VALUES = np.arange(0.005, 0.06001, 0.0005)

QUARK_MASS_GEV = {
    "u": 0.00216, "d": 0.00467, "s": 0.0934,
    "c": 1.27, "b": 4.18, "t": 172.69,
}
NC = 3.0
CMS_DIJET_FLOOR = 0.03


class ResolvedEq4p4DIS(Case6DIS):
    """v0.42d Eq4.4 kernel evaluated at gq=1, kappa=a.

    For arbitrary gq the 100-PeV ratio is exactly gq^2 times this base
    result because kappa=gq*a while the validity boundary depends only on a.
    """

    def _d2_base(self, E: float, x: float, y: float, a: float) -> float:
        s = 2.0 * M_N * E
        q2 = x * y * s
        if q2 < self.cfg.qmin_gev**2:
            return 0.0
        p = self._pdf_components(x, q2)
        c4 = a / (MZ * MZ)  # gq=1 => kappa=a
        val = (4.0 / math.pi) * c4*c4 * x * s*s * y * (1.0-y) * sum(p.values())
        val *= (MZ*MZ / (MZ*MZ + q2))**2
        return val * GEV2_TO_CM2

    def dsnp_base(self, E: float, y: float, a: float) -> float:
        s = 2.0 * M_N * E
        xmin = max(self.cfg.qmin_gev**2 / (s*y), 1e-9)
        xmax = min(1.0 - 1e-10, (1.0/a)**2 / s)
        if xmin >= xmax:
            return 0.0
        return self._integrate_logx(lambda x: self._d2_base(E, x, y, a), xmin, xmax)


def ratio_chunk(a_values):
    calc = ResolvedEq4p4DIS(Case6Config(pdf_set=PDF, nx=52, y_log_bins=22, y_linear_bins=22))
    _, ys, dy = calc.y_grid()
    sm_nu = 0.0
    sm_anu = 0.0
    for y, d in zip(ys, dy):
        cc, nc = calc.dsdy_sm(TARGET_E, y, antinu=False)
        sm_nu += (cc + nc) * d
        cc, nc = calc.dsdy_sm(TARGET_E, y, antinu=True)
        sm_anu += (cc + nc) * d

    out = []
    for a in a_values:
        snp = 0.0
        for y, d in zip(ys, dy):
            snp += calc.dsnp_base(TARGET_E, y, float(a)) * d
        r = 0.5 * (snp/max(sm_nu,1e-300) + snp/max(sm_anu,1e-300))
        out.append((float(a), float(r)))
    return out


def gamma_qq_one(M: float, gq: float, mq: float) -> float:
    if M <= 2.0*mq:
        return 0.0
    r = (mq/M)**2
    return NC*gq*gq*M/(12.0*math.pi) * (1.0+2.0*r) * math.sqrt(1.0-4.0*r)


def first_crossing(a: np.ndarray, f: np.ndarray, target: float):
    """Return log-interpolated first upward crossing and bracket index."""
    if f[0] >= target:
        return float(a[0]), 0
    idx = np.where((f[:-1] < target) & (f[1:] >= target))[0]
    if not len(idx):
        return None, None
    i = int(idx[0])
    x0, x1 = math.log(a[i]), math.log(a[i+1])
    y0, y1 = math.log(max(f[i],1e-300)), math.log(max(f[i+1],1e-300))
    yt = math.log(target)
    t = (yt-y0)/(y1-y0) if y1 != y0 else 0.5
    return float(math.exp(x0 + t*(x1-x0))), i


def main() -> int:
    agrid = np.geomspace(A_MIN, A_MAX, NA)
    nproc = min(max(1, int(os.environ.get("NUTAU_NPROC", "14"))), NA)
    chunks = [x for x in np.array_split(agrid, nproc) if len(x)]
    print(f"[1/2] CT18 100-PeV base envelope: {NA} a-points on {nproc} processes", flush=True)
    pairs = []
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        for q in ex.map(ratio_chunk, chunks):
            pairs.extend(q)
    pairs.sort()
    a = np.asarray([x[0] for x in pairs])
    f = np.asarray([x[1] for x in pairs])
    pd.DataFrame({"a_GeV_inv":a, "ratio100_at_gq1":f}).to_csv(OUT/"base_envelope.csv", index=False)

    fmax = float(np.max(f))
    imax = int(np.argmax(f))
    gq_crit = math.sqrt(TARGET_RATIO/fmax)

    # Internal reproduction check of the v0.42d reference benchmark.
    target_ref = TARGET_RATIO/(GQ_REF*GQ_REF)
    a_ref_reco, _ = first_crossing(a, f, target_ref)
    if a_ref_reco is None:
        raise RuntimeError("v0.74 failed to reproduce the known gq=0.03 5% crossing")
    ref_rel = abs(a_ref_reco-A_REF)/A_REF
    if ref_rel > 0.03:
        raise RuntimeError(f"reference crossing mismatch {ref_rel:.2%}: reco a={a_ref_reco}, v0.42d a={A_REF}")

    rows = []
    print("[2/2] width + calibrated monojet NWA envelope", flush=True)
    for gq in GQ_VALUES:
        target_base = TARGET_RATIO/(gq*gq)
        a_sol, ib = first_crossing(a, f, target_base)
        if a_sol is None:
            rows.append(dict(gq=gq, has_5pct=False, target_ratio=TARGET_RATIO,
                             capacity_ratio=gq*gq*fmax, a_GeV_inv=np.nan, kappa_GeV_inv=np.nan,
                             Lambda_vertex_GeV=np.nan, Gamma_over_M=np.nan, BR_inv=np.nan,
                             gq_equiv_dijet=np.nan, monojet_r_pred=np.nan, monojet_class="NO_5PCT"))
            continue

        kappa = gq*a_sol
        gamma_q = sum(gamma_qq_one(MZ, gq, mq) for mq in QUARK_MASS_GEV.values())
        # Physical Eq4.4 O+h.c. result fixed by v0.73.
        gamma_inv = a_sol*a_sol*MZ**3/(24.0*math.pi)
        gamma_tot = gamma_q + gamma_inv
        brinv = gamma_inv/gamma_tot
        brqq = gamma_q/gamma_tot
        rmono = R_MONOJET_REF * (gq*gq*brinv)/(GQ_REF*GQ_REF*BRINV_REF)
        if rmono < 0.8:
            mclass = "ROBUST_SURVIVAL_SCREEN"
        elif rmono < 1.0:
            mclass = "NOMINAL_SURVIVAL_MARGINAL"
        elif rmono > 1.2:
            mclass = "ROBUST_EXCLUSION_SCREEN"
        else:
            mclass = "MARGINAL"

        gqeq = gq*math.sqrt(brqq)
        rows.append(dict(
            gq=gq, has_5pct=True, target_ratio=TARGET_RATIO,
            capacity_ratio=gq*gq*fmax, a_GeV_inv=a_sol, kappa_GeV_inv=kappa,
            Lambda_vertex_GeV=1.0/a_sol, Gamma_q_GeV=gamma_q,
            Gamma_inv_GeV=gamma_inv, Gamma_total_GeV=gamma_tot,
            Gamma_over_M=gamma_tot/MZ, BR_inv=brinv, BR_qq=brqq,
            gq_equiv_dijet=gqeq, dijet_floor_pass=bool(gqeq<CMS_DIJET_FLOOR),
            vertex_scale_pass=bool(1.0/a_sol > MZ), narrow_width_pass=bool(gamma_tot/MZ < 0.10),
            monojet_r_pred=rmono, monojet_class=mclass,
            crossing_bracket_lo=float(a[ib]) if ib is not None else np.nan,
            crossing_bracket_hi=float(a[min(ib+1,len(a)-1)]) if ib is not None else np.nan,
        ))

    df = pd.DataFrame(rows)
    df.to_csv(OUT/"survival_envelope.csv", index=False)
    good_base = df[(df.has_5pct==True) & (df.dijet_floor_pass==True) &
                   (df.vertex_scale_pass==True) & (df.narrow_width_pass==True)]
    nominal = good_base[good_base.monojet_r_pred < 1.0]
    robust = good_base[good_base.monojet_r_pred < 0.8]

    status = "NO_NOMINAL_SURVIVING_5PCT_WINDOW"
    if len(robust):
        status = "ROBUST_SURVIVING_5PCT_WINDOW_FOUND"
    elif len(nominal):
        status = "MARGINAL_SURVIVING_5PCT_WINDOW_FOUND"

    best = None
    if len(good_base):
        b = good_base.loc[good_base.monojet_r_pred.idxmin()]
        best = {k:(bool(v) if isinstance(v,(np.bool_,bool)) else float(v) if isinstance(v,(np.floating,float)) else str(v))
                for k,v in b.items() if k in ["gq","a_GeV_inv","kappa_GeV_inv","Lambda_vertex_GeV","Gamma_over_M","BR_inv","gq_equiv_dijet","monojet_r_pred","monojet_class"]}

    summary = {
        "version":"v0.74_eq4p4_survival_envelope",
        "M_GeV":MZ,
        "target":"sigma_NP/sigma_SM=5% at 100 PeV",
        "physical_width":"Eq4.4 O+h.c.; Gamma_inv=a^2 M^3/(24 pi)",
        "nwa_monojet_scaling":"same M: r=r_ref*(gq^2 BRinv)/(gq_ref^2 BRinv_ref)",
        "gq_critical_grid_estimate":gq_crit,
        "base_envelope_max":fmax,
        "base_envelope_peak_a_GeV_inv":float(a[imax]),
        "reference_a_reconstructed":a_ref_reco,
        "reference_a_relative_error":ref_rel,
        "nominal_surviving_rows":int(len(nominal)),
        "robust_surviving_rows":int(len(robust)),
        "best_screen_point":best,
        "status":status,
        "limitations":[
            "NWA monojet extrapolation is calibrated at one M=50 GeV CheckMATE point; a surviving point needs one fresh detector validation.",
            "CMS dijet comparison uses the conservative global 0.03 floor, not an exact 50-GeV likelihood.",
            "This is a resolved-EFT screen with sqrt(shat)<1/a; exact partial-wave/gauge-completion constraints remain separate.",
        ],
    }
    (OUT/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")

    print("=== v0.74 Eq4.4 M=50 survival envelope ===")
    print(f"reference a reproduction: {a_ref_reco:.8g} (rel diff {ref_rel:.3%})")
    print(f"5% capacity requires gq >= ~{gq_crit:.5f} within Lambda_vertex>M")
    print(f"nominal surviving rows (r<1 + width/dijet/vertex screen): {len(nominal)}")
    print(f"robust surviving rows  (r<0.8): {len(robust)}")
    if best is not None:
        print("best screen point:")
        for k,v in best.items(): print(f"  {k}: {v}")
    print(f"status: {status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
