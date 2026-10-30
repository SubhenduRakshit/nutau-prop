#!/usr/bin/env python3
"""v0.45: explicit PKR Z' width/BR and low-mass dijet consistency gate.

Reads the converged v0.42d 5%-at-100-PeV benchmark and computes:
  * exact tree-level Z'->q qbar vector widths including quark masses;
  * neutrino widths for the two dim-5 PKR vertices in the massless limit;
  * a one-Weyl-channel and operator+h.c. convention bracket;
  * total width/mass and branching fractions;
  * the NWA-equivalent dijet coupling gq*sqrt(BRqq), to compare conservatively
    with the smallest observed CMS EXO-24-007 universal-gq limit (0.03 over
    the published 50--250 GeV range).

Neutrino-width normalization used per massless Weyl pair:
  Eq.(4.4) derivative: Gamma = a^2 M^3/(48 pi)
  Eq.(4.3) dipole:     Gamma = a^2 M^3/(12 pi)
where a=c_i/Lambda=kappa_i/gq.  The h.c.-inclusive bracket doubles this.
The factor-of-two bracket is retained explicitly rather than hidden in a flavor
or charge-conjugation convention.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "results" / "v0.42d" / "benchmark_gate.csv"
OUT = ROOT / "results" / "v0.45"
OUT.mkdir(parents=True, exist_ok=True)

# MSbar-ish/PDG masses are more than adequate here; all except top are tiny at M=50 GeV.
QUARK_MASS_GEV = {
    "u": 0.00216,
    "d": 0.00467,
    "s": 0.0934,
    "c": 1.27,
    "b": 4.18,
    "t": 172.69,
}
NC = 3.0
CMS_EXO_24_007_GLOBAL_MIN_GQ = 0.03  # observed-limit floor quoted by CMS, 50--250 GeV


def gamma_qq_one(M: float, gq: float, mq: float) -> float:
    """Vector coupling gq Z'_mu qbar gamma^mu q."""
    if M <= 2.0 * mq:
        return 0.0
    r = (mq / M) ** 2
    return NC * gq**2 * M / (12.0 * math.pi) * (1.0 + 2.0*r) * math.sqrt(1.0 - 4.0*r)


def gamma_nu_one(model: str, M: float, a: float) -> float:
    if model == "eq4p4_derivative":
        return a*a*M**3/(48.0*math.pi)
    if model == "eq4p3_dipole":
        return a*a*M**3/(12.0*math.pi)
    raise ValueError(model)


def main() -> None:
    if not IN.exists():
        raise RuntimeError(f"Missing input: {IN}")
    b = pd.read_csv(IN)
    rows = []
    qrows = []

    for _, r in b.iterrows():
        model = str(r["model"])
        if str(r["status"]) != "OPEN_CONDITIONAL":
            continue
        M = float(r["M_GeV"])
        gq = float(r["gq"])
        kappa = float(r["kappa"])
        a = kappa / gq
        Lambda_vert = 1.0 / a

        gq_parts = {}
        for q, mq in QUARK_MASS_GEV.items():
            g = gamma_qq_one(M, gq, mq)
            gq_parts[q] = g
            qrows.append(dict(model=model, quark=q, mq_GeV=mq, Gamma_qq_GeV=g, open=bool(g > 0)))
        Gamma_q = float(sum(gq_parts.values()))
        Gamma_nu_1 = gamma_nu_one(model, M, a)

        for conv, mult in [("one_Weyl_pair", 1.0), ("operator_plus_hc", 2.0)]:
            Gamma_nu = mult * Gamma_nu_1
            Gamma_tot = Gamma_q + Gamma_nu
            brq = Gamma_q / Gamma_tot
            brinv = Gamma_nu / Gamma_tot
            gq_equiv_dijet = gq * math.sqrt(brq)
            gq_equiv_inv = gq * math.sqrt(brinv)
            rows.append(dict(
                model=model,
                convention=conv,
                M_GeV=M,
                gq=gq,
                kappa_GeV_inv=kappa,
                c_over_Lambda_GeV_inv=a,
                Lambda_vertex_GeV=Lambda_vert,
                M_over_Lambda_vertex=M/Lambda_vert,
                Gamma_quarks_GeV=Gamma_q,
                Gamma_nu_GeV=Gamma_nu,
                Gamma_total_GeV=Gamma_tot,
                Gamma_over_M=Gamma_tot/M,
                BR_qq=brq,
                BR_invisible=brinv,
                gq_equiv_dijet=gq_equiv_dijet,
                gq_equiv_invisible=gq_equiv_inv,
                cms_global_min_limit_gq=CMS_EXO_24_007_GLOBAL_MIN_GQ,
                dijet_floor_pass=bool(gq_equiv_dijet < CMS_EXO_24_007_GLOBAL_MIN_GQ),
                narrow_width_pass=bool(Gamma_tot/M < 0.10),
                vertex_scale_pass=bool(M < Lambda_vert),
            ))

    R = pd.DataFrame(rows)
    Q = pd.DataFrame(qrows)
    R.to_csv(OUT / "width_br_dijet_gate.csv", index=False)
    Q.to_csv(OUT / "quark_partial_widths.csv", index=False)

    summary = []
    for model, q in R.groupby("model"):
        summary.append(dict(
            model=model,
            max_Gamma_over_M=float(q.Gamma_over_M.max()),
            min_BR_qq=float(q.BR_qq.min()),
            max_BR_invisible=float(q.BR_invisible.max()),
            max_gq_equiv_dijet=float(q.gq_equiv_dijet.max()),
            dijet_floor_pass_all=bool(q.dijet_floor_pass.all()),
            narrow_width_pass_all=bool(q.narrow_width_pass.all()),
            vertex_scale_pass_all=bool(q.vertex_scale_pass.all()),
            verdict=("WIDTH_AND_DIJET_PRESCREEN_PASS" if q.dijet_floor_pass.all() and q.narrow_width_pass.all() and q.vertex_scale_pass.all() else "FAIL_OR_RECAST_REQUIRED"),
        ))
    S = pd.DataFrame(summary)
    S.to_csv(OUT / "summary.csv", index=False)

    meta = {
        "version": "v0.45",
        "input": "results/v0.42d/benchmark_gate.csv",
        "scope": "width/branching + conservative low-mass dijet prescreen; NOT monojet certification",
        "quark_width": "Nc*gq^2*M/(12*pi)*(1+2r)*sqrt(1-4r)",
        "nu_width_per_Weyl_pair": {
            "eq4p4_derivative": "a^2*M^3/(48*pi)",
            "eq4p3_dipole": "a^2*M^3/(12*pi)",
        },
        "normalization_bracket": "one Weyl pair vs operator+h.c. = factor 2 invisible width",
        "cms_anchor": "CMS EXO-24-007 (2026) quotes observed universal-gq limits spanning 0.03--0.13 for 50--250 GeV, assuming BR(Z'->qq)=1. We compare gq*sqrt(BRqq) to the strongest possible 0.03 floor, so passing is conservative within NWA production scaling.",
        "limitations": [
            "Not an operator-specific monojet recast.",
            "No additional invisible flavor multiplicity beyond the explicit one-pair / h.c. bracket.",
            "Dijet comparison uses narrow-width production scaling and the published global minimum limit rather than the exact 50-GeV HEPData point.",
            "Gauge completion/anomaly cancellation/kinetic mixing are not certified here.",
        ],
    }
    (OUT / "metadata.json").write_text(json.dumps(meta, indent=2))

    print("\n=== v0.45 width / BR rows ===")
    print(R[["model","convention","Gamma_over_M","BR_qq","BR_invisible","gq_equiv_dijet","M_over_Lambda_vertex","dijet_floor_pass","narrow_width_pass"]].to_string(index=False))
    print("\n=== v0.45 certification summary ===")
    print(S.to_string(index=False))


if __name__ == "__main__":
    main()
