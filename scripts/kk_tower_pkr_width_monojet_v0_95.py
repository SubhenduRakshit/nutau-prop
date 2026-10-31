#!/usr/bin/env python3
"""v0.95: tower-aware PKR Eq.(4.3) width + monojet/dijet screen.

This supersedes ONLY the PKR numerical section of v0.94.  The v0.94 contact
S/T calculations remain historical/corrected screens.

Why a new normalization is required
------------------------------------
v0.94 inherited the pre-v0.86 resolved-dipole kernel.  The accepted v0.86
normalization for the explicit active->sterile distinct-state realization is

    d sigma_43 = 4 d sigma_44
    Gamma_inv,43(O+h.c.) = a^2 M^3/(6 pi).

The present calculation therefore starts from src/nutau_prop/pkr_resolved.py
and the exact v0.87 CMS dijet + calibrated monojet screens.

Dense KK tower convention
-------------------------
Let a_inc be the *inclusive-width-equivalent* Eq.(4.3) coefficient at the
on-shell mediator mass M, defined by

    Gamma_inv^tower(M) = a_inc^2 M^3/(6 pi).

For a flat delta-dimensional dense tower with m0=0 and approximately
mode-independent brane couplings, the overall KK density and compactification
radius cancel between the on-shell width normalization and UHE production.
The remaining exact continuum factor includes the massive KK phase space.

For one KK state of mass m, r=m^2/shat, the exact Eq.(4.3) hard kernel is
proportional to

    A(r,y) = r^2 - 2r + r y + 2 - 2y,     0 <= y <= 1-r,

while the on-shell vector-decay phase-space factor is

    F_D(r) = (1-r)^2 (1+2r).

After integrating the tower density, the factor multiplying the massless
single-state Eq.(4.3) differential kernel at the same a_inc is

    H_delta(shat,y;M)
      = (sqrt(shat)/M)^delta * (1-y)^(delta/2)
        * [ (1-y)/(delta+4)
            -(2-y)/(delta+2) + 2/delta ]
        / [2 D_delta],

    D_delta = 1/delta - 3/(delta+4) + 2/(delta+6).

This is not the v0.94 arbitrary-Qref ansatz: it is normalized directly to the
same inclusive invisible width that enters the monojet and dijet branching
fractions.

Collider treatment
------------------
For fixed a_inc, the tower width is by definition the accepted v0.86 width,
so the v0.87 NWA equations apply unchanged:
  * exact mass-specific CMS dijet input from results/v0.87_constraints;
  * calibrated monojet r=K_M gq^2 BR_inv;
  * Gamma/M < 0.1, gq^2 < 4pi, and the conservative cutoff sqrt(shat)<1/a_inc.

Because every KK daughter is invisible, an on-shell Z' monojet recoil depends
on the total BR_inv, not on how that invisible width is partitioned among KK
modes.  This makes the NWA width/monojet mapping tower-aware without a new
detector simulation.  It is still a screen, not a combined experimental CL.

Limitations
-----------
* m0=0 continuum tower only; finite spacing/bulk mass are not yet folded.
* R cancels only in the dense continuum.  The script reports the IceCube-2026
  R=0.17 micrometre benchmark as a density diagnostic, not as a universal
  constraint on this specific dipole model.
* The cutoff 1/a_inc is deliberately conservative; a UV completion may map the
  higher-dimensional coupling to a different perturbative scale.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import platform
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _k in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_k] = "1"

import numpy as np
import pandas as pd
from numpy.polynomial.legendre import leggauss

from nutau_prop.case6 import Case6Config, Case6DIS
from nutau_prop.constants import M_N, GEV2_TO_CM2
from nutau_prop.pkr_resolved import (
    allowed_gq_squared,
    invisible_width,
    nucleon_d2_from_xpdf,
    quark_width_coefficient,
)
from pkr_normalization_fold_v0_86 import ROOT, collider_refs

OUT = ROOT / "results" / "v0.95_kk_pkr_width_monojet"
DIJET = ROOT / "results" / "v0.87_constraints" / "exact_dijet_limits.csv"
HBARC_GEV_M = 1.973269804e-16
ICECUBE_R90_M = 0.17e-6

CALCS = {}
SMS = {}


def calc_for(nquad: int) -> Case6DIS:
    if nquad not in CALCS:
        import lhapdf
        lhapdf.setVerbosity(0)
        CALCS[nquad] = Case6DIS(Case6Config(pdf_set="CT18NLO", nx=nquad))
    return CALCS[nquad]


def yrule(ymin: float, n: int):
    nodes, weights = leggauss(n)
    for lo, hi, log in ((ymin, 0.08, True), (0.08, 1.0, False)):
        if lo >= hi:
            continue
        if log:
            aa, bb = math.log(lo), math.log(hi)
            yy = np.exp((nodes + 1.0) * (bb - aa) / 2.0 + aa)
            ww = weights * (bb - aa) / 2.0 * yy
        else:
            yy = (nodes + 1.0) * (hi - lo) / 2.0 + lo
            ww = weights * (hi - lo) / 2.0
        yield from zip(yy, ww)


def sm_pair(E: float, nquad: int):
    key = (E, nquad)
    if key in SMS:
        return SMS[key]
    calc = calc_for(nquad)
    sN = 2.0 * M_N * E
    ymin = max(calc.cfg.qmin_gev**2 / sN, 1e-12)
    values = np.zeros(2)
    for y, wy in yrule(ymin, nquad):
        for anti in (False, True):
            values[int(anti)] += wy * sum(calc.dsdy_sm(E, y, antinu=anti))
    SMS[key] = values
    return values


def decay_density_integral(delta: int) -> float:
    d = float(delta)
    return 1.0 / d - 3.0 / (d + 4.0) + 2.0 / (d + 6.0)


def tower_factor(shat: float, y: float, M: float, delta: int) -> float:
    """Exact m0=0 dense-continuum factor normalized to inclusive width at M."""
    if shat <= 0.0 or M <= 0.0 or delta < 1 or not (0.0 <= y <= 1.0):
        return 0.0
    z = 1.0 - y
    if z <= 0.0:
        return 0.0
    d = float(delta)
    D = decay_density_integral(delta)
    bracket = z / (d + 4.0) - (1.0 + z) / (d + 2.0) + 2.0 / d
    H = (math.sqrt(shat) / M) ** d * z ** (d / 2.0) * bracket / (2.0 * D)
    return max(0.0, H)


def integrate_pair(E: float, M: float, a_inc: float, delta: int, nquad: int):
    """Return massless-single and tower Eq4.3 ratios for gq=1."""
    calc = calc_for(nquad)
    sN = 2.0 * M_N * E
    cutoff = 1.0 / a_inc
    xmax = min(1.0 - 1e-10, cutoff**2 / sN)
    if xmax <= 1e-9:
        return dict(single=0.0, tower=0.0, cutoff_GeV=cutoff)
    qmin2 = calc.cfg.qmin_gev**2
    ymin = max(qmin2 / (sN * xmax), 1e-12)
    if ymin >= 1.0:
        return dict(single=0.0, tower=0.0, cutoff_GeV=cutoff)

    sig_single = 0.0
    sig_tower = 0.0
    for y, wy in yrule(ymin, nquad):
        xmin = max(qmin2 / (sN * y), 1e-9)
        if xmin >= xmax:
            continue

        def base_integrand(x):
            q2 = x * y * sN
            xf = sum(calc.pdf.xfx_proton(pid, x, q2) for pid in range(-5, 6) if pid)
            return nucleon_d2_from_xpdf(
                sN, x, y, xf, M, a_inc, 1.0, operator="4.3"
            ) * GEV2_TO_CM2

        def tower_integrand(x):
            shat = x * sN
            return base_integrand(x) * tower_factor(shat, y, M, delta)

        sig_single += wy * calc._integrate_logx(base_integrand, xmin, xmax)
        sig_tower += wy * calc._integrate_logx(tower_integrand, xmin, xmax)

    sm = sm_pair(E, nquad)
    rs = 0.5 * sig_single * (1.0 / sm[0] + 1.0 / sm[1])
    rt = 0.5 * sig_tower * (1.0 / sm[0] + 1.0 / sm[1])
    return dict(single=rs, tower=rt, cutoff_GeV=cutoff)


def fold_task(task):
    E, M, a_inc, delta, nquad = task
    r = integrate_pair(E, M, a_inc, delta, nquad)
    return dict(
        E_GeV=E,
        E_PeV=E / 1e6,
        M_GeV=M,
        a_inc_GeV_inv=a_inc,
        delta=delta,
        nquad=nquad,
        cutoff_GeV=r["cutoff_GeV"],
        ratio_single_gq1=r["single"],
        ratio_tower_gq1=r["tower"],
        tower_over_single=(r["tower"] / r["single"] if r["single"] > 0 else math.nan),
    )


def collider_allowed(M: float, a_inc: float, rcap: float, dijet: float, K: float):
    """v0.87 exact-NWA screen, with a_inc interpreted as inclusive width coefficient."""
    z = allowed_gq_squared(M, a_inc, "4.3", rcap / K, dijet)
    I = invisible_width(M, a_inc, "4.3")
    B = quark_width_coefficient(M)
    P = z * I / (B * z + I) if z > 0.0 else 0.0
    return dict(
        gq2=z,
        gq=math.sqrt(z) if z > 0 else 0.0,
        Gamma_inv_GeV=I,
        Gamma_total_GeV=B * z + I,
        Gamma_over_M=(B * z + I) / M,
        BR_inv=(I / (B * z + I) if z > 0 else 1.0),
        monojet_r=K * P,
        gq_equiv_dijet=math.sqrt(max(0.0, z - P)),
    )


def radius_diagnostics(M: float, delta: int):
    spacing = HBARC_GEV_M / ICECUBE_R90_M
    n_axis = M / spacing
    return dict(
        icecube_R90_m=ICECUBE_R90_M,
        kk_spacing_at_R90_GeV=spacing,
        kk_spacing_at_R90_eV=spacing * 1e9,
        M_over_spacing_at_R90=n_axis,
        log10_axis_modes_at_R90=math.log10(n_axis),
        log10_hypercube_mode_scale_at_R90=delta * math.log10(n_axis),
        dense_requirement_R_for_MR10_m=10.0 * HBARC_GEV_M / M,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--masses-gev", default="50,100,220")
    ap.add_argument("--energies-pev", default="100,300,1000")
    ap.add_argument("--deltas", default="1,2,3")
    ap.add_argument("--rcaps", default="0.8,1.0,1.2")
    ap.add_argument("--na", type=int, default=48)
    ap.add_argument("--nquad", type=int, default=48)
    ap.add_argument("--nproc", type=int, default=int(os.environ.get("NUTAU_NPROC", "14")))
    args = ap.parse_args()

    masses = tuple(float(x) for x in args.masses_gev.split(",") if x.strip())
    energies = tuple(float(x) * 1e6 for x in args.energies_pev.split(",") if x.strip())
    deltas = tuple(int(x) for x in args.deltas.split(",") if x.strip())
    rcaps = tuple(float(x) for x in args.rcaps.split(",") if x.strip())

    OUT.mkdir(parents=True, exist_ok=True)
    if not DIJET.exists():
        raise RuntimeError(f"Missing exact v0.87 dijet table: {DIJET}")
    dtab = pd.read_csv(DIJET)
    exact_dijet = {
        float(r.M_GeV): float(r.observed_gq_limit) for r in dtab.itertuples()
    }
    refs = collider_refs()

    for M in masses:
        if M not in exact_dijet:
            raise RuntimeError(f"M={M:g} GeV is not an exact v0.87 dijet table mass")
        if M not in refs:
            raise RuntimeError(
                f"M={M:g} GeV lacks a calibrated monojet K_M in v0.79/v0.81; "
                "use only 50,100,220 GeV in this screen"
            )

    tasks = []
    for M in masses:
        aa = np.geomspace(2e-5, 0.999 / M, args.na)
        for a in aa:
            for d in deltas:
                for E in energies:
                    tasks.append((E, M, float(a), d, args.nquad))

    start = time.time()
    workers = max(1, min(args.nproc, len(tasks)))
    print(
        f"v0.95: {len(tasks)} corrected Eq4.3 tower folds on {workers} processes; "
        f"masses={masses}, deltas={deltas}, nquad={args.nquad}",
        flush=True,
    )
    folds = []
    if workers == 1:
        for i, task in enumerate(tasks, 1):
            folds.append(fold_task(task))
            if i % 24 == 0 or i == len(tasks):
                print(f"{i}/{len(tasks)} folds complete", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for i, row in enumerate(pool.map(fold_task, tasks), 1):
                folds.append(row)
                if i % 24 == 0 or i == len(tasks):
                    print(f"{i}/{len(tasks)} folds complete", flush=True)

    F = pd.DataFrame(folds)
    F.to_csv(OUT / "tower_folds.csv", index=False)

    rows = []
    for r in F.itertuples():
        dij = exact_dijet[float(r.M_GeV)]
        K = refs[float(r.M_GeV)]["K"]
        for cap in rcaps:
            c = collider_allowed(float(r.M_GeV), float(r.a_inc_GeV_inv), cap, dij, K)
            rows.append(dict(
                M_GeV=float(r.M_GeV),
                E_PeV=float(r.E_PeV),
                delta=int(r.delta),
                r_cap=cap,
                a_inc_GeV_inv=float(r.a_inc_GeV_inv),
                cutoff_GeV=float(r.cutoff_GeV),
                dijet_gq_observed=dij,
                monojet_K=K,
                gq=c["gq"],
                Gamma_inv_GeV=c["Gamma_inv_GeV"],
                Gamma_over_M=c["Gamma_over_M"],
                BR_inv=c["BR_inv"],
                monojet_r=c["monojet_r"],
                gq_equiv_dijet=c["gq_equiv_dijet"],
                ratio_single=float(r.ratio_single_gq1) * c["gq2"],
                ratio_tower=float(r.ratio_tower_gq1) * c["gq2"],
                tower_over_single=float(r.tower_over_single),
                size_gate_5pct=bool(float(r.ratio_tower_gq1) * c["gq2"] >= 0.05),
            ))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "tower_collider_scan.csv", index=False)

    maxima = []
    for (M, E, d, cap), q in R.groupby(["M_GeV", "E_PeV", "delta", "r_cap"]):
        b = q.loc[q.ratio_tower.idxmax()].to_dict()
        b["status"] = "PASS_SIZE_AFTER_COLLIDER_SCREEN" if b["ratio_tower"] >= 0.05 else "FAIL_SIZE"
        maxima.append(b)
    X = pd.DataFrame(maxima).sort_values(
        ["r_cap", "ratio_tower"], ascending=[True, False]
    )
    X.to_csv(OUT / "maxima.csv", index=False)

    diagnostics = []
    for M in masses:
        for d in deltas:
            z = dict(M_GeV=M, delta=d)
            z.update(radius_diagnostics(M, d))
            diagnostics.append(z)
    pd.DataFrame(diagnostics).to_csv(OUT / "radius_density_diagnostics.csv", index=False)

    nominal = X[np.isclose(X.r_cap, 1.0)].copy()
    lines = [
        "# v0.95 tower-aware PKR Eq. (4.3) width + monojet/dijet screen",
        "",
        "This supersedes the **PKR numerical section only** of v0.94 because v0.94 inherited the pre-v0.86 dipole hard normalization.",
        "The scalar/tensor contact sections of v0.94 are not changed by this calculation.",
        "",
        "## Physics convention",
        "",
        "`a_inc` is defined by the full on-shell KK invisible width:",
        "",
        r"\[\Gamma_{\rm inv}^{\rm tower}(M)=a_{\rm inc}^2 M^3/(6\pi).\]",
        "",
        "Therefore the existing v0.87 monojet/dijet branching-fraction equations apply directly.",
        "The UHE fold uses the exact massive-KK continuum phase space for an m0=0 dense flat tower.",
        "The conservative historical validity convention `sqrt(shat)<1/a_inc` is retained.",
        "",
        "## Nominal collider screen (monojet r <= 1)",
        "",
        "| M GeV | d | E PeV | max tower/SM | single/SM at same point | enhancement | a_inc GeV^-1 | gq | Gamma/M | gate |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in nominal.sort_values(["M_GeV", "delta", "E_PeV"]).itertuples():
        lines.append(
            f"| {r.M_GeV:.0f} | {int(r.delta)} | {r.E_PeV:.0f} | "
            f"{100*r.ratio_tower:.4g}% | {100*r.ratio_single:.4g}% | "
            f"{r.tower_over_single:.3g} | {r.a_inc_GeV_inv:.4g} | "
            f"{r.gq:.4g} | {r.Gamma_over_M:.4g} | "
            f"{'PASS' if r.ratio_tower >= .05 else 'FAIL'} |"
        )

    passed = nominal[nominal.ratio_tower >= 0.05]
    lines += ["", "## Decision", ""]
    if len(passed):
        lines += [
            f"{len(passed)} nominal mass/dimension/energy combinations exceed the 5% size gate after the v0.87 NWA collider screen.",
            "These are the first KK-PKR points worth taking to a finite-spacing/bulk-mass calculation.",
            "Do **not** run Earth propagation yet unless the surviving point also remains viable after the explicit 5D/6D coupling map, oscillation/cosmology constraints, and a finite-mode check.",
        ]
    else:
        lines += [
            "No nominal point exceeds 5% after the corrected v0.86/v0.87 normalization and tower-aware collider screen.",
            "In that case the PKR KK rescue is closed in this dense m0=0 benchmark; do not run Earth propagation.",
        ]
    lines += [
        "",
        "## Radius diagnostic",
        "",
        "The 2026 IceCube LED analysis quotes R <~ 0.17 micrometre (90% CL) in its benchmark model.",
        "At R=0.17 micrometre the KK spacing is ~1.16 eV, so a 50--220 GeV mediator contains an enormous number of modes and the dense-tower approximation is numerically excellent.",
        "This IceCube radius bound is used only as a density benchmark here; it is not automatically a constraint on the present dipole coupling realization.",
        "",
        "## Remaining caveats",
        "",
        "- m0=0; no bulk Dirac mass yet.",
        "- Mode-independent brane couplings in the continuum approximation.",
        "- NWA collider response reused from v0.87; this is not a combined confidence-level fit.",
        "- `1/a_inc` is a conservative validity convention, not a derived higher-dimensional UV cutoff.",
        "",
    ]
    (OUT / "PKR_KK_TOWER_WIDTH_MONOJET_v0.95.md").write_text("\n".join(lines))

    metadata = dict(
        version="v0.95",
        canonical_compute_dir="~/Research/nutau-prop",
        supersedes="PKR numerical section of v0.94 only",
        normalization="v0.86 exact Eq4.3 factor-4 hard kernel and Gamma_inv=a_inc^2 M^3/(6pi)",
        collider="v0.87 exact CMS dijet table + calibrated NWA monojet K_M",
        masses_GeV=masses,
        energies_PeV=[E / 1e6 for E in energies],
        deltas=deltas,
        rcaps=rcaps,
        na=args.na,
        nquad=args.nquad,
        workers=workers,
        seconds=time.time() - start,
        host=platform.node(),
        icecube_2026_radius_benchmark_m=ICECUBE_R90_M,
        assumptions=[
            "dense flat KK continuum",
            "m0=0",
            "approximately mode-independent brane couplings",
            "all open KK daughters invisible for monojet acceptance",
            "sqrt(shat)<1/a_inc validity convention",
        ],
    )
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    print("\n=== v0.95 nominal r<=1 maxima ===")
    cols = [
        "M_GeV", "delta", "E_PeV", "ratio_tower", "ratio_single",
        "tower_over_single", "a_inc_GeV_inv", "gq", "Gamma_over_M",
        "monojet_r", "gq_equiv_dijet", "size_gate_5pct",
    ]
    print(nominal[cols].sort_values("ratio_tower", ascending=False).to_string(index=False))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
