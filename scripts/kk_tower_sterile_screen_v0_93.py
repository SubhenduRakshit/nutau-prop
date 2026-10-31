#!/usr/bin/env python3
"""
nutau-prop v0.93: collider-normalized KK-tower screen for sterile production.

Purpose
-------
Test whether a bulk sterile-neutrino KK tower can rescue sterile channels that
were too small in the single-state analysis, without making the unphysical
replacement sigma -> N_KK * sigma at fixed coupling.

The per-mode coupling is re-normalized to the existing collider/width bound.
In the continuum regime (Q_ref R >> 1), for delta flat bulk dimensions and
approximately mode-independent brane couplings,

    sigma_KK(shat) / sigma_single(shat)
        = N_eff(sqrt(shat)) / N_eff(Q_ref)
        ~ (sqrt(shat) / Q_ref)^delta,

provided both scales are below the same UV cutoff.  If the published bound is
anchored above the EFT/vertex cutoff we conservatively set Q_ref -> cutoff.
This immediately guarantees no enhancement in a strict c=1 contact domain
when the UHE fold is also truncated at that same cutoff.

Implemented parents
-------------------
1) gauge-invariant d=6 active->sterile scalar/tensor O_lNqd branch, using the
   corrected v0.32 direct-SMEFT kernels and current HighPT benchmark bounds;
2) resolved PKR Eq.(4.3) transition-dipole Z' branch, reusing the validated
   v0.41d hard kernel.  Here Q_ref is scanned because an on-shell-width anchor
   probes ~M_Z', while a hard monojet recast may effectively probe a larger
   invariant scale.

This is a SCREEN, not a constraints certification.  A >5% survivor still
requires an explicit higher-dimensional UV model and a tower-aware collider,
width, cosmology/oscillation and unitarity recast.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from nutau_prop.case6 import Case6Config, Case6DIS
from nutau_prop.constants import M_N, GEV2_TO_CM2

PEV_TO_GEV = 1.0e6
DEFAULT_ENERGIES_PEV = (100.0, 300.0, 1000.0)
DEFAULT_DELTAS = (1, 2, 3)

D6 = {
    "scalar_tau23": dict(
        lorentz="S",
        G_TeV_m2=0.41,
        pids=(3, -3, 5, -5),
        weight_factor=0.5,
        note="HighPT tau23 scalar benchmark; corrected v0.32 normalization",
    ),
    "tensor_tau23": dict(
        lorentz="T",
        G_TeV_m2=0.25,
        pids=(3, -3, 5, -5),
        weight_factor=0.5,
        note="HighPT tau23 tensor benchmark; corrected v0.32 normalization",
    ),
}

PKR = dict(
    M_GeV=50.0,
    gq=0.03,
    kappa_GeV_m1=5.882689479800913e-5,
    dipole_factor=2.0,
    note="v0.41d conditional 5%-scale Eq.(4.3) transition-dipole benchmark",
)


def integrate_y(arr, dy):
    return float(np.sum(np.asarray(arr) * np.asarray(dy)))


def sm_total(calc: Case6DIS, E_GeV: float) -> float:
    _, y, dy = calc.y_grid()
    vals = []
    for yy in y:
        cc1, nc1 = calc.dsdy_sm(E_GeV, yy, antinu=False)
        cc2, nc2 = calc.dsdy_sm(E_GeV, yy, antinu=True)
        vals.append(0.5 * (cc1 + nc1 + cc2 + nc2))
    return integrate_y(vals, dy)


def qref_effective(qref_GeV: float, cutoff_GeV: float) -> float:
    """Constraint normalization scale, clipped to the controlled UV domain."""
    if cutoff_GeV <= 0 or qref_GeV <= 0:
        raise ValueError("qref and cutoff must be positive")
    return min(qref_GeV, cutoff_GeV)


def kk_factor(shat_GeV2: float, delta: int, qref_GeV: float, cutoff_GeV: float) -> float:
    """
    Collider-normalized continuum KK multiplicity.

    For a dense flat KK tower, N_eff(Q) propto Q^delta.  The overall radius,
    mode degeneracy and wavefunction normalization cancel between the UHE
    process and the reference observable.  The UHE integration itself is
    already restricted to sqrt(shat) < cutoff_GeV.
    """
    if shat_GeV2 <= 0:
        return 0.0
    q = math.sqrt(shat_GeV2)
    if q > cutoff_GeV * (1.0 + 1e-12):
        return 0.0
    qref = qref_effective(qref_GeV, cutoff_GeV)
    return (q / qref) ** int(delta)


def d6_dsdy(calc, E_GeV, yval, bench, cutoff_TeV, delta=0, qref_GeV=1.0):
    bounds = calc._x_bounds(E_GeV, yval, lambda_tev=cutoff_TeV, cut="shat")
    if bounds is None:
        return 0.0
    xmin, xmax = bounds
    s = 2.0 * M_N * E_GeV
    G = bench["G_TeV_m2"] * 1e-6
    pids = bench["pids"]
    fac = bench["weight_factor"]
    lor = bench["lorentz"]
    cutoff_GeV = 1.0e3 * cutoff_TeV

    def f(x):
        q2 = x * yval * s
        if q2 < calc.cfg.qmin_gev**2:
            return 0.0
        shat = x * s
        xf = fac * sum(calc.pdf.xfx_isoscalar(pid, x, q2) for pid in pids)
        if lor == "S":
            d2 = (G * G * s / (64.0 * math.pi)) * yval * yval * xf
        elif lor == "T":
            d2 = (G * G * s / (2.0 * math.pi)) * (2.0 - yval) ** 2 * xf
        else:
            raise ValueError(lor)
        mult = 1.0 if delta == 0 else kk_factor(shat, delta, qref_GeV, cutoff_GeV)
        return d2 * mult * GEV2_TO_CM2

    return calc._integrate_logx(f, xmin, xmax)


def d6_sigma(calc, E_GeV, bench, cutoff_TeV, delta=0, qref_GeV=1.0):
    _, y, dy = calc.y_grid()
    vals = [d6_dsdy(calc, E_GeV, yy, bench, cutoff_TeV, delta, qref_GeV) for yy in y]
    return integrate_y(vals, dy)


def pkr_dsdy(calc, E_GeV, yval, delta=0, qref_GeV=50.0):
    M = PKR["M_GeV"]
    gq = PKR["gq"]
    kappa = PKR["kappa_GeV_m1"]
    f_dip = PKR["dipole_factor"]
    cutoff_GeV = gq / kappa
    s = 2.0 * M_N * E_GeV
    xmin = max(calc.cfg.qmin_gev**2 / (s * yval), 1e-9)
    xmax = min(1.0 - 1e-10, cutoff_GeV**2 / s)
    if xmin >= xmax:
        return 0.0

    def f(x):
        q2 = x * yval * s
        if q2 < calc.cfg.qmin_gev**2:
            return 0.0
        shat = x * s
        p = calc._pdf_components(x, q2)
        c = kappa / M**2
        d2 = (
            f_dip
            * (4.0 / math.pi)
            * c
            * c
            * x
            * s
            * s
            * yval
            * (1.0 - yval)
            * sum(p.values())
            * (M * M / (M * M + q2)) ** 2
        )
        mult = 1.0 if delta == 0 else kk_factor(shat, delta, qref_GeV, cutoff_GeV)
        return d2 * mult * GEV2_TO_CM2

    return calc._integrate_logx(f, xmin, xmax)


def pkr_sigma(calc, E_GeV, delta=0, qref_GeV=50.0):
    _, y, dy = calc.y_grid()
    vals = [pkr_dsdy(calc, E_GeV, yy, delta, qref_GeV) for yy in y]
    return integrate_y(vals, dy)


def one_job(job):
    pdf, nx, energy_pev, channel, domain, delta, qref = job
    calc = Case6DIS(Case6Config(pdf_set=pdf, nx=nx))
    E = energy_pev * PEV_TO_GEV
    sm = sm_total(calc, E)

    if channel in D6:
        b = D6[channel]
        G = b["G_TeV_m2"]
        if domain == "c1":
            cutoff_TeV = 1.0 / math.sqrt(G)
        elif domain == "4pi":
            cutoff_TeV = math.sqrt(4.0 * math.pi / G)
        else:
            raise ValueError(domain)
        cutoff_GeV = 1.0e3 * cutoff_TeV
        sig0 = d6_sigma(calc, E, b, cutoff_TeV, delta=0)
        sig = d6_sigma(calc, E, b, cutoff_TeV, delta=delta, qref_GeV=qref)
        qref_eff = qref_effective(qref, cutoff_GeV)
        note = b["note"]
    elif channel == "pkr_eq4p3_dipole":
        cutoff_GeV = PKR["gq"] / PKR["kappa_GeV_m1"]
        sig0 = pkr_sigma(calc, E, delta=0)
        sig = pkr_sigma(calc, E, delta=delta, qref_GeV=qref)
        qref_eff = qref_effective(qref, cutoff_GeV)
        note = PKR["note"]
    else:
        raise ValueError(channel)

    r0 = sig0 / sm if sm > 0 else float("nan")
    r = sig / sm if sm > 0 else float("nan")
    enh = sig / sig0 if sig0 > 0 else float("nan")
    return dict(
        pdf=pdf,
        E_PeV=energy_pev,
        channel=channel,
        domain=domain,
        delta=delta,
        qref_input_GeV=qref,
        qref_effective_GeV=qref_eff,
        cutoff_GeV=cutoff_GeV,
        sigma_SM_cm2=sm,
        sigma_single_cm2=sig0,
        sigma_KK_cm2=sig,
        r_single_over_SM=r0,
        r_KK_over_SM=r,
        enhancement_KK_over_single=enh,
        pass_5pct=bool(r >= 0.05),
        note=note,
    )


def build_jobs(pdf, nx, energies, deltas):
    jobs = []
    for channel in ("scalar_tau23", "tensor_tau23"):
        for domain in ("c1", "4pi"):
            for E in energies:
                for d in deltas:
                    for qref in (1000.0, 2000.0):
                        jobs.append((pdf, nx, E, channel, domain, d, qref))

    for E in energies:
        for d in deltas:
            for qref in (PKR["M_GeV"], 250.0, 500.0):
                jobs.append((pdf, nx, E, "pkr_eq4p3_dipole", "vertex_c1", d, qref))
    return jobs


def write_summary(rows, out):
    lines = [
        "# v0.93 collider-normalized bulk-sterile KK tower screen",
        "",
        "This screen does **not** multiply a fixed single-state coupling by an arbitrary number of KK modes.",
        "The per-mode coupling is normalized to the same collider/width reference observable, so the dense-tower factor is",
        "",
        r"\\[H(\\hat s)=N_{\\rm eff}(\\sqrt{\\hat s})/N_{\\rm eff}(Q_{\\rm ref})\\simeq(\\sqrt{\\hat s}/Q_{\\rm ref})^\\delta,\\]",
        "",
        "with both scales clipped to the same controlled UV domain.",
        "",
        "## Immediate structural result",
        "",
        "For the d=6 c=1 scalar/tensor branches, the collider reference scale is clipped to the same",
        "EFT cutoff used in the UHE fold. Therefore every accepted partonic point has H<=1: a dense KK",
        "tower cannot rescue the strict weak-coupling c=1 conclusion. Any enhancement can only appear",
        "in the looser 4pi interpretation, or in a resolved mediator where Q_ref lies below the vertex cutoff.",
        "",
        "## Numerical maxima",
        "",
        "| channel | domain | delta | max r_KK/SM | E [PeV] | Qref [GeV] | enhancement | gate |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    keys = sorted({(r["channel"], r["domain"], r["delta"]) for r in rows})
    for key in keys:
        rr = [r for r in rows if (r["channel"], r["domain"], r["delta"]) == key]
        best = max(rr, key=lambda z: z["r_KK_over_SM"])
        lines.append(
            f"| {key[0]} | {key[1]} | {key[2]} | "
            f"{100*best['r_KK_over_SM']:.4g}% | {best['E_PeV']:.0f} | "
            f"{best['qref_input_GeV']:.0f} | {best['enhancement_KK_over_single']:.3g} | "
            f"{'PASS' if best['pass_5pct'] else 'FAIL'} |"
        )

    survivors = [r for r in rows if r["pass_5pct"]]
    lines += ["", "## Decision", ""]
    if survivors:
        lines.append(
            "At least one screening point exceeds 5%. These points are **candidates only**: before Earth "
            "propagation, replace the Qref bracket by an explicit tower-aware collider/width likelihood and "
            "check bulk-neutrino cosmology/oscillation limits plus perturbative unitarity."
        )
    else:
        lines.append(
            "No collider-normalized screening point reaches 5%; do not spend transport time on this branch "
            "without changing the UV completion or constraint normalization."
        )
    lines += [
        "",
        "The electroweak photon/Z dipole portal is not folded here because the repository does not yet contain",
        "a certified massive-final-state differential kernel for it. The dipole channel here is the already-validated",
        "PKR Eq.(4.3) resolved transition dipole; adding the EW dipole is a separate v0.94 calculation.",
        "",
    ]
    (out / "KK_TOWER_STERILE_SCREEN_v0.93.md").write_text("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", default="CT18NLO")
    ap.add_argument("--nx", type=int, default=96)
    ap.add_argument("--energies-pev", default="100,300,1000")
    ap.add_argument("--deltas", default="1,2,3")
    ap.add_argument("--nproc", type=int, default=int(os.environ.get("NUTAU_NPROC", "8")))
    ap.add_argument("--outdir", default="results/v0.93_kk_tower")
    args = ap.parse_args()

    energies = tuple(float(x) for x in args.energies_pev.split(",") if x.strip())
    deltas = tuple(int(x) for x in args.deltas.split(",") if x.strip())
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    jobs = build_jobs(args.pdf, args.nx, energies, deltas)
    rows = []
    nproc = max(1, min(args.nproc, len(jobs)))
    if nproc == 1:
        rows = [one_job(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=nproc) as ex:
            futs = [ex.submit(one_job, j) for j in jobs]
            for i, fut in enumerate(as_completed(futs), 1):
                rows.append(fut.result())
                if i % 12 == 0 or i == len(futs):
                    print(f"[v0.93] completed {i}/{len(futs)} jobs", flush=True)

    rows.sort(key=lambda r: (r["channel"], r["domain"], r["delta"], r["qref_input_GeV"], r["E_PeV"]))
    fields = list(rows[0].keys())
    with open(out / "kk_tower_screen_v0_93.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    meta = dict(
        version="v0.93",
        pdf=args.pdf,
        nx=args.nx,
        energies_PeV=energies,
        deltas=deltas,
        nproc=nproc,
        tower_model="dense flat KK continuum; collider-normalized N_eff ratio; common UV cutoff",
        d6=D6,
        pkr=PKR,
        caveats=[
            "screen only; not a tower-aware collider recast",
            "continuum approximation assumes Q_ref R >> 1",
            "bulk radius cancels in the dense-tower ratio",
            "EW photon/Z dipole portal not included in this version",
        ],
    )
    (out / "metadata_v0_93.json").write_text(json.dumps(meta, indent=2))
    write_summary(rows, out)

    best = sorted(rows, key=lambda r: r["r_KK_over_SM"], reverse=True)[:12]
    print("\nTop KK-screen points:")
    for r in best:
        print(
            f"{r['channel']:22s} {r['domain']:9s} d={r['delta']} "
            f"E={r['E_PeV']:4.0f} PeV Qref={r['qref_input_GeV']:5.0f} GeV "
            f"enh={r['enhancement_KK_over_single']:8.3g} "
            f"rKK={100*r['r_KK_over_SM']:8.3g}% "
            f"{'PASS' if r['pass_5pct'] else 'FAIL'}"
        )
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
