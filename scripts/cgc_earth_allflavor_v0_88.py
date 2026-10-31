#!/usr/bin/env python3
"""v0.88 — all-flavour Earth propagation for a CGC/rcBK control model.

Physics included
----------------
* CT18NLO isoscalar SM CC+NC differential DIS baseline.
* PREM Earth profile.
* NC energy migration for all flavours.
* CC absorption for nu_e and nu_mu.
* Full nu_tau <-> tau production, continuous tau energy loss, tau decay and
  nu_tau regeneration using the validated mixed-grid sparse-exponential solver.
* Leptonic tau-decay secondary nu_e and nu_mu (and antiparticles), each with
  branching fraction 0.18 and the same massless daughter spectrum.
* Separate neutrino/antineutrino evolution; reported flavour fluxes sum both.

CGC implementation
------------------
The primary model is a literature-calibrated rcBK+DLA-2015 total-cross-section
surrogate.  CC and NC totals are multiplied by R_CGC(E), while the conditional
SM y distributions are retained.  This isolates the published saturation
cross-section effect without inventing an unpublished CGC d sigma/dy table.

The anchor points are documented in src/nutau_prop/cgc.py and in the v0.88
report.  This is therefore a propagation benchmark, not a claim of an exact
machine-readable reproduction of the 2015 differential rcBK calculation.
"""
from __future__ import annotations

import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from nutau_prop.case6 import Case6Config, Case6DIS, probability_from_dsdy
from nutau_prop.cgc import LightFlavorTransport, rcbk_dla_2015_ratio, rescale_sm_table
from nutau_prop.chirality_transport_mixed_expm import TauHelicityFlipTransportMixedExpm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.88_cgc_allflavor"
CACHE = ROOT / "cache" / "v0.88_cgc_allflavor"
OUT.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)

PDF = "CT18NLO"
NPROC = int(os.environ.get("NUTAU_NPROC", "14"))
NX = 64
YLOG = 26
YLIN = 26
TAU_REFINE = 8
MAX_STEP_KM = 48.0
GAMMAS = (2.0, 2.5)
TARGET_E = np.asarray([1e7, 3e7, 1e8, 3e8, 1e9, 3e9, 1e10], float)
EGRID = np.unique(np.r_[np.geomspace(3e4, 3e10, 67), TARGET_E])
REPORT_MUS = np.asarray([0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 1.00], float)
_GLX, _GLW = np.polynomial.legendre.leggauss(12)
QUAD_MUS = 0.5 * (_GLX + 1.0)
QUAD_W = 0.5 * _GLW
ALL_MUS = np.unique(np.r_[REPORT_MUS, QUAD_MUS])

_CALC = None


def _calc():
    global _CALC
    if _CALC is None:
        _CALC = Case6DIS(Case6Config(pdf_set=PDF, nx=NX,
                                     y_log_bins=YLOG, y_linear_bins=YLIN))
    return _CALC


def _sm_row(task):
    E, antinu = task
    c = _calc()
    _, y, dy = c.y_grid()
    cc = np.zeros_like(y); nc = np.zeros_like(y)
    for i, yy in enumerate(y):
        cc[i], nc[i] = c.dsdy_sm(float(E), float(yy), antinu=bool(antinu))
    scc, pcc = probability_from_dsdy(cc, dy)
    snc, pnc = probability_from_dsdy(nc, dy)
    return float(E), bool(antinu), scc, snc, pcc, pnc, y, dy


def build_sm_tables():
    tasks = [(float(E), anti) for anti in (False, True) for E in EGRID]
    rows = []
    with ProcessPoolExecutor(max_workers=NPROC) as ex:
        fut = [ex.submit(_sm_row, t) for t in tasks]
        for f in as_completed(fut):
            rows.append(f.result())
    tabs = {}
    for anti in (False, True):
        q = sorted((r for r in rows if r[1] == anti), key=lambda r: r[0])
        E = np.asarray([r[0] for r in q], float)
        y = np.asarray(q[0][6], float); dy = np.asarray(q[0][7], float)
        tab = dict(
            energies=E, y=y, dy=dy,
            sigma_cc=np.asarray([r[2] for r in q], float),
            sigma_nc_sm=np.asarray([r[3] for r in q], float),
            sigma_np=np.zeros(len(q), float),
            p_cc=np.asarray([r[4] for r in q], float),
            p_nc_sm=np.asarray([r[5] for r in q], float),
            p_np=np.zeros((len(q), len(y)), float),
        )
        tabs[anti] = tab
    return tabs


def save_table(path, tab):
    np.savez_compressed(path, **tab)


def load_table(path):
    z = np.load(path)
    return {k: z[k] for k in z.files}


def prepare_tables():
    pnu = CACHE / "sm_nu.npz"; pan = CACHE / "sm_antinu.npz"
    if not (pnu.exists() and pan.exists()):
        tabs = build_sm_tables()
        save_table(pnu, tabs[False]); save_table(pan, tabs[True])
    smn, sma = load_table(pnu), load_table(pan)
    ratio = rcbk_dla_2015_ratio(smn["energies"])
    cgcn = rescale_sm_table(smn, ratio)
    cgca = rescale_sm_table(sma, ratio)
    save_table(CACHE / "cgc_nu.npz", cgcn)
    save_table(CACHE / "cgc_antinu.npz", cgca)
    pd.DataFrame({"E_GeV": smn["energies"], "sigma_CGC_over_SM": ratio,
                  "sigma_cc_nu_SM_cm2": smn["sigma_cc"],
                  "sigma_nc_nu_SM_cm2": smn["sigma_nc_sm"]}).to_csv(
        OUT / "cross_section_ratio.csv", index=False)


def _transport_one(task):
    model, mu, gamma = task
    if model == "SM":
        tn = load_table(CACHE / "sm_nu.npz"); ta = load_table(CACHE / "sm_antinu.npz")
    else:
        tn = load_table(CACHE / "cgc_nu.npz"); ta = load_table(CACHE / "cgc_antinu.npz")
    E = np.asarray(tn["energies"], float)
    # Unit total incident flux per flavour; half neutrino, half antineutrino.
    f = (E / 1.0e6) ** (-float(gamma))
    n0 = 0.5 * f; a0 = 0.5 * f

    light = LightFlavorTransport(tn, ta).propagate(
        float(mu), n0, a0, max_step_km=MAX_STEP_KM)
    tau = TauHelicityFlipTransportMixedExpm(
        tn, ta, tau_refine=TAU_REFINE, secondary_branching=0.18).propagate(
            float(mu), n0, a0, mode="dirac_sink", use_np=False,
            max_step_km=MAX_STEP_KM, track_secondary_numu=True)

    # tau- -> anti-nu_l and tau+ -> nu_l. Electron and muon leptonic spectra
    # are identical in the massless-daughter approximation used by the solver.
    sec_n = np.asarray(tau["nu_mu_secondary"], float)
    sec_a = np.asarray(tau["antinu_mu_secondary"], float)

    out = dict(
        E=E,
        nu_e=np.asarray(light["nu"], float) + sec_n,
        antinu_e=np.asarray(light["antinu"], float) + sec_a,
        nu_mu=np.asarray(light["nu"], float) + sec_n,
        antinu_mu=np.asarray(light["antinu"], float) + sec_a,
        nu_tau=np.asarray(tau["nu_tau"], float),
        antinu_tau=np.asarray(tau["antinu_tau"], float),
    )
    return model, float(mu), float(gamma), out


def interpolate_log(E, y, Et):
    E = np.asarray(E, float); y = np.asarray(y, float)
    if np.all(y > 0):
        return float(np.exp(np.interp(np.log(Et), np.log(E), np.log(y))))
    return float(np.interp(np.log(Et), np.log(E), y))


def main():
    prepare_tables()
    tasks = [(m, float(mu), float(g)) for g in GAMMAS for m in ("SM", "CGC") for mu in ALL_MUS]
    rec = {}
    with ProcessPoolExecutor(max_workers=NPROC) as ex:
        fut = [ex.submit(_transport_one, t) for t in tasks]
        for f in as_completed(fut):
            model, mu, gamma, out = f.result()
            rec[(model, gamma, round(mu, 12))] = out

    flavours = ("e", "mu", "tau")
    grid_rows = []
    for gamma in GAMMAS:
        for mu in REPORT_MUS:
            sm = rec[("SM", gamma, round(float(mu), 12))]
            cg = rec[("CGC", gamma, round(float(mu), 12))]
            E = sm["E"]
            incident = (E / 1e6) ** (-gamma)
            for fl in flavours:
                s = sm[f"nu_{fl}"] + sm[f"antinu_{fl}"]
                c = cg[f"nu_{fl}"] + cg[f"antinu_{fl}"]
                for i, ee in enumerate(E):
                    grid_rows.append(dict(gamma=gamma, cosz_abs=mu, E_GeV=ee,
                                          flavour=fl, flux_SM=s[i], flux_CGC=c[i],
                                          ratio_CGC_SM=c[i]/max(s[i], 1e-300),
                                          delta_pct=100*(c[i]/max(s[i],1e-300)-1),
                                          survival_SM=s[i]/incident[i],
                                          survival_CGC=c[i]/incident[i]))
    G = pd.DataFrame(grid_rows)
    G.to_csv(OUT / "flux_grid.csv", index=False)

    # Isotropic upgoing angle average using Gauss-Legendre integration in mu.
    avg_rows = []
    for gamma in GAMMAS:
        E = rec[("SM", gamma, round(float(QUAD_MUS[0]),12))]["E"]
        incident = (E/1e6)**(-gamma)
        for fl in flavours:
            smav = np.zeros_like(E); cgav = np.zeros_like(E)
            for mu, w in zip(QUAD_MUS, QUAD_W):
                sm = rec[("SM", gamma, round(float(mu),12))]
                cg = rec[("CGC", gamma, round(float(mu),12))]
                smav += w*(sm[f"nu_{fl}"] + sm[f"antinu_{fl}"])
                cgav += w*(cg[f"nu_{fl}"] + cg[f"antinu_{fl}"])
            for i, ee in enumerate(E):
                avg_rows.append(dict(gamma=gamma,E_GeV=ee,flavour=fl,
                                     flux_SM=smav[i],flux_CGC=cgav[i],
                                     ratio_CGC_SM=cgav[i]/max(smav[i],1e-300),
                                     delta_pct=100*(cgav[i]/max(smav[i],1e-300)-1),
                                     survival_SM=smav[i]/incident[i],
                                     survival_CGC=cgav[i]/incident[i]))
    A = pd.DataFrame(avg_rows)
    A.to_csv(OUT / "upgoing_angle_average.csv", index=False)

    # Compact target-energy summary: fixed chords and angle-averaged flux.
    summary = []
    for gamma in GAMMAS:
        for Et in TARGET_E:
            Rxs = float(rcbk_dla_2015_ratio(np.asarray([Et]))[0])
            for fl in flavours:
                qa=A[(A.gamma==gamma)&(A.flavour==fl)]
                sm=interpolate_log(qa.E_GeV,qa.flux_SM,Et)
                cg=interpolate_log(qa.E_GeV,qa.flux_CGC,Et)
                summary.append(dict(gamma=gamma,location="upgoing_avg",cosz_abs=np.nan,
                                    E_GeV=Et,flavour=fl,cgc_xsec_ratio=Rxs,
                                    flux_ratio=cg/max(sm,1e-300),delta_pct=100*(cg/max(sm,1e-300)-1)))
                for mu in (0.10,0.30,0.50,0.70,1.00):
                    q=G[(G.gamma==gamma)&np.isclose(G.cosz_abs,mu)&(G.flavour==fl)]
                    sm=interpolate_log(q.E_GeV,q.flux_SM,Et)
                    cg=interpolate_log(q.E_GeV,q.flux_CGC,Et)
                    summary.append(dict(gamma=gamma,location="chord",cosz_abs=mu,
                                        E_GeV=Et,flavour=fl,cgc_xsec_ratio=Rxs,
                                        flux_ratio=cg/max(sm,1e-300),delta_pct=100*(cg/max(sm,1e-300)-1)))
    S=pd.DataFrame(summary)
    S.to_csv(OUT/"summary.csv",index=False)

    # Report only gamma=2.5, the primary benchmark, and flag opaque bins.
    primary=S[(S.gamma==2.5)&(S.E_GeV.isin([1e8,3e8,1e9]))]
    print("=== v0.88 CGC all-flavour Earth propagation ===")
    print("model: rcBK+DLA-2015 total-cross-section surrogate; SM y-shapes retained")
    print("incident: 1:1:1, nu=antinu, primary gamma=2.5 (gamma=2.0 also saved)")
    print("\nflux change [%], CGC relative to CT18/DGLAP-like SM baseline:")
    print(primary[["location","cosz_abs","E_GeV","flavour","cgc_xsec_ratio","delta_pct"]].to_string(index=False))

    p=primary[primary.location=="upgoing_avg"]
    print("\nupgoing-angle average (gamma=2.5):")
    print(p[["E_GeV","flavour","cgc_xsec_ratio","flux_ratio","delta_pct"]].to_string(index=False))

    meta={
      "version":"v0.88_cgc_allflavor",
      "baseline":"CT18NLO LO transport kernel from Case6DIS; PREM",
      "cgc":"rcBK+DLA-2015 literature-calibrated total cross-section surrogate",
      "cgc_differential_assumption":"same conditional y distributions as SM; CC and NC totals share R_CGC(E)",
      "incident":"equal e:mu:tau; equal nu:antinu; E^-2 and E^-2.5 both run",
      "tau":"mixed-grid tau production/loss/decay/regeneration; leptonic e and mu secondaries each BR=0.18",
      "glashow":"not included; report focus is >=10 PeV nucleon DIS",
      "oscillations":"neglected inside Earth at UHE",
      "nproc":NPROC,"nx":NX,"tau_refine":TAU_REFINE,"max_step_km":MAX_STEP_KM,
      "references":[
        "Albacete, Illana, Soto-Ontoso, Phys Rev D92 014027 (2015), arXiv:1505.06583",
        "Goncalves, Gratieri, Quadros, Eur Phys J C81 496 (2021)"
      ]
    }
    (OUT/"metadata.json").write_text(json.dumps(meta,indent=2)+"\n")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
