#!/usr/bin/env python3
"""v0.20.1 tau energy-grid diagnostic for the combined-generator candidate.

Diagnostics only.  No full chord and no production certification.

Purpose
-------
The v0.20 smoke showed exact discrete semigroup composition but large L1 shape
error for the surviving tau spectrum on the historical 91-point energy grid.
This script isolates that issue from PDFs and neutrino interactions by testing
the tau-only combined generator on 91, 181 and 361 log-spaced energy grids over
exactly the historical energy range.

For each grid it compares one 0.75-km constant-density slab against the analytic
CoupledTauEvolution map at 1, 10, 100 and 300 PeV.  In addition to tau L1 it
reports survival-number error, mean-logE error, a 1D Wasserstein distance in
logE, and primary/secondary daughter L1 errors.  These metrics distinguish
harmless displacement of a narrow tau packet from genuine broad numerical
diffusion.
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import time
from pathlib import Path

import numpy as np
from scipy.sparse import bmat, csr_matrix

from nutau_prop.chirality_transport_full_expm import TauHelicityFlipTransportFullExpm
from nutau_prop.tau_evolution import CoupledTauEvolution
from nutau_prop.prem import chord_profile

ROOT = Path(__file__).resolve().parents[1]
T15 = ROOT / "results" / "v0.15" / "tables"
OUT = ROOT / "results" / "v0.20.1"
NS = (91, 181, 361)
TARGETS_PEV = (1.0, 10.0, 100.0, 300.0)
STEP_KM = 0.75


def atomic_json(path, obj):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2))
    tmp.replace(path)


def historical_range():
    z = np.load(T15 / "case7_L3_nu.npz")
    E = np.asarray(z["energies"], float)
    return float(E[0]), float(E[-1])


def dummy_table(E):
    """Minimal zero-interaction table; only tau machinery is exercised."""
    E = np.asarray(E, float)
    N = len(E)
    y = np.array([0.5], float)
    p = np.ones((N, 1), float)
    z = np.zeros(N, float)
    return dict(energies=E, y=y, sigma_cc=z.copy(), sigma_nc_sm=z.copy(),
                sigma_np=z.copy(), p_cc=p.copy(), p_nc_sm=p.copy(), p_np=p.copy())


def tau_only_generator(tr, rho):
    N = tr.N
    Z = csr_matrix((N, N))
    L = float(rho) * tr._tau_loss_per_rho()
    Q, Dm, Dp, Sm, Sp = tr._tau_decay_blocks(True)
    return bmat([
        [Z,Z,Z,Z,Z,Z], [Z,Z,Z,Z,Z,Z],
        [Z,Z,Z,Z,Z,Z], [Z,Z,Z,Z,Z,Z],
        [Dm,Z,Z,Sm,L+Q,Z],
        [Z,Dp,Sp,Z,Z,L+Q],
    ], format="csr")


def rel_l1(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return float(np.sum(np.abs(a-b)) / max(np.sum(np.abs(b)), 1e-300))


def wasserstein_logE(E, a, b):
    """W1 between normalized tau shapes in lnE; number error is reported separately."""
    a = np.clip(np.asarray(a, float), 0.0, None)
    b = np.clip(np.asarray(b, float), 0.0, None)
    sa = float(np.sum(a)); sb = float(np.sum(b))
    if sa <= 0 or sb <= 0:
        return np.nan
    aa = a/sa; bb = b/sb
    x = np.log(np.asarray(E, float))
    ca = np.cumsum(aa); cb = np.cumsum(bb)
    # CDF is piecewise constant between grid points.
    return float(np.sum(np.abs(ca[:-1]-cb[:-1]) * np.diff(x)))


def mean_logE(E, a):
    a = np.clip(np.asarray(a, float), 0.0, None)
    s = float(np.sum(a))
    return np.nan if s <= 0 else float(np.sum(a*np.log(E))/s)


def worker(N):
    emin, emax = historical_range()
    E = np.geomspace(emin, emax, int(N))
    tab = dummy_table(E)
    tr = TauHelicityFlipTransportFullExpm(tab, tab)
    prof = chord_profile(0.1, max_step_km=3.0)
    rho = float(np.median(prof["rho"]))
    h = STEP_KM * 1e5
    G = tau_only_generator(tr, rho)
    exact = CoupledTauEvolution(E, beta_tau=tr.beta)
    rows = []
    t0 = time.perf_counter()
    for epev in TARGETS_PEV:
        j = int(np.argmin(np.abs(np.log(E)-np.log(epev*1e6))))
        tm0 = np.zeros(N); tm0[j] = 1.0
        z = np.zeros(N)
        x = np.concatenate([z,z,z,z,tm0,z])
        y = tr._advance(G, x, h)
        nut, anut, num, anum, tm, tp = [y[i*N:(i+1)*N] for i in range(6)]
        tm_ex, nut_ex, anum_ex, _ = exact.advance_with_kernels(
            tm0, rho, h, tr.Kdec_minus, tr.Ksec_minus)
        sn = float(np.sum(tm)); se = float(np.sum(tm_ex))
        rows.append(dict(
            N=int(N), E_requested_PeV=float(epev), E_grid_PeV=float(E[j]/1e6),
            tau_L1=rel_l1(tm, tm_ex),
            tau_number_new=sn, tau_number_exact=se,
            tau_number_relative_error=abs(sn-se)/max(abs(se),1e-300),
            tau_mean_logE_difference=mean_logE(E,tm)-mean_logE(E,tm_ex),
            tau_W1_logE=wasserstein_logE(E,tm,tm_ex),
            primary_daughter_L1=rel_l1(nut,nut_ex),
            secondary_daughter_L1=rel_l1(anum,anum_ex),
            finite=bool(all(np.all(np.isfinite(q)) for q in (tm,nut,anum))),
            negative_bins=int(sum(np.sum(q < 0) for q in (tm,nut,anum))),
        ))
    return dict(N=int(N), seconds=time.perf_counter()-t0, rho_g_cm3=rho,
                step_km=STEP_KM, rows=rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    nproc = min(int(os.environ.get("NUTAU_NPROC", "3")), len(NS))
    print(f"v0.20.1 tau energy-grid diagnostic workers={nproc}", flush=True)
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=nproc) as pool:
        rr = list(pool.imap_unordered(worker, NS, chunksize=1))
    rr = sorted(rr, key=lambda r:r["N"])
    for r in rr:
        print(f"GRID N={r['N']} sec={r['seconds']:.2f}", flush=True)
        for q in r["rows"]:
            print("  E={:g}PeV tauL1={:.6g} W1logE={:.6g} dmeanlogE={:.6g} numerr={:.6g} primary={:.6g} secondary={:.6g}".format(
                q["E_requested_PeV"], q["tau_L1"], q["tau_W1_logE"],
                q["tau_mean_logE_difference"], q["tau_number_relative_error"],
                q["primary_daughter_L1"], q["secondary_daughter_L1"]), flush=True)
    # Simple refinement ratios for physically relevant daughter errors and W1.
    trends=[]
    for epev in TARGETS_PEV:
        vals=[]
        for r in rr:
            q=next(x for x in r["rows"] if x["E_requested_PeV"]==epev)
            vals.append(q)
        tr=dict(E_PeV=epev,
                W1=[q["tau_W1_logE"] for q in vals],
                primary=[q["primary_daughter_L1"] for q in vals],
                secondary=[q["secondary_daughter_L1"] for q in vals])
        trends.append(tr)
    status=dict(version="v0.20.1-diagnostic", grids=list(NS),
                energy_range_GeV=list(historical_range()), step_km=STEP_KM,
                results=rr, trends=trends,
                interpretation_rule=(
                    "Keep the combined-generator architecture only if tau transport metrics "
                    "and regenerated-daughter errors decrease systematically with energy-grid refinement. "
                    "No full-chord production is authorized by this diagnostic."))
    atomic_json(OUT/"tau_energy_grid_status_v0_20_1.json", status)
    print("STATUS wrote", OUT/"tau_energy_grid_status_v0_20_1.json", flush=True)


if __name__ == "__main__":
    main()
