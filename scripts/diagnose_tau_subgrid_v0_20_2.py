#!/usr/bin/env python3
"""v0.20.2 mixed-grid diagnostic for the combined-generator transport.

Purpose
-------
Keep the neutrino/daughter grid fixed at the existing 91-point Case-7 table,
while refining only the tau log-energy grid by factors 1,2,4,8,16.  This tests
whether the remaining v0.20 numerical diffusion is localized to tau advection
and can be removed without rebuilding the expensive DIS tables.

The reference is a direct continuous-loss + decay quadrature, independent of
the discrete tau grid.  No full Earth chord and no production output.
"""
from __future__ import annotations

import json, math, os, time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import bmat, csr_matrix, diags, lil_matrix
from scipy.sparse.linalg import expm_multiply

from nutau_prop.constants import M_TAU, C_TAU_CM
from nutau_prop.kernels import tau_nutau_pdf, tau_secondary_numu_pdf

ROOT = Path(__file__).resolve().parents[1]
T15 = ROOT / "results" / "v0.15" / "tables"
OUT = ROOT / "results" / "v0.20.2"
BETA = 0.85e-6
RHO = 2.9
H_CM = 0.75e5
FACTORS = (1, 2, 4, 8, 16)
TARGETS_PEV = (1.0, 10.0, 100.0, 300.0)
NZ = 800
NQ = 48
BR = 0.18


def load_energy_grid():
    z = np.load(T15 / "case7_L3_nu.npz")
    return np.asarray(z["energies"], float)


def tau_grid(E, factor):
    l = np.log(np.asarray(E, float))
    out = [l[0]]
    for a, b in zip(l[:-1], l[1:]):
        out.extend(np.linspace(a, b, int(factor) + 1)[1:].tolist())
    return np.exp(np.asarray(out, float))


def deposit(row, grid, eout, weight):
    if weight == 0.0 or eout < grid[0] or eout > grid[-1]:
        return
    lg = np.log(grid); u = math.log(float(eout)); j = int(np.searchsorted(lg, u))
    if j == 0:
        row[0] += weight
    elif j >= len(grid):
        row[-1] += weight
    else:
        f = (u-lg[j-1])/(lg[j]-lg[j-1])
        row[j-1] += weight*(1.0-f); row[j] += weight*f


def z_weights(sign):
    z = np.linspace(1e-6, 1.0-1e-6, NZ)
    pol = -1.0 if sign == "minus" else +1.0
    fp = tau_nutau_pdf(z, pol)
    fp /= np.trapezoid(fp, z)
    fs = tau_secondary_numu_pdf(z, pol, BR)
    ns = np.trapezoid(fs, z)
    if ns > 0:
        fs *= BR/ns
    dz = z[1]-z[0]
    return z, fp*dz, fs*dz


def rectangular_decay_kernels(Et, Enu, sign):
    z, wp, ws = z_weights(sign)
    Kp = np.zeros((len(Et), len(Enu)), float)
    Ks = np.zeros((len(Et), len(Enu)), float)
    for j, ep in enumerate(Et):
        rp = Kp[j]; rs = Ks[j]
        for zz, a, b in zip(z, wp, ws):
            ee = ep*zz
            if a != 0.0: deposit(rp, Enu, ee, float(a))
            if b != 0.0: deposit(rs, Enu, ee, float(b))
    return csr_matrix(Kp), csr_matrix(Ks)


def tau_loss_generator(Et):
    l = np.log(Et); dl = np.diff(l)
    rate = np.empty(len(Et), float)
    rate[0] = BETA*RHO/dl[0]
    rate[1:] = BETA*RHO/dl
    Q = lil_matrix((len(Et), len(Et)), dtype=float)
    Q.setdiag(-rate)
    for j in range(1, len(Et)):
        Q[j, j-1] = rate[j]
    return Q.tocsr()


def full_tau_generator(Et, Enu, sign):
    Nt, Nn = len(Et), len(Enu)
    Znn = csr_matrix((Nn, Nn)); Znt = csr_matrix((Nn, Nt))
    Kp, Ks = rectangular_decay_kernels(Et, Enu, sign)
    d = M_TAU/(C_TAU_CM*Et)
    Dp = csr_matrix(d[:, None]*Kp.toarray())
    Ds = csr_matrix(d[:, None]*Ks.toarray())
    T = tau_loss_generator(Et) - diags(d, 0, format="csr")
    return bmat([[Znn, Znn, Znt],
                 [Znn, Znn, Znt],
                 [Dp,  Ds,  T]], format="csr")


def exact_reference(E0, Enu, Et, sign):
    pol = -1.0 if sign == "minus" else +1.0
    b = BETA*RHO
    a = M_TAU/(C_TAU_CM*E0)
    U = (a/b)*np.expm1(b*H_CM) if b != 0 else a*H_CM
    surv = math.exp(-min(float(U), 745.0))
    Eend = E0*math.exp(-b*H_CM)
    tau_proj = np.zeros(len(Et)); deposit(tau_proj, Et, Eend, surv)

    primary = np.zeros(len(Enu)); secondary = np.zeros(len(Enu))
    pdec = -np.expm1(-min(float(U), 745.0))
    if pdec > 0:
        qx, qw = leggauss(NQ)
        tt = 0.5*(qx+1.0)*pdec
        ww = 0.5*pdec*qw
        uu = -np.log1p(-tt)
        Edec = E0/(1.0 + (b/a)*uu) if b != 0 else np.full_like(uu, E0)
        z = np.linspace(1e-6, 1.0-1e-6, NZ)
        fp = tau_nutau_pdf(z, pol); fp /= np.trapezoid(fp, z)
        fs = tau_secondary_numu_pdf(z, pol, BR); ns = np.trapezoid(fs, z)
        if ns > 0: fs *= BR/ns
        dz = z[1]-z[0]
        for ed, wdec in zip(Edec, ww):
            for zz, pp, ss in zip(z, fp, fs):
                ee = float(ed*zz)
                if pp != 0: deposit(primary, Enu, ee, float(wdec*pp*dz))
                if ss != 0: deposit(secondary, Enu, ee, float(wdec*ss*dz))
    return dict(tau=tau_proj, primary=primary, secondary=secondary,
                survival=surv, Eend=Eend)


def rel_l1(a, b):
    return float(np.sum(np.abs(np.asarray(a)-np.asarray(b))) /
                 max(np.sum(np.abs(b)), 1e-300))


def run_one(Enu, factor, sign):
    t0 = time.perf_counter(); Et = tau_grid(Enu, factor)
    G = full_tau_generator(Et, Enu, sign)
    build_sec = time.perf_counter()-t0
    rows = []
    for target in TARGETS_PEV:
        jn = int(np.argmin(np.abs(np.log(Enu)-np.log(target*1e6))))
        E0 = float(Enu[jn])  # exactly a node of every nested tau grid
        jt = int(np.argmin(np.abs(np.log(Et)-np.log(E0))))
        x = np.zeros(2*len(Enu)+len(Et)); x[2*len(Enu)+jt] = 1.0
        s0 = time.perf_counter()
        y = np.asarray(expm_multiply(G.T*H_CM, x), float)
        evolve_sec = time.perf_counter()-s0
        p = y[:len(Enu)]; s = y[len(Enu):2*len(Enu)]; tau = y[2*len(Enu):]
        ref = exact_reference(E0, Enu, Et, sign)
        nt = float(np.sum(tau)); nr = float(ref["survival"])
        if nt > 0:
            prob = tau/nt
            meanlog = float(np.sum(prob*np.log(Et)))
            w1 = float(np.sum(prob*np.abs(np.log(Et)-math.log(ref["Eend"]))))
        else:
            meanlog = float("nan"); w1 = float("nan")
        row = dict(factor=int(factor), sign=sign, N_nu=len(Enu), N_tau=len(Et),
                   E_requested_PeV=target, E0_grid_PeV=E0/1e6,
                   tau_L1=rel_l1(tau, ref["tau"]),
                   tau_number_relative_error=abs(nt-nr)/max(nr,1e-300),
                   W1_logE=w1,
                   delta_mean_logE=abs(meanlog-math.log(ref["Eend"])) if np.isfinite(meanlog) else None,
                   primary_daughter_L1=rel_l1(p, ref["primary"]),
                   secondary_daughter_L1=rel_l1(s, ref["secondary"]),
                   build_seconds=build_sec, evolve_seconds=evolve_sec)
        rows.append(row)
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    Enu = load_energy_grid()
    allrows = []
    for factor in FACTORS:
        for sign in ("minus", "plus"):
            rr = run_one(Enu, factor, sign); allrows.extend(rr)
            worst_p = max(r["primary_daughter_L1"] for r in rr)
            worst_s = max(r["secondary_daughter_L1"] for r in rr)
            worst_w = max(r["W1_logE"] for r in rr)
            print(f"SUBGRID factor={factor:2d} sign={sign:5s} Ntau={rr[0]['N_tau']:4d} "
                  f"worst_primary={worst_p:.6g} worst_secondary={worst_s:.6g} "
                  f"worst_W1logE={worst_w:.6g} build_sec={rr[0]['build_seconds']:.2f}", flush=True)
            for r in rr:
                print("  E={:.0f}PeV tauL1={:.6g} W1={:.6g} numerr={:.6g} primary={:.6g} secondary={:.6g} evolve_sec={:.3f}".format(
                    r["E_requested_PeV"], r["tau_L1"], r["W1_logE"],
                    r["tau_number_relative_error"], r["primary_daughter_L1"],
                    r["secondary_daughter_L1"], r["evolve_seconds"]), flush=True)

    status = dict(version="v0.20.2-diagnostic", neutrino_grid_points=len(Enu),
                  tau_refinement_factors=list(FACTORS), rho_g_cm3=RHO,
                  slab_km=H_CM/1e5, rows=allrows,
                  interpretation="Diagnostic only. Prefer a tau-only refined subgrid if daughter and W1 errors converge without globally refining the DIS grid.")
    (OUT/"tau_subgrid_status_v0_20_2.json").write_text(json.dumps(status, indent=2))
    print("STATUS wrote", OUT/"tau_subgrid_status_v0_20_2.json", flush=True)


if __name__ == "__main__":
    main()
