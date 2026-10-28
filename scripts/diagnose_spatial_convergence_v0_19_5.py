#!/usr/bin/env python3
"""v0.19.5 targeted spatial-convergence diagnostic.

Purpose
-------
Continue the v0.19.4 Case-7 shallow-chord full A-B-A convergence study to
0.1875 and 0.09375 km, using the same PREM-interface-aligned exactly nested
construction.  Reuse the validated-by-provenance v0.19.4 0.75/0.375 km
checkpoints; do not alter any production gate.

Target: Case 7, Lambda=3 TeV, |cos z|=0.1, Majorana flip.
"""
from __future__ import annotations

import importlib.util
import json
import math
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE_SCRIPT = ROOT / "scripts" / "diagnose_transport_v0_19_4.py"
BASE_CK = ROOT / "results" / "v0.19.4" / "checkpoints"
OUT = ROOT / "results" / "v0.19.5"
CK = OUT / "checkpoints"
NEW_STEPS = (0.1875, 0.09375)
ALL_STEPS = (0.75, 0.375, 0.1875, 0.09375)
SPECIES = ("nut", "anut", "num", "anum", "tau_minus", "tau_plus")
TARGETS = (1.0, 10.0, 100.0, 220.0, 300.0)


def load_base():
    if not BASE_SCRIPT.exists():
        raise FileNotFoundError(BASE_SCRIPT)
    spec = importlib.util.spec_from_file_location("diag194", BASE_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


M = load_base()


def nested_profile_ext(mu, step_km):
    factors = {0.75: 4, 0.375: 8, 0.1875: 16, 0.09375: 32}
    step_km = float(step_km)
    if step_km not in factors:
        raise ValueError(f"unsupported nested step {step_km}")
    fac = factors[step_km]
    mu = float(mu)
    L = 2.0 * M.R_EARTH_KM * mu
    impact = M.R_EARTH_KM * np.sqrt(max(0.0, 1.0 - mu * mu))
    half = 0.5 * L
    cuts = [-half, half]
    for r in M.PREM_INTERFACES_KM:
        if r > impact:
            x = np.sqrt(max(0.0, r*r - impact*impact))
            if 1e-12 < x < half - 1e-12:
                cuts.extend((-x, x))
    cuts = np.array(sorted(set(round(float(x), 12) for x in cuts)), float)
    edges = [cuts[0]]
    for a, b in zip(cuts[:-1], cuts[1:]):
        nmacro = max(1, int(np.ceil((b-a) / 3.0)))
        macro = np.linspace(a, b, nmacro + 1)
        for u, v in zip(macro[:-1], macro[1:]):
            sub = np.linspace(u, v, fac + 1)
            edges.extend(sub[1:].tolist())
    edges = np.asarray(edges, float)
    ds = np.diff(edges)
    mids = 0.5 * (edges[:-1] + edges[1:])
    r = np.sqrt(impact*impact + mids*mids)
    rho = M.density_g_cm3(r)
    return dict(ds_km=ds, rho=np.asarray(rho, float), r_km=r, edges_km=edges)


M.nested_profile = nested_profile_ext


def tag(x):
    return (f"{float(x):.8g}").replace(".", "p")


def atomic_json(path, obj):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2))
    tmp.replace(path)


def fine_worker(step):
    CK.mkdir(parents=True, exist_ok=True)
    npz = CK / f"full_nested_step_{tag(step)}.npz"
    meta = CK / f"full_nested_step_{tag(step)}.json"
    if npz.exists() and meta.exists():
        return json.loads(meta.read_text())
    t0 = time.perf_counter()
    E, inc, out, prof = M.full_expm_profile(float(step))
    np.savez_compressed(npz, E=E, incident=inc, **out)
    rec = dict(step_km=float(step), seconds=time.perf_counter()-t0,
               n_segments=int(len(prof["ds_km"])),
               min_flux={k: float(np.min(v)) for k, v in out.items()},
               negative_bins={k: int(np.sum(v < 0)) for k, v in out.items()},
               finite=bool(all(np.all(np.isfinite(v)) for v in out.values())))
    atomic_json(meta, rec)
    return rec


def load_step(step):
    if step in (0.75, 0.375):
        p = BASE_CK / f"full_nested_step_{tag(step)}.npz"
    else:
        p = CK / f"full_nested_step_{tag(step)}.npz"
    if not p.exists():
        raise FileNotFoundError(p)
    z = np.load(p)
    return z["E"], z["incident"], {k: z[k] for k in SPECIES}


def interp(E, y, target):
    j = int(np.searchsorted(E, target))
    if j == 0 or j >= len(E):
        return np.nan
    i = j - 1
    if not (np.isfinite(y[i]) and np.isfinite(y[j])):
        return np.nan
    if y[i] > 0 and y[j] > 0:
        a = (np.log(target)-np.log(E[i]))/(np.log(E[j])-np.log(E[i]))
        return float(np.exp(np.log(y[i]) + a*(np.log(y[j])-np.log(y[i]))))
    a = (target-E[i])/(E[j]-E[i])
    return float(y[i] + a*(y[j]-y[i]))


def adjacent_metrics(coarse, fine):
    E, inc, A = load_step(coarse)
    E2, inc2, B = load_step(fine)
    if not np.allclose(E, E2) or not np.allclose(inc, inc2):
        raise RuntimeError("energy/input grids differ")
    maxima = {sp: 0.0 for sp in SPECIES}
    rows = []
    l1 = {}
    for ep in TARGETS:
        ii = interp(E, inc, ep*1e6)
        for sp in SPECIES:
            a = interp(E, A[sp], ep*1e6); b = interp(E, B[sp], ep*1e6)
            den = max(abs(b), abs(ii)*1e-12, 1e-300)
            s = (a-b)/den
            maxima[sp] = max(maxima[sp], abs(s))
            rows.append(dict(E_PeV=ep, species=sp, coarse=a, fine=b,
                             signed_relative_difference=s,
                             absolute_relative_difference=abs(s)))
    for sp in SPECIES:
        scale = max(float(np.max(np.abs(B[sp]))), 1e-300)
        mask = (inc > 1e-12*np.max(inc)) & (np.abs(B[sp]) > 1e-12*scale)
        l1[sp] = 0.0 if not np.any(mask) else float(
            np.sum(np.abs(A[sp][mask]-B[sp][mask])) /
            max(np.sum(np.abs(B[sp][mask])), 1e-300))
    return dict(coarse_step_km=coarse, fine_step_km=fine,
                maxima=maxima, L1=l1, pointwise=rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True); CK.mkdir(parents=True, exist_ok=True)
    n = min(int(os.environ.get("NUTAU_NPROC", "2")), len(NEW_STEPS))
    print(f"v0.19.5: computing new nested fine grids {NEW_STEPS} with workers={n}", flush=True)
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=n) as pool:
        meta = list(pool.imap_unordered(fine_worker, NEW_STEPS, chunksize=1))
    for m in sorted(meta, key=lambda x: x["step_km"], reverse=True):
        print(f"FINE step={m['step_km']:.6g}km sec={m['seconds']:.1f} segments={m['n_segments']} finite={m['finite']} negatives={sum(m['negative_bins'].values())}", flush=True)

    pairs = [(0.75,0.375),(0.375,0.1875),(0.1875,0.09375)]
    rr = [adjacent_metrics(a,b) for a,b in pairs]
    for r in rr:
        mx = r["maxima"]
        prim = max(mx["nut"], mx["anut"])
        sec = max(mx["num"], mx["anum"])
        tau = max(mx["tau_minus"], mx["tau_plus"])
        print(f"ADJ {r['coarse_step_km']:.6g}->{r['fine_step_km']:.6g} km primary={prim:.6g} secondary={sec:.6g} tau={tau:.6g}", flush=True)

    # Empirical order estimates from adjacent max errors, interpreted only when
    # the same asymptotic regime is approached; maxima may move between bins.
    orders = []
    for cls, sps in (("primary",("nut","anut")),("secondary",("num","anum")),("tau",("tau_minus","tau_plus"))):
        errs = [max(r["maxima"][sp] for sp in sps) for r in rr]
        p1 = math.log(errs[0]/errs[1],2.0) if errs[0] > 0 and errs[1] > 0 else None
        p2 = math.log(errs[1]/errs[2],2.0) if errs[1] > 0 and errs[2] > 0 else None
        orders.append(dict(channel=cls, errors=errs, p_075_0375_to_0375_01875=p1,
                           p_0375_01875_to_01875_009375=p2))
        print(f"ORDER {cls} p12={p1 if p1 is not None else float('nan'):.3f} p23={p2 if p2 is not None else float('nan'):.3f}", flush=True)

    status = dict(version="v0.19.5-diagnostic", target="case7_L3_cosz0.1_majorana_flip",
                  reused_steps_km=[0.75,0.375], computed_steps_km=list(NEW_STEPS),
                  adjacent=rr, empirical_orders=orders,
                  conclusion_rule="Do not run energy-grid refinement unless spatial adjacent errors plateau; diagnostics only, no production certification.")
    atomic_json(OUT/"spatial_convergence_status_v0_19_5.json", status)
    print("STATUS wrote", OUT/"spatial_convergence_status_v0_19_5.json", flush=True)


if __name__ == "__main__":
    main()
