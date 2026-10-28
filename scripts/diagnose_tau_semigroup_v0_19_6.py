#!/usr/bin/env python3
"""v0.19.6 diagnostics: tau-map semigroup defect + one finer spatial step.

Diagnostics only; does not alter any production gate.

1) Reuse v0.19.5 0.09375-km result and compute only 0.046875 km on the
   same PREM-interface-aligned nested construction.
2) Test whether the discrete tau A-map is a semigroup by comparing A(h) with
   A(h/2)A(h/2) for representative tau energies and shallow-chord densities.

Target full transport: Case 7, Lambda=3 TeV, |cos z|=0.1, Majorana flip.
"""
from __future__ import annotations

import importlib.util
import json
import math
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPT195 = ROOT / "scripts" / "diagnose_spatial_convergence_v0_19_5.py"
OUT = ROOT / "results" / "v0.19.6"
CK = OUT / "checkpoints"
OLDCK = ROOT / "results" / "v0.19.5" / "checkpoints"
STEP_FINE = 0.046875
STEP_COARSE = 0.09375
TARGETS = (1.0, 10.0, 100.0, 220.0, 300.0)
TAU_TEST_E_PEV = (1.0, 10.0, 100.0, 300.0)
TAU_TEST_STEPS_KM = (0.75, 0.375, 0.1875, 0.09375, 0.046875)


def load195():
    spec = importlib.util.spec_from_file_location("diag195", SCRIPT195)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


M5 = load195()
M = M5.M


def tag(x):
    return (f"{float(x):.8g}").replace(".", "p")


def atomic_json(path, obj):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2))
    tmp.replace(path)


def nested_profile_ext(mu, step_km):
    factors = {0.75: 4, 0.375: 8, 0.1875: 16, 0.09375: 32, 0.046875: 64}
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


# Patch the base diagnostic so full_expm_profile uses the extended nested grid.
M.nested_profile = nested_profile_ext


def compute_fine():
    CK.mkdir(parents=True, exist_ok=True)
    p = CK / f"full_nested_step_{tag(STEP_FINE)}.npz"
    meta = CK / f"full_nested_step_{tag(STEP_FINE)}.json"
    if p.exists() and meta.exists():
        return json.loads(meta.read_text())
    t0 = time.perf_counter()
    E, inc, out, prof = M.full_expm_profile(STEP_FINE)
    np.savez_compressed(p, E=E, incident=inc, **out)
    rec = dict(step_km=STEP_FINE, seconds=time.perf_counter()-t0,
               n_segments=int(len(prof["ds_km"])),
               finite=bool(all(np.all(np.isfinite(v)) for v in out.values())),
               negative_bins={k: int(np.sum(v < 0)) for k,v in out.items()},
               min_flux={k: float(np.min(v)) for k,v in out.items()})
    atomic_json(meta, rec)
    return rec


def load_step(step):
    if float(step) == STEP_COARSE:
        p = OLDCK / f"full_nested_step_{tag(step)}.npz"
    elif float(step) == STEP_FINE:
        p = CK / f"full_nested_step_{tag(step)}.npz"
    else:
        raise ValueError(step)
    z = np.load(p)
    species = ("nut","anut","num","anum","tau_minus","tau_plus")
    return z["E"], z["incident"], {k:z[k] for k in species}


def interp(E, y, target):
    j = int(np.searchsorted(E, target))
    if j == 0 or j >= len(E): return np.nan
    i = j - 1
    if not (np.isfinite(y[i]) and np.isfinite(y[j])): return np.nan
    if y[i] > 0 and y[j] > 0:
        a=(np.log(target)-np.log(E[i]))/(np.log(E[j])-np.log(E[i]))
        return float(np.exp(np.log(y[i])+a*(np.log(y[j])-np.log(y[i]))))
    a=(target-E[i])/(E[j]-E[i])
    return float(y[i]+a*(y[j]-y[i]))


def adjacent_metrics():
    E, inc, A = load_step(STEP_COARSE)
    E2, inc2, B = load_step(STEP_FINE)
    if not np.allclose(E,E2) or not np.allclose(inc,inc2):
        raise RuntimeError("energy/input grids differ")
    rows=[]; maxima={k:0.0 for k in A}; l1={}
    for ep in TARGETS:
        ii=interp(E,inc,ep*1e6)
        for sp in A:
            a=interp(E,A[sp],ep*1e6); b=interp(E,B[sp],ep*1e6)
            den=max(abs(b),abs(ii)*1e-12,1e-300)
            s=(a-b)/den
            maxima[sp]=max(maxima[sp],abs(s))
            rows.append(dict(E_PeV=ep,species=sp,coarse=a,fine=b,
                             signed_relative_difference=s,
                             absolute_relative_difference=abs(s)))
    for sp in A:
        scale=max(float(np.max(np.abs(B[sp]))),1e-300)
        mask=(inc>1e-12*np.max(inc)) & (np.abs(B[sp])>1e-12*scale)
        l1[sp]=0.0 if not np.any(mask) else float(
            np.sum(np.abs(A[sp][mask]-B[sp][mask])) /
            max(np.sum(np.abs(B[sp][mask])),1e-300))
    return dict(coarse_step_km=STEP_COARSE,fine_step_km=STEP_FINE,
                maxima=maxima,L1=l1,pointwise=rows)


def tau_semigroup_tests():
    tn,_ = M.load_table("nu"); ta,_ = M.load_table("antinu")
    tr = M.TauHelicityFlipTransportExpm(tn,ta)
    E = np.asarray(tn["energies"],float); N=len(E)
    prof = nested_profile_ext(M.MU, STEP_COARSE)
    rhos = {"median":float(np.median(prof["rho"])),
            "max":float(np.max(prof["rho"]))}
    results=[]
    for rho_name,rho in rhos.items():
        for epev in TAU_TEST_E_PEV:
            j=int(np.argmin(np.abs(np.log(E)-np.log(epev*1e6))))
            for sign in ("minus","plus"):
                for step in TAU_TEST_STEPS_KM:
                    h=step*1e5
                    n=np.zeros(N); a=np.zeros(N); mn=np.zeros(N); ma=np.zeros(N)
                    tm=np.zeros(N); tp=np.zeros(N)
                    (tm if sign=="minus" else tp)[j]=1.0
                    one=tr._tau_half(tm.copy(),tp.copy(),n.copy(),a.copy(),mn.copy(),ma.copy(),rho,h,True)
                    first=tr._tau_half(tm.copy(),tp.copy(),n.copy(),a.copy(),mn.copy(),ma.copy(),rho,0.5*h,True)
                    two=tr._tau_half(*first,rho,0.5*h,True)
                    names=("tau_minus","tau_plus","nu_tau","antinu_tau","nu_mu","antinu_mu")
                    one_d={k:np.asarray(v) for k,v in zip(names,one)}
                    two_d={k:np.asarray(v) for k,v in zip(names,two)}
                    # L1 relative to the one-step map, plus channel-wise maxima.
                    total_num=sum(float(np.sum(np.abs(two_d[k]-one_d[k]))) for k in names)
                    total_den=max(sum(float(np.sum(np.abs(one_d[k]))) for k in names),1e-300)
                    channel={}
                    for k in names:
                        den=max(float(np.sum(np.abs(one_d[k]))),1e-300)
                        channel[k]=float(np.sum(np.abs(two_d[k]-one_d[k]))/den)
                    results.append(dict(rho_label=rho_name,rho_g_cm3=rho,
                                        E_requested_PeV=epev,E_grid_PeV=float(E[j]/1e6),
                                        tau_sign=sign,step_km=step,
                                        total_L1_semigroup_defect=total_num/total_den,
                                        channel_L1_semigroup_defect=channel))
    return results


def main():
    OUT.mkdir(parents=True,exist_ok=True); CK.mkdir(parents=True,exist_ok=True)
    print("v0.19.6: one finer full-transport grid",flush=True)
    m=compute_fine()
    print(f"FINE step={m['step_km']:.6g}km sec={m['seconds']:.1f} segments={m['n_segments']} finite={m['finite']} negatives={sum(m['negative_bins'].values())}",flush=True)

    adj=adjacent_metrics(); mx=adj["maxima"]
    prim=max(mx["nut"],mx["anut"]); sec=max(mx["num"],mx["anum"]); tau=max(mx["tau_minus"],mx["tau_plus"])
    print(f"ADJ {STEP_COARSE:.6g}->{STEP_FINE:.6g} km primary={prim:.6g} secondary={sec:.6g} tau={tau:.6g}",flush=True)

    semi=tau_semigroup_tests()
    worst=sorted(semi,key=lambda r:r["total_L1_semigroup_defect"],reverse=True)[:10]
    print("TAU-SEMIGROUP worst cases",flush=True)
    for r in worst:
        print(f"rho={r['rho_label']} E={r['E_requested_PeV']:.6g}PeV sign={r['tau_sign']} step={r['step_km']:.6g}km defect={r['total_L1_semigroup_defect']:.6g}",flush=True)

    status=dict(version="v0.19.6-diagnostic",
                target="case7_L3_cosz0.1_majorana_flip",
                fine_run=m,adjacent=adj,tau_semigroup_tests=semi,
                interpretation="Diagnostics only. If the adjacent secondary error continues ~first-order while the tau A-map has a nonzero semigroup defect, redesign the discrete tau evolution rather than brute-force production with tiny spatial steps.")
    atomic_json(OUT/"tau_semigroup_status_v0_19_6.json",status)
    print("STATUS wrote",OUT/"tau_semigroup_status_v0_19_6.json",flush=True)

if __name__=="__main__":
    main()
