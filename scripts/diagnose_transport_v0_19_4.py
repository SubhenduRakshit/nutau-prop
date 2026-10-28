#!/usr/bin/env python3
"""v0.19.4 numerical diagnostics for the Case-7 Phase-B discrepancy.

Diagnostics only; this script does not certify production.

Stage 1 isolates the neutrino interaction operator B:
  old finite-max_p update -> exact discrete exp(B h)
with tau loss/decay disabled.

Stage 2 restores the full A(h/2) exp(Bh) A(h/2) evolution and studies
spatial splitting on PREM-interface-aligned, exactly nested grids
3.0 -> 1.5 -> 0.75 -> 0.375 km.

Target: Case 7, Lambda=3 TeV, |cos z|=0.1, Majorana flip.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.chirality_transport import TauHelicityFlipTransport
from nutau_prop.chirality_transport_expm import TauHelicityFlipTransportExpm
from nutau_prop.constants import R_EARTH_KM
from nutau_prop.prem import density_g_cm3

ROOT = Path(__file__).resolve().parents[1]
T15 = ROOT / "results" / "v0.15" / "tables"
UHE = ROOT / "data" / "uhe_official"
OUT = ROOT / "results" / "v0.19.4"
CK = OUT / "checkpoints"
LAM = 3.0
MU = 0.1
MODE = "majorana_flip"
MAXPS = (0.02, 0.01, 0.005, 0.0025, 0.00125)
STEPS = (3.0, 1.5, 0.75, 0.375)
TARGETS = (1.0, 10.0, 100.0, 220.0, 300.0)
PREM_INTERFACES_KM = (1221.5, 3480.0, 5701.0, 5771.0, 5971.0,
                      6151.0, 6346.6, 6356.0, 6368.0)


def tag(x):
    return (f"{float(x):.8g}").replace(".", "p")


def atomic_json(path, obj):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2))
    tmp.replace(path)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def load_table(species):
    p = T15 / f"case7_L3_{species}.npz"
    if not p.exists():
        raise FileNotFoundError(p)
    z = np.load(p)
    return {k: z[k] for k in z.files}, p


def input_flux(E):
    q = pd.read_json(UHE / "flux.json")
    ef = q["Energy [GeV]"].to_numpy(float)
    sed = q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    f = interp1d(ef, sed, bounds_error=False, fill_value=-np.inf)(E)
    return (10.0**f / E) / 6.0


def nested_profile(mu, step_km):
    """PREM-interface-aligned grids nested at factors 1,2,4,8."""
    factors = {3.0: 1, 1.5: 2, 0.75: 4, 0.375: 8}
    step_km = float(step_km)
    if step_km not in factors:
        raise ValueError(f"unsupported nested step {step_km}")
    fac = factors[step_km]
    mu = float(mu)
    L = 2.0 * R_EARTH_KM * mu
    impact = R_EARTH_KM * np.sqrt(max(0.0, 1.0 - mu * mu))
    half = 0.5 * L
    cuts = [-half, half]
    for r in PREM_INTERFACES_KM:
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
    rho = density_g_cm3(r)
    return dict(ds_km=ds, rho=np.asarray(rho, float), r_km=r, edges_km=edges)


def ip(E, y, target):
    E = np.asarray(E, float); y = np.asarray(y, float)
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


def compare(E, incident, val, ref, species):
    rows = []
    maxima = {sp: 0.0 for sp in species}
    for epev in TARGETS:
        inc = ip(E, incident, epev*1e6)
        for sp in species:
            a = ip(E, val[sp], epev*1e6); b = ip(E, ref[sp], epev*1e6)
            den = max(abs(b), abs(inc)*1e-12, 1e-300)
            signed = (a-b)/den
            maxima[sp] = max(maxima[sp], abs(signed))
            rows.append(dict(E_PeV=epev, species=sp, value=a, reference=b,
                             signed_relative_difference=signed,
                             absolute_relative_difference=abs(signed)))
    l1 = {}
    for sp in species:
        scale = max(float(np.max(np.abs(ref[sp]))), 1e-300)
        mask = (incident > 1e-12*np.max(incident)) & (np.abs(ref[sp]) > 1e-12*scale)
        l1[sp] = 0.0 if not np.any(mask) else float(
            np.sum(np.abs(val[sp][mask]-ref[sp][mask])) /
            max(np.sum(np.abs(ref[sp][mask])), 1e-300))
    return rows, maxima, l1


def interaction_expm_reference():
    tn, pnu = load_table("nu"); ta, panu = load_table("antinu")
    E = np.asarray(tn["energies"], float); inc = input_flux(E)
    tr = TauHelicityFlipTransportExpm(tn, ta); H = tr._operator(MODE, True)
    N = len(E); x = np.concatenate([inc, inc, np.zeros(4*N)])
    prof = nested_profile(MU, 3.0)
    for ds, rho in zip(prof["ds_km"], prof["rho"]):
        x = tr._advance_linear(H, x, rho, ds*1e5)
    n, a, mn, ma, tm, tp = [x[i*N:(i+1)*N] for i in range(6)]
    return E, inc, dict(nut=n, anut=a, tau_minus=tm, tau_plus=tp), pnu, panu


def interaction_old(max_p):
    tn, _ = load_table("nu"); ta, _ = load_table("antinu")
    E = np.asarray(tn["energies"], float); inc = input_flux(E)
    tr = TauHelicityFlipTransport(tn, ta)
    n = inc.copy(); a = inc.copy(); tm = np.zeros_like(n); tp = np.zeros_like(n)
    prof = nested_profile(MU, 3.0)
    for ds, rho in zip(prof["ds_km"], prof["rho"]):
        ds_cm = ds*1e5
        maxrate = max(tr.nu.matched_step_rate(rho), tr.anu.matched_step_rate(rho))
        nsub = max(1, int(np.ceil(maxrate*ds_cm/float(max_p))))
        h = ds_cm/nsub
        for _ in range(nsub):
            n0, ccn, npn = tr.nu.step(n, rho, h, use_np=True)
            a0, cca, npa = tr.anu.step(a, rho, h, use_np=True)
            n = n0 + npa; a = a0 + npn
            tm = tm + ccn @ tr.nu.Kcc; tp = tp + cca @ tr.anu.Kcc
    return E, inc, dict(nut=n, anut=a, tau_minus=tm, tau_plus=tp)


def full_expm_profile(step_km):
    tn, _ = load_table("nu"); ta, _ = load_table("antinu")
    E = np.asarray(tn["energies"], float); inc = input_flux(E)
    tr = TauHelicityFlipTransportExpm(tn, ta); H = tr._operator(MODE, True)
    N = len(E)
    n = inc.copy(); a = inc.copy(); mn = np.zeros_like(n); ma = np.zeros_like(n)
    tm = np.zeros_like(n); tp = np.zeros_like(n)
    prof = nested_profile(MU, step_km)
    for ds, rho in zip(prof["ds_km"], prof["rho"]):
        h = ds*1e5; hh = 0.5*h
        tm, tp, n, a, mn, ma = tr._tau_half(tm, tp, n, a, mn, ma, rho, hh, True)
        x = np.concatenate([n, a, mn, ma, tm, tp])
        x = tr._advance_linear(H, x, rho, h)
        n, a, mn, ma, tm, tp = [x[i*N:(i+1)*N] for i in range(6)]
        tm, tp, n, a, mn, ma = tr._tau_half(tm, tp, n, a, mn, ma, rho, hh, True)
    return E, inc, dict(nut=n, anut=a, num=mn, anum=ma,
                        tau_minus=tm, tau_plus=tp), prof


def interaction_worker(max_p):
    cp = CK / f"interaction_old_maxp_{tag(max_p)}.json"
    if cp.exists():
        return json.loads(cp.read_text())
    z = np.load(OUT / "interaction_expm_reference.npz")
    E = z["E"]; inc = z["incident"]
    ref = {k: z[k] for k in ("nut", "anut", "tau_minus", "tau_plus")}
    t0 = time.perf_counter(); E2, inc2, val = interaction_old(max_p)
    if not np.allclose(E, E2) or not np.allclose(inc, inc2):
        raise RuntimeError("interaction diagnostic grids differ")
    rows, maxima, l1 = compare(E, inc, val, ref,
                                ("nut", "anut", "tau_minus", "tau_plus"))
    out = dict(max_p=float(max_p), seconds=time.perf_counter()-t0,
               maxima=maxima, L1=l1, pointwise=rows)
    atomic_json(cp, out)
    return out


def full_worker(step):
    cp = CK / f"full_nested_step_{tag(step)}.npz"
    meta = CK / f"full_nested_step_{tag(step)}.json"
    if cp.exists() and meta.exists():
        return json.loads(meta.read_text())
    t0 = time.perf_counter(); E, inc, out, prof = full_expm_profile(step)
    np.savez_compressed(cp, E=E, incident=inc, **out)
    m = dict(step_km=float(step), seconds=time.perf_counter()-t0,
             n_segments=int(len(prof["ds_km"])),
             min_flux={k: float(np.min(v)) for k, v in out.items()},
             negative_bins={k: int(np.sum(v < 0)) for k, v in out.items()},
             finite=bool(all(np.all(np.isfinite(v)) for v in out.values())))
    atomic_json(meta, m)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nproc", type=int, default=int(os.environ.get("NUTAU_NPROC", "5")))
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); CK.mkdir(parents=True, exist_ok=True)

    print("v0.19.4 Stage 1: interaction-only old max_p -> exp(Bh)", flush=True)
    refnpz = OUT / "interaction_expm_reference.npz"
    if not refnpz.exists():
        E, inc, ref, pnu, panu = interaction_expm_reference()
        np.savez_compressed(refnpz, E=E, incident=inc, **ref)
    else:
        _, pnu = load_table("nu"); _, panu = load_table("antinu")

    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=min(a.nproc, len(MAXPS))) as pool:
        rr = list(pool.imap_unordered(interaction_worker, MAXPS, chunksize=1))
    rr = sorted(rr, key=lambda r: r["max_p"], reverse=True)
    for r in rr:
        prim = max(r["maxima"]["nut"], r["maxima"]["anut"])
        tau = max(r["maxima"]["tau_minus"], r["maxima"]["tau_plus"])
        print(f"INTERACTION max_p={r['max_p']:.5g} sec={r['seconds']:.1f} "
              f"primary={prim:.6g} tau_accumulator={tau:.6g}", flush=True)

    print("v0.19.4 Stage 2: full A-B-A on nested PREM grids", flush=True)
    with ctx.Pool(processes=min(a.nproc, len(STEPS))) as pool:
        smeta = list(pool.imap_unordered(full_worker, STEPS, chunksize=1))
    smeta = sorted(smeta, key=lambda x: x["step_km"], reverse=True)
    refz = np.load(CK / f"full_nested_step_{tag(0.375)}.npz")
    E = refz["E"]; inc = refz["incident"]
    ref = {k: refz[k] for k in ("nut", "anut", "num", "anum")}
    split = []
    for step in (3.0, 1.5, 0.75):
        z = np.load(CK / f"full_nested_step_{tag(step)}.npz")
        val = {k: z[k] for k in ("nut", "anut", "num", "anum")}
        rows, maxima, l1 = compare(E, inc, val, ref, ("nut", "anut", "num", "anum"))
        rec = dict(step_km=step, reference_step_km=0.375,
                   maxima=maxima, L1=l1, pointwise=rows)
        split.append(rec)
        prim = max(maxima["nut"], maxima["anut"])
        sec = max(maxima["num"], maxima["anum"])
        print(f"SPLIT step={step:g}km ->0.375km primary={prim:.6g} secondary={sec:.6g}", flush=True)

    provenance = dict(case="Case7 Lambda=3 TeV cosz=0.1 majorana_flip",
                      table_nu=str(pnu.relative_to(ROOT)),
                      table_antinu=str(panu.relative_to(ROOT)),
                      table_nu_sha256=sha256(pnu), table_antinu_sha256=sha256(panu),
                      flux_sha256=sha256(UHE / "flux.json"),
                      diagnostic_source_sha256=sha256(Path(__file__)))
    status = dict(version="v0.19.4-diagnostic", certification=False,
                  interaction_scan=rr, nested_split_comparisons=split,
                  full_run_metadata=smeta, provenance=provenance,
                  interpretation_note=(
                      "No automatic blame assignment. Stage 1 isolates finite-max_p interaction error. "
                      "Stage 2 measures full tau-splitting/PREM discretization on nested interface-aligned grids. "
                      "If secondary differences plateau, run the optional refined-energy-grid stage."))
    atomic_json(OUT / "diagnostic_status_v0_19_4.json", status)
    pd.DataFrame([dict(max_p=r["max_p"], seconds=r["seconds"],
                       max_primary=max(r["maxima"]["nut"], r["maxima"]["anut"]),
                       max_tau_accumulator=max(r["maxima"]["tau_minus"], r["maxima"]["tau_plus"]),
                       L1_nut=r["L1"]["nut"], L1_anut=r["L1"]["anut"],
                       L1_tau_minus=r["L1"]["tau_minus"], L1_tau_plus=r["L1"]["tau_plus"])
                  for r in rr]).to_csv(OUT / "interaction_only_convergence_v0_19_4.csv", index=False)
    print("STATUS wrote", OUT / "diagnostic_status_v0_19_4.json", flush=True)


if __name__ == "__main__":
    main()
