#!/usr/bin/env python3
"""v0.13 end-to-end Case-6 propagation + detector fold.

Purpose
-------
Close the remaining Case-6 consistency gap by combining the v0.12 propagated
nu_tau / secondary-nu_mu spectra with the same Case-6 detector vertex used in
v0.11.  Where public responses expose topology/zenith (ARCA), the fold is
channel resolved.  Where the public response is sky-averaged and not CC/NC
resolved (IceCube EHE and KM3-230213A bright-track response), the script reports
an explicit public-response envelope instead of pretending to perform an
unreleased collaboration likelihood.

Important conventions
---------------------
* Case 6 is tau-specific and NC-like.
* SM CC is unchanged.
* Direct Case-6 muon-track production is zero.
* BSM-induced changes in tau-decay secondary nu_mu are included as DELTAS
  relative to the published SM response, avoiding double counting of any
  standard tau-induced secondary contribution already present in detector MC.
* The central sky-averaged IceCube/KM3 treatment uses equal up/down solid-angle
  weight.  The low/high IceCube envelopes span the available propagated zenith
  factors and CC/NC detector-vertex ambiguity.
* Exact changed-y selection efficiency is unavailable publicly; two ARCA
  predictions are therefore written: propagation-only and propagation+vertex.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.arca import fig7_response, calibrate_fig7_response, integrate_loggrid
from nutau_prop.config import load_config
from nutau_prop.transport import NuOnlyTransport

ROOT = Path(__file__).resolve().parents[1]
V11 = ROOT / "results" / "v0.11"
V12 = ROOT / "results" / "v0.12"
OUT = ROOT / "results" / "v0.13"
TABLES = V11 / "tables"
UHE = ROOT / "data" / "uhe_official"
YEAR = 365.25 * 86400.0
LAMBDAS = (1.5, 2.0, 2.44, 3.0)


def tag(lam):
    return str(lam).replace(".", "p")


def ehlert_flux(E):
    q = pd.read_json(UHE / "flux.json")
    ef = q["Energy [GeV]"].to_numpy(float)
    sed = q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    s = interp1d(ef, sed, bounds_error=False, fill_value=-np.inf)(E)
    return 10.0**s / E**2


def load_table(lam, species):
    z = np.load(TABLES / f"case6_L{tag(lam)}_{species}.npz")
    return {k: z[k] for k in z.files}


def loginterp(E, grid, values, outside):
    E = np.asarray(E, float); grid = np.asarray(grid, float); values = np.asarray(values, float)
    out = np.full(E.shape, float(outside), dtype=float)
    m = (E >= grid.min()) & (E <= grid.max())
    if np.any(m):
        out[m] = 10.0 ** np.interp(
            np.log10(E[m]), np.log10(grid), np.log10(np.maximum(values, 1e-300))
        )
    return out


def detector_ratios(lam, E):
    """nu+antinu averaged detector NC and inclusive cross-section ratios."""
    vals_nc, vals_inc = [], []
    for species in ("nu", "antinu"):
        t = load_table(lam, species)
        scc = loginterp(E, t["energies"], t["sigma_cc"], 0.0)
        snc0 = loginterp(E, t["energies"], t["sigma_nc_sm"], 0.0)
        snp = loginterp(E, t["energies"], t["sigma_np"], 0.0)
        inside = (np.asarray(E) >= t["energies"].min()) & (np.asarray(E) <= t["energies"].max())
        rnc = np.ones_like(np.asarray(E, float)); rinc = np.ones_like(rnc)
        rnc[inside] = (snc0[inside] + snp[inside]) / np.maximum(snc0[inside], 1e-300)
        rinc[inside] = (scc[inside] + snc0[inside] + snp[inside]) / np.maximum(scc[inside] + snc0[inside], 1e-300)
        vals_nc.append(rnc); vals_inc.append(rinc)
    return 0.5 * (vals_nc[0] + vals_nc[1]), 0.5 * (vals_inc[0] + vals_inc[1])


def transport_data():
    p = V12 / "case6_earth_transport_v0_12.csv"
    if not p.exists():
        raise SystemExit(f"Missing {p}; run v0.12 first")
    d = pd.read_csv(p)
    d = d[d.spectrum == "Ehlert"].copy()
    d["sec_sm_frac"] = d.SM_secondary_numu / np.maximum(d.incident_per_logE, 1e-300)
    d["sec_bsm_frac"] = d.BSM_secondary_numu / np.maximum(d.incident_per_logE, 1e-300)
    return d


def interp_transport(d, lam, E, mu, field, horizon_value):
    q = d[np.isclose(d.Lambda_TeV, lam)]
    mus = np.sort(q.cosz_abs.unique())
    curves = []
    for m in mus:
        z = q[np.isclose(q.cosz_abs, m)].sort_values("E_GeV")
        curves.append(loginterp(E, z.E_GeV.to_numpy(float), z[field].to_numpy(float), horizon_value))
    mus2 = np.r_[0.0, mus]
    vals = np.vstack([np.full(len(E), horizon_value), *curves])
    out = np.empty(len(E), float)
    for i in range(len(E)):
        out[i] = np.interp(mu, mus2, vals[:, i])
    return out


def upgoing_average(d, lam, E, field, horizon_value, nmu=24):
    mus = (np.arange(nmu) + 0.5) / nmu
    return np.mean([interp_transport(d, lam, E, float(mu), field, horizon_value) for mu in mus], axis=0)


def upgoing_extrema(d, lam, E, field, horizon_value):
    q = d[np.isclose(d.Lambda_TeV, lam)]
    mus = np.r_[0.0, np.sort(q.cosz_abs.unique())]
    arr = np.vstack([interp_transport(d, lam, E, float(mu), field, horizon_value) for mu in mus])
    return np.min(arr, axis=0), np.max(arr, axis=0)


def normalized_ehlert_per_logE(E):
    p = E * ehlert_flux(E)
    m = np.nanmax(p)
    return p / m if m > 0 else p


def sm_view(t):
    return dict(
        energies=t["energies"], y=t["y"], dy=t["dy"],
        sigma_cc=t["sigma_cc"], sigma_nc=t["sigma_nc_sm"],
        p_cc=t["p_cc"], p_nc=t["p_nc_sm"],
    )


def primary_numu_transmission(Eout, mu, lam_for_sm=1.5):
    """Spectrum-specific SM numu+nubar transmission, used only to remove the
    Earth-propagation factor already embedded in published upgoing track/NC Aeff
    when adding a neutrino created inside Earth by tau decay.
    """
    vals = []
    for species in ("nu", "antinu"):
        t = load_table(lam_for_sm, species)
        p0 = normalized_ehlert_per_logE(t["energies"])
        tr = NuOnlyTransport(sm_view(t))
        f = tr.propagate(mu, p0, max_step_km=12.0, max_p=0.02)["nu"]
        vals.append((t["energies"], f, p0))
    grid = vals[0][0]
    comb = vals[0][1] + vals[1][1]
    inc = vals[0][2] + vals[1][2]
    trans = comb / np.maximum(inc, 1e-300)
    return loginterp(Eout, grid, trans, 1.0)


def sec_detector_delta(d, lam, E, mu):
    """Delta secondary numu flux at detector divided by the SM primary-numu
    flux at detector.  Values where primary transmission is numerically tiny
    are set conservatively to zero; their published Aeff is tiny as well.
    """
    sb = interp_transport(d, lam, E, mu, "sec_bsm_frac", 0.0)
    ss = interp_transport(d, lam, E, mu, "sec_sm_frac", 0.0)
    tm = primary_numu_transmission(E, mu)
    out = np.zeros_like(E)
    good = tm > 1e-8
    out[good] = (sb[good] - ss[good]) / tm[good]
    return out


def arca_fold(d, lam, cfg, include_vertex=True):
    E = np.logspace(5, 8, 1501)
    phi = ehlert_flux(E)
    R = fig7_response(E); sh = R["shower"]
    sf_sh, sf_tr = calibrate_fig7_response(cfg)
    expo = float(cfg["arca"]["exposure_years"]) * YEAR
    rnc, _ = detector_ratios(lam, E)
    if not include_vertex:
        rnc = np.ones_like(rnc)

    rtau = upgoing_average(d, lam, E, "nutau_ratio", 1.0)
    dsec_nc = np.mean([
        sec_detector_delta(d, lam, E, float(mu))
        for mu in (np.arange(24) + 0.5) / 24
    ], axis=0)

    # Published SM response is the baseline.  We alter only the tau contribution
    # and add the BSM-SM secondary-muon delta, avoiding double counting.
    A_sm_down = (sh["down_ecc"] + sh["down_tcc"] + sh["down_nc"]) * sf_sh
    A_bs_down = (sh["down_ecc"] + sh["down_tcc"] + sh["down_nc"] * (2.0 + rnc) / 3.0) * sf_sh

    A_sm_up = (sh["up_ecc"] + sh["up_tcc"] + sh["up_nc"]) * sf_sh
    A_bs_up = (
        sh["up_ecc"]
        + sh["up_tcc"] * rtau
        + sh["up_nc"] * (2.0 + rtau * rnc) / 3.0
        + sh["up_nc"] * dsec_nc / 3.0
    ) * sf_sh

    hemi = expo * 2.0 * math.pi * 1e4
    nsm_sh = hemi * (integrate_loggrid(E, phi * A_sm_down) + integrate_loggrid(E, phi * A_sm_up))
    nbs_sh = hemi * (integrate_loggrid(E, phi * A_bs_down) + integrate_loggrid(E, phi * A_bs_up))

    nsm_tr = 0.0; nbs_tr = 0.0
    track_rows = []
    for key, (clo, chi) in R["track_cos_bins"].items():
        mu = abs(0.5 * (clo + chi))
        dsec = sec_detector_delta(d, lam, E, mu)
        A = R["track"][key] * sf_tr
        domega = 2.0 * math.pi * (chi - clo)
        sm = expo * domega * integrate_loggrid(E, phi * A) * 1e4
        bs = expo * domega * integrate_loggrid(E, phi * A * (1.0 + dsec)) * 1e4
        nsm_tr += sm; nbs_tr += bs
        track_rows.append(dict(Lambda_TeV=lam, mode=("joint" if include_vertex else "propagation_only"),
                               track_bin=key, cos_center=0.5*(clo+chi), N_SM=sm, N_BSM=bs,
                               fractional_shift=bs/sm-1.0 if sm>0 else np.nan))

    return dict(
        Lambda_TeV=lam, mode=("joint" if include_vertex else "propagation_only"),
        N_shower_SM=nsm_sh, N_shower_BSM=nbs_sh,
        shower_fractional_shift=nbs_sh/nsm_sh-1.0,
        N_track_SM=nsm_tr, N_track_BSM=nbs_tr,
        track_fractional_shift=nbs_tr/nsm_tr-1.0,
        total_fractional_shift=(nbs_sh+nbs_tr)/(nsm_sh+nsm_tr)-1.0,
    ), track_rows


def response_areas_icecube():
    values = np.loadtxt(UHE / "icecube_ehe_effective_area.csv", delimiter=",", skiprows=1)
    E = values[:, 0]
    A = values[:, 1:] * 1e4
    if not np.allclose(A[:,0], A[:,1:].sum(axis=1), rtol=0.01):
        raise ValueError("IceCube total area does not match flavor columns")
    return E, A


def curve(E, xp, fp):
    return interp1d(xp, fp, bounds_error=False, fill_value=0.0)(E)


def counts(E, area, flux, exposure, low=None, high=None):
    m = np.ones(len(E), bool)
    if low is not None: m &= E >= low
    if high is not None: m &= E <= high
    return 4.0 * math.pi * exposure * np.trapezoid(area[m] * flux[m], E[m])


def icecube_public_scale(E, area_sm, area_bsm, phi):
    limits = np.loadtxt(UHE / "icecube_ehe_differential_limit.csv", delimiter=",", skiprows=1)
    rows = []
    emax = E.max()
    for center, observed_limit, sensitivity in limits:
        low, high = center / math.sqrt(10.0), center * math.sqrt(10.0)
        if high > emax: continue
        ref_flux = (observed_limit / (3.0 * center**2)) * (E / center) ** -1.0
        ref = counts(E, area_sm, ref_flux, 12.6*YEAR, low, high)
        val = counts(E, area_bsm, phi, 12.6*YEAR, low, high)
        rows.append(dict(center_GeV=center, lower_GeV=low, upper_GeV=high,
                         reference_limit_yield=ref, model_yield=val,
                         ratio_to_public_limit=val/ref if ref>0 else np.nan))
    finite = [r["ratio_to_public_limit"] for r in rows if np.isfinite(r["ratio_to_public_limit"]) and r["ratio_to_public_limit"]>0]
    scale = 1.0/max(finite) if finite else np.inf
    return scale, rows


def icecube_fold(d, lam):
    Eic, Araw = response_areas_icecube()
    table_max = min(load_table(lam,"nu")["energies"].max(), load_table(lam,"antinu")["energies"].max())
    E = np.logspace(5, np.log10(table_max), 5001)
    phi = ehlert_flux(E)
    Atot = curve(E, Eic, Araw[:,0]); Ae = curve(E,Eic,Araw[:,1]); Amu = curve(E,Eic,Araw[:,2]); Atau = curve(E,Eic,Araw[:,3])
    rnc, rinc = detector_ratios(lam, E)
    ravg = upgoing_average(d, lam, E, "nutau_ratio", 1.0)
    rmin, rmax = upgoing_extrema(d, lam, E, "nutau_ratio", 1.0)
    rsky = 0.5 * (1.0 + ravg)
    rlo = np.minimum(1.0, rmin)
    rhi = np.maximum(1.0, rmax)

    dsec_sky = 0.5 * np.mean([
        sec_detector_delta(d, lam, E, float(mu))
        for mu in (np.arange(24)+0.5)/24
    ], axis=0)

    # Central: 50/50 up/down solid angle, inclusive tau vertex factor, and the
    # BSM-SM secondary-muon delta.  Low/high span public zenith and CC/NC ambiguity.
    Acentral = Atot + Atau * (rsky * rinc - 1.0) + Amu * dsec_sky
    Alow = Atot + Atau * (rlo - 1.0) + np.minimum(Amu*dsec_sky, 0.0)
    Ahigh = Atot + Atau * (rhi * rnc - 1.0) + np.maximum(Amu*dsec_sky, 0.0)
    Aprop = Atot + Atau * (rsky - 1.0) + Amu * dsec_sky
    for arr in (Acentral, Alow, Ahigh, Aprop):
        arr[arr < 0] = 0.0

    expo = 12.6 * YEAR
    nsm = counts(E, Atot, phi, expo)
    out = {"Lambda_TeV":lam, "transport_grid_max_GeV":table_max,
           "IceCube_SM_within_grid":nsm,
           "IceCube_joint_central":counts(E,Acentral,phi,expo),
           "IceCube_joint_low":counts(E,Alow,phi,expo),
           "IceCube_joint_high":counts(E,Ahigh,phi,expo),
           "IceCube_propagation_only":counts(E,Aprop,phi,expo)}
    scales={}; windows=[]
    for name,area in (("central",Acentral),("low",Alow),("high",Ahigh),("propagation_only",Aprop)):
        sc, rows = icecube_public_scale(E, Atot, area, phi)
        scales[name]=sc
        for r in rows:
            r.update(Lambda_TeV=lam, response=name)
        windows.extend(rows)
        out[f"public_window_flux_scale_{name}"]=sc
    return out, windows, (E,phi,dsec_sky)


def km3_fold(lam, icecube_row, aux):
    E, phi, dsec_sky = aux
    km = pd.read_json(UHE / "km3.json")
    Akm = curve(E, km["Energy [GeV]"].to_numpy(float), km["Aeff [cm^2]"].to_numpy(float))
    T = 335.0 * 86400.0
    nsm = counts(E, Akm, phi, T)
    nb = counts(E, Akm * np.maximum(0.0, 1.0 + dsec_sky), phi, T)
    # Direct Case-6 track increment is zero; this is only the BSM-SM secondary-muon change.
    d220 = float(10.0**np.interp(np.log10(2.2e8), np.log10(E), np.log10(np.maximum(1.0+dsec_sky,1e-300)))) if E.min()<=2.2e8<=E.max() else np.nan
    out=dict(Lambda_TeV=lam, KM3_brighttrack_SM_within_grid=nsm,
             KM3_secondary_track_central=nb,
             KM3_secondary_track_fractional_shift=nb/nsm-1.0 if nsm>0 else np.nan,
             KM3_track_response_ratio_at_220PeV=d220)
    for key in ("central","low","high","propagation_only"):
        sc=icecube_row[f"public_window_flux_scale_{key}"]
        out[f"KM3_at_IceCube_scale_{key}"]=nb*sc
    return out


def omitted_icecube_sm_fraction():
    Eic,Araw=response_areas_icecube(); E=np.logspace(5,11,6001); phi=ehlert_flux(E)
    Atot=curve(E,Eic,Araw[:,0]); full=counts(E,Atot,phi,12.6*YEAR)
    maxe=load_table(1.5,"nu")["energies"].max(); inside=counts(E,Atot,phi,12.6*YEAR,high=maxe)
    return 1.0-inside/full, full, inside


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--arca-config",default="config/arca.yaml")
    a=ap.parse_args(); OUT.mkdir(parents=True,exist_ok=True)
    d=transport_data(); cfg=load_config(a.arca_config)

    omitted, icfull, icinside = omitted_icecube_sm_fraction()
    arca=[]; trackbins=[]; ice=[]; windows=[]; km3=[]
    for lam in LAMBDAS:
        for vertex in (False,True):
            row,tr=arca_fold(d,lam,cfg,include_vertex=vertex); arca.append(row); trackbins.extend(tr)
        ir,wr,aux=icecube_fold(d,lam); ice.append(ir); windows.extend(wr)
        km3.append(km3_fold(lam,ir,aux))

    A=pd.DataFrame(arca); T=pd.DataFrame(trackbins); I=pd.DataFrame(ice); W=pd.DataFrame(windows); K=pd.DataFrame(km3)
    A.to_csv(OUT/"case6_arca_joint_v0_13.csv",index=False)
    T.to_csv(OUT/"case6_arca_secondary_track_bins_v0_13.csv",index=False)
    I.to_csv(OUT/"case6_icecube_joint_v0_13.csv",index=False)
    W.to_csv(OUT/"case6_icecube_public_windows_joint_v0_13.csv",index=False)
    K.to_csv(OUT/"case6_km3_joint_v0_13.csv",index=False)

    gate={
      "version":"v0.13",
      "purpose":"end-to-end Case-6 propagation + detector fold",
      "historical_validation":"v0.12 ordered runner completed, therefore RR/DRS and v0.11 differential gates passed",
      "ARCA":"channel/topology resolved using public Fig.7 response; propagation-only and propagation+vertex both reported",
      "IceCube":"sky-averaged public-response envelope; not collaboration binned likelihood",
      "KM3NeT":"direct Case-6 track increment zero; only BSM-SM secondary-numu delta included",
      "IceCube_SM_full_Ehlert":icfull,
      "IceCube_SM_within_transport_grid":icinside,
      "IceCube_SM_fraction_above_transport_grid":omitted,
      "grid_gate_PASS":bool(omitted<0.05),
      "public_y_response_limit":"published Aeff integrates over SM y; propagation-only vs propagation+vertex is retained as a theory-response bracket",
    }
    (OUT/"case6_joint_status_v0_13.json").write_text(json.dumps(gate,indent=2))

    print("\nARCA end-to-end Case-6 fold")
    print(A.to_string(index=False,float_format=lambda x:f"{x:.5g}"))
    print("\nIceCube EHE end-to-end public-response envelope")
    print(I.to_string(index=False,float_format=lambda x:f"{x:.5g}"))
    print("\nKM3NeT bright-track response after secondary-numu delta and IceCube scaling")
    print(K.to_string(index=False,float_format=lambda x:f"{x:.5g}"))
    print("\nGRID QA:",json.dumps(gate,indent=2))
    print("WROTE",OUT)


if __name__=="__main__":
    main()
