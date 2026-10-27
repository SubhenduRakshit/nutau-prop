#!/usr/bin/env python3
"""v0.9 true-energy audit against IceCube's released differential-limit curve.

This is a public-data shape screen, not the unpublished reconstructed-bin
likelihood. The BSM effective-area rescalings are response upper proxies.
"""

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "uhe_official"
RESULTS = ROOT / "results" / "v0.9"
YEAR_SECONDS = 365.25 * 86400.0
ICECUBE_EXPOSURE = 12.6 * YEAR_SECONDS
KM3_EXPOSURE = 335.0 * 86400.0
KM3_BENCHMARK = 0.242767
KM3_CENTRAL_BENCHMARK = 0.211117
ENERGY = np.logspace(5, 11, 6001)


def response_areas():
    values = np.loadtxt(DATA / "icecube_ehe_effective_area.csv", delimiter=",", skiprows=1)
    energy = values[:, 0]
    areas = values[:, 1:] * 1e4
    if not np.allclose(areas[:, 0], areas[:, 1:].sum(axis=1), rtol=0.01):
        raise ValueError("IceCube total area does not match flavor columns")
    return energy, areas


def curve(energy, values):
    return interp1d(energy, values, bounds_error=False, fill_value=0.0)(ENERGY)


def incident_flux():
    table = pd.read_json(DATA / "flux.json")
    energy = table["Energy [GeV]"].to_numpy(float)
    sed = table["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    interpolated = interp1d(energy, sed, bounds_error=False, fill_value=-np.inf)(ENERGY)
    return 10.0**interpolated / ENERGY**2


def case6_ratio(table, cutoff, column):
    selected = table[np.isclose(table["Lambda_TeV"], cutoff)].sort_values("E_PeV")
    energy = selected["E_PeV"].to_numpy(float) * 1e6
    ratio = selected[column].to_numpy(float)
    output = np.zeros_like(ENERGY)
    valid = (ENERGY >= energy.min()) & (ENERGY <= energy.max())
    output[valid] = 10.0**np.interp(
        np.log10(ENERGY[valid]), np.log10(energy), np.log10(ratio)
    )
    return output


def counts(area, flux, exposure, low=None, high=None):
    selected = np.ones(ENERGY.size, dtype=bool)
    if low is not None:
        selected &= ENERGY >= low
    if high is not None:
        selected &= ENERGY <= high
    return 4.0 * math.pi * exposure * np.trapezoid(area[selected] * flux[selected], ENERGY[selected])


def published_windows(total_area, flux, all_response, tau_response):
    limits = np.loadtxt(DATA / "icecube_ehe_differential_limit.csv", delimiter=",", skiprows=1)
    window_rows = []
    for center, observed_limit, sensitivity in limits:
        lower, upper = center / math.sqrt(10.0), center * math.sqrt(10.0)
        if upper > 9.46e10:
            continue
        reference_flux_per_flavor = (observed_limit / (3.0 * center**2)) * (ENERGY / center) ** -1.0
        reference = counts(total_area, reference_flux_per_flavor, ICECUBE_EXPOSURE, lower, upper)
        baseline = counts(total_area, flux, ICECUBE_EXPOSURE, lower, upper)
        all_count = counts(all_response, flux, ICECUBE_EXPOSURE, lower, upper)
        tau_count = counts(tau_response, flux, ICECUBE_EXPOSURE, lower, upper)
        window_rows.append(dict(
            center_GeV=center, lower_GeV=lower, upper_GeV=upper,
            published_limit_E2Phi_all=observed_limit,
            published_sensitivity_E2Phi_all=sensitivity,
            reference_limit_yield=reference,
            Ehlert_SM_yield=baseline,
            Ehlert_allflavor_NP_proxy_yield=all_count,
            Ehlert_tauinitial_NP_proxy_yield=tau_count,
            allflavor_window_ratio=all_count / reference,
            tauinitial_window_ratio=tau_count / reference,
        ))
    return window_rows


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    energy_icecube, areas = response_areas()
    area_total = curve(energy_icecube, areas[:, 0])
    area_tau = curve(energy_icecube, areas[:, 3])
    flux = incident_flux()
    km3 = pd.read_json(DATA / "km3.json")
    area_km3 = curve(km3["Energy [GeV]"].to_numpy(float), km3["Aeff [cm^2]"].to_numpy(float))
    km3_sm = counts(area_km3, flux, KM3_EXPOSURE)
    km3_central = counts(area_km3, flux, KM3_EXPOSURE, 7.24e7, 2.57e9)
    icecube_sm = counts(area_total, flux, ICECUBE_EXPOSURE)
    if abs(km3_sm / KM3_BENCHMARK - 1.0) > 0.005:
        raise ValueError("KM3NeT full-range benchmark failed")
    if abs(km3_central / KM3_CENTRAL_BENCHMARK - 1.0) > 0.005:
        raise ValueError("KM3NeT central-range benchmark failed")
    if abs(icecube_sm / 17.674595 - 1.0) > 0.005:
        raise ValueError("IceCube m^2 to cm^2 regression failed")

    cross_sections = pd.read_csv(ROOT / "results" / "v0.8" / "case6_pdf_cut_v0_8.csv")
    summaries = []
    windows = []
    for cutoff in (1.5, 2.0, 2.44, 3.0):
        ratio = case6_ratio(cross_sections, cutoff, "r_shatcut")
        all_response = area_total * (1.0 + ratio)
        tau_response = area_total + area_tau * ratio
        cutoff_windows = published_windows(area_total, flux, all_response, tau_response)
        for row in cutoff_windows:
            row["Lambda_TeV"] = cutoff
        windows.extend(cutoff_windows)
        all_scale = min(1.0 / row["allflavor_window_ratio"] for row in cutoff_windows)
        tau_scale = min(1.0 / row["tauinitial_window_ratio"] for row in cutoff_windows)
        all_icecube = counts(all_response, flux, ICECUBE_EXPOSURE)
        tau_icecube = counts(tau_response, flux, ICECUBE_EXPOSURE)
        summaries.append(dict(
            Lambda_TeV=cutoff,
            c4_over_4pi=float(cross_sections[np.isclose(cross_sections.Lambda_TeV, cutoff)].iloc[0].c4_over_4pi),
            IceCube_SM=icecube_sm,
            IceCube_allflavor_NP_response_proxy=all_icecube,
            IceCube_tauinitial_NP_response_proxy=tau_icecube,
            KM3_brighttrack_SM=km3_sm,
            KM3_brighttrack_direct_Case6=km3_sm,
            KM3_NC_equivalent_proxy=counts(area_km3 * (1.0 + ratio), flux, KM3_EXPOSURE),
            public_window_allflavor_scale_proxy=all_scale,
            public_window_tauinitial_scale_proxy=tau_scale,
            KM3_direct_track_at_allflavor_scale_proxy=km3_sm * all_scale,
            KM3_direct_track_at_tauinitial_scale_proxy=km3_sm * tau_scale,
        ))

    pd.DataFrame(windows).to_csv(RESULTS / "icecube_public_windows_v0_9.csv", index=False)
    summary = pd.DataFrame(summaries)
    summary.to_csv(RESULTS / "case6_public_audit_v0_9.csv", index=False)
    print(f"KM3NeT SM full={km3_sm:.6f} (published {KM3_BENCHMARK:.6f}); central90={km3_central:.6f} (published {KM3_CENTRAL_BENCHMARK:.6f})")
    print(f"IceCube EHE SM response fold={icecube_sm:.6f}; release area is m^2 and the flux is per flavor")
    print(summary.to_string(index=False, float_format=lambda value: f"{value:.5g}"))
    print("PUBLIC-DATA SCREEN ONLY: differential limits use sliding E^-1 windows; no event-level EHE likelihood or BSM response matrix was released.")
    print("NC-equivalent KM3 values are not bright-track predictions; no direct Case-6 track increment was assigned.")


if __name__ == "__main__":
    main()
