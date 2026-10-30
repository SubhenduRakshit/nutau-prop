#!/usr/bin/env python3
"""v0.43 IceCube HESE public-data recast of the PKR absorption-horizon signal.

Uses the official IceCube HESE 7.5-year public event/MC release and its published
analysis binning.  Detector-systematic splines are intentionally held at nominal
values so this script has no photospline dependency; the collaboration's main
flux/background nuisance parameters are profiled.  MC statistical uncertainty is
handled with the same L_eff construction used by the public HESE likelihood.

BSM input is the converged v0.42d coupled-transport map.  The tau-flavor
astrophysical MC is reweighted by the propagated BSM/SM tau ratio.  The muon-
flavor astrophysical MC is reweighted by the change in TOTAL muon flavor flux,
primary nu_mu + tau-decay secondary nu_mu, so tracks are not assumed BSM blind.

This is a public-data recast, not the collaboration likelihood: detector nuisance
splines are fixed and the v0.42d tau ratio is nu+nubar averaged.  These limitations
are reported in the output metadata.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize
from scipy.special import gammaln

from nutau_prop.transport import NuOnlyTransport

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.43"
OUT.mkdir(parents=True, exist_ok=True)
HESE_ROOT = ROOT / "cache" / "HESE-7-year-data-release" / "HESE-7-year-data-release"
V42 = ROOT / "results" / "v0.42d"
LIVETIME = 227708167.68
PIVOT_ASTRO = 1.0e5
PIVOT_CR = 2020.0
PIVOT_PROMPT = 7887.0
GAMMA_MAP = 2.5

MODELS = ["eq4p3_dipole", "eq4p4_derivative"]
TOPOLOGIES = ["active_redistribution_proxy", "sterile_sink"]

# Official HESE nuisance set minus detector-systematic spline parameters.
PARAM_NAMES = [
    "astro_norm", "astro_gamma", "conv_norm", "prompt_norm", "muon_norm",
    "cr_delta_gamma", "kpi_ratio", "nunubar_ratio",
]
X0 = np.array([6.0, 2.87, 1.0, 1.0, 1.0, -0.05, 1.0, 1.0], float)
BOUNDS = [
    (0.0, 20.0), (1.5, 4.5), (0.0, 4.0), (0.0, 10.0),
    (0.0, 4.0), (-0.30, 0.20), (0.4, 1.6), (0.0, 2.0),
]
PRIORS = {
    "conv_norm": (1.0, 0.4),
    "muon_norm": (1.0, 0.5),
    "cr_delta_gamma": (-0.05, 0.05),
    "kpi_ratio": (1.0, 0.1),
    "nunubar_ratio": (1.0, 0.1),
}


def require_inputs() -> None:
    need = [
        V42 / "coupled_true_energy_map.csv",
        V42 / "benchmark_gate.csv",
        V42 / "tables" / "eq4p3_dipole_nu.npz",
        V42 / "tables" / "eq4p3_dipole_antinu.npz",
    ]
    missing = [str(p) for p in need if not p.exists()]
    if missing:
        raise RuntimeError("Missing v0.42d inputs:\n" + "\n".join(missing))
    if not HESE_ROOT.exists():
        raise RuntimeError(
            f"Official HESE release not found at {HESE_ROOT}. "
            "Run run_icecube_hese_public_likelihood_v0_43.sh."
        )


def import_hese_modules():
    sys.path.insert(0, str(HESE_ROOT))
    import data_loader  # type: ignore
    import binning  # type: ignore
    return data_loader, binning


def load_npz_table(path: Path) -> dict[str, np.ndarray]:
    z = np.load(path)
    return {k: z[k] for k in z.files}


def nuonly_table(tab: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {
        "energies": tab["energies"],
        "sigma_cc": tab["sigma_cc"],
        "sigma_nc": tab["sigma_nc_sm"],
        "y": tab["y"],
        "p_nc": tab["p_nc_sm"],
    }


def build_flux_modifiers(model: str, topology: str):
    """Return tau and total-muon astro flux modifiers on the v0.42d grid."""
    df = pd.read_csv(V42 / "coupled_true_energy_map.csv")
    q = df[(df.model == model) & (df.topology == topology)].copy()
    if q.empty:
        raise RuntimeError(f"No v0.42d map for {model}/{topology}")

    mus = np.sort(q.mu.unique().astype(float))
    energies = np.sort(q.E_GeV.unique().astype(float))
    loge = np.log10(energies)

    # Primary nu_mu/antinu_mu SM transport.  At HESE energies the DIS kernels
    # are flavor universal to the accuracy relevant here; use the same nu and
    # antinu SM kernels already validated for the tau calculation.
    tnu = load_npz_table(V42 / "tables" / f"{model}_nu.npz")
    tanu = load_npz_table(V42 / "tables" / f"{model}_antinu.npz")
    tr_nu = NuOnlyTransport(nuonly_table(tnu))
    tr_anu = NuOnlyTransport(nuonly_table(tanu))
    inc = (energies / 1.0e6) ** (-GAMMA_MAP)

    rtau = np.ones((len(mus), len(energies)))
    rmu = np.ones_like(rtau)
    tau_sm_grid = np.zeros_like(rtau)
    tau_bs_grid = np.zeros_like(rtau)
    mu_sm_grid = np.zeros_like(rtau)
    mu_bs_grid = np.zeros_like(rtau)

    for im, mu in enumerate(mus):
        s = q[np.isclose(q.mu, mu)].sort_values("E_GeV")
        tau_sm = s.tau_flux_SM.to_numpy(float)
        tau_bs = s.tau_flux_BSM.to_numpy(float)
        sec_sm = s.sec_numu_SM.to_numpy(float)
        sec_bs = s.sec_numu_BSM.to_numpy(float)

        # Transport the primary muon flavor through the same Earth chord.
        pnu = tr_nu.propagate(float(mu), inc, max_step_km=96.0, max_p=0.05)["nu"]
        panu = tr_anu.propagate(float(mu), inc, max_step_km=96.0, max_p=0.05)["nu"]
        primary_mu = pnu + panu
        total_mu_sm = primary_mu + sec_sm
        total_mu_bs = primary_mu + sec_bs

        rtau[im] = np.divide(tau_bs, tau_sm, out=np.ones_like(tau_sm), where=tau_sm > 0)
        rmu[im] = np.divide(total_mu_bs, total_mu_sm, out=np.ones_like(total_mu_sm), where=total_mu_sm > 0)
        tau_sm_grid[im], tau_bs_grid[im] = tau_sm, tau_bs
        mu_sm_grid[im], mu_bs_grid[im] = total_mu_sm, total_mu_bs

    itau = RegularGridInterpolator((mus, loge), rtau, bounds_error=False, fill_value=np.nan)
    imu = RegularGridInterpolator((mus, loge), rmu, bounds_error=False, fill_value=np.nan)

    diagnostics = {
        "mu_min": float(mus.min()), "mu_max": float(mus.max()),
        "E_min_GeV": float(energies.min()), "E_max_GeV": float(energies.max()),
        "max_abs_tau_modifier": float(np.max(np.abs(rtau - 1.0))),
        "max_abs_total_numu_modifier": float(np.max(np.abs(rmu - 1.0))),
    }

    def evaluate(primary_energy: np.ndarray, primary_zenith: np.ndarray, primary_type: np.ndarray):
        e = np.asarray(primary_energy, float)
        cz = np.cos(np.asarray(primary_zenith, float))
        pdg = np.asarray(primary_type, int)
        out = np.ones(len(e), float)

        # Earth-crossing events have cos(zenith)<0 in the IceCube convention.
        mu = -cz
        in_domain = (
            (cz < 0.0) & (mu >= mus.min()) & (mu <= mus.max()) &
            (e >= energies.min()) & (e <= energies.max())
        )
        if not np.any(in_domain):
            return out
        pts = np.column_stack([mu[in_domain], np.log10(e[in_domain])])
        idx = np.where(in_domain)[0]
        tau_mask = np.abs(pdg[idx]) == 16
        mu_mask = np.abs(pdg[idx]) == 14
        if np.any(tau_mask):
            vals = itau(pts[tau_mask])
            out[idx[tau_mask]] = np.where(np.isfinite(vals), vals, 1.0)
        if np.any(mu_mask):
            vals = imu(pts[mu_mask])
            out[idx[mu_mask]] = np.where(np.isfinite(vals), vals, 1.0)
        # Electron flavor and atmospheric muons are unchanged.
        return out

    return evaluate, diagnostics


def official_arrays():
    data_loader, binning = import_hese_modules()
    d = HESE_ROOT / "resources" / "data"
    mc = data_loader.load_mc([
        str(d / "HESE_mc_observable.json"),
        str(d / "HESE_mc_flux.json"),
        str(d / "HESE_mc_truth.json"),
    ])
    data = data_loader.load_data(str(d / "HESE_data.json"))
    smc, mc_slices = binning.bin_data(mc)
    sdata, data_slices = binning.bin_data(data)
    counts = np.array([len(sdata[s]) for s in data_slices], float)
    return smc, mc_slices, counts


def pars(x):
    return dict(zip(PARAM_NAMES, np.asarray(x, float)))


def event_weights(mc, x, modifier=None):
    p = pars(x)
    E = mc["primaryEnergy"].astype(float)
    pdg = mc["primaryType"].astype(int)

    astro = p["astro_norm"] * (1.0e-18 / 6.0) * (E / PIVOT_ASTRO) ** (-p["astro_gamma"])

    nu_ratio = np.where(pdg > 0, p["nunubar_ratio"], 2.0 - p["nunubar_ratio"])
    tilt_conv = p["conv_norm"] * (E / PIVOT_CR) ** (-p["cr_delta_gamma"])
    conv = (mc["pionFlux"] + p["kpi_ratio"] * mc["kaonFlux"]) * tilt_conv * nu_ratio
    conv *= mc["conventionalSelfVetoCorrection"]

    prompt = p["prompt_norm"] * (E / PIVOT_PROMPT) ** (-p["cr_delta_gamma"])
    prompt *= mc["promptFlux"] * nu_ratio * mc["promptSelfVetoCorrection"]

    # Apply BSM only to the astrophysical component.  v0.42d transport assumes
    # an equal-flavor E^-gamma incident astrophysical flux.
    if modifier is not None:
        astro = astro * modifier(E, mc["primaryZenith"], pdg)

    neutrino = (astro + conv + prompt) * mc["weightOverFluxOverLivetime"]
    muon = p["muon_norm"] * mc["muonWeightOverLivetime"]
    return LIVETIME * (neutrino + muon)


def prior_nll(x):
    p = pars(x)
    nll = 0.0
    for name, (mu, sig) in PRIORS.items():
        nll += 0.5 * ((p[name] - mu) / sig) ** 2
    return nll


def leff_loglike(k, w):
    sw = float(np.sum(w))
    sw2 = float(np.sum(w * w))
    if sw <= 0.0:
        return 0.0 if k == 0 else -1.0e100
    if sw2 <= 0.0:
        return k * math.log(sw) - sw - gammaln(k + 1.0)
    alpha = sw * sw / sw2 + 1.0
    beta = sw / sw2
    return (
        alpha * math.log(beta) + gammaln(k + alpha) - gammaln(k + 1.0)
        - (k + alpha) * math.log1p(beta) - gammaln(alpha)
    )


def leff_nll(x, mc, mc_slices, counts, modifier=None):
    w = event_weights(mc, x, modifier)
    ll = 0.0
    for k, s in zip(counts, mc_slices):
        if s.stop == s.start:
            continue
        ll += leff_loglike(float(k), w[s])
    return -ll + prior_nll(x)


def binned_expectation(mc, mc_slices, x, modifier=None):
    w = event_weights(mc, x, modifier)
    return np.array([float(np.sum(w[s])) for s in mc_slices])


def poisson_nll_from_expectation(obs, exp, x):
    exp = np.maximum(np.asarray(exp, float), 1e-300)
    obs = np.asarray(obs, float)
    return float(np.sum(exp - obs * np.log(exp) + gammaln(obs + 1.0))) + prior_nll(x)


def fit_data(mc, slices, counts, modifier=None, start=X0):
    fun = lambda x: leff_nll(x, mc, slices, counts, modifier)
    best = None
    starts = [np.asarray(start, float), X0.copy()]
    for s in starts:
        r = minimize(fun, s, method="L-BFGS-B", bounds=BOUNDS, options={"maxiter": 500, "ftol": 1e-9})
        if best is None or r.fun < best.fun:
            best = r
    return best


def fit_asimov(mc, slices, obs, modifier, start):
    def f(x):
        exp = binned_expectation(mc, slices, x, modifier)
        return poisson_nll_from_expectation(obs, exp, x)
    return minimize(f, start, method="L-BFGS-B", bounds=BOUNDS, options={"maxiter": 500, "ftol": 1e-9})


def topology_summary(mc, slices, xsm, xbsm, modifier):
    # Official ordering: cascade E-zenith bins, track E-zenith bins, double-cascade E-length bins.
    wsm = event_weights(mc, xsm, None)
    wbs = event_weights(mc, xbsm, modifier)
    nbin = len(slices)
    # Derive topology from events rather than hard-coding number of bins.
    rows = []
    for morph, label in [(0, "cascade"), (1, "track"), (2, "double_cascade")]:
        m = mc["recoMorphology"] == morph
        a = float(np.sum(wsm[m])); b = float(np.sum(wbs[m]))
        rows.append(dict(topology=label, N_SM=a, N_BSM=b, frac=(b-a)/max(a,1e-300)))
    return rows


def main():
    require_inputs()
    mc, slices, counts = official_arrays()
    print(f"HESE public sample: data={int(counts.sum())}, MC={len(mc)}, bins={len(slices)}", flush=True)

    # SM baseline validation.
    smfit = fit_data(mc, slices, counts, None)
    if not smfit.success:
        print("WARNING: SM optimizer:", smfit.message, flush=True)
    psm = pars(smfit.x)
    gamma_pull = (psm["astro_gamma"] - 2.87) / 0.20
    baseline = dict(
        n_data=int(counts.sum()), n_mc=int(len(mc)), n_bins=int(len(slices)),
        nll=float(smfit.fun), success=bool(smfit.success), message=str(smfit.message),
        published_gamma=2.87, published_gamma_sigma=0.20,
        gamma_pull_vs_published=float(gamma_pull), **{f"best_{k}":float(v) for k,v in psm.items()}
    )
    pd.DataFrame([baseline]).to_csv(OUT / "sm_baseline_validation.csv", index=False)

    # SM Asimov expectation at the public-data best fit.
    asimov = binned_expectation(mc, slices, smfit.x, None)
    sm_asimov_nll = poisson_nll_from_expectation(asimov, asimov, smfit.x)

    results, topo_rows = [], []
    for model in MODELS:
        for top in TOPOLOGIES:
            print("FIT", model, top, flush=True)
            modifier, diag = build_flux_modifiers(model, top)
            bf = fit_data(mc, slices, counts, modifier, start=smfit.x)
            pb = pars(bf.x)
            q_data_signed = 2.0 * (float(bf.fun) - float(smfit.fun))

            af = fit_asimov(mc, slices, asimov, modifier, smfit.x)
            q_asimov = max(0.0, 2.0 * (float(af.fun) - float(sm_asimov_nll)))
            z_asimov = math.sqrt(q_asimov)

            row = dict(
                model=model, topology=top,
                nll_SM=float(smfit.fun), nll_BSM=float(bf.fun),
                q_data_signed=float(q_data_signed),
                data_preference="SM" if q_data_signed > 0 else "BSM" if q_data_signed < 0 else "tie",
                q_asimov=float(q_asimov), Z_asimov=float(z_asimov),
                bsm_fit_success=bool(bf.success), asimov_fit_success=bool(af.success),
                **{f"BSM_best_{k}":float(v) for k,v in pb.items()},
                **diag,
            )
            results.append(row)
            for t in topology_summary(mc, slices, smfit.x, bf.x, modifier):
                topo_rows.append(dict(model=model, topology_assumption=top, **t))

    rdf = pd.DataFrame(results)
    rdf.to_csv(OUT / "hese_likelihood_summary.csv", index=False)
    pd.DataFrame(topo_rows).to_csv(OUT / "topology_event_summary.csv", index=False)

    meta = {
        "version": "v0.43",
        "detector": "IceCube HESE 7.5-year official public data+MC",
        "official_release_repo": "icecube/HESE-7-year-data-release",
        "likelihood": "official HESE binning + L_eff MC-statistical likelihood; detector spline nuisances fixed",
        "profiled": PARAM_NAMES,
        "BSM_source": "results/v0.42d coupled transport",
        "BSM_domain": "only Earth-crossing MC events inside the computed v0.42d mu/E grid are reweighted; all others are conservative SM controls",
        "secondary_numu": "included as total primary+tau-decay-secondary muon-flavor modifier for the astrophysical component",
        "limitations": [
            "detector-systematic splines epsilon_dom/epsilon_head_on/anisotropy are fixed rather than profiled",
            "v0.42d tau modifier is nu+nubar averaged",
            "BSM modifier outside the v0.42d computed mu/E domain is set to unity",
            "this is an independent public-data recast, not an IceCube Collaboration result",
        ],
    }
    (OUT / "metadata.json").write_text(json.dumps(meta, indent=2))

    print("\n=== v0.43 SM baseline ===")
    print(pd.DataFrame([baseline]).to_string(index=False))
    print("\n=== v0.43 HESE likelihood summary ===")
    print(rdf[["model","topology","q_data_signed","data_preference","q_asimov","Z_asimov"]].to_string(index=False))
    print("\n=== v0.43 topology event summary ===")
    print(pd.DataFrame(topo_rows).to_string(index=False))


if __name__ == "__main__":
    main()
