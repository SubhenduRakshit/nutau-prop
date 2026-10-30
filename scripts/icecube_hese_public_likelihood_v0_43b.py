#!/usr/bin/env python3
"""v0.43b statistical wrapper: re-profile the SM Asimov denominator.

Physics and event weights are defined in icecube_hese_public_likelihood_v0_43.py.
This wrapper corrects the expected-sensitivity denominator by fitting the SM
Asimov sample with the same nuisance priors used for the BSM hypotheses.
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "hese43", HERE / "icecube_hese_public_likelihood_v0_43.py"
)
m = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(m)


def main():
    m.require_inputs()
    mc, slices, counts = m.official_arrays()
    print(f"HESE public sample: data={int(counts.sum())}, MC={len(mc)}, bins={len(slices)}", flush=True)

    smfit = m.fit_data(mc, slices, counts, None)
    if not smfit.success:
        print("WARNING: SM data optimizer:", smfit.message, flush=True)
    psm = m.pars(smfit.x)
    gamma_pull = (psm["astro_gamma"] - 2.87) / 0.20
    baseline = dict(
        n_data=int(counts.sum()), n_mc=int(len(mc)), n_bins=int(len(slices)),
        nll=float(smfit.fun), success=bool(smfit.success), message=str(smfit.message),
        published_gamma=2.87, published_gamma_sigma=0.20,
        gamma_pull_vs_published=float(gamma_pull),
        **{f"best_{k}": float(v) for k, v in psm.items()},
    )
    pd.DataFrame([baseline]).to_csv(m.OUT / "sm_baseline_validation.csv", index=False)

    # Generate the SM Asimov sample at the real-data SM best fit, then re-fit it
    # under SM with the same nuisance priors.  This is the denominator for all
    # expected BSM tests.
    asimov = m.binned_expectation(mc, slices, smfit.x, None)
    sm_af = m.fit_asimov(mc, slices, asimov, None, smfit.x)
    if not sm_af.success:
        print("WARNING: SM Asimov optimizer:", sm_af.message, flush=True)
    sm_asimov_nll = float(sm_af.fun)

    results, topo_rows = [], []
    for model in m.MODELS:
        for top in m.TOPOLOGIES:
            print("FIT", model, top, flush=True)
            modifier, diag = m.build_flux_modifiers(model, top)

            bf = m.fit_data(mc, slices, counts, modifier, start=smfit.x)
            pb = m.pars(bf.x)
            q_data_signed = 2.0 * (float(bf.fun) - float(smfit.fun))

            af = m.fit_asimov(mc, slices, asimov, modifier, sm_af.x)
            q_asimov_raw = 2.0 * (float(af.fun) - sm_asimov_nll)
            q_asimov = max(0.0, q_asimov_raw)
            z_asimov = math.sqrt(q_asimov)

            results.append(dict(
                model=model, topology=top,
                nll_SM=float(smfit.fun), nll_BSM=float(bf.fun),
                q_data_signed=float(q_data_signed),
                data_preference="SM" if q_data_signed > 0 else "BSM" if q_data_signed < 0 else "tie",
                sm_asimov_nll=sm_asimov_nll,
                bsm_asimov_nll=float(af.fun),
                q_asimov_raw=float(q_asimov_raw),
                q_asimov=float(q_asimov), Z_asimov=float(z_asimov),
                bsm_fit_success=bool(bf.success), asimov_fit_success=bool(af.success),
                **{f"BSM_best_{k}": float(v) for k, v in pb.items()},
                **diag,
            ))
            for t in m.topology_summary(mc, slices, smfit.x, bf.x, modifier):
                topo_rows.append(dict(model=model, topology_assumption=top, **t))

    rdf = pd.DataFrame(results)
    rdf.to_csv(m.OUT / "hese_likelihood_summary.csv", index=False)
    pd.DataFrame(topo_rows).to_csv(m.OUT / "topology_event_summary.csv", index=False)

    meta = {
        "version": "v0.43b",
        "detector": "IceCube HESE 7.5-year official public data+MC",
        "official_release_repo": "icecube/HESE-7-year-data-release",
        "likelihood": "official HESE binning + L_eff MC-statistical likelihood for real data; Poisson Asimov likelihood; detector spline nuisances fixed",
        "asimov_denominator": "SM Asimov hypothesis is re-profiled with identical nuisance priors",
        "profiled": m.PARAM_NAMES,
        "BSM_source": "results/v0.42d coupled transport",
        "BSM_domain": "only Earth-crossing MC events inside the computed v0.42d mu/E grid are reweighted; outside region is left SM as a conservative control",
        "secondary_numu": "included through total primary+tau-decay-secondary astrophysical muon-flavor flux, affecting both tracks and muon-flavor NC showers",
        "limitations": [
            "IceCube detector-systematic splines epsilon_dom/epsilon_head_on/anisotropy are fixed rather than profiled",
            "v0.42d tau modifier is nu+nubar averaged",
            "BSM modifier outside the v0.42d computed mu/E domain is set to unity",
            "independent public-data recast; not an IceCube Collaboration result",
        ],
    }
    (m.OUT / "metadata.json").write_text(json.dumps(meta, indent=2))

    print("\n=== v0.43b SM baseline ===")
    print(pd.DataFrame([baseline]).to_string(index=False))
    print("\n=== v0.43b HESE likelihood summary ===")
    print(rdf[["model", "topology", "q_data_signed", "data_preference", "q_asimov", "Z_asimov"]].to_string(index=False))
    print("\n=== v0.43b topology event summary ===")
    print(pd.DataFrame(topo_rows).to_string(index=False))


if __name__ == "__main__":
    main()
