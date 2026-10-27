#!/usr/bin/env python3
"""Audit gate for corrected v0.15 chirality-aware contact transport.

Reuses existing v0.11/v0.15 differential tables.  No PDF integration is done.
Checks:
  1. Dirac sink monotonicity against a matched-step SM baseline;
  2. representative max_p convergence;
  3. transmission/spectral-weight diagnostics so opaque-bin ratios are not
     mistaken for detector sensitivity.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.chirality_transport import TauHelicityFlipTransport, TauMuonTransitionTransport
from nutau_prop.momentum_operators import C5_BENCHMARK, C7_BENCHMARK, lambda_4pi_tev

ROOT = Path(__file__).resolve().parents[1]
T15 = ROOT / "results" / "v0.15" / "tables"
T11 = ROOT / "results" / "v0.11" / "tables"
OUT = ROOT / "results" / "v0.15.1"
UHE = ROOT / "data" / "uhe_official"


def tag(x):
    return (f"{float(x):.6g}").replace(".", "p")


def load15(op, lam, species):
    z = np.load(T15 / f"case{op}_L{tag(lam)}_{species}.npz")
    return {k: z[k] for k in z.files}


def load6(lam, species):
    z = np.load(T11 / f"case6_L{tag(lam)}_{species}.npz")
    return {k: z[k] for k in z.files}


def incident_shape(E):
    q = pd.read_json(UHE / "flux.json")
    ef = q["Energy [GeV]"].to_numpy(float)
    sed = q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    ss = interp1d(ef, sed, bounds_error=False, fill_value=-np.inf)(E)
    # E*F(E): contribution per logarithmic energy interval, up to constants.
    w = E * (10.0 ** ss / E**2)
    m = np.nanmax(w)
    return w / m if m > 0 else w


def interp_pos(E, y, target):
    E = np.asarray(E, float) / 1e6
    y = np.asarray(y, float)
    good = np.isfinite(E) & np.isfinite(y) & (E > 0) & (y > 0)
    if not np.any(good) or target < E[good].min() or target > E[good].max():
        return np.nan
    return float(10.0 ** np.interp(np.log10(target), np.log10(E[good]), np.log10(y[good])))


def rel_excess(bs, sm, scale):
    return float(np.max(bs - sm) / max(float(scale), 1e-300))


def convergence(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    scale = max(float(np.max(np.abs(b))), 1e-300)
    mask = np.abs(b) > 1e-10 * scale
    if not np.any(mask):
        return 0.0
    return float(np.max(np.abs(a[mask] - b[mask]) / np.maximum(np.abs(b[mask]), 1e-300)))


def case5_run(lam, mu, mode, max_p):
    tn = load15(5, lam, "nu"); ta = load15(5, lam, "antinu")
    E = tn["energies"]; p0 = incident_shape(E)
    tr = TauMuonTransitionTransport(tn, ta)
    sm = tr.propagate(mu, p0, p0, p0, p0, mode=mode, use_np=False, max_p=max_p)
    bs = tr.propagate(mu, p0, p0, p0, p0, mode=mode, use_np=True, max_p=max_p)
    smtau = sm["nu_tau"] + sm["antinu_tau"]
    bstau = bs["nu_tau"] + bs["antinu_tau"]
    smmu = sm["nu_mu"] + sm["antinu_mu"]
    bsmu = bs["nu_mu"] + bs["antinu_mu"]
    return E, p0, smtau, bstau, smmu, bsmu


def caseflip_run(op, lam, mu, mode, max_p):
    if op == 6:
        tn = load6(lam, "nu"); ta = load6(lam, "antinu")
    else:
        tn = load15(op, lam, "nu"); ta = load15(op, lam, "antinu")
    E = tn["energies"]; p0 = incident_shape(E)
    tr = TauHelicityFlipTransport(tn, ta)
    sm = tr.propagate(mu, p0, p0, mode=mode, use_np=False, max_p=max_p)
    bs = tr.propagate(mu, p0, p0, mode=mode, use_np=True, max_p=max_p)
    smtau = sm["nu_tau"] + sm["antinu_tau"]
    bstau = bs["nu_tau"] + bs["antinu_tau"]
    smsec = sm["nu_mu_secondary"] + sm["antinu_mu_secondary"]
    bssec = bs["nu_mu_secondary"] + bs["antinu_mu_secondary"]
    return E, p0, smtau, bstau, smsec, bssec


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    selected_E = (1.0, 10.0, 100.0, 220.0, 300.0)
    rows = []
    gates = []

    # Case 5: endpoint is the highest-priority correction target.
    lam5 = lambda_4pi_tev(C5_BENCHMARK)
    for mu in (0.2, 0.5, 1.0):
        for mode in ("dirac_sink", "majorana_taumu"):
            E, p0, smt, bst, smm, bsm = case5_run(lam5, mu, mode, 0.02)
            E2, p02, smt2, bst2, smm2, bsm2 = case5_run(lam5, mu, mode, 0.01)
            assert np.allclose(E, E2)
            if mode == "dirac_sink":
                ex_tau = rel_excess(bst, smt, np.max(2*p0))
                ex_mu = rel_excess(bsm, smm, np.max(2*p0))
                gates.append(dict(operator=5, Lambda_TeV=lam5, cosz_abs=mu, mode=mode,
                                  gate="dirac_monotonicity", excess=max(ex_tau, ex_mu),
                                  PASS=max(ex_tau, ex_mu) <= 2e-10))
            gates.append(dict(operator=5, Lambda_TeV=lam5, cosz_abs=mu, mode=mode,
                              gate="max_p_convergence_tau", excess=convergence(bst, bst2),
                              PASS=convergence(bst, bst2) < 0.02))
            for ep in selected_E:
                inc_tau = interp_pos(E, 2*p0, ep)
                inc_pair = interp_pos(E, 4*p0, ep)
                smtau = interp_pos(E, smt, ep); bstau = interp_pos(E, bst, ep)
                smpair = interp_pos(E, smt+smm, ep); bspair = interp_pos(E, bst+bsm, ep)
                rows.append(dict(operator=5, Lambda_TeV=lam5, cosz_abs=mu, mode=mode, E_PeV=ep,
                                 incident_tau_weight=inc_tau,
                                 sm_tau_weight=smtau, bsm_tau_weight=bstau,
                                 sm_tau_transmission=smtau/inc_tau if inc_tau>0 else np.nan,
                                 bsm_tau_transmission=bstau/inc_tau if inc_tau>0 else np.nan,
                                 incident_active_pair_weight=inc_pair,
                                 sm_active_pair_weight=smpair, bsm_active_pair_weight=bspair,
                                 sm_active_pair_transmission=smpair/inc_pair if inc_pair>0 else np.nan,
                                 bsm_active_pair_transmission=bspair/inc_pair if inc_pair>0 else np.nan))

    # Case 7 endpoint, including corrected propagation of secondary muons.
    lam7 = lambda_4pi_tev(C7_BENCHMARK)
    for mu in (0.2, 0.5, 1.0):
        for mode in ("dirac_sink", "majorana_flip"):
            E, p0, smt, bst, sms, bss = caseflip_run(7, lam7, mu, mode, 0.02)
            E2, p02, smt2, bst2, sms2, bss2 = caseflip_run(7, lam7, mu, mode, 0.01)
            if mode == "dirac_sink":
                ex_tau = rel_excess(bst, smt, np.max(2*p0))
                ex_sec = rel_excess(bss, sms, np.max(2*p0))
                gates.append(dict(operator=7, Lambda_TeV=lam7, cosz_abs=mu, mode=mode,
                                  gate="dirac_monotonicity", excess=max(ex_tau, ex_sec),
                                  PASS=max(ex_tau, ex_sec) <= 2e-10))
            gates.append(dict(operator=7, Lambda_TeV=lam7, cosz_abs=mu, mode=mode,
                              gate="max_p_convergence_tau", excess=convergence(bst, bst2),
                              PASS=convergence(bst, bst2) < 0.02))
            for ep in selected_E:
                inc = interp_pos(E, 2*p0, ep)
                smtau = interp_pos(E, smt, ep); bstau = interp_pos(E, bst, ep)
                smsec = interp_pos(E, sms, ep); bssec = interp_pos(E, bss, ep)
                rows.append(dict(operator=7, Lambda_TeV=lam7, cosz_abs=mu, mode=mode, E_PeV=ep,
                                 incident_tau_weight=inc,
                                 sm_tau_weight=smtau, bsm_tau_weight=bstau,
                                 sm_tau_transmission=smtau/inc if inc>0 else np.nan,
                                 bsm_tau_transmission=bstau/inc if inc>0 else np.nan,
                                 sm_secondary_mu_weight=smsec,
                                 bsm_secondary_mu_weight=bssec))

    # Case-6 chirality robustness: representative 1.5 and 2.44 TeV points.
    for lam6 in (1.5, 2.44):
        for mu in (0.5, 1.0):
            for mode in ("dirac_sink", "majorana_flip"):
                E, p0, smt, bst, sms, bss = caseflip_run(6, lam6, mu, mode, 0.02)
                E2, p02, smt2, bst2, sms2, bss2 = caseflip_run(6, lam6, mu, mode, 0.01)
                if mode == "dirac_sink":
                    ex_tau = rel_excess(bst, smt, np.max(2*p0))
                    ex_sec = rel_excess(bss, sms, np.max(2*p0))
                    gates.append(dict(operator=6, Lambda_TeV=lam6, cosz_abs=mu, mode=mode,
                                      gate="dirac_monotonicity", excess=max(ex_tau, ex_sec),
                                      PASS=max(ex_tau, ex_sec) <= 2e-10))
                gates.append(dict(operator=6, Lambda_TeV=lam6, cosz_abs=mu, mode=mode,
                                  gate="max_p_convergence_tau", excess=convergence(bst, bst2),
                                  PASS=convergence(bst, bst2) < 0.02))

    R = pd.DataFrame(rows)
    G = pd.DataFrame(gates)
    R.to_csv(OUT / "contact_transport_transmission_v0_15_1.csv", index=False)
    G.to_csv(OUT / "contact_transport_validation_gates_v0_15_1.csv", index=False)
    status = {
        "version": "v0.15.1",
        "all_gates_PASS": bool(G.PASS.all()),
        "n_failed": int((~G.PASS).sum()),
        "warning": "Passing this numerical gate does not establish detector-level exclusion. Detector weighting remains required.",
    }
    (OUT / "contact_transport_validation_status_v0_15_1.json").write_text(json.dumps(status, indent=2))
    print("GATES")
    print(G.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print("\nTRANSMISSION (|cos z|=0.5 and 1 only)")
    print(R[R.cosz_abs.isin([0.5,1.0])].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print("\nSTATUS", json.dumps(status, indent=2))
    print("WROTE", OUT)


if __name__ == "__main__":
    main()
