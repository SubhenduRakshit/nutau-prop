#!/usr/bin/env python3
"""v0.18 production rerun with the validated finite-beta_tau transport.

Order and scope
---------------
1. Re-run momentum-enhanced contact Cases 5, 6, 7 with the v0.17 repaired
   transport.  Existing PDF/dsigma/dy tables are reused; no PDF integration.
2. Re-run only the v0.16 resolved-mediator transport from its saved tables.
3. Write physical Ehlert spectral weights, SM and BSM transmitted spectra,
   transmissions, and BSM/SM ratios.  Detector folding is intentionally NOT
   done here.

Flux convention
---------------
The KM3NeT Ehlert flux JSON is treated as an all-flavor (nu+antinu) diffuse
flux.  After oscillations, the transport input is split equally among the six
active species: nu_e, antinu_e, nu_mu, antinu_mu, nu_tau, antinu_tau.  Thus
one species receives 1/6 of E*F_all(E).  The untouched electron sector is not
propagated here but the all-flavor incident normalization is recorded so the
later detector fold can add unaffected channels consistently.

Numerics
--------
The v0.17 validation demonstrated max_p=0.02 vs 0.01 differences at the
~1e-6--2e-5 level and matter-step differences below ~6e-4 for representative
Case-6 tests.  Production therefore uses max_p=0.02 and max_step_km=12.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.chirality_transport import TauHelicityFlipTransport, TauMuonTransitionTransport
from nutau_prop.momentum_operators import C5_BENCHMARK, C7_BENCHMARK, lambda_4pi_tev

ROOT = Path(__file__).resolve().parents[1]
T11 = ROOT / "results" / "v0.11" / "tables"
T15 = ROOT / "results" / "v0.15" / "tables"
T16 = ROOT / "results" / "v0.16" / "tables"
OUT = ROOT / "results" / "v0.18"
UHE = ROOT / "data" / "uhe_official"

MAX_P = 0.02
MAX_STEP_KM = 12.0
SM_RATIO_TRANSMISSION_FLOOR = 1e-12


def tag(x):
    return (f"{float(x):.6g}").replace(".", "p")


def _load_npz(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    z = np.load(path)
    return {k: z[k] for k in z.files}


def load_contact(op, lam, species):
    if op == 6:
        # v0.11 used str(float) tags, e.g. 2p0, while other versions use 6g.
        candidates = [
            T11 / f"case6_L{str(float(lam)).replace('.', 'p')}_{species}.npz",
            T11 / f"case6_L{tag(lam)}_{species}.npz",
        ]
    else:
        candidates = [T15 / f"case{op}_L{tag(lam)}_{species}.npz"]
    for p in candidates:
        if p.exists():
            return _load_npz(p)
    raise FileNotFoundError(" or ".join(str(x) for x in candidates))


def load_resolved(mass, scale, species):
    return _load_npz(T16 / f"resolved_M{tag(mass)}_S{tag(scale)}_{species}.npz")


def ehlert_flux_per_logE(E):
    q = pd.read_json(UHE / "flux.json")
    ef = q["Energy [GeV]"].to_numpy(float)
    sed = q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    ss = interp1d(ef, sed, bounds_error=False, fill_value=-np.inf)(E)
    # F(E) = 10^sed / E^2, hence dPhi/dlnE = E F(E) = 10^sed/E.
    allflavor = 10.0**ss / E
    return allflavor, allflavor / 6.0


def safe_ratio(num, den, incident):
    raw = num / max(den, 1e-300)
    transmission = den / max(incident, 1e-300)
    masked = raw if transmission >= SM_RATIO_TRANSMISSION_FLOOR else np.nan
    return raw, masked, transmission


def base_record(kind, label, mu, mode, E, allf, relevant_incident):
    return dict(
        model_kind=kind,
        model_label=label,
        cosz_abs=float(mu),
        mode=mode,
        E_GeV=float(E),
        E_PeV=float(E / 1e6),
        incident_allflavor_per_logE=float(allf),
        incident_relevant_per_logE=float(relevant_incident),
        flux_unit="cm^-2 s^-1 sr^-1 per dlnE",
    )


def run_case5(task):
    lam, mu, mode = task
    tn = load_contact(5, lam, "nu"); ta = load_contact(5, lam, "antinu")
    E = np.asarray(tn["energies"], float)
    allf, one = ehlert_flux_per_logE(E)
    tr = TauMuonTransitionTransport(tn, ta, beta_tau=0.85e-6)
    # Equal six-species split: four relevant initial species (tau/mu, nu/antinu).
    sm = tr.propagate(mu, one, one, one, one, mode=mode, use_np=False,
                      max_step_km=MAX_STEP_KM, max_p=MAX_P)
    bs = tr.propagate(mu, one, one, one, one, mode=mode, use_np=True,
                      max_step_km=MAX_STEP_KM, max_p=MAX_P)
    smtau = sm["nu_tau"] + sm["antinu_tau"]
    bstau = bs["nu_tau"] + bs["antinu_tau"]
    smmu = sm["nu_mu"] + sm["antinu_mu"]
    bsmu = bs["nu_mu"] + bs["antinu_mu"]
    smrel = smtau + smmu; bsrel = bstau + bsmu
    rows = []
    for i, e in enumerate(E):
        inc_tau = 2.0 * one[i]; inc_mu = 2.0 * one[i]; inc_rel = 4.0 * one[i]
        rr, rm, trsm = safe_ratio(bsrel[i], smrel[i], inc_rel)
        rt, rtm, ttsm = safe_ratio(bstau[i], smtau[i], inc_tau)
        rmu, rmum, tmsm = safe_ratio(bsmu[i], smmu[i], inc_mu)
        r = base_record("contact", f"case5_L{lam:.6g}", mu, mode, e, allf[i], inc_rel)
        r.update(operator=5, Lambda_UV_TeV=float(lam),
                 SM_tau_per_logE=float(smtau[i]), BSM_tau_per_logE=float(bstau[i]),
                 SM_mu_per_logE=float(smmu[i]), BSM_mu_per_logE=float(bsmu[i]),
                 SM_relevant_per_logE=float(smrel[i]), BSM_relevant_per_logE=float(bsrel[i]),
                 SM_relevant_transmission=float(trsm),
                 BSM_relevant_transmission=float(bsrel[i]/max(inc_rel,1e-300)),
                 relevant_ratio_raw=float(rr), relevant_ratio=float(rm) if np.isfinite(rm) else np.nan,
                 tau_ratio_raw=float(rt), tau_ratio=float(rtm) if np.isfinite(rtm) else np.nan,
                 mu_ratio_raw=float(rmu), mu_ratio=float(rmum) if np.isfinite(rmum) else np.nan,
                 SM_tau_transmission=float(ttsm), SM_mu_transmission=float(tmsm))
        rows.append(r)
    return rows


def run_tauonly_contact(task):
    op, lam, mu, mode = task
    tn = load_contact(op, lam, "nu"); ta = load_contact(op, lam, "antinu")
    E = np.asarray(tn["energies"], float)
    allf, one = ehlert_flux_per_logE(E)
    tr = TauHelicityFlipTransport(tn, ta, beta_tau=0.85e-6)
    sm = tr.propagate(mu, one, one, mode=mode, use_np=False,
                      max_step_km=MAX_STEP_KM, max_p=MAX_P, track_secondary_numu=True)
    bs = tr.propagate(mu, one, one, mode=mode, use_np=True,
                      max_step_km=MAX_STEP_KM, max_p=MAX_P, track_secondary_numu=True)
    smtau = sm["nu_tau"] + sm["antinu_tau"]
    bstau = bs["nu_tau"] + bs["antinu_tau"]
    smsec = sm["nu_mu_secondary"] + sm["antinu_mu_secondary"]
    bssec = bs["nu_mu_secondary"] + bs["antinu_mu_secondary"]
    smrel = smtau + smsec; bsrel = bstau + bssec
    rows = []
    for i, e in enumerate(E):
        inc_tau = 2.0 * one[i]
        rr, rm, trsm = safe_ratio(bsrel[i], smrel[i], inc_tau)
        rt, rtm, ttsm = safe_ratio(bstau[i], smtau[i], inc_tau)
        r = base_record("contact", f"case{op}_L{lam:.6g}", mu, mode, e, allf[i], inc_tau)
        r.update(operator=int(op), Lambda_UV_TeV=float(lam),
                 SM_tau_per_logE=float(smtau[i]), BSM_tau_per_logE=float(bstau[i]),
                 SM_secondary_mu_per_logE=float(smsec[i]), BSM_secondary_mu_per_logE=float(bssec[i]),
                 SM_relevant_per_logE=float(smrel[i]), BSM_relevant_per_logE=float(bsrel[i]),
                 SM_relevant_transmission=float(trsm),
                 BSM_relevant_transmission=float(bsrel[i]/max(inc_tau,1e-300)),
                 relevant_ratio_raw=float(rr), relevant_ratio=float(rm) if np.isfinite(rm) else np.nan,
                 tau_ratio_raw=float(rt), tau_ratio=float(rtm) if np.isfinite(rtm) else np.nan,
                 SM_tau_transmission=float(ttsm))
        rows.append(r)
    return rows


def run_resolved(task):
    mass, scale, mu, mode = task
    tn = load_resolved(mass, scale, "nu"); ta = load_resolved(mass, scale, "antinu")
    E = np.asarray(tn["energies"], float)
    allf, one = ehlert_flux_per_logE(E)
    tr = TauHelicityFlipTransport(tn, ta, beta_tau=0.85e-6)
    sm = tr.propagate(mu, one, one, mode=mode, use_np=False,
                      max_step_km=MAX_STEP_KM, max_p=MAX_P, track_secondary_numu=True)
    bs = tr.propagate(mu, one, one, mode=mode, use_np=True,
                      max_step_km=MAX_STEP_KM, max_p=MAX_P, track_secondary_numu=True)
    smtau = sm["nu_tau"] + sm["antinu_tau"]
    bstau = bs["nu_tau"] + bs["antinu_tau"]
    smsec = sm["nu_mu_secondary"] + sm["antinu_mu_secondary"]
    bssec = bs["nu_mu_secondary"] + bs["antinu_mu_secondary"]
    smrel = smtau + smsec; bsrel = bstau + bssec
    rows = []
    for i, e in enumerate(E):
        inc_tau = 2.0 * one[i]
        rr, rm, trsm = safe_ratio(bsrel[i], smrel[i], inc_tau)
        rt, rtm, ttsm = safe_ratio(bstau[i], smtau[i], inc_tau)
        r = base_record("resolved_mediator", f"resolved_M{mass:.6g}_S{scale:.6g}",
                        mu, mode, e, allf[i], inc_tau)
        r.update(M_Zprime_TeV=float(mass), lowq_scale=float(scale),
                 SM_tau_per_logE=float(smtau[i]), BSM_tau_per_logE=float(bstau[i]),
                 SM_secondary_mu_per_logE=float(smsec[i]), BSM_secondary_mu_per_logE=float(bssec[i]),
                 SM_relevant_per_logE=float(smrel[i]), BSM_relevant_per_logE=float(bsrel[i]),
                 SM_relevant_transmission=float(trsm),
                 BSM_relevant_transmission=float(bsrel[i]/max(inc_tau,1e-300)),
                 relevant_ratio_raw=float(rr), relevant_ratio=float(rm) if np.isfinite(rm) else np.nan,
                 tau_ratio_raw=float(rt), tau_ratio=float(rtm) if np.isfinite(rtm) else np.nan,
                 SM_tau_transmission=float(ttsm))
        rows.append(r)
    return rows


def loginterp(q, target, col):
    x = q.E_PeV.to_numpy(float); y = q[col].to_numpy(float)
    good = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if not np.any(good) or target < x[good].min() or target > x[good].max():
        return np.nan
    return float(10.0**np.interp(np.log10(target), np.log10(x[good]), np.log10(y[good])))


def summarize(D):
    rows = []
    groups = ["model_kind", "model_label", "cosz_abs", "mode"]
    for keys, q in D.groupby(groups, dropna=False):
        kind, label, mu, mode = keys; q = q.sort_values("E_PeV")
        s = dict(model_kind=kind, model_label=label, cosz_abs=mu, mode=mode)
        energies = (0.3, 1., 10., 100., 220., 300.)
        if kind == "resolved_mediator":
            energies = energies + (1000., 10000.)
        for e in energies:
            for col in ("incident_allflavor_per_logE", "SM_relevant_per_logE",
                        "BSM_relevant_per_logE", "SM_relevant_transmission",
                        "BSM_relevant_transmission", "relevant_ratio"):
                s[f"{col}_{e:g}PeV"] = loginterp(q, e, col)
        rows.append(s)
    return pd.DataFrame(rows).sort_values(groups)


def require_inputs():
    required = [UHE / "flux.json", T15, T11, T16]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise SystemExit("Missing required inputs: " + ", ".join(missing))
    # Require the revised v0.17.1 refinement assessment to have passed.
    status_path = ROOT / "results" / "v0.17" / "tau_refinement_reassessment_status_v0_17_1.json"
    if not status_path.exists():
        raise SystemExit(f"Missing {status_path}; run scripts/reassess_tau_transport_v0_17_1.py first")
    status = json.loads(status_path.read_text())
    if not status.get("all_refinement_gates_PASS", False):
        raise SystemExit("v0.17.1 refinement reassessment did not pass")


def parallel_map(func, tasks, workers):
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=min(max(1, workers), len(tasks))) as pool:
        chunks = list(pool.imap_unordered(func, tasks, chunksize=1))
    return [r for chunk in chunks for r in chunk]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=int(os.environ.get("NUTAU_NPROC", "12")))
    ap.add_argument("--cosz", nargs="+", type=float, default=[0.05,0.1,0.2,0.5,0.8,1.0])
    ap.add_argument("--lowq-scale", type=float, default=1.0)
    a = ap.parse_args()
    require_inputs(); OUT.mkdir(parents=True, exist_ok=True)

    l5 = [3.0, lambda_4pi_tev(C5_BENCHMARK)]
    l6 = [1.5, 2.0, 2.44, 3.0]
    l7 = [3.0, lambda_4pi_tev(C7_BENCHMARK)]
    masses = [0.1,0.2,0.5,1.0,2.0,3.0]

    tasks5 = [(lam,mu,mode) for lam in l5 for mu in a.cosz
              for mode in ("dirac_sink","majorana_taumu")]
    tasks67 = [(6,lam,mu,mode) for lam in l6 for mu in a.cosz
               for mode in ("dirac_sink","majorana_flip")]
    tasks67 += [(7,lam,mu,mode) for lam in l7 for mu in a.cosz
                for mode in ("dirac_sink","majorana_flip")]
    tasksm = [(m,a.lowq_scale,mu,mode) for m in masses for mu in a.cosz
              for mode in ("dirac_sink","majorana_flip")]

    print(f"[1/3] Case 5 corrected transport: {len(tasks5)} tasks, workers={a.workers}")
    rows5 = parallel_map(run_case5, tasks5, a.workers)
    print(f"[2/3] Cases 6/7 corrected transport: {len(tasks67)} tasks")
    rows67 = parallel_map(run_tauonly_contact, tasks67, a.workers)
    contacts = pd.DataFrame(rows5 + rows67)
    contacts.to_csv(OUT / "contact_transport_corrected_v0_18.csv", index=False)
    Sc = summarize(contacts)
    Sc.to_csv(OUT / "contact_transport_corrected_summary_v0_18.csv", index=False)

    print(f"[3/3] Resolved mediator corrected transport: {len(tasksm)} tasks")
    med = pd.DataFrame(parallel_map(run_resolved, tasksm, a.workers))
    med.to_csv(OUT / "resolved_mediator_transport_corrected_v0_18.csv", index=False)
    Sm = summarize(med)
    Sm.to_csv(OUT / "resolved_mediator_transport_corrected_summary_v0_18.csv", index=False)

    status = {
        "version": "v0.18",
        "transport_core": "v0.17 exact finite-beta_tau loss+decay with Strang-style splitting",
        "max_p": MAX_P,
        "max_step_km": MAX_STEP_KM,
        "cosz_abs": list(map(float,a.cosz)),
        "flux_convention": "Ehlert all-flavor nu+antinu flux split equally among six active species",
        "ratio_mask": f"masked relevant_ratio is reported only when SM relevant transmission >= {SM_RATIO_TRANSMISSION_FLOOR:g}; raw ratio is retained separately",
        "detector_weighting_done": False,
        "interpretation": "Propagation rerun only. Do not declare exclusion/viability until detector-response weighting is performed.",
    }
    (OUT / "corrected_transport_status_v0_18.json").write_text(json.dumps(status, indent=2))

    cols = [c for c in Sc.columns if c in ("model_kind","model_label","cosz_abs","mode")
            or c.endswith("_1PeV") or c.endswith("_100PeV") or c.endswith("_220PeV")]
    print("\nCONTACT SUMMARY (selected columns)")
    print(Sc[cols].to_string(index=False, float_format=lambda x:f"{x:.4g}"))
    colsm = [c for c in Sm.columns if c in ("model_kind","model_label","cosz_abs","mode")
             or c.endswith("_100PeV") or c.endswith("_220PeV") or c.endswith("_1000PeV")]
    print("\nRESOLVED MEDIATOR SUMMARY (selected columns)")
    print(Sm[colsm].to_string(index=False, float_format=lambda x:f"{x:.4g}"))
    print("\nWROTE", OUT)


if __name__ == "__main__":
    main()
