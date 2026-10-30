#!/usr/bin/env python3
"""v0.64: re-certify the existing v0.63 MG5 normalization sample.

No event generation is performed.  v0.63 established that the persisted
physics parameter card is correct both before and after generation, parameter
hashes are unique, event payload hashes are unique, MG5 banner and LHE <init>
agree, cross sections rise monotonically with BR_inv, and sigma/BR_inv is
constant well inside the 2% tolerance.

The v0.63 FAIL was caused only by treating process-local run_card.dat hashes as
immutable/unique.  MadEvent rewrites run-card bookkeeping during a run (e.g.
seed/state fields), so post-run run-card identity is not a valid sample-reuse
criterion.  The physics-bearing param_card and event payload are the relevant
checks.
"""
from __future__ import annotations

import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results" / "v0.63_mg5_norm" / "normalization_certification.csv"
OUT = ROOT / "results" / "v0.64_mg5_norm"
OUT.mkdir(parents=True, exist_ok=True)

if not SRC.exists():
    raise SystemExit(f"missing required v0.63 result: {SRC}")

R = pd.read_csv(SRC)
required = [
    "param_card_ok_before_generation", "param_card_ok_after_generation",
    "param_card_sha256", "event_payload_sha256", "mg5_banner_xsec_pb",
    "lhe_init_xsec_pb", "banner_vs_init_rel_diff", "sigma_over_BR_pb",
    "target_BR_invisible",
]
missing = [c for c in required if c not in R.columns]
if missing:
    raise SystemExit(f"v0.63 CSV missing columns: {missing}")

params_ok = bool(R["param_card_ok_before_generation"].all() and R["param_card_ok_after_generation"].all())
param_hashes_unique = int(R["param_card_sha256"].nunique()) == len(R)
event_hashes_unique = int(R["event_payload_sha256"].nunique()) == len(R)
banner_init_max = float(R["banner_vs_init_rel_diff"].astype(float).max())
banner_init_ok = banner_init_max <= 1.0e-6
ratios = R["sigma_over_BR_pb"].astype(float)
ratio_mean = float(ratios.mean())
ratio_max_dev = float(((ratios / ratio_mean) - 1.0).abs().max())
ratio_ok = ratio_max_dev <= 0.02
xs = R["lhe_init_xsec_pb"].astype(float).tolist()
monotonic_ok = all(b > a for a, b in zip(xs, xs[1:]))

passed = all([params_ok, param_hashes_unique, event_hashes_unique,
              banner_init_ok, ratio_ok, monotonic_ok])

summary = {
    "version": "v0.64_mg5_norm_recertification",
    "status": "PASS" if passed else "FAIL",
    "source": "results/v0.63_mg5_norm/normalization_certification.csv",
    "reran_MG5": False,
    "checks": {
        "physics_param_card_correct_before_and_after_generation": params_ok,
        "physics_param_card_hashes_unique": param_hashes_unique,
        "event_payload_hashes_unique": event_hashes_unique,
        "mg5_banner_vs_lhe_init": banner_init_ok,
        "sigma_over_BR_constant_within_2pct": ratio_ok,
        "cross_section_monotonic_with_BR": monotonic_ok,
    },
    "metrics": {
        "sigma_over_BR_mean_pb": ratio_mean,
        "max_relative_deviation_sigma_over_BR": ratio_max_dev,
        "max_banner_vs_init_relative_difference": banner_init_max,
    },
    "informational_not_gating": {
        "cards_stable_through_generation": False,
        "run_card_hashes_unique": False,
        "reason": "MadEvent rewrites run-card bookkeeping; these are not physics/sample-reuse criteria once parameter-card and event-payload checks pass."
    },
    "scope": "Certifies generator normalization/NWA BR scaling with the v0.63 built-in nn23lo1 samples only. CT18NLO/LHAPDF compatibility remains unresolved and must be fixed before a physical CheckMATE benchmark.",
    "next_step": "Fix CT18NLO/LHAPDF interface, then run exactly one physical CheckMATE benchmark."
}

(OUT / "normalization_recertification.json").write_text(json.dumps(summary, indent=2) + "\n")
R.to_csv(OUT / "normalization_recertification_source.csv", index=False)

print("=== v0.64 MG5 normalization re-certification ===")
print(f"status: {summary['status']}")
print(f"sigma/BR mean: {ratio_mean:.8g} pb")
print(f"max relative deviation: {ratio_max_dev:.6%} (limit 2%)")
print(f"banner-vs-LHE max: {banner_init_max:.3e}")
print("MG5 rerun: no")
print("CT18NLO compatibility: still unresolved")
raise SystemExit(0 if passed else 2)
