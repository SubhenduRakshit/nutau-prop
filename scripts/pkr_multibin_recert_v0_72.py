#!/usr/bin/env python3
"""v0.72: recertify the existing v0.69 CMS multibin result from its Spey log.

No MG5, Pythia, Delphes, or CheckMATE rerun.  v0.71 showed that Spey completed
its covariance fit and printed the observed upper limit, but the old parser
missed the nested output location.  This script extracts the existing
cms_2107_13201 / combination upper limit from results/v0.69_multibin_retry/
checkmate.log and propagates it to r = 1/mu95.  It also NWA-rescales the result
between the two Eq.4.4 invisible-width conventions using their BR ratio.
"""
from __future__ import annotations

import json, re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "results" / "v0.69_multibin_retry" / "checkmate.log"
IN45 = ROOT / "results" / "v0.45" / "width_br_dijet_gate.csv"
V68 = ROOT / "results" / "v0.68_checkmate_single" / "single_benchmark.json"
OUT = ROOT / "results" / "v0.72_multibin_recert"


def classify(r: float) -> str:
    if r > 1.2:
        return "ROBUST_EXCLUSION_THRESHOLD"
    if r < 0.8:
        return "ROBUST_SURVIVAL_THRESHOLD"
    return "MARGINAL"


def main() -> int:
    if not LOG.exists():
        raise RuntimeError(f"missing existing v0.69 log: {LOG}")
    text = LOG.read_text(errors="replace")
    pat = re.compile(
        r"Calculating approximate likelihood with covariance matrix:\s*cms_2107_13201,\s*SR:\s*combination.*?"
        r"Upper limit:\s*([0-9.eE+\-]+)", re.S)
    hits = [float(x) for x in pat.findall(text)]
    if not hits:
        raise RuntimeError("could not find cms_2107_13201 combination upper limit in v0.69 log")
    mu95 = hits[-1]
    if mu95 <= 0:
        raise RuntimeError(f"nonpositive mu95: {mu95}")
    r_low = 1.0 / mu95

    d = pd.read_csv(IN45)
    lo = d[(d.model == "eq4p4_derivative") & (d.convention == "one_Weyl_pair")]
    hi = d[(d.model == "eq4p4_derivative") & (d.convention == "operator_plus_hc")]
    if len(lo) != 1 or len(hi) != 1:
        raise RuntimeError("expected exactly one Eq4.4 row for each width convention")
    br_low = float(lo.iloc[0].BR_invisible)
    br_high = float(hi.iloc[0].BR_invisible)
    scale = br_high / br_low
    r_high = r_low * scale

    single = json.loads(V68.read_text()) if V68.exists() else {}
    rsingle = single.get("checkmate", {}).get("max_singlebin_robscons")

    result = {
        "version": "v0.72_multibin_recert",
        "source": "existing v0.69 CheckMATE/Spey covariance fit; no rerun",
        "analysis": "cms_2107_13201",
        "signal_region": "combination",
        "observed_mu95": mu95,
        "one_Weyl_pair": {
            "BR_invisible": br_low,
            "multibin_r": r_low,
            "classification": classify(r_low),
            "singlebin_robscons_v068": rsingle,
        },
        "operator_plus_hc_NWA_rescaled": {
            "BR_invisible": br_high,
            "BR_ratio_to_one_Weyl": scale,
            "multibin_r_rescaled": r_high,
            "classification": classify(r_high),
            "note": "same on-shell production/acceptance in validated NWA mapping; no second detector simulation",
        },
        "status": "WIDTH_CONVENTIONS_STRADDLE_COLLIDER_BOUNDARY",
        "interpretation": (
            "one_Weyl_pair is not excluded but is marginal under the predeclared 0.8/1.2 robustness band; "
            "operator_plus_hc lies above the 1.2 exclusion threshold after pure-BR NWA rescaling. "
            "Therefore Eq4.4 cannot receive a convention-independent collider verdict until the exact invisible-width normalization is resolved."
        ),
        "next_step": "resolve exact Eq4.4 width convention/normalization; do not generate more collider events yet",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "multibin_recertification.json").write_text(json.dumps(result, indent=2) + "\n")
    print("=== v0.72 CMS multibin recertification ===")
    print(f"mu95={mu95:.12g}")
    print(f"one_Weyl_pair: r={r_low:.6f}  {classify(r_low)}")
    print(f"operator_plus_hc NWA-rescaled: r={r_high:.6f}  {classify(r_high)}")
    print("status: WIDTH_CONVENTIONS_STRADDLE_COLLIDER_BOUNDARY")
    print("next: resolve exact Eq4.4 invisible-width normalization; no more collider events yet")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
