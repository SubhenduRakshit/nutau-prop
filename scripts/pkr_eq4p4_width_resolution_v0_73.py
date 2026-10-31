#!/usr/bin/env python3
from __future__ import annotations

import json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.73_eq4p4_width_resolution"

M = 50.0
gq = 0.03
kappa = 1.068919824203e-4
Gamma_q = 0.017903862718532193
BR_one = 0.3701925696633259
mu95 = 1.1784367170979357

a = kappa / gq
Gamma_one = a*a*M**3/(48.0*math.pi)
Gamma_inv = 2.0*Gamma_one
Gamma_tot = Gamma_q + Gamma_inv
BR_inv = Gamma_inv / Gamma_tot
r_one = 1.0 / mu95
r_hc = r_one * BR_inv / BR_one

res = {
    "version": "v0.73_eq4p4_width_resolution",
    "model_scope": "PKR Eq4.4 active-to-sterile/off-diagonal Hermitian completion",
    "M_GeV": M,
    "gq": gq,
    "kappa_GeV_inv": kappa,
    "a_GeV_inv": a,
    "Gamma_q_GeV": Gamma_q,
    "Gamma_one_Weyl_GeV": Gamma_one,
    "Gamma_inv_operator_plus_hc_GeV": Gamma_inv,
    "Gamma_total_GeV": Gamma_tot,
    "Gamma_over_M": Gamma_tot/M,
    "BR_inv": BR_inv,
    "mu95_v072": mu95,
    "r_one_v072": r_one,
    "r_operator_plus_hc": r_hc,
    "thresholds": {"robust_survival": 0.8, "robust_exclusion": 1.2},
    "width_convention": "operator_plus_hc",
    "collider_status": "ROBUST_EXCLUSION_THRESHOLD" if r_hc > 1.2 else "NOT_ROBUST_EXCLUSION",
    "note": "No new MG5/Pythia/CheckMATE run. Uses certified v0.72 likelihood and NWA BR rescaling."
}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "width_resolution.json").write_text(json.dumps(res, indent=2) + "\n")

print("=== v0.73 PKR Eq4.4 width resolution ===")
print(f"a={a:.12g} GeV^-1")
print(f"Gamma_one={Gamma_one:.12g} GeV")
print(f"Gamma_inv(O+h.c.)={Gamma_inv:.12g} GeV")
print(f"Gamma_tot/M={Gamma_tot/M:.8g}")
print(f"BR_inv={BR_inv:.12g}")
print(f"r_one={r_one:.6f}")
print(f"r_operator_plus_hc={r_hc:.6f}")
print(f"status: {res['collider_status']}")
