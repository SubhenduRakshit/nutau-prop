#!/usr/bin/env python3
"""v0.70: diagnose why cms_2107_13201 multibin limit was skipped.

Reads only the existing v0.69 CheckMATE output. No MG5, Pythia, Delphes,
or likelihood rerun. For the CMS multibin `combination`, report signal yield s,
CheckMATE signal uncertainty ds, and the exact Spey-wrapper gate
s - 1.64*min(s,ds). The multibin covariance likelihood is skipped when the
maximum of this quantity is <= 0.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = "eq4p4_derivative__one_Weyl_pair"
CMROOT = ROOT / "cache" / "v0.69_multibin_retry" / "checkmate" / TAG
OUT = ROOT / "results" / "v0.70_multibin_mcstat"


def find_total_results() -> Path:
    # CheckMATE can insert a Name/output subdirectory below OutputDirectory.
    # v0.69's successful parser already searches recursively, so do the same.
    candidates = list(CMROOT.rglob("total_results.txt"))
    candidates = [p for p in candidates if p.parent.name == "evaluation"]
    if not candidates:
        raise RuntimeError(f"no evaluation/total_results.txt found recursively below {CMROOT}")
    if len(candidates) > 1:
        candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]


def find_varj() -> Path:
    candidates = list((ROOT / "cache" / "v0.51").rglob("cms_2107_13201_var.j"))
    if not candidates:
        raise RuntimeError("cannot find cms_2107_13201_var.j under cache/v0.51")
    return sorted(candidates, key=lambda p: ("checkmate-install" not in str(p), len(str(p))))[0]


def main() -> int:
    total = find_total_results()
    varj = find_varj()
    meta = json.loads(varj.read_text())
    names = meta["mb_signal_regions"]["combination"]

    rows = {}
    for line in total.read_text(errors="replace").splitlines()[1:]:
        p = line.split()
        if len(p) < 7 or p[0] != "cms_2107_13201":
            continue
        rows[p[1]] = dict(o=float(p[2]), b=float(p[3]), db=float(p[4]), s=float(p[5]), ds=float(p[6]))

    audit = []
    missing = []
    for name in names:
        if name not in rows:
            missing.append(name)
            continue
        r = rows[name]
        ds_eff = max(min(r["s"], r["ds"]), 1.0e-4)
        margin = r["s"] - 1.64 * ds_eff
        audit.append({"sr": name, **r, "ds_eff": ds_eff, "s_over_ds_eff": r["s"] / ds_eff if ds_eff else None,
                      "s_minus_1p64ds": margin})

    if not audit:
        raise RuntimeError("no CMS multibin SRs found in total_results.txt")
    best = max(audit, key=lambda x: x["s_minus_1p64ds"])
    npos = sum(x["s_minus_1p64ds"] > 0 for x in audit)
    criterion_pass = best["s_minus_1p64ds"] > 0

    OUT.mkdir(parents=True, exist_ok=True)
    result = {
        "version": "v0.70_multibin_mcstat",
        "source": str(total),
        "analysis": "cms_2107_13201",
        "combination_bins_expected": len(names),
        "combination_bins_found": len(audit),
        "missing_bins": missing,
        "spey_mcstat_gate_pass": criterion_pass,
        "positive_margin_bins": npos,
        "best_bin": best,
        "rows": audit,
    }
    (OUT / "mcstat_diagnostic.json").write_text(json.dumps(result, indent=2) + "\n")

    print("=== v0.70 CMS multibin MC-stat diagnostic ===")
    print(f"source: {total}")
    print(f"bins: {len(audit)}/{len(names)}; positive s-1.64ds bins: {npos}")
    print(f"best SR: {best['sr']}")
    print(f"  s={best['s']:.8g}  ds={best['ds']:.8g}  ds_eff={best['ds_eff']:.8g}")
    print(f"  s/ds_eff={best['s_over_ds_eff']:.4f}")
    print(f"  s-1.64*ds_eff={best['s_minus_1p64ds']:.8g}")
    print(f"Spey covariance gate: {'PASS' if criterion_pass else 'FAIL_MC_STAT'}")
    if not criterion_pass:
        # Pure-stat scaling estimate: ds/s ~ 1/sqrt(N), hence required event factor ~ (1.64/(s/ds))^2.
        z = best['s_over_ds_eff']
        fac = (1.64 / z) ** 2 if z and z > 0 else None
        if fac is not None:
            print(f"naive pure-stat event multiplier to cross gate: {fac:.2f}x")
            print(f"suggested safety multiplier: {max(2.0, 1.5*fac):.2f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
