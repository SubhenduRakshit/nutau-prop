#!/usr/bin/env python3
"""v0.75 companion: uncertainty/shape audit for the single Eq4.4 best-point validation.

No generation, showering, detector simulation, or likelihood rerun. Reads the
v0.75 CheckMATE output and the existing v0.69 reference output.

Checks:
  * selected-bin signal shape after dividing out the total generated cross section;
  * CheckMATE signal-error scale ds/s in the CMS multibin combination;
  * approximate MC-only component obtained by subtracting the configured 10% XS
    uncertainty in quadrature (diagnostic only);
  * an approximate MC scale band on r; this is NOT a replacement for a full
    nuisance-parameter rerun and is explicitly labelled as such.

The physics decision is deliberately coarse: 0.8 < r < 1.2 is a marginal
validation-only region. No significance is attached to the fourth decimal place.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.75_eq4p4_bestpoint_checkmate"
NEW_JSON = OUT / "bestpoint_checkmate.json"
REF_JSON = ROOT / "results" / "v0.68_checkmate_single" / "single_benchmark.json"
NEW_CM = ROOT / "cache" / "v0.75_eq4p4_bestpoint_checkmate" / "checkmate"
REF_CM = ROOT / "cache" / "v0.69_multibin_retry" / "checkmate"
ANALYSIS = "cms_2107_13201"
XSEC_REL = 0.10


def find_total(root: Path) -> Path:
    hits = list(root.rglob("total_results.txt"))
    if not hits:
        raise RuntimeError(f"no total_results.txt under {root}")
    return sorted(hits, key=lambda p: p.stat().st_mtime)[-1]


def parse_total(path: Path) -> dict[str, dict[str, float]]:
    rows: dict[str, dict[str, float]] = {}
    lines = path.read_text(errors="replace").splitlines()
    for line in lines[1:]:
        p = line.split()
        if len(p) < 7 or p[0] != ANALYSIS:
            continue
        try:
            rows[p[1]] = {"s": float(p[5]), "ds": float(p[6])}
        except ValueError:
            continue
    if not rows:
        raise RuntimeError(f"no {ANALYSIS} rows in {path}")
    return rows


def weighted_quantile(values, weights, q):
    values = np.asarray(values, float); weights = np.asarray(weights, float)
    order = np.argsort(values); values = values[order]; weights = weights[order]
    c = np.cumsum(weights) / max(np.sum(weights), 1e-300)
    return float(np.interp(q, c, values))


def main() -> int:
    if not NEW_JSON.exists() or not REF_JSON.exists():
        raise RuntimeError("missing v0.75 or v0.68 JSON result")
    new = json.loads(NEW_JSON.read_text())
    ref = json.loads(REF_JSON.read_text())
    r = new.get("checkmate", {}).get("multibin_r")
    if r is None:
        raise RuntimeError("v0.75 multibin_r unavailable; uncertainty audit cannot certify a boundary")
    r = float(r)
    xnew = float(new["generator"]["cross_section_pb"])
    xref = float(ref["cross_section_pb"])

    pnew = find_total(NEW_CM); pref = find_total(REF_CM)
    nrows = parse_total(pnew); rrows = parse_total(pref)
    common = sorted(set(nrows) & set(rrows))
    if not common:
        raise RuntimeError("no common CMS signal regions between v0.75 and v0.69")

    # Restrict shape comparison to bins with non-negligible generated signal in both samples.
    shape = []
    ds_rel = []
    mc_rel = []
    sig_weights = []
    for sr in common:
        sn, dsn = nrows[sr]["s"], nrows[sr]["ds"]
        sr0, dsr0 = rrows[sr]["s"], rrows[sr]["ds"]
        if sn <= 0 or sr0 <= 0:
            continue
        ds_rel_i = dsn / sn
        # CheckMATE input used XSectErr=10%. If ds is quadrature-combined, this
        # removes that normalization term. We label the result diagnostic.
        mc_abs = math.sqrt(max(dsn*dsn - (XSEC_REL*sn)**2, 0.0))
        mc_rel_i = mc_abs / sn
        ds_rel.append(ds_rel_i); mc_rel.append(mc_rel_i); sig_weights.append(sn)
        if sn >= 1.0 and sr0 >= 1.0:
            shape.append((sr, (sn/xnew)/(sr0/xref), sn))

    if not ds_rel:
        raise RuntimeError("no populated bins for uncertainty audit")
    w = np.asarray(sig_weights, float)
    ds_arr = np.asarray(ds_rel, float); mc_arr = np.asarray(mc_rel, float)
    ds_med = weighted_quantile(ds_arr, w, 0.5)
    ds_p84 = weighted_quantile(ds_arr, w, 0.84)
    mc_med = weighted_quantile(mc_arr, w, 0.5)
    mc_p84 = weighted_quantile(mc_arr, w, 0.84)

    # Effective independent-bin MC error on the summed selected signal.
    mc_abs_all = []
    s_all = []
    for sr in common:
        s = nrows[sr]["s"]; ds = nrows[sr]["ds"]
        if s <= 0: continue
        mc_abs_all.append(math.sqrt(max(ds*ds - (XSEC_REL*s)**2, 0.0)))
        s_all.append(s)
    mc_scale_eff = math.sqrt(sum(x*x for x in mc_abs_all)) / max(sum(s_all), 1e-300)

    shape_mean = shape_rms = shape_max = None
    if shape:
        ratios = np.array([x[1] for x in shape], float)
        sw = np.array([x[2] for x in shape], float)
        shape_mean = float(np.sum(sw*ratios)/np.sum(sw))
        rel = ratios/shape_mean - 1.0
        shape_rms = float(math.sqrt(np.sum(sw*rel*rel)/np.sum(sw)))
        shape_max = float(np.max(np.abs(rel)))

    r_mc_lo = r * (1.0 - mc_scale_eff)
    r_mc_hi = r * (1.0 + mc_scale_eff)
    if 0.8 < r < 1.2:
        status = "MARGINAL_BOUNDARY_UNRESOLVED_VALIDATION_ONLY"
    elif r >= 1.2:
        status = "ROBUST_EXCLUSION_THRESHOLD"
    else:
        status = "ROBUST_SURVIVAL_THRESHOLD"

    result = {
        "version": "v0.75_uncertainty_companion",
        "source_total_results": str(pnew),
        "reference_total_results": str(pref),
        "common_signal_regions": len(common),
        "multibin_r_nominal": r,
        "decision_status": status,
        "configured_xsec_relative_uncertainty": XSEC_REL,
        "selected_bin_total_ds_over_s_signal_weighted_median": ds_med,
        "selected_bin_total_ds_over_s_signal_weighted_p84": ds_p84,
        "inferred_mc_only_ds_over_s_signal_weighted_median": mc_med,
        "inferred_mc_only_ds_over_s_signal_weighted_p84": mc_p84,
        "inferred_independent_bin_mc_scale_on_summed_signal": mc_scale_eff,
        "approx_r_mc_scale_band": [r_mc_lo, r_mc_hi],
        "shape_after_xsec_scaling": {
            "bins_with_s_ge_1_in_both_samples": len(shape),
            "weighted_mean_ratio": shape_mean,
            "weighted_rms_fractional_deviation_about_mean": shape_rms,
            "max_fractional_deviation_about_mean": shape_max,
        },
        "caveats": [
            "The MC-only component assumes CheckMATE combines the configured 10% XS error and MC error in quadrature.",
            "The r MC band is an overall-scale diagnostic, not a full nuisance-parameter likelihood propagation.",
            "Collider shower/matching/detector and surrogate-operator modelling uncertainties are not reduced by increasing event count.",
            "Do not interpret sub-percent movement around r=1 as a physical exclusion/survival decision."
        ],
    }
    (OUT / "uncertainty_shape_audit.json").write_text(json.dumps(result, indent=2) + "\n")

    print("=== v0.75 uncertainty / shape audit ===")
    print(f"multibin r: {r:.4f}  status: {status}")
    print(f"selected-bin ds/s weighted median: {ds_med:.2%}  p84: {ds_p84:.2%}")
    print(f"inferred MC-only ds/s weighted median: {mc_med:.2%}  p84: {mc_p84:.2%}")
    print(f"effective independent-bin MC scale: {mc_scale_eff:.2%}")
    print(f"approx r MC-scale band: [{r_mc_lo:.3f}, {r_mc_hi:.3f}]")
    if shape_mean is not None:
        print(f"shape/xsec ratio mean: {shape_mean:.4f}; weighted RMS about mean: {shape_rms:.2%}; max: {shape_max:.2%}")
    print("interpretation: validation only inside 0.8<r<1.2; no fourth-decimal boundary claim")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
