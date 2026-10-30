#!/usr/bin/env python3
"""v0.48b numerical-stability wrapper for the v0.48 CMS monojet recast.

v0.48 failed on current NumPy/SciPy because the rounded public covariance can
lead the legacy simplified-likelihood implementation to an invalid slogdet,
after which the gq root solve has no finite bracket.  This wrapper preserves
v0.48 unchanged and patches only the numerical layer:

  * enforce exact symmetry of the published covariance;
  * diagonalise it and lift only numerically non-positive eigenvalues to a
    tiny relative floor;
  * record the size of that repair and refuse a large alteration;
  * require a finite positive mu95;
  * use an adaptive gq bracket and return NaN (with diagnostics) if the target
    yield is outside the physical NWA scaling range.

All physics, signal templates, morphing, likelihood data and output logic are
those of v0.48.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import pkr_monojet_template_recast_v0_48 as core

OUT = ROOT / "results" / "v0.48b"
OUT.mkdir(parents=True, exist_ok=True)
core.OUT = OUT

COV_REL_FLOOR = 1.0e-12
MAX_REL_FROBENIUS_REPAIR = 1.0e-5
_cov_diag = {}
_root_diag = []

_orig_load = core.load_cms_likelihood
_orig_ul = core.ul_on_mu


def robust_load_cms_likelihood():
    sr, nobs, nbg, C0 = _orig_load()
    C0 = np.asarray(C0, dtype=float)
    if not np.isfinite(C0).all():
        bad = np.argwhere(~np.isfinite(C0))
        raise RuntimeError(f"CMS covariance contains non-finite entries at {bad[:10].tolist()}")

    Csym = 0.5 * (C0 + C0.T)
    sym_rel = float(np.linalg.norm(C0 - C0.T, ord='fro') / max(np.linalg.norm(Csym, ord='fro'), 1e-300))
    w, V = np.linalg.eigh(Csym)
    wmin0 = float(w.min())
    wmax0 = float(w.max())
    floor = max(COV_REL_FLOOR * max(wmax0, 1.0), np.finfo(float).eps)
    wfix = np.maximum(w, floor)
    C = (V * wfix) @ V.T
    C = 0.5 * (C + C.T)

    rel_repair = float(np.linalg.norm(C - C0, ord='fro') / max(np.linalg.norm(C0, ord='fro'), 1e-300))
    sign0, logdet0 = np.linalg.slogdet(C0)
    sign1, logdet1 = np.linalg.slogdet(C)
    cond1 = float(np.linalg.cond(C))

    _cov_diag.update({
        'finite_original': True,
        'symmetry_relative_error': sym_rel,
        'min_eigenvalue_original': wmin0,
        'max_eigenvalue_original': wmax0,
        'negative_eigenvalues_original': int(np.sum(w < 0.0)),
        'eigenvalue_floor': float(floor),
        'relative_frobenius_repair': rel_repair,
        'slogdet_original_sign': float(sign0),
        'slogdet_original_logabsdet': float(logdet0) if np.isfinite(logdet0) else None,
        'slogdet_repaired_sign': float(sign1),
        'slogdet_repaired_logabsdet': float(logdet1) if np.isfinite(logdet1) else None,
        'condition_number_repaired': cond1,
    })
    (OUT / 'covariance_diagnostics.json').write_text(json.dumps(_cov_diag, indent=2))

    if rel_repair > MAX_REL_FROBENIUS_REPAIR:
        raise RuntimeError(
            f"Covariance PSD repair is too large ({rel_repair:.3e}); "
            "refusing to alter the public likelihood."
        )
    if sign1 <= 0 or not np.isfinite(logdet1):
        raise RuntimeError("Repaired covariance is not positive definite.")

    print("\n=== v0.48b covariance diagnostics ===")
    print(f"min eig before={wmin0:.6e}  max eig={wmax0:.6e}  negatives={int(np.sum(w<0))}")
    print(f"relative repair={rel_repair:.3e}  repaired cond={cond1:.3e}")
    return sr, nobs, nbg, C


def robust_ul_on_mu(ns, nobs, nbg, C):
    mu = _orig_ul(ns, nobs, nbg, C)
    if not np.isfinite(mu) or mu <= 0.0:
        raise RuntimeError(f"Simplified-likelihood returned invalid mu95={mu!r} after covariance repair")
    return float(mu)


def robust_solve_gq_from_mu(mu_ul, m, mdm_tpl):
    if not np.isfinite(mu_ul) or mu_ul <= 0.0:
        return float('nan')
    _, _, brref = core.vector_widths_over_m(m, mdm_tpl, core.GQ_REF, core.GDM_REF)
    target_mdm = m / 3.0

    def scale(g):
        _, _, br = core.vector_widths_over_m(m, target_mdm, g, core.GDM_REF)
        return (g * g * br) / (core.GQ_REF * core.GQ_REF * brref)

    # The NWA signal scaling is not guaranteed to grow forever because BR_inv
    # decreases at large gq.  Find a genuine sign-changing interval on a grid.
    grid = np.geomspace(1.0e-6, 3.0, 500)
    vals = np.asarray([scale(float(g)) - mu_ul for g in grid])
    roots = []
    for i in range(len(grid) - 1):
        if not (np.isfinite(vals[i]) and np.isfinite(vals[i + 1])):
            continue
        if vals[i] == 0.0:
            roots.append(float(grid[i]))
        elif vals[i] * vals[i + 1] < 0.0:
            roots.append(float(brentq(lambda z: scale(z) - mu_ul, float(grid[i]), float(grid[i + 1]))))
    max_scale = float(max(scale(float(g)) for g in grid))
    if not roots:
        _root_diag.append(dict(m_GeV=float(m), mu95=float(mu_ul), max_physical_scale=max_scale, root=None))
        return float('nan')
    # Physical exclusion curve uses the smaller-coupling branch.
    root = min(roots)
    _root_diag.append(dict(m_GeV=float(m), mu95=float(mu_ul), max_physical_scale=max_scale, root=root))
    return root


core.load_cms_likelihood = robust_load_cms_likelihood
core.ul_on_mu = robust_ul_on_mu
core.solve_gq_from_mu = robust_solve_gq_from_mu


def main():
    core.main()
    (OUT / 'root_diagnostics.json').write_text(json.dumps(_root_diag, indent=2))
    mp = OUT / 'metadata.json'
    if mp.exists():
        meta = json.loads(mp.read_text())
        meta['version'] = 'v0.48b'
        meta['supersedes_numerically'] = 'v0.48 failed before producing a physics result'
        meta['covariance_stabilization'] = {
            'relative_eigenvalue_floor': COV_REL_FLOOR,
            'max_allowed_relative_frobenius_repair': MAX_REL_FROBENIUS_REPAIR,
            'diagnostics_file': 'covariance_diagnostics.json',
        }
        meta['root_handling'] = 'adaptive scan for a physical sign-changing NWA gq branch; diagnostics recorded'
        mp.write_text(json.dumps(meta, indent=2))
    print("\n=== v0.48b numerical certification ===")
    print("covariance repair accepted; see covariance_diagnostics.json")


if __name__ == '__main__':
    main()
