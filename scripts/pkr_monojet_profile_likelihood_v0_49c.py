#!/usr/bin/env python3
"""v0.49c: stationarity-controlled nuisance profiling for the CMS monojet recast.

v0.49b exposed that the previous fallback could stop on a function-change
criterion while ||grad NLL||_inf was still too large.  This wrapper preserves
the v0.49 statistical model and physics, but fixes the numerical optimizer:

* damped Newton with a proper Armijo directional-decrease condition;
* bounded least_squares on grad_theta NLL = 0 as the fallback;
* exact analytic Hessian/Jacobian;
* explicit final stationarity requirement before any physics result is used.

v0.49 and v0.49b remain preserved as failed numerical audit versions.
"""
from __future__ import annotations

import json, sys
from pathlib import Path
import numpy as np
from scipy.linalg import cho_solve
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import pkr_monojet_profile_likelihood_v0_49 as core

OUT = ROOT / "results" / "v0.49c"
OUT.mkdir(parents=True, exist_ok=True)
core.OUT = OUT

STRICT_TOL = 1.0e-8
FINAL_TOL = 1.0e-6
MAX_NEWTON = 160
ARMIJO_C1 = 1.0e-4
_profile_diag = []


def robust_profile_theta(self, mu, s, n=None, theta0=None, tol=STRICT_TOL, maxiter=MAX_NEWTON):
    s = np.asarray(s, float)
    n = self.n if n is None else np.asarray(n, float)
    base = self.b + float(mu) * s
    Cinv = cho_solve(self.cf, np.eye(self.nbin), check_finite=False)

    # Keep lambda_i=base_i+theta_i strictly positive.
    eps = np.maximum(1.0e-10, 1.0e-12 * np.maximum(base, 1.0))
    lower = -base + eps
    theta = np.zeros(self.nbin) if theta0 is None else np.asarray(theta0, float).copy()
    theta = np.maximum(theta, lower + eps)

    def nll(th):
        lam = base + th
        if np.any(lam <= 0.0) or not np.isfinite(lam).all():
            return np.inf
        return float(np.sum(lam - n * np.log(lam)) + 0.5 * th @ (Cinv @ th))

    def gradfun(th):
        lam = base + th
        if np.any(lam <= 0.0) or not np.isfinite(lam).all():
            return np.full_like(th, 1.0e30)
        return 1.0 - n / lam + Cinv @ th

    def jacfun(th):
        lam = base + th
        if np.any(lam <= 0.0):
            return np.eye(self.nbin) * 1.0e30
        return np.diag(n / (lam * lam)) + Cinv

    f = nll(theta)
    last_g = np.inf
    # Convex Newton solve with a scale-independent Armijo condition.
    for it in range(maxiter):
        grad = gradfun(theta)
        ginf = float(np.max(np.abs(grad)))
        last_g = ginf
        if ginf < tol:
            _profile_diag.append(dict(mu=float(mu), method='newton_armijo', niter=it+1,
                                      grad_inf=ginf, accepted='strict'))
            return theta, f
        H = jacfun(theta)
        try:
            step = np.linalg.solve(H, -grad)
        except np.linalg.LinAlgError:
            step = -np.linalg.pinv(H, rcond=1e-12) @ grad
        directional = float(grad @ step)
        if not np.isfinite(directional) or directional >= 0.0:
            break

        alpha = 1.0
        accepted = False
        while alpha > 1.0e-12:
            cand = theta + alpha * step
            if np.any(cand <= lower):
                alpha *= 0.5
                continue
            fc = nll(cand)
            if np.isfinite(fc) and fc <= f + ARMIJO_C1 * alpha * directional:
                theta, f = cand, fc
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            break

    # Target the actual stationarity equations rather than function decrease.
    x0 = np.maximum(theta, lower + eps)
    res = least_squares(
        gradfun, x0, jac=jacfun, bounds=(lower, np.full(self.nbin, np.inf)),
        method='trf', x_scale='jac', ftol=1.0e-13, xtol=1.0e-13, gtol=1.0e-13,
        max_nfev=6000, verbose=0,
    )
    th = np.asarray(res.x, float)

    # One exact Newton polishing step is often useful after least_squares.
    for _ in range(8):
        gg = gradfun(th)
        if float(np.max(np.abs(gg))) < FINAL_TOL:
            break
        H = jacfun(th)
        try:
            step = np.linalg.solve(H, -gg)
        except np.linalg.LinAlgError:
            break
        alpha = 1.0
        ff = nll(th)
        dd = float(gg @ step)
        while alpha > 1e-12:
            cand = th + alpha * step
            if np.any(cand <= lower):
                alpha *= 0.5
                continue
            fc = nll(cand)
            if np.isfinite(fc) and fc <= ff + ARMIJO_C1 * alpha * dd:
                th = cand
                break
            alpha *= 0.5
        else:
            break

    ff = nll(th)
    gg = gradfun(th)
    ginf = float(np.max(np.abs(gg)))
    _profile_diag.append(dict(
        mu=float(mu), method='least_squares_stationarity', niter=int(getattr(res, 'nfev', -1)),
        grad_inf=ginf, lsq_cost=float(res.cost), lsq_optimality=float(res.optimality),
        success=bool(res.success), status=int(res.status), message=str(res.message),
        accepted=bool(np.isfinite(ff) and ginf < FINAL_TOL), newton_grad_before=float(last_g),
    ))
    if not np.isfinite(ff) or ginf >= FINAL_TOL:
        raise RuntimeError(
            f'theta stationarity solve failed at mu={mu}, grad_inf={ginf:.3e}, '
            f'least_squares={res.message}'
        )
    return th, ff


core.ProfileLikelihood.profile_theta = robust_profile_theta


def main():
    try:
        core.main()
    finally:
        (OUT / 'profile_convergence_diagnostics.json').write_text(json.dumps(_profile_diag, indent=2))

    mp = OUT / 'metadata.json'
    if mp.exists():
        meta = json.loads(mp.read_text())
        meta['version'] = 'v0.49c'
        meta['supersedes_numerically'] = [
            'v0.49 stalled on an over-strict Newton stopping rule',
            'v0.49b L-BFGS-B stopped on relative function reduction with grad_inf=3.779e-5',
        ]
        meta['profile_solver'] = 'Armijo damped Newton + bounded least_squares on grad_theta NLL = 0 + Newton polish'
        meta['profile_tolerances'] = {'strict_newton': STRICT_TOL, 'final_stationarity': FINAL_TOL}
        meta['profile_diagnostics'] = 'profile_convergence_diagnostics.json'
        mp.write_text(json.dumps(meta, indent=2))

    print('\n=== v0.49c numerical certification ===')
    if _profile_diag:
        vals = [d['grad_inf'] for d in _profile_diag if d.get('accepted')]
        print(f"profiles={len(_profile_diag)}  accepted={len(vals)}  max accepted grad_inf={max(vals) if vals else float('nan'):.3e}")


if __name__ == '__main__':
    main()
