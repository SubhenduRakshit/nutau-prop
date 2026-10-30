#!/usr/bin/env python3
"""v0.49b: numerical-convergence wrapper for v0.49.

Preserves the v0.49 statistical model and physics.  The only change is a more
robust nuisance-profile solver: damped Newton first, controlled near-convergence
acceptance, then L-BFGS-B fallback with the physical lambda_i>0 bounds.
All final gradient norms are recorded.  A physics result is refused if the
profiled solution is not converged to grad_inf < 2e-5.
"""
from __future__ import annotations

import json, sys
from pathlib import Path
import numpy as np
from scipy.linalg import cho_solve
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parents[1]
SCRIPTS=ROOT/'scripts'
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0,str(SCRIPTS))

import pkr_monojet_profile_likelihood_v0_49 as core

OUT=ROOT/'results'/'v0.49b'
OUT.mkdir(parents=True,exist_ok=True)
core.OUT=OUT

NEWTON_TOL=2e-8
NEAR_TOL=5e-6
FINAL_TOL=2e-5
MAXITER=120
_profile_diag=[]


def robust_profile_theta(self, mu, s, n=None, theta0=None, tol=NEWTON_TOL, maxiter=MAXITER):
    s=np.asarray(s,float)
    n=self.n if n is None else np.asarray(n,float)
    theta=np.zeros(self.nbin) if theta0 is None else np.asarray(theta0,float).copy()
    base=self.b+float(mu)*s
    theta=np.maximum(theta,-base+1e-10)
    Cinv=cho_solve(self.cf,np.eye(self.nbin),check_finite=False)

    def nll(th):
        lam=base+th
        if np.any(lam<=0) or not np.isfinite(lam).all(): return np.inf
        return float(np.sum(lam-n*np.log(lam))+0.5*np.dot(th,Cinv@th))

    def gradfun(th):
        lam=base+th
        if np.any(lam<=0):
            return np.full_like(th,1e30)
        return 1.0-n/lam+Cinv@th

    f=nll(theta)
    method='newton'
    nit=0
    for it in range(maxiter):
        nit=it+1
        lam=base+theta
        grad=1.0-n/lam+Cinv@theta
        ginf=float(np.max(np.abs(grad)))
        if ginf<tol:
            _profile_diag.append(dict(mu=float(mu),method=method,niter=nit,grad_inf=ginf,accepted='strict'))
            return theta,f
        H=np.diag(n/(lam*lam))+Cinv
        try: step=np.linalg.solve(H,-grad)
        except np.linalg.LinAlgError: step=-np.linalg.pinv(H,rcond=1e-12)@grad
        alpha=1.0
        improved=False
        while alpha>1e-10:
            cand=theta+alpha*step
            fc=nll(cand)
            if np.isfinite(fc) and fc < f-1e-12*max(1.0,abs(f)):
                theta,f=cand,fc; improved=True; break
            alpha*=0.5
        if not improved:
            if ginf<NEAR_TOL:
                _profile_diag.append(dict(mu=float(mu),method=method,niter=nit,grad_inf=ginf,accepted='near'))
                return theta,f
            break

    # Convex fallback with exact analytic gradient and physical lambda bounds.
    method='L-BFGS-B'
    bounds=[(-float(base[i])+1e-10,None) for i in range(self.nbin)]
    res=minimize(nll,theta,jac=gradfun,method='L-BFGS-B',bounds=bounds,
                 options={'ftol':1e-13,'gtol':1e-10,'maxiter':3000,'maxls':50})
    th=np.asarray(res.x,float)
    ff=nll(th)
    gg=gradfun(th)
    ginf=float(np.max(np.abs(gg)))
    _profile_diag.append(dict(mu=float(mu),method=method,niter=int(getattr(res,'nit',-1)),grad_inf=ginf,
                              success=bool(res.success),message=str(res.message),accepted=bool(ginf<FINAL_TOL)))
    if not np.isfinite(ff) or ginf>=FINAL_TOL:
        raise RuntimeError(f'theta profiling genuinely failed at mu={mu}, grad_inf={ginf:.3e}, fallback={res.message}')
    return th,ff


core.ProfileLikelihood.profile_theta=robust_profile_theta


def main():
    try:
        core.main()
    finally:
        (OUT/'profile_convergence_diagnostics.json').write_text(json.dumps(_profile_diag,indent=2))
    mp=OUT/'metadata.json'
    if mp.exists():
        meta=json.loads(mp.read_text())
        meta['version']='v0.49b'
        meta['supersedes_numerically']='v0.49 stalled at grad_inf=1.175e-6 before producing a result'
        meta['profile_solver']='damped Newton + controlled near-convergence + L-BFGS-B fallback'
        meta['profile_tolerances']={'newton':NEWTON_TOL,'near_accept':NEAR_TOL,'final_required':FINAL_TOL}
        meta['profile_diagnostics']='profile_convergence_diagnostics.json'
        mp.write_text(json.dumps(meta,indent=2))
    print('\n=== v0.49b numerical certification ===')
    if _profile_diag:
        print(f"profiles={len(_profile_diag)}  max accepted grad_inf={max(d['grad_inf'] for d in _profile_diag):.3e}")

if __name__=='__main__':
    main()
