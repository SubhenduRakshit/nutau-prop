#!/usr/bin/env python3
"""v0.49: independent CMS EXO-20-004 monojet profile-likelihood recast.

Purpose
-------
v0.48/v0.48b exposed a failure of the legacy LLP-CDFO UpperLimitComputer on
this current NumPy/SciPy stack: its returned mu95~5.6 is incompatible with the
published low-mass spin-1 templates and therefore cannot be used for physics.

This script keeps the same public CMS inputs and the same 100->50 GeV signal
morphing, but replaces the legacy statistical layer by an explicit profile
likelihood:

  n_i ~ Poisson(mu*s_i + b_i + theta_i)
  theta ~ N(0, C)

The 66 correlated Gaussian nuisances are profiled with a damped Newton solve
using a Cholesky factorisation of the positive-definite public covariance.
Observed upper limits are obtained with the asymptotic CLs construction.  The
100 and 250 GeV public vector templates are mandatory validation points before
any 50 GeV PKR verdict is accepted.
"""
from __future__ import annotations

import json, math, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import pkr_monojet_template_recast_v0_48 as core

OUT = ROOT / "results" / "v0.49"
OUT.mkdir(parents=True, exist_ok=True)
core.OUT = OUT

VALIDATION_RATIO_MAX = 2.0   # monojet-only may be weaker than published monojet+mono-V
MORPH_FLOOR = 0.15


def covariance_pd(C0):
    C0 = np.asarray(C0, float)
    Csym = 0.5 * (C0 + C0.T)
    w, V = np.linalg.eigh(Csym)
    wmin, wmax = float(w.min()), float(w.max())
    floor = max(1e-12 * max(wmax, 1.0), np.finfo(float).eps)
    wfix = np.maximum(w, floor)
    C = (V * wfix) @ V.T
    C = 0.5 * (C + C.T)
    rel = float(np.linalg.norm(C-C0, ord='fro') / max(np.linalg.norm(C0, ord='fro'),1e-300))
    if rel > 1e-5:
        raise RuntimeError(f"Covariance PSD repair too large: {rel:.3e}")
    cf = cho_factor(C, lower=True, check_finite=True)
    return C, cf, dict(min_eigenvalue=wmin,max_eigenvalue=wmax,
                       negative_eigenvalues=int(np.sum(w<0)),
                       relative_frobenius_repair=rel,
                       condition_number=float(np.linalg.cond(C)))


class ProfileLikelihood:
    def __init__(self, nobs, bkg, cov):
        self.n = np.asarray(nobs, float)
        self.b = np.asarray(bkg, float)
        self.C, self.cf, self.covdiag = covariance_pd(cov)
        self.nbin = len(self.n)

    def _cinv(self, x):
        return cho_solve(self.cf, x, check_finite=False)

    def profile_theta(self, mu, s, n=None, theta0=None, tol=2e-9, maxiter=100):
        s = np.asarray(s,float)
        n = self.n if n is None else np.asarray(n,float)
        theta = np.zeros(self.nbin) if theta0 is None else np.asarray(theta0,float).copy()
        base = self.b + float(mu)*s
        # ensure a feasible starting point
        lam = base + theta
        if np.any(lam <= 0):
            theta = np.maximum(theta, -base + 1e-9)

        def nll(th):
            l = base + th
            if np.any(l <= 0) or not np.isfinite(l).all():
                return np.inf
            # constants independent of mu/theta omitted
            return float(np.sum(l - n*np.log(l)) + 0.5*np.dot(th,self._cinv(th)))

        f = nll(theta)
        for _ in range(maxiter):
            lam = base + theta
            invCtheta = self._cinv(theta)
            grad = 1.0 - n/lam + invCtheta
            if np.max(np.abs(grad)) < tol:
                return theta, f
            H = np.diag(n/(lam*lam))
            # H += C^{-1}, constructed columnwise via Cholesky solves
            H += cho_solve(self.cf, np.eye(self.nbin), check_finite=False)
            try:
                step = np.linalg.solve(H, -grad)
            except np.linalg.LinAlgError:
                step = -np.linalg.pinv(H, rcond=1e-12) @ grad
            # feasibility + Armijo-like backtracking
            alpha = 1.0
            while alpha > 1e-8:
                cand = theta + alpha*step
                fc = nll(cand)
                if np.isfinite(fc) and fc < f:
                    theta, f = cand, fc
                    break
                alpha *= 0.5
            else:
                # If Newton stalls extremely close to optimum, accept small gradient.
                if np.max(np.abs(grad)) < 1e-6:
                    return theta, f
                raise RuntimeError(f"theta profiling stalled at mu={mu}, grad_inf={np.max(np.abs(grad)):.3e}")
        raise RuntimeError(f"theta profiling failed to converge at mu={mu}")

    def profiled_nll(self, mu, s, n=None):
        _, f = self.profile_theta(mu,s,n=n)
        return f

    def muhat(self, s, n=None):
        n = self.n if n is None else np.asarray(n,float)
        f0 = self.profiled_nll(0.0,s,n=n)
        # Scale estimate for an upper search bound from Fisher information.
        var = self.b + np.diag(self.C)
        info = float(np.sum(np.asarray(s,float)**2 / np.maximum(var,1.0)))
        sig = 1.0/math.sqrt(max(info,1e-18))
        upper = max(5.0*sig, 1.0)
        # expand if minimum is sitting on the upper boundary
        best = None
        for _ in range(8):
            res = minimize_scalar(lambda x:self.profiled_nll(float(x),s,n=n),bounds=(0.0,upper),method='bounded',options={'xatol':1e-6})
            best = res
            if res.x < 0.95*upper:
                break
            upper *= 4.0
        mu = max(0.0,float(best.x))
        fm = float(best.fun)
        if f0 <= fm + 1e-8:
            return 0.0, f0
        return mu, fm

    def cls(self, mu, s):
        mu = float(mu)
        muhat, nllhat = self.muhat(s)
        nllmu = self.profiled_nll(mu,s)
        # one-sided q_mu
        qmu = 0.0 if muhat > mu else max(0.0,2.0*(nllmu-nllhat))

        # background Asimov dataset; profile exactly, do not use Gaussian shortcut
        nA = self.b.copy()
        muhatA, nllhatA = self.muhat(s,n=nA)
        nllmuA = self.profiled_nll(mu,s,n=nA)
        qA = 0.0 if muhatA > mu else max(0.0,2.0*(nllmuA-nllhatA))

        sqmu, sqA = math.sqrt(qmu), math.sqrt(qA)
        # Same asymptotic CLs mapping used by the LLP-CDFO statistics module.
        if qA >= qmu:
            CLsb = 1.0 - norm.cdf(sqmu)
            CLb = norm.cdf(sqA - sqmu)
        else:
            if qA == 0.0:
                CLsb, CLb = 1.0, 1.0
            else:
                CLsb = 1.0 - norm.cdf((qmu+qA)/(2.0*sqA))
                CLb = 1.0 - norm.cdf((qmu-qA)/(2.0*sqA))
        return CLsb/max(CLb,1e-300), dict(muhat=muhat,qmu=qmu,qA=qA)

    def upper_limit(self, s, alpha=0.05):
        s=np.asarray(s,float)
        if not np.any(s>0):
            raise RuntimeError('zero signal')
        def f(mu): return self.cls(mu,s)[0]-alpha
        lo=0.0; flo=f(lo)
        hi=1e-4; fhi=f(hi)
        while fhi > 0.0 and hi < 1e6:
            hi *= 2.0
            fhi=f(hi)
        if not (np.isfinite(flo) and np.isfinite(fhi) and flo*fhi <= 0):
            raise RuntimeError(f'could not bracket CLs=0.05: f0={flo}, hi={hi}, fhi={fhi}')
        mu=brentq(f,lo,hi,xtol=1e-6,rtol=1e-5,maxiter=100)
        clsval,diag=self.cls(mu,s)
        return float(mu),float(clsval),diag


def gq_from_mu(mu_ul,m,mdm_tpl):
    _,_,brref=core.vector_widths_over_m(m,mdm_tpl,core.GQ_REF,core.GDM_REF)
    target=m/3.0
    def scale(g):
        _,_,br=core.vector_widths_over_m(m,target,g,core.GDM_REF)
        return g*g*br/(core.GQ_REF**2*brref)
    grid=np.geomspace(1e-6,3.0,1000)
    vals=np.array([scale(float(g))-mu_ul for g in grid])
    roots=[]
    for a,b,fa,fb in zip(grid[:-1],grid[1:],vals[:-1],vals[1:]):
        if fa==0: roots.append(float(a))
        elif fa*fb<0: roots.append(float(brentq(lambda z:scale(z)-mu_ul,float(a),float(b))))
    return min(roots) if roots else float('nan')


def main():
    core.require_inputs()
    blocks=core.parse_blocks(core.SIGFILE)
    sr,nobs,nbg,C0=core.load_cms_likelihood()
    like=ProfileLikelihood(nobs,nbg,C0)
    (OUT/'covariance_diagnostics.json').write_text(json.dumps(like.covdiag,indent=2))

    pub=core.published_gq_observed()
    p100=core.choose_point(blocks,100.0,100.0/3.0)
    p250=core.choose_point(blocks,250.0,250.0/3.0)
    ns100=core.vector_from_template(core.template_lookup(blocks,p100),sr)
    ns250=core.vector_from_template(core.template_lookup(blocks,p250),sr)

    validation=[]
    for m,p,ns in [(100.0,p100,ns100),(250.0,p250,ns250)]:
        mu95,cls,diag=like.upper_limit(ns)
        g95=gq_from_mu(mu95,m,p[1])
        gpub=float(pub.get(m,np.nan))
        ratio=g95/gpub if np.isfinite(g95) and np.isfinite(gpub) else np.nan
        validation.append(dict(mmed=m,template_mdm=p[1],signal_sum=float(ns.sum()),signal_first_bin=float(ns[0]),mu95=mu95,CLs_at_limit=cls,monojet_recast_gq95=g95,published_combined_gq95=gpub,ratio=ratio,muhat=diag['muhat'],qmu=diag['qmu'],qA=diag['qA']))
    V=pd.DataFrame(validation)
    V.to_csv(OUT/'public_likelihood_validation.csv',index=False)
    if not np.isfinite(V['monojet_recast_gq95']).all():
        raise RuntimeError('public-template validation produced non-finite gq95')
    if (V['ratio'] > VALIDATION_RATIO_MAX).any() or (V['ratio'] < 0.5).any():
        raise RuntimeError('independent likelihood fails public coupling-scale validation; refusing PKR verdict')

    # Same v0.48 morphing, with closure uncertainty.
    pdf=core.LHAPDFNucleon(core.PDF_SET,0)
    r250,mrows250=core.morph_ratios(pdf,250.0,100.0)
    pred250=core.apply_bin_ratio(ns100,sr,r250)
    _,_,br100=core.vector_widths_over_m(100.,p100[1],core.GQ_REF,core.GDM_REF)
    _,_,br250=core.vector_widths_over_m(250.,p250[1],core.GQ_REF,core.GDM_REF)
    pred250*=br250/max(br100,1e-300)
    nz=ns250>0
    frac=(pred250[nz]-ns250[nz])/ns250[nz]
    weighted_total=pred250.sum()/max(ns250.sum(),1e-300)-1.0
    rms=float(np.sqrt(np.mean(np.clip(frac,-5,5)**2)))
    morph_unc=max(MORPH_FLOOR,abs(weighted_total),min(rms,0.50))

    r50,mrows50=core.morph_ratios(pdf,50.0,100.0)
    pd.DataFrame(mrows250+mrows50).to_csv(OUT/'lo_mass_morphing.csv',index=False)
    base50=core.apply_bin_ratio(ns100,sr,r50)

    d45=pd.read_csv(core.IN45)
    rows=[]
    for _,x in d45.iterrows():
        scale=(float(x.gq)**2*float(x.BR_invisible))/(core.GQ_REF**2*br100)
        ns=base50*scale
        # Compute exclusion ratio directly: r=1/mu95 for the PKR signal hypothesis.
        mu95,cls,diag=like.upper_limit(ns)
        r=1.0/mu95
        rlo=r*(1.0-morph_unc); rhi=r*(1.0+morph_unc)
        if rlo>1.0: verdict='ROBUSTLY_EXCLUDED_MONOJET_ONLY'
        elif rhi<1.0: verdict='SURVIVES_MONOJET_ONLY_WITH_MORPH_UNCERTAINTY'
        else: verdict='MARGINAL_FULL_50GEV_MC_REQUIRED'
        rows.append(dict(model=x.model,convention=x.convention,gq=float(x.gq),BR_invisible=float(x.BR_invisible),signal_sum=float(ns.sum()),mu95=mu95,r_obs=r,morph_uncertainty=morph_unc,r_low=rlo,r_high=rhi,muhat=diag['muhat'],qmu=diag['qmu'],qA=diag['qA'],verdict=verdict))
    R=pd.DataFrame(rows); R.to_csv(OUT/'pkr_50gev_monojet_gate.csv',index=False)

    summary=[]
    for model,q in R.groupby('model'):
        if np.all(q.r_low>1): st='EXCLUDED_FOR_BOTH_WIDTH_CONVENTIONS'
        elif np.all(q.r_high<1): st='SURVIVES_MONOJET_ONLY_BUT_MONOV_AND_FULL_MC_REMAIN'
        else: st='WIDTH_OR_MORPH_UNCERTAINTY_DECIDES_FULL_MC_REQUIRED'
        summary.append(dict(model=model,min_r_low=float(q.r_low.min()),max_r_high=float(q.r_high.max()),status=st))
    S=pd.DataFrame(summary); S.to_csv(OUT/'summary.csv',index=False)

    meta=dict(version='v0.49',analysis='CMS EXO-20-004 independent monojet profile likelihood',statistical_model='Poisson counts + correlated Gaussian background nuisance; direct nuisance profiling; asymptotic CLs',legacy_v048_status='FAILED_VALIDATION_DISCARDED',validation_requirement='100 and 250 GeV monojet-only coupling scale must be within factor 2 of published combined gq limits',morph_uncertainty=morph_unc,limitations=['50 GeV obtained by LO CT18NLO mass morph of public 100 GeV template','mono-V omitted','published gq curve is combined monojet+mono-V, so validation is scale-level not exact equality'])
    (OUT/'metadata.json').write_text(json.dumps(meta,indent=2))

    print('\n=== v0.49 covariance diagnostics ===')
    print(json.dumps(like.covdiag,indent=2))
    print('\n=== v0.49 public likelihood validation ===')
    print(V.to_string(index=False))
    print('\n=== v0.49 100->250 morph closure ===')
    print(f'total fractional error={weighted_total:+.4f}  RMS bin error={rms:.4f}  adopted morph uncertainty={morph_unc:.4f}')
    print('\n=== v0.49 PKR 50 GeV monojet gate ===')
    print(R[['model','convention','r_obs','r_low','r_high','verdict']].to_string(index=False))
    print('\n=== v0.49 summary ===')
    print(S.to_string(index=False))

if __name__=='__main__':
    main()
