#!/usr/bin/env python3
"""Independent finite-beta_tau validation for repaired v0.17 transport.

No PDF integration is performed. Existing v0.11 Case-6 tables are reused.
The gates are intentionally physics/numerics gates, not detector claims.
"""
from __future__ import annotations

import argparse, json, multiprocessing as mp, os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad

from nutau_prop.constants import M_TAU, C_TAU_CM
from nutau_prop.tau_evolution import CoupledTauEvolution
from nutau_prop.transport_case6 import TauTransportFlavorAware
from nutau_prop.util import powerlaw_per_logE

ROOT = Path(__file__).resolve().parents[1]
T11 = ROOT / "results" / "v0.11" / "tables"
OUT = ROOT / "results" / "v0.17"


def tag(x):
    return (f"{float(x):g}").replace(".", "p")


def load_case6(lam, species):
    z = np.load(T11 / f"case6_L{tag(lam)}_{species}.npz")
    return {k: z[k] for k in z.files}


def sm_view(t):
    return dict(energies=t["energies"], y=t["y"], dy=t["dy"],
                sigma_cc=t["sigma_cc"], sigma_nc=t["sigma_nc_sm"],
                p_cc=t["p_cc"], p_nc=t["p_nc_sm"])


def bsm_view(t):
    return dict(energies=t["energies"], y=t["y"], dy=t["dy"],
                sigma_cc=t["sigma_cc"], sigma_nc=t["sigma_nc"],
                p_cc=t["p_cc"], p_nc=t["p_nc"])


def logmean(E, w):
    w = np.asarray(w, float); s = float(np.sum(w))
    return np.nan if s <= 0 else float(np.exp(np.sum(w*np.log(E))/s))


def primitive_tests():
    E = np.geomspace(1e4, 1e10, 241)
    evo = CoupledTauEvolution(E, beta_tau=0.85e-6, nquad=32)
    rows, gates = [], []

    for E0 in (1e6, 1e7, 1e8, 1e9):
        j = int(np.argmin(np.abs(np.log(E/E0)))); E0g = float(E[j])
        for rho in (1.0, 3.0, 8.0):
            for ds_km in (0.05, 0.5, 5.0, 20.0):
                ds = ds_km*1e5; tau = np.zeros_like(E); tau[j] = 1.0
                tout, dpar = evo.advance(tau, rho, ds)
                S, Eend, U = evo.analytic_parent(E0g, rho, ds)
                sout = float(np.sum(tout)); pdec = float(np.sum(dpar))
                got_end = logmean(E, tout)
                b = evo.beta*rho; a = M_TAU/(C_TAU_CM*E0g)

                if pdec > 1e-14:
                    def f(t):
                        if b == 0: return np.log(E0g)
                        u = -np.log1p(-t)
                        return np.log(E0g/(1.0+(b/a)*u))
                    ref_log = quad(f, 0.0, pdec, epsabs=1e-11,
                                   epsrel=1e-11, limit=200)[0]/pdec
                    got_log = float(np.sum(dpar*np.log(E))/max(pdec,1e-300))
                    decay_logerr = abs(got_log-ref_log)
                else:
                    decay_logerr = 0.0

                surv_err = abs(sout-S)
                dec_err = abs(pdec-(1.0-S))
                end_logerr = (abs(np.log(got_end/Eend))
                              if S>1e-12 and np.isfinite(got_end) else 0.0)
                rows.append(dict(E0_GeV=E0g, rho=rho, ds_km=ds_km, U=U,
                                 survival_exact=S, survival_discrete=sout,
                                 decay_exact=1-S, decay_discrete=pdec,
                                 survivor_logE_error=end_logerr,
                                 decay_mean_logE_error=decay_logerr))
                for name,val,tol in (
                    ("survival",surv_err,5e-12),
                    ("probability",dec_err,5e-12),
                    ("survivor_energy",end_logerr,2e-12),
                    ("decay_energy_moment",decay_logerr,3e-4),
                ):
                    gates.append(dict(group="primitive",test=name,value=val,
                                      tolerance=tol,PASS=val<tol))

    # Exact beta_tau=0 regression.
    evo0 = CoupledTauEvolution(E, beta_tau=0.0, nquad=8)
    for E0 in (1e6,1e8):
        j=int(np.argmin(np.abs(np.log(E/E0)))); E0g=float(E[j])
        tau=np.zeros_like(E); tau[j]=1.0
        for ds_km in (0.05,5.0):
            ds=ds_km*1e5; tout,dpar=evo0.advance(tau,5.0,ds)
            S=np.exp(-(M_TAU/(C_TAU_CM*E0g))*ds)
            err=max(abs(np.sum(tout)-S),abs(np.sum(dpar)-(1-S)),
                    abs(logmean(E,tout)/E0g-1) if np.sum(tout)>1e-12 else 0,
                    abs(logmean(E,dpar)/E0g-1) if np.sum(dpar)>1e-12 else 0)
            gates.append(dict(group="beta0",test="analytic_regression",value=err,
                              tolerance=5e-12,PASS=err<5e-12))

    # Semigroup/refinement check. The 3e-3 tolerance allows only the known
    # log-grid interpolation diffusion when the surviving population is tiny.
    for E0 in (1e7,1e8,1e9):
        j=int(np.argmin(np.abs(np.log(E/E0)))); tau=np.zeros_like(E); tau[j]=1.0
        for rho,ds_km in ((3.0,5.0),(8.0,20.0)):
            ds=ds_km*1e5
            ft,fd=evo.advance(tau,rho,ds)
            ht,d1=evo.advance(tau,rho,0.5*ds)
            tt,d2=evo.advance(ht,rho,0.5*ds)
            pop=max(abs(np.sum(ft)-np.sum(tt)),abs(np.sum(fd)-np.sum(d1+d2)))
            eerr=0.0
            if np.sum(ft)>1e-8 and np.sum(tt)>1e-8:
                eerr=abs(np.log(logmean(E,ft)/logmean(E,tt)))
            val=max(pop,eerr)
            gates.append(dict(group="semigroup",test="one_vs_two_half",value=val,
                              tolerance=3e-3,PASS=val<3e-3))

    return pd.DataFrame(rows), pd.DataFrame(gates)


def meaningful_rel_l1(a,b,incident):
    a=np.asarray(a,float); b=np.asarray(b,float); incident=np.asarray(incident,float)
    scale=max(float(np.max(np.abs(b))),1e-300)
    iscale=max(float(np.max(np.abs(incident))),1e-300)
    mask=(np.abs(b)>1e-10*scale)&(incident>1e-12*iscale)
    if not np.any(mask): return 0.0
    return float(np.sum(np.abs(a[mask]-b[mask]))/max(np.sum(np.abs(b[mask])),1e-300))


def _full_case(task):
    lam,mu,bsm=task; results={}; p0_ref=None
    for max_p in (0.04,0.02,0.01):
        tsum=None; msum=None
        for species in ("nu","antinu"):
            raw=load_case6(lam,species); base=sm_view(raw)
            primary=bsm_view(raw) if bsm else base; E=raw["energies"]
            p0=powerlaw_per_logE(E,2.0); p0/=max(float(np.max(p0)),1e-300)
            P=-1.0 if species=="nu" else +1.0
            tr=TauTransportFlavorAware(primary,base,beta_tau=0.85e-6,
                                       tau_polarization=P,secondary_branching=0.18)
            out=tr.propagate(mu,p0,max_step_km=12.0,max_p=max_p,
                             track_secondary_numu=True)
            tsum=out["nu"].copy() if tsum is None else tsum+out["nu"]
            msum=out["numu"].copy() if msum is None else msum+out["numu"]
            p0_ref=2*p0
        results[max_p]=(tsum,msum)

    tsum6=None; msum6=None
    for species in ("nu","antinu"):
        raw=load_case6(lam,species); base=sm_view(raw)
        primary=bsm_view(raw) if bsm else base; E=raw["energies"]
        p0=powerlaw_per_logE(E,2.0); p0/=max(float(np.max(p0)),1e-300)
        P=-1.0 if species=="nu" else +1.0
        tr=TauTransportFlavorAware(primary,base,beta_tau=0.85e-6,
                                   tau_polarization=P,secondary_branching=0.18)
        out=tr.propagate(mu,p0,max_step_km=6.0,max_p=0.01,
                         track_secondary_numu=True)
        tsum6=out["nu"].copy() if tsum6 is None else tsum6+out["nu"]
        msum6=out["numu"].copy() if msum6 is None else msum6+out["numu"]

    d40=meaningful_rel_l1(results[0.04][0],results[0.02][0],p0_ref)
    d20=meaningful_rel_l1(results[0.02][0],results[0.01][0],p0_ref)
    dm20=meaningful_rel_l1(results[0.02][1],results[0.01][1],p0_ref)
    ds=meaningful_rel_l1(results[0.01][0],tsum6,p0_ref)
    dsm=meaningful_rel_l1(results[0.01][1],msum6,p0_ref)
    arrays=[x for pair in results.values() for x in pair]+[tsum6,msum6]
    return dict(Lambda_TeV=lam,cosz_abs=mu,bsm=bool(bsm),
                d40_20_tau=d40,d20_10_tau=d20,d20_10_mu=dm20,
                matter_step_tau=ds,matter_step_mu=dsm,
                min_flux=min(float(np.min(x)) for x in arrays),
                finite=all(np.all(np.isfinite(x)) for x in arrays))


def full_transport_tests(workers):
    missing=[T11/f"case6_L{tag(l)}_{s}.npz" for l in (1.5,2.44)
             for s in ("nu","antinu") if not (T11/f"case6_L{tag(l)}_{s}.npz").exists()]
    if missing: raise SystemExit("Missing v0.11 tables: "+", ".join(map(str,missing)))
    tasks=[(l,m,b) for l in (1.5,2.44) for m in (0.5,1.0) for b in (False,True)]
    ctx=mp.get_context("spawn")
    with ctx.Pool(processes=min(max(1,workers),len(tasks))) as pool:
        D=pd.DataFrame(list(pool.imap_unordered(_full_case,tasks,chunksize=1)))
    D=D.sort_values(["Lambda_TeV","cosz_abs","bsm"])
    gates=[]
    for _,r in D.iterrows():
        checks=(("max_p_tau",r.d20_10_tau,0.01),
                ("max_p_secondary_mu",r.d20_10_mu,0.015),
                ("matter_step_tau",r.matter_step_tau,0.01),
                ("matter_step_secondary_mu",r.matter_step_mu,0.015))
        for name,val,tol in checks:
            gates.append(dict(group="full_transport",test=name,value=val,
                              tolerance=tol,PASS=val<tol))
        gates.append(dict(group="full_transport",test="positivity_finite",
                          value=min(r.min_flux,0.0),tolerance=0.0,
                          PASS=bool(r.finite and r.min_flux>=-1e-14)))
        ratio=r.d20_10_tau/max(r.d40_20_tau,1e-300)
        gates.append(dict(group="full_transport",test="refinement_trend_tau",
                          value=ratio,tolerance=1.05,
                          PASS=bool(r.d20_10_tau<=1.05*r.d40_20_tau+1e-8)))
    return D,pd.DataFrame(gates)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--workers",type=int,default=int(os.environ.get("NUTAU_NPROC","12")))
    a=ap.parse_args(); OUT.mkdir(parents=True,exist_ok=True)
    prim,g1=primitive_tests(); full,g2=full_transport_tests(a.workers)
    G=pd.concat([g1,g2],ignore_index=True)
    prim.to_csv(OUT/"tau_primitive_validation_v0_17.csv",index=False)
    full.to_csv(OUT/"tau_full_transport_convergence_v0_17.csv",index=False)
    G.to_csv(OUT/"tau_validation_gates_v0_17.csv",index=False)
    status=dict(version="v0.17",
                finite_beta_tau_repair="exact constant-density loss+decay plus Strang-style splitting",
                all_gates_PASS=bool(G.PASS.all()),n_gates=int(len(G)),
                n_failed=int((~G.PASS).sum()),
                note="Historical beta=0 RR/DRS gates are run separately by run_tau_transport_validation_v0_17.sh.")
    (OUT/"tau_validation_status_v0_17.json").write_text(json.dumps(status,indent=2))
    print("GATES"); print(G.to_string(index=False,float_format=lambda x:f"{x:.5g}"))
    print("\nFULL TRANSPORT"); print(full.to_string(index=False,float_format=lambda x:f"{x:.5g}"))
    print("\nSTATUS",json.dumps(status,indent=2)); print("WROTE",OUT)
    if not status["all_gates_PASS"]: raise SystemExit("v0.17 finite-beta_tau validation FAILED")


if __name__=="__main__": main()
