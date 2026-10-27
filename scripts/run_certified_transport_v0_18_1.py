#!/usr/bin/env python3
"""Certified v0.18.1 transport production.

Fixes the v0.18 audit findings:
- requires composite numerical certification, including the actual chirality solvers;
- enforces table sidecar provenance and records SHA256 checksums;
- preserves nu and antinu spectra separately for later detector folding;
- never defines transmission outside the incident-flux support;
- summary interpolation never bridges invalid/masked gaps;
- computes each matched SM solution once per model/zenith and reuses it across modes.

Detector weighting is intentionally deferred.
"""
from __future__ import annotations

import argparse, hashlib, json, multiprocessing as mp, os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from nutau_prop.chirality_transport import TauHelicityFlipTransport, TauMuonTransitionTransport
from nutau_prop.momentum_operators import C5_BENCHMARK,C7_BENCHMARK,lambda_4pi_tev

ROOT=Path(__file__).resolve().parents[1]
T11=ROOT/"results"/"v0.11"/"tables"; T15=ROOT/"results"/"v0.15"/"tables"; T16=ROOT/"results"/"v0.16"/"tables"
UHE=ROOT/"data"/"uhe_official"; OUT=ROOT/"results"/"v0.18.1"
MAX_P=0.02; MAX_STEP_KM=12.0; SM_TRANSMISSION_FLOOR=1e-12; EXPECTED_PDF="CT10nlo"


def tag(x): return (f"{float(x):.6g}").replace(".","p")
def sha256(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()

def _npz(p):
    z=np.load(p); return {k:z[k] for k in z.files}

def _meta(p):
    j=p.with_suffix(".json")
    if not j.exists(): raise SystemExit(f"Missing provenance sidecar {j}")
    return json.loads(j.read_text())

def close(a,b,tol=1e-8): return abs(float(a)-float(b))<=tol*max(1.0,abs(float(a)),abs(float(b)))


def contact_path(op,par,s):
    if op==6:
        cands=(T11/f"case6_L{str(float(par)).replace('.', 'p')}_{s}.npz",T11/f"case6_L{tag(par)}_{s}.npz")
    else: cands=(T15/f"case{op}_L{tag(par)}_{s}.npz",)
    for p in cands:
        if p.exists(): return p
    raise SystemExit("Missing table: "+" or ".join(map(str,cands)))
def resolved_path(m,scale,s): return T16/f"resolved_M{tag(m)}_S{tag(scale)}_{s}.npz"


def verify_table(p,kind,par,s,scale=1.0):
    if not p.exists(): raise SystemExit(f"Missing table {p}")
    m=_meta(p)
    if m.get("pdf")!=EXPECTED_PDF: raise SystemExit(f"Unexpected PDF in {p}: {m.get('pdf')}")
    if kind in ("case5","case7"):
        op=int(kind[-1])
        if m.get("version")!="v0.15" or int(m.get("operator",-1))!=op or m.get("cut")!="shat" or not close(m.get("Lambda_UV_TeV"),par):
            raise SystemExit(f"Provenance mismatch for {p}: {m}")
    elif kind=="case6":
        if m.get("version")!="v0.11" or m.get("cut")!="shat" or not close(m.get("Lambda_TeV"),par):
            raise SystemExit(f"Provenance mismatch for {p}: {m}")
    elif kind=="resolved":
        if m.get("version")!="v0.16" or not close(m.get("M_Zprime_TeV"),par) or not close(m.get("lowQ_C4_GeV^-3"),8.6e-10*scale):
            raise SystemExit(f"Resolved provenance mismatch for {p}: {m}")
        if m.get("vertex_cut_TeV") is not None:
            raise SystemExit(f"Refusing resolved table with vertex cutoff: {p} metadata={m}")
    return {"file":str(p.relative_to(ROOT)),"species":s,"sha256":sha256(p),"metadata":m}


def preflight(l5,l7,case6_lams,masses,scale):
    cert=ROOT/"results"/"v0.18.1"/"transport_production_certification_v0_18_1.json"
    if not cert.exists(): raise SystemExit(f"Missing {cert}; run certification first")
    c=json.loads(cert.read_text())
    if not c.get("production_certified",False): raise SystemExit("v0.18.1 production is not certified")
    items=[]
    for op,pars in ((5,l5),(7,l7),(6,case6_lams)):
        for par in pars:
            for s in ("nu","antinu"): items.append(verify_table(contact_path(op,par,s),f"case{op}",par,s,scale))
    for m in masses:
        for s in ("nu","antinu"): items.append(verify_table(resolved_path(m,scale,s),"resolved",m,s,scale))
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"input_provenance_manifest_v0_18_1.json").write_text(json.dumps({"version":"v0.18.1","tables":items},indent=2))


def load_contact(op,par,s): return _npz(contact_path(op,par,s))
def load_resolved(m,scale,s): return _npz(resolved_path(m,scale,s))


def ehlert(E):
    q=pd.read_json(UHE/"flux.json")
    ef=q["Energy [GeV]"].to_numpy(float); sed=q["log10(E^2 F(E) [GeV.cm^-2.s^-1.sr^-1])"].to_numpy(float)
    support=(E>=ef.min())&(E<=ef.max())
    ss=interp1d(ef,sed,bounds_error=False,fill_value=np.nan)(E)
    allf=np.zeros_like(E,float); allf[support]=10.0**ss[support]/E[support]
    return allf,allf/6.0,support,float(ef.min()),float(ef.max())


def ratio_fields(num,den,inc,support):
    raw=(num/den) if den>0 else np.nan
    trans=(den/inc) if (support and inc>0) else np.nan
    masked=raw if (support and np.isfinite(trans) and trans>=SM_TRANSMISSION_FLOOR) else np.nan
    return raw,masked,trans


def base(kind,label,mu,mode,E,allf,support):
    return dict(model_kind=kind,model_label=label,cosz_abs=float(mu),mode=mode,E_GeV=float(E),E_PeV=float(E/1e6),
                in_incident_flux_support=bool(support),incident_allflavor_per_logE=float(allf),
                flux_unit="cm^-2 s^-1 sr^-1 per dlnE")


def case5_pair(task):
    lam,mu=task; tn=load_contact(5,lam,"nu"); ta=load_contact(5,lam,"antinu"); E=tn["energies"]
    allf,one,supp,_,_=ehlert(E); tr=TauMuonTransitionTransport(tn,ta)
    sm=tr.propagate(mu,one,one,one,one,mode="dirac_sink",use_np=False,max_p=MAX_P,max_step_km=MAX_STEP_KM)
    rows=[]
    for mode in ("dirac_sink","majorana_taumu"):
        bs=tr.propagate(mu,one,one,one,one,mode=mode,use_np=True,max_p=MAX_P,max_step_km=MAX_STEP_KM)
        for i,e in enumerate(E):
            smrel=sm["nu_tau"][i]+sm["antinu_tau"][i]+sm["nu_mu"][i]+sm["antinu_mu"][i]
            bsrel=bs["nu_tau"][i]+bs["antinu_tau"][i]+bs["nu_mu"][i]+bs["antinu_mu"][i]
            inc=4*one[i]; rr,rm,trsm=ratio_fields(bsrel,smrel,inc,bool(supp[i]))
            r=base("contact",f"case5_L{lam:.6g}",mu,mode,e,allf[i],supp[i])
            r.update(operator=5,Lambda_UV_TeV=float(lam),incident_relevant_per_logE=float(inc),
                SM_nu_tau=float(sm["nu_tau"][i]),SM_antinu_tau=float(sm["antinu_tau"][i]),SM_nu_mu=float(sm["nu_mu"][i]),SM_antinu_mu=float(sm["antinu_mu"][i]),
                BSM_nu_tau=float(bs["nu_tau"][i]),BSM_antinu_tau=float(bs["antinu_tau"][i]),BSM_nu_mu=float(bs["nu_mu"][i]),BSM_antinu_mu=float(bs["antinu_mu"][i]),
                SM_relevant_per_logE=float(smrel),BSM_relevant_per_logE=float(bsrel),SM_relevant_transmission=trsm,
                BSM_relevant_transmission=(bsrel/inc if supp[i] and inc>0 else np.nan),relevant_ratio_raw=rr,relevant_ratio=rm)
            rows.append(r)
    return rows


def tauonly_pair(task):
    kind,par,mu,scale=task
    if kind=="case6": tn,ta=load_contact(6,par,"nu"),load_contact(6,par,"antinu")
    elif kind=="case7": tn,ta=load_contact(7,par,"nu"),load_contact(7,par,"antinu")
    else: tn,ta=load_resolved(par,scale,"nu"),load_resolved(par,scale,"antinu")
    E=tn["energies"]; allf,one,supp,_,_=ehlert(E); tr=TauHelicityFlipTransport(tn,ta)
    sm=tr.propagate(mu,one,one,mode="dirac_sink",use_np=False,max_p=MAX_P,max_step_km=MAX_STEP_KM,track_secondary_numu=True)
    rows=[]
    for mode in ("dirac_sink","majorana_flip"):
        bs=tr.propagate(mu,one,one,mode=mode,use_np=True,max_p=MAX_P,max_step_km=MAX_STEP_KM,track_secondary_numu=True)
        for i,e in enumerate(E):
            smrel=sm["nu_tau"][i]+sm["antinu_tau"][i]+sm["nu_mu_secondary"][i]+sm["antinu_mu_secondary"][i]
            bsrel=bs["nu_tau"][i]+bs["antinu_tau"][i]+bs["nu_mu_secondary"][i]+bs["antinu_mu_secondary"][i]
            inc=2*one[i]; rr,rm,trsm=ratio_fields(bsrel,smrel,inc,bool(supp[i]))
            label=(f"{kind}_L{par:.6g}" if kind!="resolved" else f"resolved_M{par:.6g}_S{scale:.6g}")
            r=base("resolved_mediator" if kind=="resolved" else "contact",label,mu,mode,e,allf[i],supp[i])
            if kind=="resolved": r.update(M_Zprime_TeV=float(par),lowq_scale=float(scale))
            else: r.update(operator=int(kind[-1]),Lambda_UV_TeV=float(par))
            r.update(incident_relevant_per_logE=float(inc),
                SM_nu_tau=float(sm["nu_tau"][i]),SM_antinu_tau=float(sm["antinu_tau"][i]),
                SM_nu_mu_secondary=float(sm["nu_mu_secondary"][i]),SM_antinu_mu_secondary=float(sm["antinu_mu_secondary"][i]),
                BSM_nu_tau=float(bs["nu_tau"][i]),BSM_antinu_tau=float(bs["antinu_tau"][i]),
                BSM_nu_mu_secondary=float(bs["nu_mu_secondary"][i]),BSM_antinu_mu_secondary=float(bs["antinu_mu_secondary"][i]),
                SM_relevant_per_logE=float(smrel),BSM_relevant_per_logE=float(bsrel),SM_relevant_transmission=trsm,
                BSM_relevant_transmission=(bsrel/inc if supp[i] and inc>0 else np.nan),relevant_ratio_raw=rr,relevant_ratio=rm)
            rows.append(r)
    return rows


def local_interp(q,target,col,require_support=False):
    q=q.sort_values("E_PeV"); x=q.E_PeV.to_numpy(float); y=q[col].to_numpy(float); sup=q.in_incident_flux_support.to_numpy(bool)
    j=int(np.searchsorted(x,target))
    if j==0 or j>=len(x): return np.nan
    i=j-1
    if require_support and not (sup[i] and sup[j]): return np.nan
    if not (np.isfinite(y[i]) and np.isfinite(y[j])): return np.nan
    if y[i]>0 and y[j]>0:
        t=(np.log(target)-np.log(x[i]))/(np.log(x[j])-np.log(x[i])); return float(np.exp(np.log(y[i])+t*(np.log(y[j])-np.log(y[i]))))
    t=(target-x[i])/(x[j]-x[i]); return float(y[i]+t*(y[j]-y[i]))


def summarize(D):
    rows=[]; groups=["model_kind","model_label","cosz_abs","mode"]
    species_cols=[c for c in D.columns if c.startswith("SM_nu") or c.startswith("SM_antinu") or c.startswith("BSM_nu") or c.startswith("BSM_antinu")]
    for keys,q in D.groupby(groups,dropna=False):
        kind,label,mu,mode=keys; s=dict(model_kind=kind,model_label=label,cosz_abs=mu,mode=mode)
        energies=(0.3,1.,10.,100.,220.,300.)+((1000.,10000.) if kind=="resolved_mediator" else ())
        for ep in energies:
            for col in ["incident_allflavor_per_logE","SM_relevant_per_logE","BSM_relevant_per_logE"]+species_cols:
                if col in q: s[f"{col}_{ep:g}PeV"]=local_interp(q,ep,col,require_support=False)
            for col in ("SM_relevant_transmission","BSM_relevant_transmission","relevant_ratio"):
                s[f"{col}_{ep:g}PeV"]=local_interp(q,ep,col,require_support=True)
        rows.append(s)
    return pd.DataFrame(rows).sort_values(groups)


def pmap(fn,tasks,n):
    ctx=mp.get_context("spawn")
    with ctx.Pool(processes=min(max(1,n),len(tasks))) as pool: z=list(pool.imap_unordered(fn,tasks,chunksize=1))
    return [r for a in z for r in a]


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--workers",type=int,default=int(os.environ.get("NUTAU_NPROC","12")))
    ap.add_argument("--cosz",nargs="+",type=float,default=[0.05,0.1,0.2,0.5,0.8,1.0]); ap.add_argument("--lowq-scale",type=float,default=1.0)
    a=ap.parse_args(); l5=[3.0,lambda_4pi_tev(C5_BENCHMARK)]; l7=[3.0,lambda_4pi_tev(C7_BENCHMARK)]; l6=[1.5,2.0,2.44,3.0]; masses=[0.1,0.2,0.5,1.0,2.0,3.0]
    preflight(l5,l7,l6,masses,a.lowq_scale)
    c5=pmap(case5_pair,[(l,m) for l in l5 for m in a.cosz],a.workers)
    other=[]
    tasks=[("case6",l,m,a.lowq_scale) for l in l6 for m in a.cosz]+[("case7",l,m,a.lowq_scale) for l in l7 for m in a.cosz]
    other=pmap(tauonly_pair,tasks,a.workers)
    med=pmap(tauonly_pair,[("resolved",x,m,a.lowq_scale) for x in masses for m in a.cosz],a.workers)
    C=pd.DataFrame(c5+other); M=pd.DataFrame(med)
    C.to_csv(OUT/"contact_transport_certified_v0_18_1.csv",index=False); M.to_csv(OUT/"resolved_transport_certified_v0_18_1.csv",index=False)
    CS=summarize(C); MS=summarize(M); CS.to_csv(OUT/"contact_transport_certified_summary_v0_18_1.csv",index=False); MS.to_csv(OUT/"resolved_transport_certified_summary_v0_18_1.csv",index=False)
    q=pd.read_json(UHE/"flux.json"); support=[float(q["Energy [GeV]"].min()),float(q["Energy [GeV]"].max())]
    status={"version":"v0.18.1","numerically_certified":True,"detector_weighting_done":False,
            "incident_flux_support_GeV":support,"transmission_outside_support":"undefined/NaN",
            "nu_antinu_exported_separately":True,"SM_solution_reused_across_modes":True,
            "caution":["dirac_sink is an absorbing-state bracket, not a complete RH-state transport model","resolved mediator does not validate the dimension-five nu-Zprime vertex or laboratory constraints"]}
    (OUT/"production_status_v0_18_1.json").write_text(json.dumps(status,indent=2))
    print("CONTACT SUMMARY\n",CS.to_string(index=False,float_format=lambda x:f"{x:.4g}")); print("\nRESOLVED SUMMARY\n",MS.to_string(index=False,float_format=lambda x:f"{x:.4g}")); print("WROTE",OUT)
if __name__=="__main__": main()
