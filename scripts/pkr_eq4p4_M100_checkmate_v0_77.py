#!/usr/bin/env python3
"""v0.77 — exactly one new-mass Eq4.4 CheckMATE validation at M=100 GeV.

Reads the v0.76 best robust-dijet-screen point at 100 GeV, uses the physical
O+h.c. width fixed in v0.73, generates one CT18NLO sample, and runs the same
CMS EXO-20-004 multibin CheckMATE analysis used for v0.75.

This is a validation benchmark, not a scan.  We do NOT transfer the 50-GeV
monojet r-value across masses.  The previously certified MG5/CT18/NWA machinery
is reused; at this new mass we require stable persisted parameters, exact
banner/LHE agreement, LHAID 14400, narrow width, and the v0.76 dijet screen.
"""
from __future__ import annotations

import gzip, importlib.util, json, math, os, re, shutil, subprocess, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "cache" / "v0.77_eq4p4_M100_checkmate"
OUT = ROOT / "results" / "v0.77_eq4p4_M100_checkmate"
SUM76 = ROOT / "results" / "v0.76_eq4p4_mass_prescreen" / "mass_summary.csv"
TAG = "eq4p4_derivative__M100_best"
ANALYSIS = os.environ.get("PKR_CM_ANALYSIS", "cms_2107_13201")
NCORE = int(os.environ.get("PKR_MG5_NCORE", "14"))

P67 = ROOT / "scripts" / "pkr_ct18_smoke_v0_67.py"
s67 = importlib.util.spec_from_file_location("v067", P67)
if s67 is None or s67.loader is None: raise RuntimeError(f"cannot import {P67}")
v067 = importlib.util.module_from_spec(s67); s67.loader.exec_module(v067)
base, cfg = v067.v065.base, v067.v065.cfg

P51 = ROOT / "scripts" / "pkr_run_checkmate_v0_51.py"
s51 = importlib.util.spec_from_file_location("cm51", P51)
if s51 is None or s51.loader is None: raise RuntimeError(f"cannot import {P51}")
cm51 = importlib.util.module_from_spec(s51); s51.loader.exec_module(cm51)

P69 = ROOT / "scripts" / "pkr_checkmate_multibin_retry_v0_69.py"
s69 = importlib.util.spec_from_file_location("v069", P69)
if s69 is None or s69.loader is None: raise RuntimeError(f"cannot import {P69}")
v069 = importlib.util.module_from_spec(s69); s69.loader.exec_module(v069)

QUARK_MASS = {"u":0.00216,"d":0.00467,"s":0.0934,"c":1.27,"b":4.18,"t":172.69}


def gamma_qq(M,gq,mq):
    if M <= 2*mq: return 0.0
    r=(mq/M)**2
    return 3*gq*gq*M/(12*math.pi)*(1+2*r)*math.sqrt(1-4*r)


def setopt(text,key,value):
    pat=re.compile(rf"(?m)^\s*#?\s*{re.escape(key)}\s*=.*$")
    line=f"{key} = {value}"
    return pat.sub(line,text,count=1) if pat.search(text) else text.rstrip()+"\n"+line+"\n"


def configure_multicore(procdir):
    p=procdir/"Cards"/"me5_configuration.txt"
    t=p.read_text(errors="replace")
    t=setopt(t,"run_mode","2" if NCORE>1 else "0")
    t=setopt(t,"nb_core",str(NCORE)); p.write_text(t)
    print(f"MG5 multicore: {NCORE} workers")


def run(cmd,log=None,allow_nonzero=False):
    p=subprocess.run([str(x) for x in cmd],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if log: Path(log).write_text(p.stdout)
    if p.returncode and not allow_nonzero:
        raise RuntimeError("command failed rc=%d\n%s"%(p.returncode,"\n".join(p.stdout.splitlines()[-60:])))
    return p.returncode,p.stdout


def plain_lhe(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True)
    if src.suffix==".gz":
        with gzip.open(src,"rb") as fi,dst.open("wb") as fo: shutil.copyfileobj(fi,fo)
    else: shutil.copy2(src,dst)
    return dst


def write_cm_card(hepmc,xsec):
    card=WORK/"cards"/f"checkmate_{TAG}.in"; cmout=WORK/"checkmate"/TAG
    card.parent.mkdir(parents=True,exist_ok=True)
    card.write_text(f"""[Parameters]
Backend: numpy
Name: pkr_v077_{TAG}
Analyses: {ANALYSIS}
OutputDirectory: {cmout}
OutputExists: Overwrite
SkipParamCheck: True
EffTab: True
Multibin: detailed
Expected: True
Uplim: True
Model: simple
InvisiblePIDs: 5000521
RandomSeed: 7701
TotalResultFileColumns: analysis, sr, o, b, db, s, ds, s95obs, s95exp, robscons, rexpcons
BestPerAnalysisResultFileColumns: analysis, sr, o, b, db, s, ds, s95obs, s95exp, robscons, rexpcons

[signal]
Events: {hepmc}
XSect: {xsec:.12g} PB
XSectErr: 10 %%
""")
    return card,cmout


def main():
    OUT.mkdir(parents=True,exist_ok=True); WORK.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(SUM76)
    q=d[(d.M_GeV==100.0) & (d.status=="ROBUST_DIJET_SCREEN_WINDOW")]
    if len(q)!=1: raise RuntimeError(f"expected one M100 summary row, found {len(q)}")
    b=q.iloc[0]
    M=float(b.M_GeV); gq=float(b.gq); a=float(b.a_GeV_inv)
    gamma_inv=a*a*M**3/(24*math.pi)
    gamma_q=sum(gamma_qq(M,gq,mq) for mq in QUARK_MASS.values())
    gamma_tot=gamma_q+gamma_inv; br=gamma_inv/gamma_tot
    gchi=cfg.gchi_for_width(M,cfg.MCHI,gamma_inv)
    gqeq=gq*math.sqrt(1-br)

    proc=WORK/"mg5"/TAG
    if proc.exists(): shutil.rmtree(proc)
    build=WORK/"cards"/f"build_{TAG}.txt"; build.parent.mkdir(parents=True,exist_ok=True)
    base.build_process(proc,build,OUT/"mg5_build.log")
    param=proc/"Cards"/"param_card.dat"; rcard=proc/"Cards"/"run_card.dat"
    base.edit_param_card(param,M,cfg.MCHI,gq,gchi,gamma_tot)
    base.configure_run_card(rcard,seed=7701)
    ok0,*_=base.certify_persisted_param_card(param,M,cfg.MCHI,gq,gchi,gamma_tot)
    if not ok0: raise RuntimeError("parameter card failed before generation")
    v067.patch_lhapdf_lambda_fallback_v067(proc); configure_multicore(proc)
    base.generate_events(proc,OUT/"mg5_generate.log")
    ok1,actual_M,actual_mchi,actual_gchi,actual_width,actual_gqs=base.certify_persisted_param_card(param,M,cfg.MCHI,gq,gchi,gamma_tot)
    if not ok1: raise RuntimeError("parameter card changed during generation")

    lhe_src=cfg.find_lhe(proc); banner=cfg.find_banner(proc)
    xsec=cfg.parse_lhe_init_xsec_pb(lhe_src); xb=cfg.parse_banner_xsec_pb(banner)
    banner_rel=abs(xb-xsec)/max(abs(xsec),1e-30)
    has14400=bool(re.search(r"(?im)^\s*14400\s*=\s*lhaid\b",banner.read_text(errors="replace")))
    if banner_rel>1e-10: raise RuntimeError(f"banner/LHE mismatch {banner_rel:.3e}")
    if not has14400: raise RuntimeError("LHAID 14400 absent from banner")

    lhe=plain_lhe(lhe_src,WORK/"lhe"/f"{TAG}.lhe")
    shower=Path(os.environ.get("PKR_SHOWER",ROOT/"cache"/"v0.51"/"bin"/"pythia_lhe_to_hepmc_v0_51"))
    checkmate=Path(os.environ.get("PKR_CHECKMATE",ROOT/"cache"/"v0.51"/"checkmate-install"/"bin"/"CheckMATE"))
    hepmc=WORK/"hepmc"/f"{TAG}.hepmc"; hepmc.parent.mkdir(parents=True,exist_ok=True)
    run([shower,lhe,hepmc,"7701"],OUT/"pythia.log")
    if not hepmc.exists() or hepmc.stat().st_size<1000: raise RuntimeError("no usable HepMC")

    for mod in ("spey","spey_pyhf","pyhf"): __import__(mod)
    v069.patch_optional_jax()
    card,cmout=write_cm_card(hepmc,xsec)
    if cmout.exists(): shutil.rmtree(cmout)
    rc,stdout=run([checkmate,card],OUT/"checkmate.log",allow_nonzero=True)
    mu95,rmb,cls=cm51.parse_multibin(cmout); rsingle=cm51.parse_singlebin(cmout)
    if mu95 is None:
        vals=re.findall(r"Upper limit:\s*([0-9.eE+-]+)",stdout)
        if vals:
            mu95=float(vals[-1]); rmb=1.0/mu95 if mu95>0 else None
    if rmb is None and rsingle is None: raise RuntimeError(f"no usable CheckMATE result rc={rc}")

    if rmb is None: status="MULTIBIN_UNRESOLVED"
    elif rmb>=1.2: status="ROBUST_EXCLUSION_THRESHOLD"
    elif rmb<=0.8: status="ROBUST_SURVIVAL_THRESHOLD"
    else: status="MARGINAL_BOUNDARY_UNRESOLVED_VALIDATION_ONLY"

    result={
      "version":"v0.77_eq4p4_M100_checkmate",
      "benchmark":{"M_GeV":M,"gq":gq,"a_GeV_inv":a,"kappa_GeV_inv":gq*a,
                   "Lambda_vertex_GeV":1/a,"Gamma_inv_GeV":gamma_inv,
                   "Gamma_total_GeV":gamma_tot,"Gamma_over_M":gamma_tot/M,
                   "BR_invisible":br,"gq_equiv_dijet":gqeq,
                   "cms_gq95_digitized":float(b.cms_gq95_digitized)},
      "generator":{"cross_section_pb":xsec,"banner_vs_lhe_rel":banner_rel,
                   "lhaid14400":has14400,"nevents_requested":int(os.environ.get("PKR_NORM_NEVENTS","100000")),
                   "mg5_workers":NCORE,"param_card_before":bool(ok0),"param_card_after":bool(ok1)},
      "checkmate":{"analysis":ANALYSIS,"returncode":rc,"multibin_mu95":mu95,
                   "multibin_r":rmb,"multibin_CLs_mu1":cls,"max_singlebin_robscons":rsingle,
                   "status":status},
      "interpretation":"This is the first exact monojet validation at M=100 GeV. Do not use the 50-GeV r-value as a cross-mass prediction. Values 0.8<r<1.2 remain boundary-unresolved. The dijet number is still a digitized screen and needs exact mass-specific confirmation if monojet survives."
    }
    (OUT/"M100_checkmate.json").write_text(json.dumps(result,indent=2)+"\n")
    print("=== v0.77 exact Eq4.4 M=100 CheckMATE ===")
    print(f"gq={gq:.6g}  a={a:.8g} GeV^-1  BRinv={br:.6f}  gq*sqrt(BRqq)={gqeq:.6f}")
    print(f"Gamma/M={gamma_tot/M:.6g}  CT18 xsec={xsec:.8g} pb  banner/LHE={banner_rel:.2e}")
    print(f"CheckMATE rc={rc}  mu95={mu95}  multibin_r={rmb}  singlebin={rsingle}")
    print(f"status: {status}")
    return 0

if __name__=="__main__": sys.exit(main())
