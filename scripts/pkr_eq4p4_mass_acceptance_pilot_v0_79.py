#!/usr/bin/env python3
"""v0.79 — Cross-mass production x acceptance pilot for PKR Eq.4.4.

Purpose
-------
Correct the v0.76 cross-mass ranking mistake.  For each mass use the v0.78
representative point (the valid >=5%-UHE point with the smallest gq^2 BRinv at
that mass), and compare the actual CMS monojet selected-bin vector

  s_b = L * [gq^2 BRinv] * sigma_prod^(gq=1)(M) * A_b(M).

50 GeV reuses the existing v0.75 100k detector sample.  Only 100 and 220 GeV
are newly simulated, with a small pilot sample (default 20k events each).
This is a ranking diagnostic only: no exclusion/survival conclusion is allowed
from the pilot.

Outputs include:
  * sigma/(gq^2 BRinv), the mass-dependent unit production factor;
  * per-bin K_b = s_b/(gq^2 BRinv), which contains production x acceptance;
  * the representative expected s_b vector;
  * a diagonal background-weighted ranking score sqrt(sum s_b^2/(b+db^2));
  * a pilot MC-scale diagnostic and any CheckMATE multibin r as secondary info.

The diagonal score ignores bin correlations and is used only to decide whether
one mass is clearly preferable for a later full-statistics validation.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.79_eq4p4_mass_acceptance_pilot"
WORK = ROOT / "cache" / "v0.79_eq4p4_mass_acceptance_pilot"
SUM78 = ROOT / "results" / "v0.78_eq4p4_full_envelope" / "mass_summary.csv"
REF75 = ROOT / "results" / "v0.75_eq4p4_bestpoint_checkmate" / "bestpoint_checkmate.json"
CM75 = ROOT / "cache" / "v0.75_eq4p4_bestpoint_checkmate" / "checkmate"
ANALYSIS = "cms_2107_13201"
NCORE = int(os.environ.get("PKR_MG5_NCORE", "14"))
XSEC_REL = 0.10
NEW_MASSES = (100.0, 220.0)
QUARK_MASS = {"u":0.00216,"d":0.00467,"s":0.0934,"c":1.27,"b":4.18,"t":172.69}

# Certified MG5/CT18 machinery.
P67 = ROOT / "scripts" / "pkr_ct18_smoke_v0_67.py"
s67 = importlib.util.spec_from_file_location("v067", P67)
if s67 is None or s67.loader is None:
    raise RuntimeError(f"cannot import {P67}")
v067 = importlib.util.module_from_spec(s67); s67.loader.exec_module(v067)
base, cfg = v067.v065.base, v067.v065.cfg

# CheckMATE parser + Spey/JAX compatibility patch.
P51 = ROOT / "scripts" / "pkr_run_checkmate_v0_51.py"
s51 = importlib.util.spec_from_file_location("cm51", P51)
if s51 is None or s51.loader is None:
    raise RuntimeError(f"cannot import {P51}")
cm51 = importlib.util.module_from_spec(s51); s51.loader.exec_module(cm51)
P69 = ROOT / "scripts" / "pkr_checkmate_multibin_retry_v0_69.py"
s69 = importlib.util.spec_from_file_location("v069", P69)
if s69 is None or s69.loader is None:
    raise RuntimeError(f"cannot import {P69}")
v069 = importlib.util.module_from_spec(s69); s69.loader.exec_module(v069)


def gamma_qq(M: float, gq: float, mq: float) -> float:
    if M <= 2*mq:
        return 0.0
    r=(mq/M)**2
    return 3*gq*gq*M/(12*math.pi)*(1+2*r)*math.sqrt(1-4*r)


def setopt(text: str, key: str, value: str) -> str:
    pat=re.compile(rf"(?m)^\s*#?\s*{re.escape(key)}\s*=.*$")
    line=f"{key} = {value}"
    return pat.sub(line,text,count=1) if pat.search(text) else text.rstrip()+"\n"+line+"\n"


def configure_multicore(procdir: Path) -> None:
    p=procdir/"Cards"/"me5_configuration.txt"
    t=p.read_text(errors="replace")
    t=setopt(t,"run_mode","2" if NCORE>1 else "0")
    t=setopt(t,"nb_core",str(NCORE))
    p.write_text(t)


def run(cmd, log: Path | None=None, allow_nonzero=False):
    p=subprocess.run([str(x) for x in cmd],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if log is not None:
        log.parent.mkdir(parents=True,exist_ok=True)
        log.write_text(p.stdout)
    if p.returncode and not allow_nonzero:
        raise RuntimeError("command failed rc=%d\n%s"%(p.returncode,"\n".join(p.stdout.splitlines()[-60:])))
    return p.returncode,p.stdout


def plain_lhe(src: Path,dst: Path) -> Path:
    dst.parent.mkdir(parents=True,exist_ok=True)
    if src.suffix==".gz":
        with gzip.open(src,"rb") as fi,dst.open("wb") as fo:
            shutil.copyfileobj(fi,fo)
    else:
        shutil.copy2(src,dst)
    return dst


def find_total(root: Path) -> Path:
    hits=list(root.rglob("total_results.txt"))
    if not hits:
        raise RuntimeError(f"no total_results.txt under {root}")
    return sorted(hits,key=lambda p:p.stat().st_mtime)[-1]


def parse_total(path: Path) -> dict[str,dict[str,float]]:
    rows={}
    lines=path.read_text(errors="replace").splitlines()
    for line in lines[1:]:
        p=line.split()
        if len(p)<7 or p[0]!=ANALYSIS:
            continue
        try:
            rows[p[1]]={"o":float(p[2]),"b":float(p[3]),"db":float(p[4]),"s":float(p[5]),"ds":float(p[6])}
        except ValueError:
            pass
    if not rows:
        raise RuntimeError(f"no {ANALYSIS} rows in {path}")
    return rows


def write_card(tag: str, hepmc: Path, xsec: float, seed: int) -> tuple[Path,Path]:
    card=WORK/"cards"/f"checkmate_{tag}.in"
    cmout=WORK/"checkmate"/tag
    card.parent.mkdir(parents=True,exist_ok=True)
    card.write_text(f"""[Parameters]
Backend: numpy
Name: pkr_v079_{tag}
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
RandomSeed: {seed}
TotalResultFileColumns: analysis, sr, o, b, db, s, ds, s95obs, s95exp, robscons, rexpcons
BestPerAnalysisResultFileColumns: analysis, sr, o, b, db, s, ds, s95obs, s95exp, robscons, rexpcons

[signal]
Events: {hepmc}
XSect: {xsec:.12g} PB
XSectErr: 10 %%
""")
    return card,cmout


def representative_points() -> dict[float,dict[str,float]]:
    if not SUM78.exists():
        raise RuntimeError(f"missing {SUM78}; rerun corrected v0.78 first")
    s=pd.read_csv(SUM78)
    out={}
    for M in (50.0,100.0,220.0):
        q=s[np.isclose(s.M_GeV,M)]
        if len(q)!=1 or not str(q.iloc[0].status).startswith("FULL_A_ROBUST"):
            raise RuntimeError(f"mass {M} has no v0.78 robust-dijet representative")
        r=q.iloc[0]
        out[M]={"gq":float(r.representative_gq),"a":float(r.representative_a),"BR_inv":float(r.representative_BR_inv)}
    return out


def generate_mass(M: float, pt: dict[str,float], seed: int):
    gq=pt["gq"]; a=pt["a"]
    gamma_inv=a*a*M**3/(24*math.pi)
    gamma_q=sum(gamma_qq(M,gq,mq) for mq in QUARK_MASS.values())
    gamma_tot=gamma_q+gamma_inv
    br=gamma_inv/gamma_tot
    if abs(br-pt["BR_inv"])/max(br,1e-30)>0.01:
        raise RuntimeError(f"M={M}: v0.78 BR mismatch: recalculated {br}, stored {pt['BR_inv']}")
    gchi=cfg.gchi_for_width(M,cfg.MCHI,gamma_inv)
    tag=f"M{int(M)}"
    proc=WORK/"mg5"/tag
    if proc.exists(): shutil.rmtree(proc)
    build=WORK/"cards"/f"build_{tag}.txt"; build.parent.mkdir(parents=True,exist_ok=True)
    base.build_process(proc,build,OUT/f"mg5_build_{tag}.log")
    param=proc/"Cards"/"param_card.dat"; rcard=proc/"Cards"/"run_card.dat"
    base.edit_param_card(param,M,cfg.MCHI,gq,gchi,gamma_tot)
    base.configure_run_card(rcard,seed=seed)
    ok0,*_=base.certify_persisted_param_card(param,M,cfg.MCHI,gq,gchi,gamma_tot)
    if not ok0: raise RuntimeError(f"M={M}: parameter card failed before generation")
    v067.patch_lhapdf_lambda_fallback_v067(proc)
    configure_multicore(proc)
    base.generate_events(proc,OUT/f"mg5_generate_{tag}.log")
    ok1,actual_M,actual_mchi,actual_gchi,actual_width,actual_gqs=base.certify_persisted_param_card(param,M,cfg.MCHI,gq,gchi,gamma_tot)
    if not ok1: raise RuntimeError(f"M={M}: parameter card changed during generation")
    lhe_src=cfg.find_lhe(proc); banner=cfg.find_banner(proc)
    xsec=cfg.parse_lhe_init_xsec_pb(lhe_src); xb=cfg.parse_banner_xsec_pb(banner)
    banner_rel=abs(xb-xsec)/max(abs(xsec),1e-30)
    if banner_rel>1e-8: raise RuntimeError(f"M={M}: banner/LHE xsec mismatch {banner_rel}")
    lhe=plain_lhe(lhe_src,WORK/"lhe"/f"{tag}.lhe")
    shower=Path(os.environ.get("PKR_SHOWER",ROOT/"cache"/"v0.51"/"bin"/"pythia_lhe_to_hepmc_v0_51"))
    checkmate=Path(os.environ.get("PKR_CHECKMATE",ROOT/"cache"/"v0.51"/"checkmate-install"/"bin"/"CheckMATE"))
    hepmc=WORK/"hepmc"/f"{tag}.hepmc"; hepmc.parent.mkdir(parents=True,exist_ok=True)
    run([shower,lhe,hepmc,str(seed)],OUT/f"pythia_{tag}.log")
    for mod in ("spey","spey_pyhf","pyhf"): __import__(mod)
    v069.patch_optional_jax()
    card,cmout=write_card(tag,hepmc,xsec,seed)
    if cmout.exists(): shutil.rmtree(cmout)
    rc,stdout=run([checkmate,card],OUT/f"checkmate_{tag}.log",allow_nonzero=True)
    mu95,rmb,cls=cm51.parse_multibin(cmout)
    if mu95 is None:
        vals=re.findall(r"Upper limit:\s*([0-9.eE+-]+)",stdout)
        if vals:
            mu95=float(vals[-1]); rmb=1/mu95 if mu95>0 else None
    total=find_total(cmout)
    rows=parse_total(total)
    return {"M":M,"gq":gq,"a":a,"BR_inv":br,"Gamma_over_M":gamma_tot/M,"P":gq*gq*br,
            "xsec_pb":xsec,"sigma_unit_pb":xsec/(gq*gq*br),"banner_rel":banner_rel,
            "cm_rc":rc,"mu95":mu95,"rmb":rmb,"rows":rows,"total":str(total)}


def load_M50(pt: dict[str,float]):
    if not REF75.exists(): raise RuntimeError(f"missing {REF75}")
    j=json.loads(REF75.read_text())
    b=j["benchmark"]; gq=float(b["gq"]); br=float(b["BR_invisible"])
    Psample=gq*gq*br
    xsec=float(j["generator"]["cross_section_pb"])
    rows=parse_total(find_total(CM75))
    # Acceptance/unit-production vector is independent of the representative
    # coupling point in the narrow-width surrogate. Rescale later to v0.78 P.
    Prep=pt["gq"]**2*pt["BR_inv"]
    return {"M":50.0,"gq":pt["gq"],"a":pt["a"],"BR_inv":pt["BR_inv"],"P":Prep,
            "xsec_pb":xsec*(Prep/Psample),"sigma_unit_pb":xsec/Psample,"banner_rel":float(j["generator"]["banner_vs_lhe_rel"]),
            "cm_rc":int(j["checkmate"]["returncode"]),"mu95":j["checkmate"]["multibin_mu95"],"rmb":j["checkmate"]["multibin_r"],
            "rows_sample":rows,"sample_P":Psample,"rows":None,"total":str(find_total(CM75))}


def summarize(entry):
    P=entry["P"]
    if entry["M"]==50.0:
        scale=P/entry["sample_P"]
        rows={sr:{**v,"s":v["s"]*scale,"ds":v["ds"]*scale} for sr,v in entry["rows_sample"].items()}
    else:
        rows=entry["rows"]
    sum_s=sum(max(v["s"],0.0) for v in rows.values())
    mc_abs=[]; ds_rel=[]
    z2=0.0; bins=[]
    for sr,v in sorted(rows.items()):
        s=max(v["s"],0.0); ds=max(v["ds"],0.0); b=max(v["b"],0.0); db=max(v["db"],0.0)
        if s>0:
            mca=math.sqrt(max(ds*ds-(XSEC_REL*s)**2,0.0)); mc_abs.append(mca); ds_rel.append(ds/s)
        denom=max(b+db*db,1e-12)
        z2 += s*s/denom
        bins.append(dict(M_GeV=entry["M"],sr=sr,s=s,ds=ds,b=b,db=db,
                         coupling_factor=P,K_b=s/max(P,1e-300)))
    mc_scale=math.sqrt(sum(x*x for x in mc_abs))/max(sum_s,1e-300) if mc_abs else math.nan
    score=math.sqrt(z2)
    return dict(M_GeV=entry["M"],gq=entry["gq"],a_GeV_inv=entry["a"],BR_inv=entry["BR_inv"],
                coupling_factor=P,sigma_unit_pb=entry["sigma_unit_pb"],representative_xsec_pb=entry["xsec_pb"],
                sum_selected_signal=sum_s,diag_background_weighted_score=score,
                pilot_mc_scale=mc_scale,score_lo=score*(1-mc_scale),score_hi=score*(1+mc_scale),
                pilot_multibin_r=entry.get("rmb"),pilot_mu95=entry.get("mu95"),source_total_results=entry["total"]), bins


def main():
    OUT.mkdir(parents=True,exist_ok=True); WORK.mkdir(parents=True,exist_ok=True)
    pts=representative_points()
    print("=== v0.79 production x acceptance pilot ===",flush=True)
    print("50 GeV: reuse v0.75; 100/220 GeV: fresh pilot samples",flush=True)
    entries=[load_M50(pts[50.0])]
    for M,seed in ((100.0,7901),(220.0,7902)):
        print(f"running M={M:.0f} GeV pilot",flush=True)
        entries.append(generate_mass(M,pts[M],seed))
    summaries=[]; allbins=[]
    for e in entries:
        s,b=summarize(e); summaries.append(s); allbins.extend(b)
    S=pd.DataFrame(summaries).sort_values("diag_background_weighted_score")
    B=pd.DataFrame(allbins)
    S.to_csv(OUT/"mass_ranking.csv",index=False); B.to_csv(OUT/"bin_yields_and_Kb.csv",index=False)

    clear=False; candidate=None
    if len(S)>=2:
        a=S.iloc[0]; b=S.iloc[1]
        clear=bool(float(a.score_hi)<float(b.score_lo))
        if clear: candidate=float(a.M_GeV)
    meta={"version":"v0.79_eq4p4_mass_acceptance_pilot","pilot_events_per_new_mass":int(os.environ.get("PKR_NORM_NEVENTS","20000")),
          "ranking_metric":"sqrt(sum_b s_b^2/(b_b+db_b^2)); diagonal-background diagnostic only",
          "clear_rank_requires":"best score_hi < second-best score_lo using pilot MC-scale diagnostic",
          "clear_rank":clear,"candidate_M_GeV":candidate,
          "rule":"No exclusion claim from pilot. If no clear rank, do not authorize a full-statistics mass run."}
    (OUT/"summary.json").write_text(json.dumps(meta,indent=2)+"\n")
    print(S[["M_GeV","coupling_factor","sigma_unit_pb","representative_xsec_pb","sum_selected_signal","diag_background_weighted_score","pilot_mc_scale","score_lo","score_hi","pilot_multibin_r"]].to_string(index=False))
    print(f"clear cross-mass rank: {clear}")
    print(f"full-statistics candidate: {candidate if candidate is not None else 'NONE'}")
    print("interpretation: pilot ranks production x acceptance only; no collider exclusion/survival claim")
    return 0

if __name__=="__main__":
    raise SystemExit(main())
