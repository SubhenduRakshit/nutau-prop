#!/usr/bin/env python3
"""v0.81 — corrected deterministic same-seed reweight of the M=50 Eq4.4 sample.

Reuses the v0.75 100k HepMC sample and *the same CheckMATE/Delphes random seed
7501*.  Only the signal cross section is rescaled to the v0.78 representative
coupling factor.  No MG5, no Pythia, no new events.

The run is accepted only if the detector-level signal vector agrees with the
v0.75 vector times one common scale (within numerical tolerance).  This closes
the v0.80 seed bug before interpreting the multibin likelihood.
"""
from __future__ import annotations

import importlib.util, json, os, re, shutil, subprocess
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.81_eq4p4_M50_same_seed_reweight"
WORK = ROOT / "cache" / "v0.81_eq4p4_M50_same_seed_reweight"
SUM78 = ROOT / "results" / "v0.78_eq4p4_full_envelope" / "mass_summary.csv"
REF75 = ROOT / "results" / "v0.75_eq4p4_bestpoint_checkmate" / "bestpoint_checkmate.json"
HEPMC75 = ROOT / "cache" / "v0.75_eq4p4_bestpoint_checkmate" / "hepmc" / "eq4p4_derivative__v074_best.hepmc"
CM75 = ROOT / "cache" / "v0.75_eq4p4_bestpoint_checkmate" / "checkmate"
ANALYSIS = "cms_2107_13201"
SEED = 7501
TOL = 2e-6

P51 = ROOT / "scripts" / "pkr_run_checkmate_v0_51.py"
s51 = importlib.util.spec_from_file_location("cm51", P51)
if s51 is None or s51.loader is None: raise RuntimeError(f"cannot import {P51}")
cm51 = importlib.util.module_from_spec(s51); s51.loader.exec_module(cm51)
P69 = ROOT / "scripts" / "pkr_checkmate_multibin_retry_v0_69.py"
s69 = importlib.util.spec_from_file_location("v069", P69)
if s69 is None or s69.loader is None: raise RuntimeError(f"cannot import {P69}")
v069 = importlib.util.module_from_spec(s69); s69.loader.exec_module(v069)


def find_total(root: Path) -> Path:
    hits = list(root.rglob("evaluation/total_results.txt"))
    if not hits: raise RuntimeError(f"no total_results.txt under {root}")
    return sorted(hits, key=lambda p: p.stat().st_mtime)[-1]


def parse_total(path: Path):
    out = {}
    for line in path.read_text(errors="replace").splitlines()[1:]:
        p=line.split()
        if len(p)<7 or p[0] != ANALYSIS: continue
        try: out[p[1]]=(float(p[5]),float(p[6]))
        except ValueError: pass
    if not out: raise RuntimeError(f"no {ANALYSIS} bins in {path}")
    return out


def main() -> int:
    OUT.mkdir(parents=True,exist_ok=True); WORK.mkdir(parents=True,exist_ok=True)
    if not SUM78.exists() or not REF75.exists() or not HEPMC75.exists():
        raise RuntimeError("missing v0.78 summary, v0.75 JSON, or v0.75 HepMC")

    q=pd.read_csv(SUM78); q=q[q.M_GeV==50.0]
    if len(q)!=1: raise RuntimeError("expected one M=50 row in v0.78 summary")
    row=q.iloc[0]
    P_rep=float(row.representative_monojet_coupling_factor)
    gq_rep=float(row.representative_gq); br_rep=float(row.representative_BR_inv); a_rep=float(row.representative_a)

    ref=json.loads(REF75.read_text())
    gq_ref=float(ref["benchmark"]["gq"]); br_ref=float(ref["benchmark"]["BR_invisible"])
    P_ref=gq_ref*gq_ref*br_ref; xsec_ref=float(ref["generator"]["cross_section_pb"])
    scale=P_rep/P_ref; xsec_rep=xsec_ref*scale

    for mod in ("spey","spey_pyhf","pyhf"): __import__(mod)
    v069.patch_optional_jax()

    cmout=WORK/"checkmate"/"M50_rep_same_seed"
    card=WORK/"checkmate_M50_rep_same_seed.in"
    card.write_text(f"""[Parameters]\nBackend: numpy\nName: pkr_v081_M50_rep_same_seed\nAnalyses: {ANALYSIS}\nOutputDirectory: {cmout}\nOutputExists: Overwrite\nSkipParamCheck: True\nEffTab: True\nMultibin: detailed\nExpected: True\nUplim: True\nModel: simple\nInvisiblePIDs: 5000521\nRandomSeed: {SEED}\nTotalResultFileColumns: analysis, sr, o, b, db, s, ds, s95obs, s95exp, robscons, rexpcons\nBestPerAnalysisResultFileColumns: analysis, sr, o, b, db, s, ds, s95obs, s95exp, robscons, rexpcons\n\n[signal]\nEvents: {HEPMC75}\nXSect: {xsec_rep:.12g} PB\nXSectErr: 10 %%\n""")
    if cmout.exists(): shutil.rmtree(cmout)
    checkmate=Path(os.environ.get("PKR_CHECKMATE",ROOT/"cache"/"v0.51"/"checkmate-install"/"bin"/"CheckMATE"))
    p=subprocess.run([str(checkmate),str(card)],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (OUT/"checkmate.log").write_text(p.stdout)

    # Hard fixed-shape certification against v0.75 detector output.
    old=parse_total(find_total(CM75)); new=parse_total(find_total(cmout))
    common=sorted(set(old)&set(new))
    if not common: raise RuntimeError("no common detector bins")
    sdev=[]; dsdev=[]
    for sr in common:
        so,dso=old[sr]; sn,dsn=new[sr]
        if abs(so)>1e-12: sdev.append(abs(sn/(scale*so)-1.0))
        if abs(dso)>1e-12: dsdev.append(abs(dsn/(scale*dso)-1.0))
    max_sdev=max(sdev) if sdev else 0.0
    max_dsdev=max(dsdev) if dsdev else 0.0
    shape_pass=(max_sdev<TOL and max_dsdev<TOL)
    if not shape_pass:
        raise RuntimeError(f"same-seed fixed-shape certification failed: max s dev={max_sdev:.3e}, ds dev={max_dsdev:.3e}")

    mu95,rmb,cls=cm51.parse_multibin(cmout); rsingle=cm51.parse_singlebin(cmout)
    if mu95 is None:
        vals=re.findall(r"Upper limit:\s*([0-9.eE+-]+)",p.stdout)
        if vals:
            mu95=float(vals[-1]); rmb=1.0/mu95 if mu95>0 else None
    if rmb is None: raise RuntimeError(f"multibin unavailable; rc={p.returncode}")

    expected_linear=float(ref["checkmate"]["multibin_r"])*scale
    linear_rel=abs(rmb-expected_linear)/max(abs(expected_linear),1e-30)
    if rmb<0.8: status="ROBUST_SURVIVAL_THRESHOLD"
    elif rmb>1.2: status="ROBUST_EXCLUSION_THRESHOLD"
    else: status="MARGINAL_BOUNDARY_UNRESOLVED"

    out={
      "version":"v0.81_eq4p4_M50_same_seed_reweight",
      "source":"same v0.75 HepMC and same Delphes seed 7501; normalization only",
      "invalid_predecessor":"v0.80 changed seed to 8001 and is diagnostic only",
      "representative":{"M_GeV":50.0,"gq":gq_rep,"a_GeV_inv":a_rep,"BR_inv":br_rep,"coupling_factor_gq2BR":P_rep},
      "scale":scale,"xsec_ref_pb":xsec_ref,"xsec_rep_pb":xsec_rep,
      "fixed_shape_check":{"n_common_bins":len(common),"max_signal_scale_deviation":max_sdev,"max_ds_scale_deviation":max_dsdev,"tolerance":TOL,"pass":shape_pass},
      "checkmate":{"returncode":p.returncode,"mu95":mu95,"multibin_r":rmb,"singlebin_robscons":rsingle,"status":status},
      "linear_diagnostic":{"r_expected_from_v075_times_scale":expected_linear,"relative_difference":linear_rel},
      "interpretation":"Only v0.81, not v0.80, is the deterministic same-shape normalization reweight."
    }
    (OUT/"M50_same_seed_reweight.json").write_text(json.dumps(out,indent=2)+"\n")

    print("=== v0.81 Eq4.4 M50 deterministic same-seed reweight ===")
    print(f"P_rep/P_v075 = {scale:.6f}")
    print(f"fixed-shape check: PASS; max s dev={max_sdev:.3e}, max ds dev={max_dsdev:.3e}")
    print(f"linear diagnostic r~{expected_linear:.6f}")
    print(f"mu95={mu95}  multibin_r={rmb}  singlebin={rsingle}")
    print(f"status: {status}")
    print("MG5 rerun: no; Pythia rerun: no; new events: no; Delphes seed identical to v0.75")
    return 0

if __name__=="__main__": raise SystemExit(main())
