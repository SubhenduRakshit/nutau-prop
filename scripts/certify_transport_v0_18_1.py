#!/usr/bin/env python3
"""Aggregate certification gate for v0.18.1 production transport.

Production is allowed only if:
  * historical RR and DRS beta_tau=0 gates pass;
  * all substantive v0.17 primitive/full-transport gates pass, excluding only
    the superseded refinement_trend_tau heuristic;
  * the v0.17.1 refinement reassessment passes;
  * the v0.18.1 actual chirality-solver convergence gate passes.
"""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
H=ROOT/"results"/"historical"
V17=ROOT/"results"/"v0.17"
V181=ROOT/"results"/"v0.18.1"
OUT=V181


def readj(p):
    if not p.exists(): raise SystemExit(f"Missing certification input: {p}")
    return json.loads(p.read_text())

def as_bool(x):
    if isinstance(x,bool): return x
    if isinstance(x,(int,float)): return bool(x)
    s=str(x).strip().lower()
    if s in ("true","1","yes"): return True
    if s in ("false","0","no"): return False
    raise SystemExit(f"Unrecognized PASS value {x!r}")

checks={}
rr=readj(H/"rr_gate.json"); drs=readj(H/"drs_gate.json")
checks["historical_RR_PASS"]=as_bool(rr.get("PASS",False))
checks["historical_DRS_PASS"]=as_bool(drs.get("PASS",False))

p=V17/"tau_validation_gates_v0_17.csv"
if not p.exists(): raise SystemExit(f"Missing certification input: {p}")
G=pd.read_csv(p)
if "PASS" not in G or "test" not in G: raise SystemExit("Malformed v0.17 gate CSV")
# refinement_trend_tau was replaced by the explicit v0.17.1 floor-aware reassessment.
sub=G[G.test!="refinement_trend_tau"].copy()
pvals=sub.PASS.map(as_bool)
checks["v0.17_all_substantive_gates_PASS"]=bool(pvals.all())
checks["v0.17_substantive_gate_count"]=int(len(sub))
checks["v0.17_substantive_failed"]=int((~pvals).sum())

r=readj(V17/"tau_refinement_reassessment_status_v0_17_1.json")
checks["v0.17.1_refinement_PASS"]=as_bool(r.get("all_refinement_gates_PASS",False))

c=readj(V181/"chirality_validation_status_v0_18_1.json")
checks["v0.18.1_chirality_solver_PASS"]=as_bool(c.get("all_gates_PASS",False))

allpass=(checks["historical_RR_PASS"] and checks["historical_DRS_PASS"] and
         checks["v0.17_all_substantive_gates_PASS"] and checks["v0.17.1_refinement_PASS"] and
         checks["v0.18.1_chirality_solver_PASS"])
status={"version":"v0.18.1","production_certified":bool(allpass),"checks":checks,
        "scope":"numerical transport certification only; not detector or UV-model validation"}
OUT.mkdir(parents=True,exist_ok=True)
(OUT/"transport_production_certification_v0_18_1.json").write_text(json.dumps(status,indent=2))
print(json.dumps(status,indent=2))
if not allpass: raise SystemExit("v0.18.1 transport production certification FAILED")
