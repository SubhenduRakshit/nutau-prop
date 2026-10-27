#!/usr/bin/env python3
import json
from pathlib import Path
def load(p):
    p=Path(p)
    return json.loads(p.read_text()) if p.exists() else {"PASS":False,"missing":str(p)}
rr=load("results/historical/rr_gate.json")
drs=load("results/historical/drs_gate.json")
d={"rakshit_reya":rr,"dutta_reno_sarcevic":drs,
   "PASS":bool(rr.get("PASS")) and bool(drs.get("PASS"))}
Path("results/historical_gate.json").write_text(json.dumps(d,indent=2))
print(json.dumps(d,indent=2))
if not d["PASS"]:
    raise SystemExit("HISTORICAL VALIDATION GATE FAILED: BSM scan remains blocked")
print("HISTORICAL VALIDATION GATE: PASS")
