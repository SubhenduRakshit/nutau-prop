#!/usr/bin/env python3
"""
Transport-step convergence check using an already-built SM/BSM table.
Runs the same benchmark at successively smaller max interaction probabilities.
"""
import argparse
import numpy as np
from nutau_prop.config import load_config
from nutau_prop.tables import load_table
from nutau_prop.transport import TauTransport
from nutau_prop.util import powerlaw_per_logE

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",default="config/default.yaml")
    ap.add_argument("--table",default="cache/sm.npz")
    ap.add_argument("--cosz",type=float,default=1.0)
    args=ap.parse_args()
    c=load_config(args.config); t=load_table(args.table)
    tr=TauTransport(t,float(c["physics"]["tau_beta_cm2_per_g"]),
                    float(c["physics"]["tau_polarization"]))
    p0=powerlaw_per_logE(t["energies"],float(c["physics"]["spectral_index"]))
    outs=[]
    for mp in [0.05,0.025,0.0125]:
        r=tr.propagate(args.cosz,p0,max_step_km=float(c["transport"]["max_step_km"]),max_p=mp)
        outs.append(r["nu"])
    for a,b,label in [(outs[0],outs[1],"0.05->0.025"),(outs[1],outs[2],"0.025->0.0125")]:
        rel=np.max(np.abs(a-b)/np.maximum(np.abs(b),1e-300))
        print(label,rel)

if __name__=="__main__":
    main()
