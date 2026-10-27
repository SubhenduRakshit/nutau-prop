#!/usr/bin/env python3
import argparse
import numpy as np
from nutau_prop.config import load_config
from nutau_prop.pdfs import LHAPDFNucleon
from nutau_prop.cross_sections import DISCalculator,DISConfig
from nutau_prop.bsm import BSMModel

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",default="config/default.yaml")
    args=ap.parse_args()
    c=load_config(args.config); ic=c["integration"]
    pdf=LHAPDFNucleon(c["pdf"]["set"],c["pdf"]["member"])
    calc=DISCalculator(pdf,DISConfig(float(ic["q2_min_gev2"]),float(ic["w2_min_gev2"]),
                                     int(ic["nx"]),float(c["physics"]["sin2thetaW"])))
    print("# E_GeV sigma_CC_cm2 sigma_NC_cm2 NC_over_CC")
    for E in [1e4,1e5,1e6,1e7]:
        cc=calc.spectrum(E,"CC",model=BSMModel(),
             y_min=float(ic["y_min"]),y_split=float(ic["y_split"]),
             nlog=int(ic["y_log_bins"]),nlin=int(ic["y_linear_bins"]))["sigma"]
        nc=calc.spectrum(E,"NC",model=BSMModel(),
             y_min=float(ic["y_min"]),y_split=float(ic["y_split"]),
             nlog=int(ic["y_log_bins"]),nlin=int(ic["y_linear_bins"]))["sigma"]
        print(f"{E:.6e} {cc:.6e} {nc:.6e} {nc/cc:.5f}")
    print("\nReference check: modern calculations find NC/CC ~0.41 near 1 PeV.")
    print("LO code need not match NNLO normalization exactly; use ratio and trends as validation.")

if __name__=="__main__":
    main()
