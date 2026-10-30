#!/usr/bin/env python3
"""Reuse-first v0.42 discovery pass: no repeated kappa root scan.
The kappas are the established v0.41d M=50 GeV, gq=0.03 conditional 5% benchmarks.
"""
from __future__ import annotations
import importlib.util
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('pkr_v042_core',HERE/'pkr_absorption_horizon_v0_42.py')
mod=importlib.util.module_from_spec(spec); assert spec.loader is not None; spec.loader.exec_module(mod)
mod.OUT=mod.ROOT/'results'/'v0.42r'; mod.OUT.mkdir(parents=True,exist_ok=True)
mod.TABLEDIR=mod.OUT/'tables'; mod.TABLEDIR.mkdir(parents=True,exist_ok=True)
mod.MUS=np.array([0.12,0.25,0.45,0.70])
mod.EGRID=np.geomspace(1e4,3e8,25)
mod.TAU_REFINE=2
mod.MAX_STEP_KM=192.0
# v0.41d benchmark_gate.csv
K44=1.0579695766626163e-4
K43=5.882689479800913e-5
def fixed_solve(factor):
    k=K44 if factor==1.0 else K43
    return k,k,mod.ratio100(k,factor)
mod.solve_kappa=fixed_solve
if __name__=='__main__': mod.main()
