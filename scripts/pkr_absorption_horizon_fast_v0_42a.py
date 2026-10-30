#!/usr/bin/env python3
"""Fast horizon-focused screen using the corrected v0.42 physics implementation.
Runs only representative absorption-horizon chords. Any surviving best cell is
then checked with the finer v0.42 convergence settings before certification.
"""
from __future__ import annotations
import importlib.util
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
CORE=HERE/'pkr_absorption_horizon_v0_42.py'
spec=importlib.util.spec_from_file_location('pkr_v042_core',CORE)
mod=importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(mod)

mod.OUT=mod.ROOT/'results'/'v0.42a'
mod.OUT.mkdir(parents=True,exist_ok=True)
mod.TABLEDIR=mod.OUT/'tables'
mod.TABLEDIR.mkdir(parents=True,exist_ok=True)
mod.MUS=np.array([0.12,0.25,0.45,0.70])
mod.EGRID=np.geomspace(1e4,3e8,41)
mod.TAU_REFINE=4
mod.MAX_STEP_KM=96.0

if __name__=='__main__':
    mod.main()
