#!/usr/bin/env python3
"""Coarse discovery pass only; identifies horizon cells for later convergence."""
from __future__ import annotations
import importlib.util
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('pkr_v042_core',HERE/'pkr_absorption_horizon_v0_42.py')
mod=importlib.util.module_from_spec(spec); assert spec.loader is not None; spec.loader.exec_module(mod)
mod.OUT=mod.ROOT/'results'/'v0.42q'; mod.OUT.mkdir(parents=True,exist_ok=True)
mod.TABLEDIR=mod.OUT/'tables'; mod.TABLEDIR.mkdir(parents=True,exist_ok=True)
mod.MUS=np.array([0.12,0.25,0.45,0.70])
mod.EGRID=np.geomspace(1e4,3e8,31)
mod.TAU_REFINE=2
mod.MAX_STEP_KM=192.0
if __name__=='__main__': mod.main()
