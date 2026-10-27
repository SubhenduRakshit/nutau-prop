#!/usr/bin/env python3
"""Cached launcher for the v0.13 joint Case-6 fold.

The core script repeatedly asks for the same spectrum-specific SM primary-numu
Earth transmission while folding several detector channels.  Cache those
transmission curves by (mu, energy-grid) so the post-processing run stays fast.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
CORE = HERE / "run_case6_joint_event_fold_v0_13.py"
spec = importlib.util.spec_from_file_location("case6_joint_core_v013", CORE)
mod = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(mod)

_original = mod.primary_numu_transmission
_cache = {}


def cached_primary_numu_transmission(Eout, mu, lam_for_sm=1.5):
    Eout = np.asarray(Eout, float)
    key = (round(float(mu), 10), float(lam_for_sm), len(Eout), float(Eout[0]), float(Eout[-1]))
    if key not in _cache:
        _cache[key] = _original(Eout, mu, lam_for_sm=lam_for_sm)
    return _cache[key]


mod.primary_numu_transmission = cached_primary_numu_transmission

if __name__ == "__main__":
    mod.main()
