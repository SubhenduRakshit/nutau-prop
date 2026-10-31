#!/usr/bin/env python3
"""v0.88b: numerical-grid correction for v0.88 CGC propagation.

The v0.88 physics implementation is unchanged.  This wrapper removes
floating-point near-duplicate energy nodes introduced when explicit target
energies coincide with values already present in np.geomspace.  Those nodes
produced zero dlogE after tau-grid refinement.  A fresh cache/output directory
prevents reuse of the failed v0.88 tables.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import cgc_earth_allflavor_v0_88 as base


def merge_near_duplicates(values, rtol=1.0e-12):
    vals = np.sort(np.asarray(values, dtype=float))
    keep = [float(vals[0])]
    for x in vals[1:]:
        if (float(x) - keep[-1]) / keep[-1] > rtol:
            keep.append(float(x))
    out = np.asarray(keep, dtype=float)
    if np.any(np.diff(np.log(out)) <= 0.0):
        raise RuntimeError("energy grid is not strictly increasing after de-duplication")
    return out


base.EGRID = merge_near_duplicates(
    np.r_[np.geomspace(3.0e4, 3.0e10, 67), base.TARGET_E]
)
base.OUT = base.ROOT / "results" / "v0.88b_cgc_allflavor"
base.CACHE = base.ROOT / "cache" / "v0.88b_cgc_allflavor"
base.OUT.mkdir(parents=True, exist_ok=True)
base.CACHE.mkdir(parents=True, exist_ok=True)

if __name__ == "__main__":
    print(f"v0.88b grid: {len(base.EGRID)} strictly increasing energy nodes; "
          f"min dlogE={np.min(np.diff(np.log(base.EGRID))):.6e}")
    raise SystemExit(base.main())
