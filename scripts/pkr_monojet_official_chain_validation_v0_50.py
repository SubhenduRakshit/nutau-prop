#!/usr/bin/env python3
"""v0.50: validate the unmodified LLP-CDFO CMS-EXO-20-004 limit chain.

This deliberately does not evaluate the PKR benchmark.  It tests whether the
public recast authors' own processed CMS DataFrame reproduces its stored robs /
rexp values when the repository's unmodified cms_exo_20_004_UpperLimits.py is
run on the present Mac/Python stack.

Interpretation
--------------
PASS_EXACT/GOOD:
  official public chain is numerically reproducible here.  Then the failure of
  v0.48--v0.49c came from applying public DMsimp templates outside the validated
  recast pipeline/statistical mapping; proceed to genuine 50-GeV MC + recast.
FAIL_RUNTIME_OR_REPRODUCIBILITY:
  the public chain itself is incompatible with the current runtime; create a
  compatibility environment before doing any 50-GeV physics.
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MIRROR = ROOT / "cache" / "LLP-CDFO"
OUT = ROOT / "results" / "v0.50"
OUT.mkdir(parents=True, exist_ok=True)

SRC = MIRROR / "results_dataFrames" / "pp2BB1j_cms_latest.pcl"
ANA_DIR = MIRROR / "CMS-EXO-20-004"
LIMIT_SCRIPT = ANA_DIR / "cms_exo_20_004_UpperLimits.py"
WORK_IN = ROOT / "cache" / "v0.50_official_input.pcl"
WORK_OUT = ROOT / "cache" / "v0.50_official_recomputed.pcl"


def finite_rel(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if not np.any(m):
        return np.array([], float)
    den = np.maximum(np.abs(a[m]), 1e-15)
    return np.abs(b[m]-a[m]) / den


def main():
    missing = [p for p in (SRC, LIMIT_SCRIPT) if not p.exists()]
    if missing:
        raise RuntimeError("Missing official recast inputs:\n" + "\n".join(map(str, missing)))

    original = pd.read_pickle(SRC)
    cols = original.columns.tolist()
    has_robs = "robs" in cols
    has_rexp = "rexp" in cols

    shutil.copy2(SRC, WORK_IN)
    if WORK_OUT.exists():
        WORK_OUT.unlink()

    cmd = [sys.executable, str(LIMIT_SCRIPT), "-f", str(WORK_IN), "-o", str(WORK_OUT)]
    proc = subprocess.run(cmd, cwd=str(ANA_DIR), text=True, capture_output=True)
    (OUT / "official_upperlimits_stdout.log").write_text(proc.stdout)
    (OUT / "official_upperlimits_stderr.log").write_text(proc.stderr)

    rows = []
    status = "FAIL_RUNTIME_OR_REPRODUCIBILITY"
    metrics = {
        "returncode": int(proc.returncode),
        "input_rows": int(len(original)),
        "has_stored_robs": bool(has_robs),
        "has_stored_rexp": bool(has_rexp),
    }

    if proc.returncode == 0 and WORK_OUT.exists():
        new = pd.read_pickle(WORK_OUT)
        if len(new) != len(original):
            raise RuntimeError(f"Official script changed row count: {len(original)} -> {len(new)}")

        for field in ("robs", "rexp"):
            if field not in original.columns or field not in new.columns:
                continue
            old = pd.to_numeric(original[field], errors="coerce").to_numpy(float)
            now = pd.to_numeric(new[field], errors="coerce").to_numpy(float)
            m = np.isfinite(old) & np.isfinite(now)
            rel = finite_rel(old, now)
            absd = np.abs(now[m]-old[m]) if np.any(m) else np.array([], float)
            metrics[f"{field}_ncomp"] = int(np.sum(m))
            metrics[f"{field}_max_abs"] = float(np.max(absd)) if len(absd) else None
            metrics[f"{field}_max_rel"] = float(np.max(rel)) if len(rel) else None
            metrics[f"{field}_median_rel"] = float(np.median(rel)) if len(rel) else None
            for i in np.where(m)[0][:20]:
                rows.append({"row": int(i), "quantity": field,
                             "stored": float(old[i]), "recomputed": float(now[i]),
                             "abs_diff": float(abs(now[i]-old[i])),
                             "rel_diff": float(abs(now[i]-old[i]) / max(abs(old[i]),1e-15))})

        comparisons = [v for k,v in metrics.items() if k.endswith("_max_rel") and v is not None]
        ncomp = sum(int(metrics.get(f"{q}_ncomp",0)) for q in ("robs","rexp"))
        if ncomp == 0:
            status = "NO_STORED_LIMITS_TO_COMPARE"
        elif comparisons and max(comparisons) < 1e-8:
            status = "PASS_EXACT_REPRODUCTION"
        elif comparisons and max(comparisons) < 0.02:
            status = "PASS_GOOD_REPRODUCTION"
        else:
            status = "FAIL_REPRODUCIBILITY"

    pd.DataFrame(rows).to_csv(OUT / "official_chain_comparison.csv", index=False)
    summary = pd.DataFrame([{"status": status, **metrics}])
    summary.to_csv(OUT / "summary.csv", index=False)
    meta = {
        "version": "v0.50",
        "purpose": "unmodified LLP-CDFO CMS-EXO-20-004 chain validation on authors' processed DataFrame",
        "source_dataframe": "results_dataFrames/pp2BB1j_cms_latest.pcl",
        "official_limit_script": "CMS-EXO-20-004/cms_exo_20_004_UpperLimits.py",
        "python": sys.version,
        "status": status,
        "important": "No PKR exclusion/survival verdict is produced by v0.50.",
        "previous_status": "v0.48--v0.49c discarded for failed public-template validation",
    }
    (OUT / "metadata.json").write_text(json.dumps(meta, indent=2))

    print("\n=== v0.50 official CMS recast-chain validation ===")
    print(summary.to_string(index=False))
    if rows:
        print("\n=== v0.50 first stored-vs-recomputed checks ===")
        print(pd.DataFrame(rows).head(12).to_string(index=False))
    print("\nOfficial script return code:", proc.returncode)
    if proc.returncode != 0:
        print("See results/v0.50/official_upperlimits_stderr.log")


if __name__ == "__main__":
    main()
