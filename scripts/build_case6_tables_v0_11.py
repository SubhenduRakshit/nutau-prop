#!/usr/bin/env python3
"""Build versioned Case-6 differential NC tables.

The expensive PDF integration is parallelised over energy/species.  Each
worker computes the SM CC/NC dσ/dy once and all requested Case-6 cutoffs.
Outputs are transport-ready NPZ files and an audit CSV.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path

import numpy as np
import pandas as pd

from nutau_prop.case6 import Case6Config, Case6DIS, mean_y, probability_from_dsdy

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "v0.11"
TABLES = OUT / "tables"
_CALC = None
_LAMBDAS = None
_CUT = None


def _init_worker(cfg_dict, lambdas, cut):
    global _CALC, _LAMBDAS, _CUT
    _CALC = Case6DIS(Case6Config(**cfg_dict))
    _LAMBDAS = tuple(float(x) for x in lambdas)
    _CUT = cut


def _worker(task):
    E, antinu = task
    s = _CALC.spectra(E, _LAMBDAS, antinu=antinu, cut=_CUT)
    return {
        "E": float(E), "antinu": bool(antinu), "y": s["y"], "dy": s["dy"],
        "sm_cc": s["sm_cc"], "sm_nc": s["sm_nc"], "np": s["np"],
    }


def tag(lam):
    return str(lam).replace(".", "p")


def build_table(records, lam, antinu):
    q = sorted((r for r in records if r["antinu"] == antinu), key=lambda r: r["E"])
    energies = np.asarray([r["E"] for r in q], float)
    y = np.asarray(q[0]["y"], float)
    dy = np.asarray(q[0]["dy"], float)
    sigma_cc, sigma_nc_sm, sigma_np = [], [], []
    p_cc, p_nc_sm, p_nc, p_np = [], [], [], []
    mean_np, mean_nc = [], []
    for r in q:
        cc = np.asarray(r["sm_cc"], float)
        ncsm = np.asarray(r["sm_nc"], float)
        npds = np.asarray(r["np"][float(lam)], float)
        scc, pcc = probability_from_dsdy(cc, dy)
        sncsm, pncsm = probability_from_dsdy(ncsm, dy)
        snp, pnp = probability_from_dsdy(npds, dy)
        snc, pnc = probability_from_dsdy(ncsm + npds, dy)
        sigma_cc.append(scc); sigma_nc_sm.append(sncsm); sigma_np.append(snp)
        p_cc.append(pcc); p_nc_sm.append(pncsm); p_nc.append(pnc); p_np.append(pnp)
        mean_np.append(mean_y(npds, y, dy))
        mean_nc.append(mean_y(ncsm + npds, y, dy))
    return dict(
        energies=energies, y=y, dy=dy,
        sigma_cc=np.asarray(sigma_cc),
        sigma_nc=np.asarray(sigma_nc_sm) + np.asarray(sigma_np),
        sigma_nc_sm=np.asarray(sigma_nc_sm),
        sigma_np=np.asarray(sigma_np),
        p_cc=np.asarray(p_cc), p_nc_sm=np.asarray(p_nc_sm),
        p_nc=np.asarray(p_nc), p_np=np.asarray(p_np),
        mean_y_np=np.asarray(mean_np), mean_y_nc_total=np.asarray(mean_nc),
    )


def interp_log(x, xp, fp):
    x = np.asarray(x, float); xp = np.asarray(xp, float); fp = np.asarray(fp, float)
    return 10.0 ** np.interp(np.log10(x), np.log10(xp), np.log10(np.maximum(fp, 1e-300)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=int(os.environ.get("NUTAU_NPROC", "14")))
    ap.add_argument("--nenergy", type=int, default=91)
    ap.add_argument("--emin-gev", type=float, default=1.0e4)
    ap.add_argument("--emax-gev", type=float, default=2.0e10)
    ap.add_argument("--lambdas", nargs="+", type=float, default=[1.5, 2.0, 2.44, 3.0])
    ap.add_argument("--cut", choices=["shat", "Q"], default="shat")
    ap.add_argument("--nx", type=int, default=96)
    ap.add_argument("--y-log-bins", type=int, default=36)
    ap.add_argument("--y-linear-bins", type=int, default=36)
    a = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True); TABLES.mkdir(parents=True, exist_ok=True)
    energies = np.geomspace(a.emin_gev, a.emax_gev, a.nenergy)
    cfg = Case6Config(nx=a.nx, y_log_bins=a.y_log_bins, y_linear_bins=a.y_linear_bins)
    cfg_dict = dict(cfg.__dict__)
    tasks = [(float(E), antinu) for antinu in (False, True) for E in energies]
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=max(1, a.workers), initializer=_init_worker,
                  initargs=(cfg_dict, a.lambdas, a.cut)) as pool:
        records = list(pool.imap_unordered(_worker, tasks, chunksize=1))

    audit = []
    v08_path = ROOT / "results" / "v0.8" / "case6_pdf_cut_v0_8.csv"
    v08 = pd.read_csv(v08_path) if v08_path.exists() else None
    for lam in a.lambdas:
        for antinu in (False, True):
            table = build_table(records, lam, antinu)
            species = "antinu" if antinu else "nu"
            path = TABLES / f"case6_L{tag(lam)}_{species}.npz"
            np.savez_compressed(path, **table)
            meta = {
                "version": "v0.11", "operator": "Case-6 dimension-7 momentum-growing NC-like",
                "C4_GeV^-3": cfg.c4_gev_m3, "Lambda_TeV": lam,
                "cut": a.cut, "pdf": cfg.pdf_set, "antinu": antinu,
                "note": "CC is SM; Case-6 is added only to NC. Both matched SM p_nc_sm and combined p_nc are stored.",
            }
            path.with_suffix(".json").write_text(json.dumps(meta, indent=2))
            sm_total = table["sigma_cc"] + table["sigma_nc_sm"]
            ratio = table["sigma_np"] / np.maximum(sm_total, 1e-300)
            for E, r, my in zip(table["energies"], ratio, table["mean_y_np"]):
                row = dict(E_GeV=E, E_PeV=E/1e6, Lambda_TeV=lam, species=species,
                           sigma_np_over_sm_total=r, mean_y_np=my,
                           sigma_np_cm2=float(interp_log([E], table["energies"], table["sigma_np"])[0]))
                if v08 is not None and not antinu:
                    q = v08[np.isclose(v08.Lambda_TeV, lam)].sort_values("E_PeV")
                    if len(q) and E/1e6 >= q.E_PeV.min() and E/1e6 <= q.E_PeV.max():
                        ref = interp_log([E/1e6], q.E_PeV, q.r_shatcut)[0]
                        row["v08_r_shatcut"] = ref
                        row["fractional_difference_vs_v08"] = r/ref - 1.0
                audit.append(row)

    A = pd.DataFrame(audit)
    A.to_csv(OUT / "case6_differential_audit_v0_11.csv", index=False)
    check = A[(A.species == "nu") & A.v08_r_shatcut.notna()] if "v08_r_shatcut" in A else pd.DataFrame()
    maxdiff = float(np.nanmax(np.abs(check.fractional_difference_vs_v08))) if len(check) else float("nan")
    summary = {
        "version": "v0.11", "workers": a.workers, "cut": a.cut,
        "energy_grid_GeV": [a.emin_gev, a.emax_gev, a.nenergy],
        "lambdas_TeV": a.lambdas,
        "max_abs_fractional_difference_vs_v08_nu": maxdiff,
        "validation_note": "Differences reflect a true dσ/dy discretisation; investigate if >10% in the UHE overlap.",
    }
    (OUT / "case6_differential_build_summary_v0_11.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    for epev in (1.0, 10.0, 100.0, 220.0, 300.0):
        print(f"\nE={epev:g} PeV")
        for lam in a.lambdas:
            q=A[(A.species=="nu") & np.isclose(A.Lambda_TeV,lam)].sort_values("E_PeV")
            r=interp_log([epev],q.E_PeV,q.sigma_np_over_sm_total)[0]
            my=float(np.interp(np.log10(epev),np.log10(q.E_PeV),q.mean_y_np))
            print(f"  Lambda={lam:4.2f} TeV  r_NP/SM={r:.4g}  <y>_NP={my:.3f}")
    print("WROTE", OUT)


if __name__ == "__main__":
    main()
