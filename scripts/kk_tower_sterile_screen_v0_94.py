#!/usr/bin/env python3
"""nutau-prop v0.94: corrected KK-tower sterile screen.

Correction to v0.93
-------------------
v0.93 correctly evaluated H(shat)=(sqrt(shat)/Qref)^delta, but its summary
incorrectly stated that the d=6 c=1 contact rows necessarily had H<=1.
That is true only for the conservative choice Qref=cutoff.  If an external
collider reference scale lies below the cutoff, H>1 is possible and must be
reported as a sensitivity scan, not as a certified collider-normalized result.

This version keeps the validated v0.93 hard kernels and separates:
  * contact_cutoff_conservative: Qref = EFT cutoff (guarantees H<=1);
  * contact_external_qref: explicit Qref sensitivity scan;
  * pkr_width_anchor: Qref = M_Z' = 50 GeV;
  * pkr_hard_qref: 250/500 GeV sensitivity scans.

All >5% rows are SIZE gates only.  None is called allowed until a tower-aware
collider/width + bulk-neutrino constraint recast is performed.
"""
from __future__ import annotations

import argparse, csv, json, math, os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import kk_tower_sterile_screen_v0_93 as v93


def contact_cutoff_GeV(channel: str, domain: str) -> float:
    G = v93.D6[channel]["G_TeV_m2"]
    if domain == "c1":
        lam_TeV = 1.0 / math.sqrt(G)
    elif domain == "4pi":
        lam_TeV = math.sqrt(4.0 * math.pi / G)
    else:
        raise ValueError(domain)
    return 1.0e3 * lam_TeV


def continuum_diag(qref_eff, cutoff, delta, spacing):
    if spacing <= 0:
        return dict(kk_spacing_GeV=0.0, qref_over_spacing=float("nan"),
                    cutoff_over_spacing=float("nan"), dense_ref_assumption=None)
    return dict(
        kk_spacing_GeV=spacing,
        qref_over_spacing=qref_eff / spacing,
        cutoff_over_spacing=cutoff / spacing,
        dense_ref_assumption=bool(qref_eff / spacing >= 10.0),
    )


def one_job(job):
    (pdf, nx, E_PeV, channel, domain, delta, scheme, qref_input, spacing) = job
    calc = v93.Case6DIS(v93.Case6Config(pdf_set=pdf, nx=nx))
    E = E_PeV * v93.PEV_TO_GEV
    sm = v93.sm_total(calc, E)

    if channel in v93.D6:
        b = v93.D6[channel]
        cutoff = contact_cutoff_GeV(channel, domain)
        cutoff_TeV = cutoff / 1.0e3
        if scheme == "contact_cutoff_conservative":
            qref = cutoff
            interpretation = "conservative contact normalization; H<=1 pointwise"
        elif scheme == "contact_external_qref":
            qref = qref_input
            interpretation = "external-Qref sensitivity only; requires explicit tower-aware collider mapping"
        else:
            raise ValueError(scheme)
        sig0 = v93.d6_sigma(calc, E, b, cutoff_TeV, delta=0)
        sig = v93.d6_sigma(calc, E, b, cutoff_TeV, delta=delta, qref_GeV=qref)
    elif channel == "pkr_eq4p3_dipole":
        cutoff = v93.PKR["gq"] / v93.PKR["kappa_GeV_m1"]
        qref = qref_input
        if scheme == "pkr_width_anchor":
            interpretation = "50-GeV on-shell width anchor; tower-aware total width still required"
        elif scheme == "pkr_hard_qref":
            interpretation = "hard-scale sensitivity scan; not a collider recast"
        else:
            raise ValueError(scheme)
        sig0 = v93.pkr_sigma(calc, E, delta=0)
        sig = v93.pkr_sigma(calc, E, delta=delta, qref_GeV=qref)
    else:
        raise ValueError(channel)

    qeff = v93.qref_effective(qref, cutoff)
    r0 = sig0 / sm
    r = sig / sm
    enh = sig / sig0 if sig0 > 0 else float("nan")
    row = dict(
        pdf=pdf, E_PeV=E_PeV, channel=channel, domain=domain, delta=delta,
        normalization_scheme=scheme, qref_input_GeV=qref,
        qref_effective_GeV=qeff, cutoff_GeV=cutoff,
        sigma_SM_cm2=sm, sigma_single_cm2=sig0, sigma_KK_cm2=sig,
        r_single_over_SM=r0, r_KK_over_SM=r,
        enhancement_KK_over_single=enh, size_gate_5pct=bool(r >= 0.05),
        interpretation=interpretation,
    )
    row.update(continuum_diag(qeff, cutoff, delta, spacing))
    return row


def build_jobs(pdf, nx, energies, deltas, contact_qrefs, pkr_qrefs, spacing):
    jobs=[]
    for ch in ("scalar_tau23", "tensor_tau23"):
        for dom in ("c1", "4pi"):
            for E in energies:
                for d in deltas:
                    jobs.append((pdf,nx,E,ch,dom,d,"contact_cutoff_conservative",0.0,spacing))
                    for q in contact_qrefs:
                        jobs.append((pdf,nx,E,ch,dom,d,"contact_external_qref",q,spacing))
    for E in energies:
        for d in deltas:
            for q in pkr_qrefs:
                scheme = "pkr_width_anchor" if abs(q-v93.PKR["M_GeV"]) < 1e-12 else "pkr_hard_qref"
                jobs.append((pdf,nx,E,"pkr_eq4p3_dipole","vertex_c1",d,scheme,q,spacing))
    return jobs


def best_rows(rows):
    keys=sorted({(r["channel"],r["domain"],r["normalization_scheme"],r["delta"],r["qref_input_GeV"]) for r in rows})
    out=[]
    for k in keys:
        rr=[r for r in rows if (r["channel"],r["domain"],r["normalization_scheme"],r["delta"],r["qref_input_GeV"])==k]
        out.append(max(rr,key=lambda z:z["r_KK_over_SM"]))
    return out


def write_summary(rows,out):
    best=best_rows(rows)
    conservative=[r for r in best if r["normalization_scheme"]=="contact_cutoff_conservative"]
    external=[r for r in best if r["normalization_scheme"]=="contact_external_qref"]
    pkr=[r for r in best if r["channel"]=="pkr_eq4p3_dipole"]
    lines=[
        "# v0.94 corrected bulk-sterile KK tower screen","",
        "## Correction to v0.93","",
        "v0.93's numerical folds are retained, but its statement that all c=1 contact rows must have H<=1 was too strong.",
        "That statement holds only when the constraint normalization is conservatively placed at the EFT cutoff.",
        "For an external reference scale below the cutoff, H can exceed unity; those rows are therefore sensitivity scans, not certified bounds.","",
        "The continuum tower factor remains",
        r"\[H(\hat s)=N_{\rm eff}(\sqrt{\hat s})/N_{\rm eff}(Q_{\rm ref})\simeq(\sqrt{\hat s}/Q_{\rm ref})^\delta.\]","",
        "Radius/mode degeneracy cancel only in the dense-tower approximation with approximately mode-independent brane couplings and negligible final-state masses.","",
        "## A. Conservative contact normalization: Qref = cutoff","",
        "| channel | domain | d | max rKK/SM | E PeV | enhancement | gate |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for r in conservative:
        lines.append(f"| {r['channel']} | {r['domain']} | {r['delta']} | {100*r['r_KK_over_SM']:.4g}% | {r['E_PeV']:.0f} | {r['enhancement_KK_over_single']:.3g} | {'PASS' if r['size_gate_5pct'] else 'FAIL'} |")
    lines += ["", "These rows are the only contact-operator rows for which H<=1 is guaranteed pointwise.","",
              "## B. Contact external-Qref sensitivity","",
              "| channel | domain | d | Qref GeV | max rKK/SM | E PeV | enhancement | gate |",
              "|---|---|---:|---:|---:|---:|---:|---|"]
    for r in external:
        lines.append(f"| {r['channel']} | {r['domain']} | {r['delta']} | {r['qref_input_GeV']:.0f} | {100*r['r_KK_over_SM']:.4g}% | {r['E_PeV']:.0f} | {r['enhancement_KK_over_single']:.3g} | {'PASS' if r['size_gate_5pct'] else 'FAIL'} |")
    lines += ["", "These are not collider-certified until Qref is derived from an explicit tower-aware likelihood.","",
              "## C. Resolved PKR transition-dipole","",
              "| scheme | d | Qref GeV | max rKK/SM | E PeV | enhancement | gate |",
              "|---|---:|---:|---:|---:|---:|---|"]
    for r in pkr:
        lines.append(f"| {r['normalization_scheme']} | {r['delta']} | {r['qref_input_GeV']:.0f} | {100*r['r_KK_over_SM']:.4g}% | {r['E_PeV']:.0f} | {r['enhancement_KK_over_single']:.3g} | {'PASS' if r['size_gate_5pct'] else 'FAIL'} |")
    lines += ["", "The Qref=50 GeV rows are width-anchor screens, not allowed points. A KK tower also enhances the mediator invisible width, so the same tower must be included in the width/monojet fit.","",
              "## Decision","",
              "Do not run full Earth propagation yet. First perform the tower-aware PKR width + monojet recast and map the bulk-neutrino parameterization (R, mode couplings, masses) onto existing constraints. If a recast-surviving point remains above 5%, then propagate it.",""]
    (out/"KK_TOWER_STERILE_SCREEN_v0.94.md").write_text("\n".join(lines))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--pdf",default="CT18NLO")
    ap.add_argument("--nx",type=int,default=96)
    ap.add_argument("--energies-pev",default="100,300,1000")
    ap.add_argument("--deltas",default="1,2,3")
    ap.add_argument("--contact-qrefs-gev",default="1000,2000")
    ap.add_argument("--pkr-qrefs-gev",default="50,250,500")
    ap.add_argument("--kk-spacing-gev",type=float,default=0.0,
                    help="optional spacing diagnostic; 0 leaves dense-tower assumption unquantified")
    ap.add_argument("--nproc",type=int,default=int(os.environ.get("NUTAU_NPROC","14")))
    ap.add_argument("--outdir",default="results/v0.94_kk_tower_corrected")
    a=ap.parse_args()
    energies=tuple(float(x) for x in a.energies_pev.split(',') if x.strip())
    deltas=tuple(int(x) for x in a.deltas.split(',') if x.strip())
    cq=tuple(float(x) for x in a.contact_qrefs_gev.split(',') if x.strip())
    pq=tuple(float(x) for x in a.pkr_qrefs_gev.split(',') if x.strip())
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    jobs=build_jobs(a.pdf,a.nx,energies,deltas,cq,pq,a.kk_spacing_gev)
    rows=[]; nproc=max(1,min(a.nproc,len(jobs)))
    if nproc==1:
        rows=[one_job(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=nproc) as ex:
            fs=[ex.submit(one_job,j) for j in jobs]
            for i,f in enumerate(as_completed(fs),1):
                rows.append(f.result())
                if i%12==0 or i==len(fs): print(f"[v0.94] completed {i}/{len(fs)} jobs",flush=True)
    rows.sort(key=lambda r:(r['channel'],r['domain'],r['normalization_scheme'],r['delta'],r['qref_input_GeV'],r['E_PeV']))
    with open(out/'kk_tower_screen_v0_94.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
    meta=dict(version='v0.94',correction='v0.93 summary overclaimed H<=1 for c1 external-Qref rows',
              pdf=a.pdf,nx=a.nx,energies_PeV=energies,deltas=deltas,nproc=nproc,
              contact_qrefs_GeV=cq,pkr_qrefs_GeV=pq,kk_spacing_GeV=a.kk_spacing_gev,
              continuum_assumption='N_eff(Q) proportional to Q^delta; mode masses/phase space neglected in this screen',
              canonical_compute_dir='~/Research/nutau-prop')
    (out/'metadata_v0_94.json').write_text(json.dumps(meta,indent=2))
    write_summary(rows,out)
    top=sorted(rows,key=lambda r:r['r_KK_over_SM'],reverse=True)[:15]
    print('\nTop corrected KK-screen points:')
    for r in top:
        print(f"{r['channel']:24s} {r['normalization_scheme']:27s} d={r['delta']} E={r['E_PeV']:4.0f} PeV Qref={r['qref_input_GeV']:7.1f} GeV enh={r['enhancement_KK_over_single']:8.3g} rKK={100*r['r_KK_over_SM']:8.3g}% {'PASS' if r['size_gate_5pct'] else 'FAIL'}")
    print('wrote',out)

if __name__=='__main__': main()
