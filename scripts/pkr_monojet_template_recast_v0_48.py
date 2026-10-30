#!/usr/bin/env python3
"""v0.48 CMS EXO-20-004 monojet public-template recast for PKR M=50 GeV.

HEPData itself may return HTTP 403 to scripted clients, so this calculation uses
an independent public GitHub mirror of the exact HEPData v2 CSV archive
(andlessa/LLP-CDFO).  No HEPData-server access is required.

Method
------
1. Read the official CMS 66-bin monojet simplified likelihood (data, background,
   covariance) and spin-1 signal templates from the mirrored HEPData archive.
2. Validate the likelihood normalization at the published vector-mediator
   coupling points mmed=100 and 250 GeV.
3. Morph the official 100-GeV vector recoil template to 50 GeV using the exact
   LO q qbar -> V g and q g -> V q mass dependence folded with CT18NLO PDFs.
   Only the 50/100 mass ratio is used, so overall QCD/EW normalization cancels.
4. Rescale the morphed template by gq^2 BR_inv for each v0.45 PKR width
   convention and evaluate the CMS simplified likelihood.

This is stronger than the v0.46 NWA anchor screen but still not a full detector
simulation at 50 GeV.  The 100->250 public-template closure test is used to
assign a conservative morphing uncertainty; a 15% floor is retained because
50 GeV is an extrapolation below the published coupling grid.
"""
from __future__ import annotations

import json, math, re, sys
from pathlib import Path
from numpy.polynomial.legendre import leggauss

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from nutau_prop.pdfs import LHAPDFNucleon

ROOT = Path(__file__).resolve().parents[1]
MIRROR = ROOT / "cache" / "LLP-CDFO"
CMS = MIRROR / "CMS-EXO-20-004" / "CMS_data" / "HEPData-ins1894408-v2-csv"
STATROOT = MIRROR
OUT = ROOT / "results" / "v0.48"
OUT.mkdir(parents=True, exist_ok=True)
IN45 = ROOT / "results" / "v0.45" / "width_br_dijet_gate.csv"

SIGFILE = CMS / "Signaltemplates,DMsimp,spin-1,Monojet.csv"
YIELDFILE = CMS / "Simplifiedlikelihood:Yields(Monojet).csv"
COVFILE = CMS / "Simplifiedlikelihood:covariancematrix(Monojet).csv"
GQFILE = CMS / "Couplinglimitsong_{q},vectormediator.csv"

MET_EDGES = np.array([250,280,310,340,370,400,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250,1400], float)
SQRTS = 13000.0
GQ_REF = 0.25
GDM_REF = 1.0
PDF_SET = "CT18NLO"
MORPH_FLOOR = 0.15


def require_inputs():
    miss=[p for p in [IN45,SIGFILE,YIELDFILE,COVFILE,GQFILE] if not p.exists()]
    if miss:
        raise RuntimeError("Missing public-mirror inputs:\n"+"\n".join(str(p) for p in miss)+"\nRun run_pkr_monojet_template_recast_v0_48.sh")


def vector_widths_over_m(m: float, mdm: float, gq: float, gdm: float=1.0):
    # Five effectively massless quarks at 50--250 GeV; top closed.
    gqq = 5.0*3.0*gq*gq/(12.0*math.pi)
    r=(mdm/m)**2
    if 4*r >= 1:
        ginv=0.0
    else:
        ginv=gdm*gdm/(12.0*math.pi)*(1+2*r)*math.sqrt(1-4*r)
    br=ginv/max(ginv+gqq,1e-300)
    return gqq,ginv,br


def parse_blocks(path: Path):
    lines=path.read_text(errors='replace').splitlines(); blocks=[]; meta={}; i=0
    while i < len(lines):
        ln=lines[i].strip()
        if ln.startswith('#:'):
            body=ln[2:].strip(); parts=[p.strip() for p in body.split(',')]
            key=parts[0]; val=next((p for p in reversed(parts[1:]) if p!=''), '')
            if 'Coupling type' in key: meta['coupling']=val
            elif '$m_{med}$' in key: meta['mmed']=float(val)
            elif '$m_{DM}$' in key: meta['mdm']=float(val)
            elif '$g_{q}$' in key: meta['gq']=float(val)
            elif '$g_{DM}$' in key: meta['gdm']=float(val)
            elif 'Data-taking period' in key: meta['period']=int(float(val))
            i+=1; continue
        if ln.startswith('Recoil (GeV),'):
            rows=[]; i+=1
            while i < len(lines):
                s=lines[i].strip()
                if not s or s.startswith('#:'): break
                p=[x.strip() for x in s.split(',')]
                try: rows.append(tuple(float(x) for x in p[:6]))
                except Exception: pass
                i+=1
            if rows and meta.get('coupling')=='Vector':
                blocks.append(dict(meta=dict(meta), rows=np.asarray(rows,float)))
            continue
        i+=1
    return blocks


def available_points(blocks):
    return sorted(set((b['meta']['mmed'],b['meta']['mdm']) for b in blocks if all(k in b['meta'] for k in ('mmed','mdm'))))


def choose_point(blocks, mmed: float, target_mdm: float):
    pts=available_points(blocks)
    cand=[p for p in pts if abs(p[0]-mmed)<1e-8 and p[1] < 0.5*mmed]
    if not cand: raise RuntimeError(f'No on-shell vector template at mmed={mmed}; available={pts[:20]}...')
    return min(cand,key=lambda p:abs(p[1]-target_mdm))


def template_lookup(blocks, point):
    m,md=point; out={}
    for b in blocks:
        z=b['meta']
        if z.get('coupling')=='Vector' and abs(z.get('mmed',-1)-m)<1e-8 and abs(z.get('mdm',-1)-md)<1e-8:
            per=z.get('period')
            if per is not None:
                out[int(per)]={(float(r[1]),float(r[2])):float(r[3]) for r in b['rows']}
    if not out: raise RuntimeError(f'No blocks for point {point}')
    return out


def load_cms_likelihood():
    # Follow the public validated recast exactly.
    bg=np.genfromtxt(YIELDFILE,delimiter=',',skip_header=5,skip_footer=67,names=True,dtype=None,encoding=None)
    dat=np.genfromtxt(YIELDFILE,delimiter=',',skip_header=73,names=True,dtype=None,encoding=None)
    cov=np.genfromtxt(COVFILE,delimiter=',',skip_header=5,names=True,dtype=None,encoding='UTF-8')
    sr=np.asarray(cov['Second_Bin'][:66])
    C=np.zeros((len(sr),len(sr)))
    for p in cov:
        ii=np.where(sr==p['First_Bin'])[0]; jj=np.where(sr==p['Second_Bin'])[0]
        if len(ii) and len(jj): C[ii[0],jj[0]]=float(p['Covariance'])
    nobs=np.array([int(dat[dat['Bin']==s]['Data_yield'][0]) for s in sr])
    nbg=np.array([float(bg[bg['Bin']==s]['Background_yield'][0]) for s in sr])
    return sr,nobs,nbg,C


def vector_from_template(tpl, sr):
    ns=[]
    for s in sr:
        txt=s.decode() if isinstance(s,bytes) else str(s)
        bits=txt.split('_')
        if len(bits)<3: raise RuntimeError(f'Unexpected SR label {txt}')
        per=int(bits[1]); ib=int(bits[2].replace('bin',''))
        key=(float(MET_EDGES[ib]),float(MET_EDGES[ib+1]))
        if per not in tpl or key not in tpl[per]: raise RuntimeError(f'Missing template {per} {key}')
        ns.append(tpl[per][key])
    return np.asarray(ns,float)


def ul_on_mu(ns,nobs,nbg,C):
    sys.path.insert(0,str(STATROOT))
    from statisticalTools.simplifiedLikelihoods import Data,UpperLimitComputer
    d=Data(observed=nobs,backgrounds=nbg,covariance=C,nsignal=np.asarray(ns,float),deltas_rel=0.0)
    return float(UpperLimitComputer().getUpperLimitOnMu(d))


def published_gq_observed():
    # Parse only the first (Observed) block.
    vals={}; started=False
    for ln in GQFILE.read_text(errors='replace').splitlines():
        if ln.startswith('#: Quantile,Observed'): started=True; continue
        if started and ln.startswith('#: Quantile,') and 'Observed' not in ln: break
        if started and re.match(r'^\s*\d',ln):
            p=ln.split(',')
            try: vals[float(p[0])]=float(p[1])
            except Exception: pass
    return vals


def solve_gq_from_mu(mu_ul,m,mdm_tpl):
    _,_,brref=vector_widths_over_m(m,mdm_tpl,GQ_REF,GDM_REF)
    target_mdm=m/3.0
    def scale(g):
        _,_,br=vector_widths_over_m(m,target_mdm,g,GDM_REF)
        return (g*g*br)/(GQ_REF*GQ_REF*brref)
    return brentq(lambda g:scale(g)-mu_ul,1e-5,GQ_REF)


def lo_weight_bin(M,plo,phi,pdf,npt=8,nyj=10,nyv=12):
    # Overall constants cancel in mass ratios. Formula uses x*f(x,Q), as returned by LHAPDF.
    xp,wp=leggauss(npt); xj,wj=leggauss(nyj); xv,wv=leggauss(nyv)
    pnodes=.5*(phi-plo)*xp+.5*(phi+plo); pweights=.5*(phi-plo)*wp
    yjnodes=2.4*xj; yjweights=2.4*wj
    yvnodes=5.0*xv; yvweights=5.0*wv
    total=0.0
    for pt,wpt in zip(pnodes,pweights):
        mt=math.sqrt(M*M+pt*pt); Q=mt; q2=Q*Q
        try: alphas=float(pdf.pdf.alphasQ(Q))
        except Exception: alphas=0.1
        for yj,wyj in zip(yjnodes,yjweights):
            ejp=math.exp(yj); ejm=1.0/ejp
            for yv,wyv in zip(yvnodes,yvweights):
                evp=math.exp(yv); evm=1.0/evp
                x1=(mt*evp+pt*ejp)/SQRTS; x2=(mt*evm+pt*ejm)/SQRTS
                if x1<=1e-7 or x2<=1e-7 or x1>=1 or x2>=1: continue
                sh=x1*x2*SQRTS*SQRTS
                th=M*M-x1*SQRTS*mt*evm
                uh=M*M-x2*SQRTS*mt*evp
                if th>=0 or uh>=0 or sh<=M*M: continue
                def xf(pid,x):
                    try:return max(float(pdf.xfx_proton(pid,x,q2)),0.0)
                    except Exception:return 0.0
                g1=xf(21,x1); g2=xf(21,x2)
                ann=0.0; q1=0.0; q2s=0.0
                for q in (1,2,3,4,5):
                    a1=xf(q,x1); b1=xf(-q,x1); a2=xf(q,x2); b2=xf(-q,x2)
                    ann += a1*b2+b1*a2
                    q1 += a1+b1; q2s += a2+b2
                Fann=(th/uh+uh/th+2*M*M*sh/(th*uh))
                Fqg=-(sh/uh+uh/sh+2*M*M*th/(sh*uh))
                Fgq=-(sh/th+th/sh+2*M*M*uh/(sh*th))
                hard=(8.0/9.0)*ann*Fann+(1.0/3.0)*(q1*g2*Fqg+g1*q2s*Fgq)
                if hard<=0 or not np.isfinite(hard): continue
                # d sigma / d pT^2 dyV dyj; integrate d pT => factor 2 pT.
                total += wpt*wyj*wyv*(2*pt)*alphas*hard/(sh*sh)
    return float(total)


def morph_ratios(pdf,Mnew,Mref=100.0):
    rows=[]; ratios=[]
    for lo,hi in zip(MET_EDGES[:-1],MET_EDGES[1:]):
        a=lo_weight_bin(Mnew,float(lo),float(hi),pdf)
        b=lo_weight_bin(Mref,float(lo),float(hi),pdf)
        r=a/max(b,1e-300); ratios.append(r)
        rows.append(dict(M_new_GeV=Mnew,M_ref_GeV=Mref,pt_lo=lo,pt_hi=hi,I_new=a,I_ref=b,ratio=r))
    return np.asarray(ratios,float),rows


def apply_bin_ratio(ns,sr,ratios):
    out=np.array(ns,float,copy=True)
    for i,s in enumerate(sr):
        txt=s.decode() if isinstance(s,bytes) else str(s); ib=int(txt.split('_')[2].replace('bin',''))
        out[i]*=ratios[ib]
    return out


def main():
    require_inputs(); OUT.mkdir(parents=True,exist_ok=True)
    blocks=parse_blocks(SIGFILE); sr,nobs,nbg,C=load_cms_likelihood(); pub=published_gq_observed()
    p100=choose_point(blocks,100.0,100.0/3.0); p250=choose_point(blocks,250.0,250.0/3.0)
    t100=template_lookup(blocks,p100); t250=template_lookup(blocks,p250)
    ns100=vector_from_template(t100,sr); ns250=vector_from_template(t250,sr)

    # Public-likelihood validation at published mass points.
    mu100=ul_on_mu(ns100,nobs,nbg,C); g100=solve_gq_from_mu(mu100,100.0,p100[1])
    mu250=ul_on_mu(ns250,nobs,nbg,C); g250=solve_gq_from_mu(mu250,250.0,p250[1])
    val=pd.DataFrame([
      dict(mmed=100.,template_mdm=p100[1],monojet_recast_gq95=g100,published_combined_gq95=pub.get(100.,np.nan),ratio=g100/pub.get(100.,np.nan)),
      dict(mmed=250.,template_mdm=p250[1],monojet_recast_gq95=g250,published_combined_gq95=pub.get(250.,np.nan),ratio=g250/pub.get(250.,np.nan)),
    ])
    val.to_csv(OUT/'public_likelihood_validation.csv',index=False)

    pdf=LHAPDFNucleon(PDF_SET,0)
    r250,mrows250=morph_ratios(pdf,250.0,100.0)
    pred250=apply_bin_ratio(ns100,sr,r250)
    # Correct for different reference invisible BRs in the official templates.
    _,_,br100=vector_widths_over_m(100.,p100[1],GQ_REF,GDM_REF); _,_,br250=vector_widths_over_m(250.,p250[1],GQ_REF,GDM_REF)
    pred250*=br250/max(br100,1e-300)
    nz=ns250>0
    frac=(pred250[nz]-ns250[nz])/ns250[nz]
    weighted_total=pred250.sum()/max(ns250.sum(),1e-300)-1.0
    rms=float(np.sqrt(np.mean(np.clip(frac,-5,5)**2))) if np.any(nz) else np.inf
    morph_unc=max(MORPH_FLOOR,abs(weighted_total),min(rms,0.50))

    r50,mrows50=morph_ratios(pdf,50.0,100.0)
    pd.DataFrame(mrows250+mrows50).to_csv(OUT/'lo_mass_morphing.csv',index=False)
    base50=apply_bin_ratio(ns100,sr,r50)

    d45=pd.read_csv(IN45); results=[]
    for _,x in d45.iterrows():
        # ns100 corresponds to gq=.25 and its template BR. Production scales as gq^2.
        scale=(float(x.gq)**2*float(x.BR_invisible))/(GQ_REF**2*br100)
        ns=base50*scale
        mu=ul_on_mu(ns,nobs,nbg,C); robs=1.0/max(mu,1e-300)
        rlo=robs*(1.0-morph_unc); rhi=robs*(1.0+morph_unc)
        if rlo>1: verdict='ROBUSTLY_EXCLUDED_MONOJET_ONLY'
        elif rhi<1: verdict='SURVIVES_MONOJET_ONLY_WITH_MORPH_UNCERTAINTY'
        else: verdict='MARGINAL_FULL_50GEV_MC_REQUIRED'
        results.append(dict(model=x.model,convention=x.convention,M_GeV=50.,gq=float(x.gq),BR_invisible=float(x.BR_invisible),mu95=float(mu),r_obs=float(robs),morph_uncertainty=float(morph_unc),r_low=float(rlo),r_high=float(rhi),verdict=verdict))
    R=pd.DataFrame(results); R.to_csv(OUT/'pkr_50gev_monojet_gate.csv',index=False)
    summary=[]
    for model,q in R.groupby('model'):
        vs=set(q.verdict)
        if vs=={'ROBUSTLY_EXCLUDED_MONOJET_ONLY'}: st='EXCLUDED_FOR_BOTH_WIDTH_CONVENTIONS'
        elif all(v=='SURVIVES_MONOJET_ONLY_WITH_MORPH_UNCERTAINTY' for v in vs): st='SURVIVES_MONOJET_ONLY_BUT_MONOV_AND_FULL_MC_REMAIN'
        else: st='WIDTH_OR_MORPH_UNCERTAINTY_DECIDES_FULL_MC_REQUIRED'
        summary.append(dict(model=model,min_r_low=float(q.r_low.min()),max_r_high=float(q.r_high.max()),status=st))
    S=pd.DataFrame(summary); S.to_csv(OUT/'summary.csv',index=False)

    meta=dict(version='v0.48',analysis='CMS EXO-20-004 monojet simplified likelihood',public_data_source='HEPData v2 mirrored in andlessa/LLP-CDFO',hepdata_server_403='v0.47 direct-server attempt retained as failed audit trail',signal_reference_100_GeV={'template_mdm':p100[1]},validation_250_GeV={'template_mdm':p250[1],'total_morph_fractional_error':weighted_total,'rms_bin_fractional_error':rms},morph_uncertainty=morph_unc,pdf=PDF_SET,channels='monojet only; mono-V omitted conservatively',method='official 100 GeV signal template morphed to 50 GeV with LO V+j CT18NLO mass ratio; exact CMS 66-bin simplified likelihood',limitations=['50 GeV is below published coupling grid','LO mass morphing, not generated 50-GeV detector sample','mono-V channel omitted','morphing uncertainty has 15% floor'])
    (OUT/'metadata.json').write_text(json.dumps(meta,indent=2))

    print('\n=== v0.48 public likelihood validation ==='); print(val.to_string(index=False))
    print('\n=== v0.48 100->250 morph closure ==='); print(f'total fractional error={weighted_total:+.4f}  RMS bin error={rms:.4f}  adopted morph uncertainty={morph_unc:.4f}')
    print('\n=== v0.48 PKR 50 GeV monojet gate ==='); print(R[['model','convention','r_obs','r_low','r_high','verdict']].to_string(index=False))
    print('\n=== v0.48 summary ==='); print(S.to_string(index=False))

if __name__=='__main__': main()
