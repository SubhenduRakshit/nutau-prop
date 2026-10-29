#!/usr/bin/env python3
"""v0.25 certification fold for active d=8 tau-neutrino twist-2 operators.

Operators (Alioli/Murphy/HighPT convention, no 1/2 in <->D):
  O_luD22[3322] = (L3bar gamma_mu i<->D_nu L3)(cRbar gamma^mu i<->D^nu cR)
  O_ldD22[3322] = (L3bar gamma_mu i<->D_nu L3)(sRbar gamma^mu i<->D^nu sR)
with L_NP = C/Lambda^4 O.

Computes dynamic-Q CT18NLO neutrino/antineutrino DIS including:
- exact pure d=8 term,
- SM Z interference with the positivity-allowed sign C<=0,
- SM CC+NC reference,
- partonic cuts shat < xi Lambda^2,
- fixed-helicity partial-wave unitarity diagnostics,
- endpoint-support fractions.

Requires: LHAPDF with CT18NLO installed.
No Earth propagation is performed.
"""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
import lhapdf

MN=0.939
GF=1.1663787e-5
MW=80.379
MZ=91.1876
SW2=0.23122
GEV2_TO_PB=0.389379338e9
QMIN=1.3

# Strict fixed-helicity elastic partial-wave bounds.
# q:     a1 = 7 C rho^2/(96 pi), a2=C rho^2/(160 pi)
# anti-q:a0 = 3 C rho^2/(16 pi), a1=C rho^2/(48 pi)
# with rho = shat/Lambda^2.  The anti-q J=0 channel is binding.
CUNIT_STRICT=8*math.pi/3       # |Re a0| <= 1/2 at rho=1
CUNIT_LOOSE=16*math.pi/3       # |a0| <= 1 at rho=1


def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument('--pdf',default='CT18NLO')
    p.add_argument('--member',type=int,default=0)
    p.add_argument('--lambda-tev',type=float,default=4.0)
    p.add_argument('--C',type=float,default=None,
                   help='dimensionless WC; positivity requires C<=0. Default: strict unitarity saturation')
    p.add_argument('--unitarity-mode',choices=['to_lambda','to_cut'],default='to_lambda')
    p.add_argument('--xis',default='1,0.5,0.25')
    p.add_argument('--energies-pev',default='1,10,100,300')
    p.add_argument('--nx',type=int,default=260)
    p.add_argument('--ny',type=int,default=40)
    p.add_argument('--outdir',default='results/v0.25')
    return p.parse_args()


def xfx(pdf,pid,x,q2):
    try:
        return max(0.0,float(pdf.xfxQ2(int(pid),float(x),float(q2))))
    except Exception:
        return 0.0


def iso_xf(pdf,x,q2):
    up,dn=xfx(pdf,2,x,q2),xfx(pdf,1,x,q2)
    ub,db=xfx(pdf,-2,x,q2),xfx(pdf,-1,x,q2)
    return {
       2:0.5*(up+dn), 1:0.5*(up+dn),
      -2:0.5*(ub+db),-1:0.5*(ub+db),
       3:xfx(pdf,3,x,q2),-3:xfx(pdf,-3,x,q2),
       4:xfx(pdf,4,x,q2),-4:xfx(pdf,-4,x,q2),
       5:xfx(pdf,5,x,q2),-5:xfx(pdf,-5,x,q2),
    }


def sm_d2_pb(E,x,y,pdf,antinu=False):
    s=2*MN*E; q2=x*y*s
    if q2<QMIN**2: return 0.0
    p=iso_xf(pdf,x,q2)
    propW=MW**4/(MW**2+q2)**2
    propZ=MZ**4/(MZ**2+q2)**2
    pref=2*GF**2*MN*E/math.pi
    omy2=(1-y)**2
    u,d,ss,c,b=p[2],p[1],p[3],p[4],p[5]
    ub,db,sb,cb,bb=p[-2],p[-1],p[-3],p[-4],p[-5]
    if not antinu: cc=(d+ss+b)+omy2*(ub+cb)
    else: cc=(db+sb+bb)+omy2*(u+c)
    cc*=pref*propW
    gLu,gRu=0.5-(2/3)*SW2,-(2/3)*SW2
    gLd,gRd=-0.5+(1/3)*SW2,(1/3)*SW2
    nc=0.0
    for q,qb,gL,gR in ((u,ub,gLu,gRu),(d,db,gLd,gRd),(ss,sb,gLd,gRd),(c,cb,gLu,gRu),(b,bb,gLd,gRd)):
        if not antinu: nc+=(gL*gL+gR*gR*omy2)*q+(gR*gR+gL*gL*omy2)*qb
        else: nc+=(gR*gR+gL*gL*omy2)*q+(gL*gL+gR*gR*omy2)*qb
    return (cc+pref*propZ*nc)*GEV2_TO_PB


def sm_nc_d2_pb(E,x,y,pdf,antinu=False):
    s=2*MN*E; q2=x*y*s
    if q2<QMIN**2: return 0.0
    p=iso_xf(pdf,x,q2)
    propZ=MZ**4/(MZ**2+q2)**2
    pref=2*GF**2*MN*E/math.pi
    omy2=(1-y)**2
    gLu,gRu=0.5-(2/3)*SW2,-(2/3)*SW2
    gLd,gRd=-0.5+(1/3)*SW2,(1/3)*SW2
    nc=0.0
    for q,qb,gL,gR in ((p[2],p[-2],gLu,gRu),(p[1],p[-1],gLd,gRd),(p[3],p[-3],gLd,gRd),(p[4],p[-4],gLu,gRu),(p[5],p[-5],gLd,gRd)):
        if not antinu: nc+=(gL*gL+gR*gR*omy2)*q+(gR*gR+gL*gL*omy2)*qb
        else: nc+=(gR*gR+gL*gL*omy2)*q+(gL*gL+gR*gR*omy2)*qb
    return pref*propZ*nc*GEV2_TO_PB


def np_d2_pb(E,x,y,pdf,pid,C,lam,xi,antinu=False):
    """Return (pure, interference) d2sigma/dxdy in pb.

    Formula uses an unpolarized f PDF (initial q spin averaged) and a fixed
    left-handed neutrino.  For antineutrinos the q/qbar angular weights swap.
    """
    s=2*MN*E; sh=x*s
    if sh>xi*lam*lam: return (0.0,0.0)
    q2=x*y*s
    if q2<QMIN**2: return (0.0,0.0)
    p=iso_xf(pdf,x,q2)
    qxf=p[pid]; qbxf=p[-pid]
    if pid==4: gR=-(2/3)*SW2
    elif pid==3: gR=(1/3)*SW2
    else: raise ValueError(pid)
    om=(1-y); A=(2-y)
    if not antinu:
        Fq=A*A*om*om; Fqb=A*A
        Iq=A*om*om; Iqb=A
    else:
        Fq=A*A; Fqb=A*A*om*om
        Iq=A; Iqb=A*om*om
    G8=C/lam**4
    pure=(G8*G8*s**3*x*x/(8*math.pi))*(Fq*qxf+Fqb*qbxf)
    prop_amp=MZ*MZ/(MZ*MZ+q2)
    inter=(math.sqrt(2)*GF*gR*prop_amp*G8*s*s*x/(2*math.pi))*(Iq*qxf+Iqb*qbxf)
    return pure*GEV2_TO_PB, inter*GEV2_TO_PB


def integrate(E,pdf,pid,C,lam,xi,nx,ny):
    xgrid=np.geomspace(1e-8,0.999999,nx)
    yn,yw=leggauss(ny)
    accum={k:np.zeros(nx) for k in ['sm','smnc','pure','int','smbar','smncbar','purebar','intbar']}
    for ix,x in enumerate(xgrid):
        s=2*MN*E
        ymin=QMIN**2/(x*s)
        if ymin>=1: continue
        ys=0.5*(yn+1)*(1-ymin)+ymin
        ws=0.5*(1-ymin)*yw
        for y,w in zip(ys,ws):
            accum['sm'][ix]+=w*sm_d2_pb(E,x,y,pdf,False)
            accum['smbar'][ix]+=w*sm_d2_pb(E,x,y,pdf,True)
            accum['smnc'][ix]+=w*sm_nc_d2_pb(E,x,y,pdf,False)
            accum['smncbar'][ix]+=w*sm_nc_d2_pb(E,x,y,pdf,True)
            pu,ii=np_d2_pb(E,x,y,pdf,pid,C,lam,xi,False)
            pb,ib=np_d2_pb(E,x,y,pdf,pid,C,lam,xi,True)
            accum['pure'][ix]+=w*pu; accum['int'][ix]+=w*ii
            accum['purebar'][ix]+=w*pb; accum['intbar'][ix]+=w*ib
    vals={k:float(np.trapezoid(v,xgrid)) for k,v in accum.items()}
    sm=0.5*(vals['sm']+vals['smbar']); smnc=0.5*(vals['smnc']+vals['smncbar'])
    pure=0.5*(vals['pure']+vals['purebar']); inter=0.5*(vals['int']+vals['intbar'])

    support_bins=[0.0,0.25,0.5,0.75,1.0000001]
    rho_contrib=np.zeros(4)
    for ix in range(len(xgrid)-1):
        xa,xb=xgrid[ix],xgrid[ix+1]; xm=math.sqrt(xa*xb)
        rho=xm*(2*MN*E)/(lam*lam)
        if rho>xi: continue
        ibin=np.searchsorted(support_bins,rho,side='right')-1
        if not (0<=ibin<4): continue
        ymin=QMIN**2/(xm*2*MN*E)
        if ymin>=1: continue
        ys=0.5*(yn+1)*(1-ymin)+ymin; ws=0.5*(1-ymin)*yw
        val=0.0
        for y,w in zip(ys,ws):
            pu,_=np_d2_pb(E,xm,y,pdf,pid,C,lam,xi,False)
            pb,_=np_d2_pb(E,xm,y,pdf,pid,C,lam,xi,True)
            val+=w*0.5*(pu+pb)
        rho_contrib[ibin]+=val*(xb-xa)
    stot=rho_contrib.sum()
    frac_gt25=(rho_contrib[1:].sum()/stot) if stot else 0
    frac_gt50=(rho_contrib[2:].sum()/stot) if stot else 0
    return dict(sm_pb=sm,sm_nc_pb=smnc,np_pure_pb=pure,np_int_pb=inter,
                np_total_pb=pure+inter,
                r_nc=(pure+inter)/smnc if smnc else float('nan'),
                r_total=(pure+inter)/sm if sm else float('nan'),
                pure_r_nc=pure/smnc if smnc else float('nan'),
                int_r_nc=inter/smnc if smnc else float('nan'),
                support_gt_rho025=frac_gt25,support_gt_rho050=frac_gt50)


def main():
    a=parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    pdf=lhapdf.mkPDF(a.pdf,a.member)
    xis=[float(z) for z in a.xis.split(',')]
    Es=[float(z) for z in a.energies_pev.split(',')]
    lam=a.lambda_tev*1e3
    rows=[]
    for xi in xis:
        if a.C is None:
            cabs=CUNIT_STRICT if a.unitarity_mode=='to_lambda' else CUNIT_STRICT/(xi*xi)
            C=-cabs
        else:
            C=a.C
            if C>0: raise SystemExit('Positivity requires C<=0 in the Alioli convention.')
        aq_a1=7*C*xi*xi/(96*math.pi)
        aq_a2=C*xi*xi/(160*math.pi)
        ab_a0=3*C*xi*xi/(16*math.pi)
        ab_a1=C*xi*xi/(48*math.pi)
        for name,pid in [('charm',4),('strange',3)]:
            for epev in Es:
                res=integrate(epev*1e6,pdf,pid,C,lam,xi,a.nx,a.ny)
                rows.append(dict(operator=name,E_PeV=epev,Lambda_TeV=a.lambda_tev,xi=xi,C=C,
                    G8_TeV_m4=C/(a.lambda_tev**4),
                    a1_q=aq_a1,a2_q=aq_a2,a0_qbar=ab_a0,a1_qbar=ab_a1,**res))
    csvp=out/'d8_uhe_dynamicQ_v0_25.csv'
    with csvp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    meta={
      'pdf':a.pdf,'member':a.member,'Lambda_TeV':a.lambda_tev,'unitarity_mode':a.unitarity_mode,
      'strict_C_bound_at_shat_eq_Lambda2':CUNIT_STRICT,
      'loose_C_bound_at_shat_eq_Lambda2':CUNIT_LOOSE,
      'positivity':'C<=0 for O8,lu partial2 and O8,ld partial2 (Alioli convention)',
      'cross_section_convention':'unpolarized quark PDF; fixed LH neutrino',
      'note':'Earth propagation not run.'
    }
    (out/'d8_uhe_dynamicQ_v0_25.json').write_text(json.dumps(meta,indent=2))
    print(f'WROTE {csvp}')
    for r in rows:
        if r['E_PeV'] in (10.0,100.0,300.0):
            print(f"{r['operator']:7s} xi={r['xi']:.2g} E={r['E_PeV']:g}PeV C={r['C']:.4g} "
                  f"dNC={100*r['r_nc']:+.3g}% dTot={100*r['r_total']:+.3g}% "
                  f"[int={100*r['int_r_nc']:+.3g}% pure={100*r['pure_r_nc']:.3g}%] "
                  f"support(rho>.25,.5)=({r['support_gt_rho025']:.3f},{r['support_gt_rho050']:.3f})")

if __name__=='__main__': main()
