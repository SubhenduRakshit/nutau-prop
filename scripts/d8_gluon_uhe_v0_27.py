#!/usr/bin/env python3
"""Exact UHE screen for Q_l2G2D[33].

Operator:
  Q = i (L3bar gamma^mu <->D^nu L3) G^A_{mu rho} G^{A rho}_{nu}
  L_NP = C/Lambda^4 Q

Partonic result (fixed LH neutrino, unpolarized initial gluon):
  d shat_sigma/dy = C^2 shat^3/(8 pi Lambda^8)
                     * (1-y) * [1+(1-y)^2]

Strict Jacob-Wick partial waves:
  a_1/2 = - C rho^2/(16 pi)
  a_3/2 = - C rho^2/(32 pi), rho=shat/Lambda^2.
Thus |Re a_J|<=1/2 gives |C| rho^2 <= 8 pi.
Published positivity in the Murphy convention gives C<=0.

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
CUNIT_STRICT=8*math.pi


def parse_args():
    p=argparse.ArgumentParser()
    p.add_argument('--pdf',default='CT18NLO')
    p.add_argument('--member',type=int,default=0)
    p.add_argument('--lambda-tev',type=float,default=4.0)
    p.add_argument('--C',type=float,default=None)
    p.add_argument('--unitarity-mode',choices=['to_lambda','to_cut'],default='to_lambda')
    p.add_argument('--xis',default='1,0.5,0.25')
    p.add_argument('--energies-pev',default='1,10,100,300')
    p.add_argument('--nx',type=int,default=300)
    p.add_argument('--ny',type=int,default=48)
    p.add_argument('--outdir',default='results/v0.27')
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


def np_gluon_d2_pb(E,x,y,pdf,C,lam,xi):
    """d2sigma/dxdy in pb for the pure Q_l2G2D contribution."""
    s=2*MN*E; sh=x*s
    if sh>xi*lam*lam: return 0.0
    q2=x*y*s
    if q2<QMIN**2: return 0.0
    gxf=xfx(pdf,21,x,q2)  # x*g(x,Q2)
    om=1-y
    F=om*(1+om*om)
    G8=C/lam**4
    # convolution: g(x)=xg/x and shat^3=(xs)^3 -> s^3*x^2*[xg]
    val=(G8*G8*s**3*x*x/(8*math.pi))*F*gxf
    return val*GEV2_TO_PB


def integrate(E,pdf,C,lam,xi,nx,ny):
    xgrid=np.geomspace(1e-8,0.999999,nx)
    yn,yw=leggauss(ny)
    acc={k:np.zeros(nx) for k in ['sm','smnc','smbar','smncbar','np']}
    for ix,x in enumerate(xgrid):
        s=2*MN*E
        ymin=QMIN**2/(x*s)
        if ymin>=1: continue
        ys=0.5*(yn+1)*(1-ymin)+ymin
        ws=0.5*(1-ymin)*yw
        for y,w in zip(ys,ws):
            acc['sm'][ix]+=w*sm_d2_pb(E,x,y,pdf,False)
            acc['smbar'][ix]+=w*sm_d2_pb(E,x,y,pdf,True)
            acc['smnc'][ix]+=w*sm_nc_d2_pb(E,x,y,pdf,False)
            acc['smncbar'][ix]+=w*sm_nc_d2_pb(E,x,y,pdf,True)
            acc['np'][ix]+=w*np_gluon_d2_pb(E,x,y,pdf,C,lam,xi)
    vals={k:float(np.trapezoid(v,xgrid)) for k,v in acc.items()}
    sm=0.5*(vals['sm']+vals['smbar'])
    smnc=0.5*(vals['smnc']+vals['smncbar'])
    npv=vals['np']

    bins=[0.0,0.25,0.5,0.75,1.0000001]
    rho_contrib=np.zeros(4)
    for ix in range(len(xgrid)-1):
        xa,xb=xgrid[ix],xgrid[ix+1]; xm=math.sqrt(xa*xb)
        rho=xm*(2*MN*E)/(lam*lam)
        if rho>xi: continue
        ib=np.searchsorted(bins,rho,side='right')-1
        if not (0<=ib<4): continue
        ymin=QMIN**2/(xm*2*MN*E)
        if ymin>=1: continue
        ys=0.5*(yn+1)*(1-ymin)+ymin
        ws=0.5*(1-ymin)*yw
        val=sum(w*np_gluon_d2_pb(E,xm,y,pdf,C,lam,xi) for y,w in zip(ys,ws))
        rho_contrib[ib]+=val*(xb-xa)
    stot=rho_contrib.sum()
    return dict(sm_pb=sm,sm_nc_pb=smnc,np_pure_pb=npv,
                r_nc=npv/smnc if smnc else float('nan'),
                r_total=npv/sm if sm else float('nan'),
                support_gt_rho025=(rho_contrib[1:].sum()/stot if stot else 0.0),
                support_gt_rho050=(rho_contrib[2:].sum()/stot if stot else 0.0))


def main():
    a=parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    pdf=lhapdf.mkPDF(a.pdf,a.member)
    lam=a.lambda_tev*1e3
    xis=[float(z) for z in a.xis.split(',')]
    Es=[float(z) for z in a.energies_pev.split(',')]
    rows=[]
    for xi in xis:
        if a.C is None:
            cabs=CUNIT_STRICT if a.unitarity_mode=='to_lambda' else CUNIT_STRICT/(xi*xi)
            C=-cabs
        else:
            C=a.C
            if C>0: raise SystemExit('Positivity requires C<=0 for Q_l2G2D in the Murphy convention.')
        a12=-C*xi*xi/(16*math.pi)
        a32=-C*xi*xi/(32*math.pi)
        for epev in Es:
            res=integrate(epev*1e6,pdf,C,lam,xi,a.nx,a.ny)
            rows.append(dict(operator='l2G2D_33',E_PeV=epev,Lambda_TeV=a.lambda_tev,
                xi=xi,C=C,G8_TeV_m4=C/(a.lambda_tev**4),a_half=a12,a_threehalf=a32,**res))
    csvp=out/'d8_gluon_uhe_v0_27.csv'
    with csvp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    meta={
      'operator':'Q_l2G2D[33] = i(L3bar gamma^mu <->D^nu L3) G^A_{mu rho} G^{A rho}_{nu}',
      'pdf':a.pdf,'member':a.member,'Lambda_TeV':a.lambda_tev,
      'unitarity_mode':a.unitarity_mode,'strict_C_bound_at_shat_eq_Lambda2':CUNIT_STRICT,
      'positivity':'C<=0 (Murphy basis; Chala-Li PRD109 065015 Appendix B)',
      'partonic':'dsigma/dy=C^2 shat^3/(8pi Lambda^8)*(1-y)*(1+(1-y)^2)',
      'interference':'no tree-level SM elastic nu-g amplitude; omitted in first screen',
      'note':'Earth propagation not run.'
    }
    (out/'d8_gluon_uhe_v0_27.json').write_text(json.dumps(meta,indent=2))
    print(f'WROTE {csvp}')
    for r in rows:
        print(f"xi={r['xi']:.2g} E={r['E_PeV']:g}PeV C={r['C']:.5g} "
              f"dNC={100*r['r_nc']:.4g}% dTot={100*r['r_total']:.4g}% "
              f"support(rho>.25,.5)=({r['support_gt_rho025']:.3f},{r['support_gt_rho050']:.3f})")

if __name__=='__main__':
    main()
