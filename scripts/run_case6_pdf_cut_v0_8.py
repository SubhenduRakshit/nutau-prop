#!/usr/bin/env python3
"""v0.7: Case-6 dimension-7 neutrino-quark DIS with explicit EFT validity cuts.

Primary physics output is r_shatcut, using sqrt(shat)<Lambda.  r_Qcut is retained
as a more aggressive comparison.  CT10nlo is used to match the 2019 analysis.
"""
import csv, math
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
import lhapdf

PDFSET = "CT10nlo"
C4 = 8.6e-10                 # GeV^-3, collider-normalised benchmark
MN = 0.939                    # GeV
GF = 1.1663787e-5             # GeV^-2
MW, MZ = 80.379, 91.1876      # GeV
SW2 = 0.23122
GEV2_TO_CM2 = 0.389379338e-27
QMIN = 1.3                    # GeV
LAM_TEV = (1.0, 1.5, 2.0, 2.44, 3.0, 4.0)
EPEV_FULL = np.geomspace(0.1, 20000.0, 60)
EPEV_SHOW = (1.0, 10.0, 100.0, 300.0)
NX, NY = 220, 28

pdf = lhapdf.mkPDF(PDFSET, 0)
yn, yw = leggauss(NY)
xgrid = np.geomspace(1e-8, 0.999999, NX)

def xfq(pid, x, Q):
    try:
        return max(0.0, float(pdf.xfxQ(pid, x, Q)))
    except Exception:
        return 0.0

def pdfs_iso_x(x, Q):
    up, dn = xfq(2,x,Q), xfq(1,x,Q)
    ub, db = xfq(-2,x,Q), xfq(-1,x,Q)
    uiso = diso = 0.5*(up+dn)
    ubiso = dbiso = 0.5*(ub+db)
    return {2:uiso,-2:ubiso,1:diso,-1:dbiso,
            3:xfq(3,x,Q),-3:xfq(-3,x,Q),
            4:xfq(4,x,Q),-4:xfq(-4,x,Q),
            5:xfq(5,x,Q),-5:xfq(-5,x,Q)}

def ybounds(x, s, lam=None, cut="none"):
    ymin = QMIN**2/(x*s)
    if ymin >= 1.0: return None
    ymax = 1.0
    if cut == "Q":
        ymax = min(ymax, lam**2/(x*s))
    elif cut == "shat":
        if x*s >= lam**2: return None
    if ymax <= ymin: return None
    ys = 0.5*(yn+1)*(ymax-ymin)+ymin
    ws = 0.5*(ymax-ymin)*yw
    return ys, ws

def sigma_np(E, lam_TeV=None, cut="none"):
    # Case 6: |M|^2 = 64 C4^2 shat^3 y(1-y)
    # d2sigma/dxdy = |M|^2 f_j/(16*pi*x*s)
    s = 2*MN*E
    lam = None if lam_TeV is None else 1e3*lam_TeV
    vals = np.zeros_like(xgrid)
    for ix,x in enumerate(xgrid):
        bnd = ybounds(x,s,lam,cut)
        if bnd is None: continue
        ys,ws = bnd
        subtotal = 0.0
        for y,w in zip(ys,ws):
            Q = math.sqrt(x*y*s)
            sum_xf = sum(pdfs_iso_x(x,Q).values())
            d2 = (4/math.pi)*C4**2*x*s**2*y*(1-y)*sum_xf
            subtotal += w*d2
        vals[ix] = subtotal
    return float(np.trapezoid(vals,xgrid))*GEV2_TO_CM2

def sm_d2(E,x,y,antinu=False):
    s = 2*MN*E; Q2=x*y*s; Q=math.sqrt(Q2); p=pdfs_iso_x(x,Q)
    propW=MW**4/(MW**2+Q2)**2; propZ=MZ**4/(MZ**2+Q2)**2
    pref=2*GF**2*MN*E/math.pi; omy2=(1-y)**2
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
    return cc+pref*propZ*nc

def sigma_sm(E):
    s=2*MN*E; vn=np.zeros_like(xgrid); va=np.zeros_like(xgrid)
    for ix,x in enumerate(xgrid):
        bnd=ybounds(x,s)
        if bnd is None: continue
        ys,ws=bnd
        vn[ix]=sum(w*sm_d2(E,x,y,False) for y,w in zip(ys,ws))
        va[ix]=sum(w*sm_d2(E,x,y,True) for y,w in zip(ys,ws))
    sig=0.5*(np.trapezoid(vn,xgrid)+np.trapezoid(va,xgrid))
    return float(sig)*GEV2_TO_CM2

def main():
    out=Path("results/v0.8"); out.mkdir(parents=True,exist_ok=True); rows=[]
    for epev in EPEV_FULL:
        E=epev*1e6; sm=sigma_sm(E); nocut=sigma_np(E)
        for lam in LAM_TEV:
            qcut=sigma_np(E,lam,"Q"); scut=sigma_np(E,lam,"shat"); c4uv=C4*(lam*1e3)**3
            rows.append(dict(E_PeV=epev,Lambda_TeV=lam,c4_dimensionless=c4uv,c4_over_4pi=c4uv/(4*math.pi),
                sigma_SM_cm2=sm,sigma_NP_nocut_cm2=nocut,sigma_NP_Qcut_cm2=qcut,sigma_NP_shatcut_cm2=scut,
                r_Qcut=qcut/sm,r_shatcut=scut/sm))
    csvpath=out/"case6_pdf_cut_v0_8.csv"
    with csvpath.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    for lam in LAM_TEV:
        rq=[]; rs=[]
        for epev in EPEV_SHOW:
            E=epev*1e6; sm=sigma_sm(E); rq.append(sigma_np(E,lam,"Q")/sm); rs.append(sigma_np(E,lam,"shat")/sm)
        c4uv=C4*(lam*1e3)**3
        print(f"L={lam:4.2f}TeV c4={c4uv:6.2f} c4/4pi={c4uv/(4*math.pi):5.2f} "
              f"rQ(1,10,100,300PeV)={','.join(f'{z:.3g}' for z in rq)} rS={','.join(f'{z:.3g}' for z in rs)}")
    vv=[]
    for epev in (0.1,1.0):
        E=epev*1e6; vv.append(sigma_np(E)/sigma_sm(E))
    print("VALIDATION no-cut rNP/SM at 0.1,1 PeV="+",".join(f"{z:.3g}" for z in vv))
    print("WROTE",csvpath)
if __name__=="__main__": main()
