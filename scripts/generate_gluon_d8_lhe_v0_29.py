#!/usr/bin/env python3
"""Generate unweighted LHE events for gg -> tau-_L tau+_R from Q_l2G2D[33].

This is an independent open cross-check of the FeynRules/UFO implementation.
Canonical Murphy normalization:
  L = (C_G/Lambda^4) i(bar L3 gamma^mu <->D^nu L3) G^A_{mu rho} G^{A rho}_{nu}

Massless hard-process result, averaged over both initial-gluon helicities/colors,
for the single chiral final state tau-_L tau+_R:
  |M|^2_bar = (C_G/Lambda^4)^2 * shat^4/64 * (1-z^4)
  d shat_sigma/dz = (C_G/Lambda^4)^2 * shat^3/(2048*pi) * (1-z^4)
  shat_sigma = (C_G/Lambda^4)^2 * shat^3/(1280*pi)

The physical tau mass is retained in the event four-vectors; its effect in the
hard matrix element is negligible in the high-mass fiducial region used here.
"""
from __future__ import annotations
import argparse, math
from pathlib import Path
import numpy as np
import lhapdf

GEV2_TO_PB = 0.389379338e9
MTAU = 1.77686
SQRTS = 13000.0
S = SQRTS**2


def args():
    p=argparse.ArgumentParser()
    p.add_argument('--pdf',default='NNPDF31_lo_as_0130')
    p.add_argument('--member',type=int,default=0)
    p.add_argument('--C',type=float,default=-8*math.pi)
    p.add_argument('--lambda-gev',type=float,default=4000.)
    p.add_argument('--xi',type=float,default=1.0)
    p.add_argument('--mmin',type=float,default=100.)
    p.add_argument('--events',type=int,default=50000)
    p.add_argument('--prescan',type=int,default=200000)
    p.add_argument('--safety',type=float,default=2.0)
    p.add_argument('--seed',type=int,default=290529)
    p.add_argument('--out',default='results/v0.29/gluon_d8_C8pi.lhe')
    return p.parse_args()


def fgluon(pdf,x,q2):
    if not (0 < x < 1): return 0.0
    try: return max(0.0,float(pdf.xfxQ2(21,x,q2))/x)
    except Exception: return 0.0


def sigmahat_pb(M,C,lam):
    g8=C/lam**4
    return g8*g8*M**6/(1280*math.pi)*GEV2_TO_PB


def proposal_weight(pdf,u,Y,C,lam,du):
    M=math.exp(u); ymax=math.log(SQRTS/M)
    if ymax <= 0 or abs(Y)>ymax: return 0.0
    x1=M/SQRTS*math.exp(Y); x2=M/SQRTS*math.exp(-Y)
    f1=fgluon(pdf,x1,M*M); f2=fgluon(pdf,x2,M*M)
    # dx1 dx2 = 2 M^2/S du dY.  Proposal: u uniform and Y uniform in [-ymax,ymax].
    return du*(2*ymax)*f1*f2*(2*M*M/S)*sigmahat_pb(M,C,lam)


def sample_proposal(rng,umin,umax):
    u=rng.uniform(umin,umax); M=math.exp(u); ymax=math.log(SQRTS/M)
    Y=rng.uniform(-ymax,ymax)
    return u,Y


def sample_z(rng):
    # normalized p(z)=5/8*(1-z^4), -1<z<1
    while True:
        z=rng.uniform(-1,1)
        if rng.random() < 1-z**4: return z


def boost_z(E,pz,Y):
    c=math.cosh(Y); s=math.sinh(Y)
    return c*E+s*pz, s*E+c*pz


def event_lines(M,Y,z):
    E=M/2; p=math.sqrt(max(0,E*E-MTAU*MTAU)); st=math.sqrt(max(0,1-z*z))
    px=p*st; pz=p*z
    E1,pz1=boost_z(E,pz,Y); E2,pz2=boost_z(E,-pz,Y)
    x1=M/SQRTS*math.exp(Y); x2=M/SQRTS*math.exp(-Y)
    # Closed color loop for a color-singlet gg state. SPINUP = 2h: tau- L=-1, tau+ R=+1.
    return [
      f" 21 -1 0 0 501 502 0.0 0.0 {x1*SQRTS/2:.12e} {x1*SQRTS/2:.12e} 0.0 0.0 9.0",
      f" 21 -1 0 0 502 501 0.0 0.0 {-x2*SQRTS/2:.12e} {x2*SQRTS/2:.12e} 0.0 0.0 9.0",
      f" 15  1 1 2 0 0 {px:.12e} 0.0 {pz1:.12e} {E1:.12e} {MTAU:.12e} 0.0 -1.0",
      f"-15  1 1 2 0 0 {-px:.12e} 0.0 {pz2:.12e} {E2:.12e} {MTAU:.12e} 0.0 1.0",
    ]


def main():
    a=args(); rng=np.random.default_rng(a.seed); pdf=lhapdf.mkPDF(a.pdf,a.member)
    mmax=min(SQRTS*0.999, math.sqrt(a.xi)*a.lambda_gev)
    if mmax<=a.mmin: raise SystemExit('mmax <= mmin')
    umin,umax=math.log(a.mmin),math.log(mmax); du=umax-umin

    ws=np.empty(a.prescan)
    for i in range(a.prescan):
        u,Y=sample_proposal(rng,umin,umax)
        ws[i]=proposal_weight(pdf,u,Y,a.C,a.lambda_gev,du)
    xs=float(ws.mean()); err=float(ws.std(ddof=1)/math.sqrt(a.prescan)); wmax=float(ws.max()*a.safety)
    if xs<=0 or wmax<=0: raise SystemExit('zero cross section')
    print(f'PDF={a.pdf} C={a.C:.10g} Lambda={a.lambda_gev:g} xi={a.xi:g}')
    print(f'pp cross section (M>{a.mmin:g} GeV, M<{mmax:g} GeV) = {xs:.6g} +- {err:.2g} pb')
    print(f'prescan max={ws.max():.6g} pb; rejection ceiling={wmax:.6g} pb; estimated acceptance={xs/wmax:.3g}')

    accepted=[]; trials=0
    while len(accepted)<a.events:
        u,Y=sample_proposal(rng,umin,umax); w=proposal_weight(pdf,u,Y,a.C,a.lambda_gev,du); trials+=1
        if w>wmax:
            raise SystemExit(f'weight ceiling exceeded ({w} > {wmax}); increase --safety and restart')
        if rng.random()*wmax <= w:
            M=math.exp(u); accepted.append((M,Y,sample_z(rng)))
    print(f'accepted {a.events}/{trials} = {a.events/trials:.4f}')

    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w') as f:
        f.write('<LesHouchesEvents version="3.0">\n<header>\n')
        f.write('# Q_l2G2D[33] analytic generator v0.29\n')
        f.write(f'# C_G={a.C:.16g} Lambda_GeV={a.lambda_gev:.16g} xi={a.xi:.16g} PDF={a.pdf}\n')
        f.write('</header>\n<init>\n')
        # IDWTUP=3: unweighted positive events; XSECUP in pb.
        f.write(f' 2212 2212 6500.0 6500.0 0 0 0 0 3 1\n {xs:.12e} {err:.12e} {xs:.12e} 1\n')
        f.write('</init>\n')
        for M,Y,z in accepted:
            f.write('<event>\n')
            f.write(f' 4 1 1.0 {M:.12e} 0.0 0.0\n')
            f.write('\n'.join(event_lines(M,Y,z))+'\n')
            f.write('</event>\n')
        f.write('</LesHouchesEvents>\n')
    print('WROTE',out)

if __name__=='__main__': main()
