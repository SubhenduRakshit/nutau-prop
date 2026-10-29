#!/usr/bin/env python3
"""Auxiliary v0.26 pure-NP cross-check for active d=8 tau-neutrino twist-2 operators.

This is deliberately secondary to scripts/check_d8_analytic_v0_25.py and
scripts/d8_uhe_cert_v0_25.py, which are used by run_d8_cert_v0_26.sh.

Important correction: the binding fixed-helicity elastic partial wave is the
antiquark J=0 channel, not the quark J=1 channel.  The DIS formula averages the
initial quark spin; unitarity must be applied to the interacting fixed-helicity
amplitude.  For the RH quark current the corresponding antiquark has helicity
-1/2, so nu_L qbar has lambda=0 and admits J=0.
"""
from __future__ import annotations
import csv, json, math
from pathlib import Path
import numpy as np
from nutau_prop.case6 import Case6Config, Case6DIS
from nutau_prop.constants import M_N, GEV2_TO_CM2

OUT=Path("results/v0.26")
PDFSET="CT18NLO"
LAMBDA_TEV=4.0
ENERGIES_PEV=(1.0,10.0,100.0,300.0)
XIS=(1.0,0.5,0.25)
FLAVORS={"c_R":4,"s_R":3}

# Fixed-helicity partial waves at rho=shat/Lambda^2=1.
# q: M=2 C Fq, Fq=(3+c)(1+c)/4=(7/4)d^1_11+(1/4)d^2_11
A1_Q_PER_C=7.0/(96.0*math.pi)
A2_Q_PER_C=1.0/(160.0*math.pi)
# qbar: M=C(3+c)=3 C P0 + C P1
A0_QBAR_PER_C=3.0/(16.0*math.pi)
A1_QBAR_PER_C=1.0/(48.0*math.pi)
C_MAX_STRICT=0.5/A0_QBAR_PER_C   # 8 pi / 3
C_MAX_LOOSE=1.0/A0_QBAR_PER_C    # 16 pi / 3


def angular_q(y): return (2.0-y)**2*(1.0-y)**2
def angular_qbar(y): return (2.0-y)**2


def analytic_self_checks():
    iq=31.0/30.0; iqb=7.0/3.0
    return {
      "Iq":iq,"Iqbar":iqb,"Iqbar_over_Iq":iqb/iq,
      "sea_gain_over_valence_q_only":(iq+iqb)/iq,
      "a1_q_per_C":A1_Q_PER_C,"a2_q_per_C":A2_Q_PER_C,
      "a0_qbar_per_C":A0_QBAR_PER_C,"a1_qbar_per_C":A1_QBAR_PER_C,
      "Cmax_strict_binding_qbar_J0":C_MAX_STRICT,
      "Cmax_loose_binding_qbar_J0":C_MAX_LOOSE,
      "binding_channel":"nu_L + anti-q_L, J=0"
    }


class D8DIS:
    def __init__(self):
        self.base=Case6DIS(Case6Config(pdf_set=PDFSET,nx=128))
    def _xbounds(self,E,y,xi,lamtev):
        s=2*M_N*E
        xmin=max(self.base.cfg.qmin_gev**2/(s*y),1e-9)
        xmax=min(1-1e-10,xi*(1e3*lamtev)**2/s)
        return None if xmin>=xmax else (xmin,xmax)
    def _d2(self,E,x,y,pid,C,lamtev,antinu=False):
        s=2*M_N*E; sh=x*s; q2=x*y*s
        if q2<self.base.cfg.qmin_gev**2: return 0.0
        xq=max(0.0,self.base.pdf.xfx_isoscalar(pid,x,q2))
        xb=max(0.0,self.base.pdf.xfx_isoscalar(-pid,x,q2))
        w=(angular_q(y)*xq+angular_qbar(y)*xb) if not antinu else (angular_qbar(y)*xq+angular_q(y)*xb)
        lam=1e3*lamtev
        return (C*C*sh**3/(8*math.pi*lam**8))*(w/x)*GEV2_TO_CM2
    def sigma_np(self,E,pid,C,lamtev,xi,antinu=False):
        _,ys,dys=self.base.y_grid(); total=0.0
        for y,dy in zip(ys,dys):
            b=self._xbounds(E,y,xi,lamtev)
            if b is None: continue
            total += self.base._integrate_logx(lambda x:self._d2(E,x,y,pid,C,lamtev,antinu),*b)*dy
        return total
    def sigma_sm(self,E,antinu=False):
        _,ys,dys=self.base.y_grid(); cc=nc=0.0
        for y,dy in zip(ys,dys):
            a,b=self.base.dsdy_sm(E,y,antinu=antinu); cc+=a*dy; nc+=b*dy
        return cc,nc


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    dis=D8DIS(); rows=[]; C=C_MAX_STRICT
    for ep in ENERGIES_PEV:
        E=ep*1e6; sn=dis.sigma_sm(E,False); sa=dis.sigma_sm(E,True)
        cc=.5*(sn[0]+sa[0]); nc=.5*(sn[1]+sa[1]); tot=cc+nc
        for xi in XIS:
            for name,pid in FLAVORS.items():
                npv=.5*(dis.sigma_np(E,pid,C,LAMBDA_TEV,xi,False)+dis.sigma_np(E,pid,C,LAMBDA_TEV,xi,True))
                rows.append(dict(E_PeV=ep,flavor=name,Lambda_TeV=LAMBDA_TEV,xi=xi,C_strict_unitarity=C,
                  sigma_SM_CC_cm2=cc,sigma_SM_NC_cm2=nc,sigma_SM_total_cm2=tot,sigma_NP_cm2=npv,
                  r_NP_over_NC=npv/nc,r_NP_over_total=npv/tot,
                  Cmax_if_unitarity_only_to_xi_domain=C_MAX_STRICT/(xi*xi)))
    p=OUT/"d8_dynamicQ_pureNP_aux_v0_26.csv"
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    status={"version":"v0.26-aux","analytic":analytic_self_checks(),
      "note":"Pure-NP auxiliary check only. Canonical v0.26 runner uses d8_uhe_cert_v0_25.py with interference and positivity."}
    (OUT/"d8_aux_status_v0_26.json").write_text(json.dumps(status,indent=2)+"\n")
    print(f"binding=anti-q J=0 Cmax_strict={C_MAX_STRICT:.9f}")
    print("WROTE",p)

if __name__=="__main__": main()
