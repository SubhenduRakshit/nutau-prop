#!/usr/bin/env python3
"""Analytic two-state check: stable sterile state is not an irreversible sink.

This is a fixed-energy limiting system, NOT an Earth flux prediction. Its
purpose is to catch an omitted inverse transition before a long transport run.
"""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm

def solution(S,T,t):
    delta=np.hypot(S,2*T)
    if delta==0:return np.array([1.,0.])
    slow=(S+2*T-delta)/2
    fast=(S+2*T+delta)/2
    a=.5*(1-S/delta)*np.exp(-slow*t)+.5*(1+S/delta)*np.exp(-fast*t)
    n=T/delta*(np.exp(-slow*t)-np.exp(-fast*t))
    return np.array([a,n])

checks=[]
for S,T,t in [(1.,0.,3.),(0.,1.,2.),(1.,.01,2.),(1.,.00625,40.)]:
    G=np.array([[-S-T,T],[T,-T]])
    exact=solution(S,T,t); direct=expm(G*t)@np.array([1.,0.])
    assert np.allclose(exact,direct,rtol=1e-10,atol=1e-14)
    assert np.all(exact>=-1e-14)
    if S==0: assert np.isclose(exact.sum(),1.)
    checks.append(dict(S=S,T=T,column=t,active=float(exact[0]),sterile=float(exact[1]),
                       pure_sink_active=float(np.exp(-(S+T)*t)),SM_active=float(np.exp(-S*t))))
out=Path(__file__).resolve().parents[1]/'results/v0.87_constraints'
out.mkdir(parents=True,exist_ok=True)
result=dict(status='PASS',scope='exact two-state fixed-energy toy; no claim of a full Earth calculation',
    generator='[[-S-T,T],[T,-T]]',
    implication='For a stable N with the same Hermitian transition, inverse scattering must be included or bounded.',
    checks=checks)
(out/'inverse_transition_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
