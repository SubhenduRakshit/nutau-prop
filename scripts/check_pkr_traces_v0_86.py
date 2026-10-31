#!/usr/bin/env python3
"""Independent exact traces: no use of the DIS implementation being audited."""
import json
from pathlib import Path
import sympy as S

I = S.I
zero = S.zeros(2)
eye = S.eye(2)
pauli = [S.Matrix([[0, 1], [1, 0]]), S.Matrix([[0, -I], [I, 0]]), S.diag(1, -1)]
g = [S.diag(1, 1, -1, -1)] + [zero.row_join(s).col_join((-s).row_join(zero)) for s in pauli]
metric = [1, -1, -1, -1]
g5 = I*g[0]*g[1]*g[2]*g[3]
PL = (S.eye(4)-g5)/2
sigma = [[I*(g[m]*g[n]-g[n]*g[m])/2 for n in range(4)] for m in range(4)]
C = I*g[2]*g[0]
assert C.T == -C
assert (C*PL).T == -C*PL
for m in range(4):
    for n in range(4):
        assert (C*sigma[m][n]*PL).T == C*sigma[m][n]*PL

def slash(v):
    return sum((metric[m]*v[m]*g[m] for m in range(4)), S.zeros(4))

def dot(v, w):
    return sum(metric[m]*v[m]*w[m] for m in range(4))

def adj(v):
    return g[0]*v.conjugate().T*g[0]

def vertices(K, q):
    d = [K[m]*PL for m in range(4)]
    t = [2*sum((metric[n]*q[n]*sigma[n][m] for n in range(4)), S.zeros(4))*PL for m in range(4)]
    return d, t

def lepton(kout, kin, v):
    return [[S.simplify(S.trace(slash(kout)*v[m]*slash(kin)*adj(v[n]))) for n in range(4)] for m in range(4)]

checks = []
r = S.Symbol('r', positive=True)
for cosine, sine in [(S.Rational(3,5),S.Rational(4,5)), (S.Rational(-3,5),S.Rational(4,5)), (0,1),
                     ((1-r*r)/(1+r*r),2*r/(1+r*r))]:
    k, p = S.Matrix([1,0,0,1]), S.Matrix([1,0,0,-1])
    kp, pp = S.Matrix([1,sine,0,cosine]), S.Matrix([1,-sine,0,-cosine])
    q, K = k-kp, k+kp
    s, t, u = dot(k+p,k+p), dot(q,q), dot(k-pp,k-pp)
    H = [[S.trace(slash(pp)*g[m]*slash(p)*g[n])/2 for n in range(4)] for m in range(4)]
    vals = []
    for v in vertices(K,q):
        L = lepton(kp,k,v)
        vals.append(S.simplify(sum(metric[m]*metric[n]*L[m][n]*H[m][n] for m in range(4) for n in range(4))))
        ward = S.simplify(sum(metric[m]*metric[n]*q[m]*q[n]*L[m][n] for m in range(4) for n in range(4)))
        assert ward == 0, ward
    assert S.simplify(vals[0]-4*s*t*u) == 0, vals
    assert S.simplify(vals[1]-4*vals[0]) == 0, vals
    checks.append(dict(cosine=str(cosine), derivative=str(vals[0]), tensor=str(vals[1]), ratio=4, ward=0))

mass = S.Symbol('M', positive=True)
k, kp = S.Matrix([mass/2,0,0,mass/2]), S.Matrix([mass/2,0,0,-mass/2])
decay = []
for v in vertices(k-kp,k+kp):
    L = lepton(k,kp,v)
    decay.append(S.simplify(sum(L[m][m] for m in [1,2,3])/3))
assert S.simplify(decay[0]-mass**4/3) == 0
assert S.simplify(decay[1]-4*mass**4/3) == 0
# Two distinct conjugate channels; each massless two-body width is |M|^2/(16 pi M).
widths = [S.simplify(2*x/(16*S.pi*mass)) for x in decay]
y = S.Symbol('y', nonnegative=True)
assert S.integrate(y*(1-y),(y,0,1)) == S.Rational(1,6)
result = dict(version='v0.86', status='PASS', scattering=checks,
              charge_conjugated_flavour_selection='derivative and tensor bilinears antisymmetric in flavour; diagonal vanishes',
              decay_spin_averaged_amplitude_squared_over_a2=list(map(str,decay)),
              total_conjugate_width_over_a2=list(map(str,widths)),
              scope='one distinct-state chiral transition and conjugate; massless; no initial neutrino spin average')
out = Path(__file__).resolve().parents[1]/'results/v0.86_normalization'
out.mkdir(parents=True,exist_ok=True)
(out/'exact_traces.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
