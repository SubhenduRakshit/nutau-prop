"""Analytic projection of tau-decay PDFs onto the log-energy neutrino grid.

For the implemented DRS/RR decay PDFs, each smooth z interval is polynomial.
Within a daughter-grid interval, the linear log-E deposition basis is affine in
ln z.  Matrix elements are therefore evaluated from closed-form primitives of
z^n and z^n ln z, with no inner z quadrature.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.sparse import csr_matrix

from .kernels import MTAU, MPI, MRHO, MA1


def _physical_bounds():
    b = [0.0, 0.3,
         1.0-(MA1/MTAU)**2,
         1.0-(MRHO/MTAU)**2,
         1.0-(MPI/MTAU)**2,
         1.0]
    return sorted(set(max(0.0, min(1.0, float(x))) for x in b))


def _primary_coeffs(pol, zmid):
    """Ascending polynomial coefficients c_n for sum c_n z^n before clip."""
    P = float(pol)
    c = np.zeros(4, float)
    g0 = np.array([5/3, 0.0, -3.0, 4/3], float)
    g1 = np.array([1/3, 0.0, -3.0, 8/3], float)
    c += 0.36*(g0 + P*g1)

    for mass, br, vector in ((MPI,0.12,False),(MRHO,0.26,True),(MA1,0.13,True)):
        r = (mass/MTAU)**2
        support = 1.0-r
        if zmid <= support:
            if vector:
                q = (1.0-2.0*r)/(1.0+2.0*r)
                c[0] += br*(1.0/(1.0-r) + P*q)
                c[1] += br*P*(-2.0*q/(1.0-r))
            else:
                c[0] += br*(1.0/(1.0-r) + P/(1.0-r))
                c[1] += br*P*(-2.0/(1.0-r)**2)
    if zmid <= 0.3:
        c[0] += 0.13/0.3
    return c


def _secondary_coeffs(pol, branching):
    P = float(pol); B = float(branching)
    return B*np.array([2.0-2.0*P, 12.0*P,
                       -6.0-18.0*P, 4.0+8.0*P], float)


def _poly_int(c, a, b):
    return sum(float(c[n])*(b**(n+1)-a**(n+1))/(n+1)
               for n in range(len(c)))


def _poly_log_int(c, a, b):
    def F(n, z):
        if z <= 0.0:
            return 0.0
        m = n+1.0
        return z**m*(math.log(z)/m - 1.0/(m*m))
    return sum(float(c[n])*(F(n,b)-F(n,a)) for n in range(len(c)))


def _positive_segments(pol, secondary=False, branching=0.18):
    """Return normalized positive polynomial pieces [(a,b,c),...].

    Numerical polynomial roots are used only to respect possible clip(...,0,None)
    boundaries; integration on each retained segment is analytic.
    """
    bounds = _physical_bounds() if not secondary else [0.0, 1.0]
    pieces = []
    for a,b in zip(bounds[:-1],bounds[1:]):
        c = (_secondary_coeffs(pol, branching) if secondary
             else _primary_coeffs(pol, 0.5*(a+b)))
        roots = np.polynomial.Polynomial(c).roots()
        cuts = [a,b]
        for r in roots:
            if abs(float(np.imag(r))) < 1e-11:
                x = float(np.real(r))
                if a+1e-13 < x < b-1e-13:
                    cuts.append(x)
        cuts = sorted(set(cuts))
        for u,v in zip(cuts[:-1],cuts[1:]):
            mid = 0.5*(u+v)
            val = float(np.polynomial.polynomial.polyval(mid,c))
            if val > 0.0:
                pieces.append([u,v,c.copy()])

    target = float(branching) if secondary else 1.0
    total = sum(_poly_int(c,a,b) for a,b,c in pieces)
    if not np.isfinite(total) or total <= 0.0:
        raise RuntimeError("analytic decay PDF has non-positive normalization")
    scale = target/total
    return [(a,b,c*scale) for a,b,c in pieces]


def projected_decay_row(E_tau, E_nu, polarization, secondary=False,
                        branching=0.18):
    """Exact projected daughter row on the represented neutrino log-E grid."""
    E = np.asarray(E_nu,float)
    if np.any(np.diff(E) <= 0):
        raise ValueError("E_nu must be strictly increasing")
    Ep = float(E_tau)
    row = np.zeros(len(E),float)
    if Ep <= 0.0 or len(E) < 2:
        return row
    lnE = np.log(E); lnEp = math.log(Ep)
    pieces = _positive_segments(polarization,secondary,branching)

    for i in range(len(E)-1):
        za = max(0.0, E[i]/Ep)
        zb = min(1.0, E[i+1]/Ep)
        if zb <= za:
            continue
        dl = lnE[i+1]-lnE[i]
        for pa,pb,c in pieces:
            a = max(za,pa); b = min(zb,pb)
            if b <= a:
                continue
            I0 = _poly_int(c,a,b)
            IL = _poly_log_int(c,a,b)
            left = ((lnE[i+1]-lnEp)/dl)*I0 - IL/dl
            right = ((lnEp-lnE[i])/dl)*I0 + IL/dl
            row[i] += left
            row[i+1] += right
    # Roundoff only; exact integrals of nonnegative basis functions are >=0.
    tiny = 5e-15*max(1.0,float(np.max(np.abs(row))))
    row[(row < 0.0) & (row > -tiny)] = 0.0
    return row


def projected_decay_kernel(E_tau, E_nu, polarization, secondary=False,
                           branching=0.18):
    K = np.vstack([projected_decay_row(Ep,E_nu,polarization,secondary,branching)
                   for Ep in np.asarray(E_tau,float)])
    return csr_matrix(K)
