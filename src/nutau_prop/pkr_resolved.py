"""v0.86 explicit one-channel chiral PKR-like vertices, with matching widths.

This does not change the historical Case6/ATZ implementation. Distinct
massless external states and two conjugate invisible decay channels are
specified, rather than inferred from a charge-conjugated flavour sum.
"""
from __future__ import annotations
import math

QUARK_MASSES = (0.00216, 0.00467, 0.0934, 1.27, 4.18, 172.69)

def operator_factor(operator: str) -> float:
    if operator == '4.4':
        return 1.0
    if operator == '4.3':
        return 4.0
    raise ValueError('operator must be 4.3 or 4.4')

def partonic_dsdy(shat, y, mass, a, gq, operator='4.4'):
    """GeV^-2; quark spin averaged, active neutrino helicity not averaged."""
    if shat <= 0 or mass <= 0 or not 0 <= y <= 1:
        raise ValueError('invalid kinematics')
    return operator_factor(operator)*(a*gq)**2*shat**2*y*(1-y)/(4*math.pi*(mass**2+shat*y)**2)

def nucleon_d2_from_xpdf(s_n, x, y, sum_xpdf, mass, a, gq, operator='4.4'):
    """d^2sigma/dxdy in GeV^-2; input is sum of x*f, not f."""
    return partonic_dsdy(x*s_n,y,mass,a,gq,operator)*sum_xpdf/x

def invisible_width(mass, a, operator='4.4'):
    """O+h.c., two distinct conjugate massless final states, in GeV."""
    return operator_factor(operator)*a*a*mass**3/(24*math.pi)

def quark_width_coefficient(mass):
    """Gamma_qq/gq^2 for universal vector coupling, including thresholds."""
    return sum(3*mass/(12*math.pi)*(1+2*(mq/mass)**2)*math.sqrt(1-4*(mq/mass)**2)
               for mq in QUARK_MASSES if mass > 2*mq)

def allowed_gq_squared(mass, a, operator, pmax, dijet_gq=None, width_fraction=0.1):
    """Analytically solve monotone NWA screens in z=gq^2, not a grid scan.

    P=z*I/(B*z+I) <= pmax. Optional D=z-P <= dijet_gq^2.
    These are screening constraints, not a new experimental likelihood.
    """
    I = invisible_width(mass,a,operator)
    B = quark_width_coefficient(mass)
    if I <= 0 or I >= width_fraction*mass:
        return 0.0
    bounds = [4*math.pi, (width_fraction*mass-I)/B]
    if I > B*pmax:
        bounds.append(pmax*I/(I-B*pmax))
    if dijet_gq is not None:
        D = dijet_gq**2
        bounds.append((D+math.sqrt(D*D+4*D*I/B))/2)
    return max(0.0,min(bounds))
