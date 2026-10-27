"""Momentum-growing contact-operator kernels beyond Case 6.

v0.14 analytic/numerical prescreen for Pandey-Karmakar-Rakshit Cases 5 and 7.

Conventions are matched to the Case-6 module and CT10nlo implementation.
The coefficients are the 2019 monojet benchmark limits:
  Case 5: C5 = 1.8e-10 GeV^-3
  Case 7: C7 = 1.6e-10 GeV^-3

Case 5 is the derivative tensor/vector operator.  Its DIS matrix-element
magnitude is |M|^2 = 32 C5^2 shat^3 y(1-y), i.e. half the Case-6
normalisation with the same y shape.  The external-state/chirality issue is
*not* hidden here: this module only computes the inclusive hard kernel.

Case 7 is the neutrino-gluon Rayleigh operator.  In the coefficient convention
L = C7 (nu-bilinear) G_{mu nu} G^{mu nu}, the DIS matrix-element magnitude is
|M|^2 = 16 C7^2 shat^3 y^3.  This follows from the standard low-energy EFT
normalisation after absorbing alpha_s/(12 pi) into C7.

The hard sqrt(shat)<Lambda_UV restriction is an EFT-domain diagnostic only,
not a physical turnover.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from .case6 import Case6Config, Case6DIS
from .constants import M_N, GEV2_TO_CM2

C5_BENCHMARK = 1.8e-10  # GeV^-3
C7_BENCHMARK = 1.6e-10  # GeV^-3


def lambda_4pi_tev(C_gev_m3: float) -> float:
    """Scale where c=C Lambda^3 reaches 4 pi."""
    return (4.0 * math.pi / float(C_gev_m3)) ** (1.0 / 3.0) / 1.0e3


@dataclass(frozen=True)
class MomentumOperatorConfig:
    base: Case6Config = Case6Config()
    c5_gev_m3: float = C5_BENCHMARK
    c7_gev_m3: float = C7_BENCHMARK


class MomentumOperatorDIS:
    """Reuse the validated Case-6 PDF/SM machinery for Cases 5 and 7."""

    def __init__(self, config: MomentumOperatorConfig | None = None):
        self.cfg = config or MomentumOperatorConfig()
        self.base = Case6DIS(self.cfg.base)

    def y_grid(self):
        return self.base.y_grid()

    def dsdy_sm(self, E, y, antinu=False):
        return self.base.dsdy_sm(E, y, antinu=antinu)

    def _d2_case5(self, E, x, y):
        s = 2.0 * M_N * E
        q2 = x * y * s
        if q2 < self.cfg.base.qmin_gev**2:
            return 0.0
        p = self.base._pdf_components(x, q2)
        sum_xf = sum(p.values())
        # |M5|^2 = 32 C5^2 shat^3 y(1-y)
        d2 = (
            (2.0 / math.pi)
            * self.cfg.c5_gev_m3**2
            * x
            * s**2
            * y
            * (1.0 - y)
            * sum_xf
        )
        return d2 * GEV2_TO_CM2

    def _d2_case7(self, E, x, y):
        s = 2.0 * M_N * E
        q2 = x * y * s
        if q2 < self.cfg.base.qmin_gev**2:
            return 0.0
        xg = self.base.pdf.xfx_isoscalar(21, x, q2)
        # |M7|^2 = 16 C7^2 shat^3 y^3.
        # With LHAPDF returning x*g, d2 sigma/dxdy = C7^2/pi*x*s^2*y^3*(xg).
        d2 = (
            (1.0 / math.pi)
            * self.cfg.c7_gev_m3**2
            * x
            * s**2
            * y**3
            * xg
        )
        return d2 * GEV2_TO_CM2

    def dsdy_np(self, operator, E, y, lambda_tev, cut="shat"):
        bounds = self.base._x_bounds(E, y, lambda_tev=lambda_tev, cut=cut)
        if bounds is None:
            return 0.0
        xmin, xmax = bounds
        if int(operator) == 5:
            f = lambda x: self._d2_case5(E, x, y)
        elif int(operator) == 7:
            f = lambda x: self._d2_case7(E, x, y)
        else:
            raise ValueError("operator must be 5 or 7")
        return self.base._integrate_logx(f, xmin, xmax)

    def spectra(self, operator, E, lambda_tev, cut="shat"):
        edges, y, dy = self.y_grid()
        npds = np.array([self.dsdy_np(operator, E, yy, lambda_tev, cut=cut) for yy in y])
        return dict(edges=edges, y=y, dy=dy, np=npds)
