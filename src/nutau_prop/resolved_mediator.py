"""Resolved derivative-Z' completion of the Case-6 momentum-growing kernel.

Physics definition
------------------
L ⊃ (c2/Lambda) (bar(nu^c) i<->partial_mu nu) Z'^mu
    + g_q (bar q gamma_mu q) Z'^mu.

Define kappa = c2*g_q/Lambda [GeV^-1].  At Q^2 << M_Z'^2 the contact
coefficient is
    C4_lowQ = kappa / M_Z'^2  [GeV^-3].

Thus, when matched to a chosen low-Q C4, the exact resolved kernel is the
Case-6 contact kernel times
    [M_Z'^2 / (M_Z'^2 + Q^2)]^2.
This is a physical propagator softening, unlike the hard sqrt(shat)<Lambda_UV
EFT-domain cut used for the contact operator.

Important caveat: the neutrino-Z' derivative vertex is itself dimension five.
The resolved propagator removes the contact-mediator approximation but does not
by itself provide the UV completion of that vertex.  An optional vertex
sqrt(shat) cut is therefore available only as a diagnostic; the primary v0.16
shape scan uses no hard cut so that the propagator effect is isolated.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from .case6 import Case6Config, Case6DIS
from .constants import M_N, GEV2_TO_CM2

C4_CASE6_BENCHMARK = 8.6e-10  # GeV^-3


@dataclass(frozen=True)
class ResolvedMediatorConfig:
    base: Case6Config = Case6Config()
    c4_lowq_gev_m3: float = C4_CASE6_BENCHMARK


class ResolvedMediatorDIS:
    """Exact-PDF derivative-Z' neutrino-quark DIS kernel."""

    def __init__(self, config: ResolvedMediatorConfig | None = None):
        self.cfg = config or ResolvedMediatorConfig()
        self.base = Case6DIS(self.cfg.base)

    def y_grid(self):
        return self.base.y_grid()

    def dsdy_sm(self, E, y, antinu=False):
        return self.base.dsdy_sm(E, y, antinu=antinu)

    @staticmethod
    def kappa_from_lowq_c4(c4_gev_m3: float, mediator_mass_tev: float) -> float:
        """kappa=c2*gq/Lambda [GeV^-1] required to match C4 at low Q."""
        m = 1.0e3 * float(mediator_mass_tev)
        return float(c4_gev_m3) * m * m

    @staticmethod
    def propagator_factor(q2_gev2: float, mediator_mass_tev: float) -> float:
        m2 = (1.0e3 * float(mediator_mass_tev)) ** 2
        return (m2 / (m2 + float(q2_gev2))) ** 2

    def _d2_resolved(self, E, x, y, mediator_mass_tev, lowq_scale=1.0):
        s = 2.0 * M_N * E
        q2 = x * y * s
        if q2 < self.cfg.base.qmin_gev**2:
            return 0.0
        p = self.base._pdf_components(x, q2)
        sum_xf = sum(p.values())
        c4 = self.cfg.c4_lowq_gev_m3 * float(lowq_scale)
        prop = self.propagator_factor(q2, mediator_mass_tev)
        # Same Case-6 Lorentz/y structure, now with the exact t-channel propagator.
        d2_gev_m2 = (
            (4.0 / math.pi)
            * c4**2
            * x
            * s**2
            * y
            * (1.0 - y)
            * sum_xf
            * prop
        )
        return d2_gev_m2 * GEV2_TO_CM2

    def _x_bounds(self, E, y, vertex_cut_tev=None):
        # Primary resolved-mediator scan: no hard contact-EFT cutoff.
        if vertex_cut_tev is None:
            return self.base._x_bounds(E, y, lambda_tev=None, cut="none")
        # Optional diagnostic for the still-effective dimension-5 neutrino-Z' vertex.
        return self.base._x_bounds(E, y, lambda_tev=vertex_cut_tev, cut="shat")

    def dsdy_np(self, E, y, mediator_mass_tev, lowq_scale=1.0, vertex_cut_tev=None):
        bounds = self._x_bounds(E, y, vertex_cut_tev=vertex_cut_tev)
        if bounds is None:
            return 0.0
        xmin, xmax = bounds
        return self.base._integrate_logx(
            lambda x: self._d2_resolved(E, x, y, mediator_mass_tev, lowq_scale),
            xmin,
            xmax,
        )

    def spectra(self, E, mediator_mass_tev, lowq_scale=1.0, antinu=False,
                vertex_cut_tev=None):
        edges, y, dy = self.y_grid()
        sm_cc = np.zeros_like(y)
        sm_nc = np.zeros_like(y)
        npds = np.zeros_like(y)
        for i, yy in enumerate(y):
            sm_cc[i], sm_nc[i] = self.dsdy_sm(E, yy, antinu=antinu)
            npds[i] = self.dsdy_np(
                E, yy, mediator_mass_tev, lowq_scale=lowq_scale,
                vertex_cut_tev=vertex_cut_tev,
            )
        return dict(edges=edges, y=y, dy=dy, sm_cc=sm_cc, sm_nc=sm_nc, np=npds)
