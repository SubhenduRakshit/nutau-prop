"""Differential DIS kernel for the momentum-growing Case-6 operator.

This module keeps the v0.7/v0.8 normalization and CT10nlo convention but
returns d sigma / dy, so the same interaction kernel can be used in both
Earth transport and detector-level NC shower reweighting.

Primary EFT-domain prescription: sqrt(shat) < Lambda_UV.
The Q < Lambda_UV option is retained only as an aggressive diagnostic.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

from .constants import GF, M_N, M_W, M_Z, GEV2_TO_CM2
from .pdfs import LHAPDFNucleon

C4_BENCHMARK = 8.6e-10  # GeV^-3, Pandey-Karmakar-Rakshit Case 6 benchmark
SIN2W = 0.23122
UP_TYPES = (2, 4)
DOWN_TYPES = (1, 3, 5)


@dataclass(frozen=True)
class Case6Config:
    pdf_set: str = "CT10nlo"
    member: int = 0
    c4_gev_m3: float = C4_BENCHMARK
    qmin_gev: float = 1.3
    nx: int = 96
    y_min: float = 1.0e-5
    y_split: float = 0.08
    y_log_bins: int = 36
    y_linear_bins: int = 36
    w2_min_gev2: float = 0.0  # zero reproduces the v0.8 high-energy setup


class Case6DIS:
    def __init__(self, config: Case6Config | None = None):
        self.cfg = config or Case6Config()
        self.pdf = LHAPDFNucleon(self.cfg.pdf_set, self.cfg.member)
        self.nodes, self.weights = leggauss(self.cfg.nx)

    def y_grid(self):
        e1 = np.geomspace(self.cfg.y_min, self.cfg.y_split, self.cfg.y_log_bins + 1)
        e2 = np.linspace(self.cfg.y_split, 1.0 - 1.0e-6, self.cfg.y_linear_bins + 1)
        edges = np.unique(np.r_[e1, e2[1:]])
        y = 0.5 * (edges[:-1] + edges[1:])
        dy = np.diff(edges)
        return edges, y, dy

    def _x_bounds(self, E, y, lambda_tev=None, cut="none"):
        s = 2.0 * M_N * E
        xmin = self.cfg.qmin_gev**2 / (s * y)
        xmax = 1.0 - 1.0e-10
        if self.cfg.w2_min_gev2 > 0:
            xmax = min(
                xmax,
                1.0 - max(0.0, self.cfg.w2_min_gev2 - M_N * M_N) / (s * y),
            )
        if lambda_tev is not None:
            lam2 = (1.0e3 * float(lambda_tev)) ** 2
            if cut == "shat":
                xmax = min(xmax, lam2 / s)
            elif cut == "Q":
                xmax = min(xmax, lam2 / (s * y))
            elif cut != "none":
                raise ValueError("cut must be none, shat, or Q")
        xmin = max(xmin, 1.0e-9)
        if xmin >= xmax:
            return None
        return xmin, xmax

    def _integrate_logx(self, func, xmin, xmax):
        a, b = math.log(xmin), math.log(xmax)
        t = 0.5 * (self.nodes + 1.0) * (b - a) + a
        x = np.exp(t)
        jac = x * 0.5 * (b - a)
        values = np.asarray([func(xx) for xx in x], dtype=float)
        return float(np.sum(self.weights * values * jac))

    def _pdf_components(self, x, q2):
        p = {pid: self.pdf.xfx_isoscalar(pid, x, q2) for pid in range(-5, 6) if pid != 0}
        return p

    @staticmethod
    def _nc_couplings(pid):
        if abs(pid) in UP_TYPES:
            return 0.5 - (2.0 / 3.0) * SIN2W, -(2.0 / 3.0) * SIN2W
        return -0.5 + (1.0 / 3.0) * SIN2W, (1.0 / 3.0) * SIN2W

    def _d2_sm_components(self, E, x, y, antinu=False):
        s = 2.0 * M_N * E
        q2 = x * y * s
        if q2 < self.cfg.qmin_gev**2:
            return 0.0, 0.0
        p = self._pdf_components(x, q2)
        omy2 = (1.0 - y) ** 2
        pref = 2.0 * GF * GF * M_N * E / math.pi
        prop_w = M_W**4 / (M_W**2 + q2) ** 2
        prop_z = M_Z**4 / (M_Z**2 + q2) ** 2

        down = sum(p[k] for k in DOWN_TYPES)
        up = sum(p[k] for k in UP_TYPES)
        adown = sum(p[-k] for k in DOWN_TYPES)
        aup = sum(p[-k] for k in UP_TYPES)
        if not antinu:
            cc_weight = down + aup * omy2
        else:
            cc_weight = adown + up * omy2
        cc = pref * prop_w * cc_weight

        nc_weight = 0.0
        for pid in (1, 2, 3, 4, 5):
            q, qb = p[pid], p[-pid]
            gL, gR = self._nc_couplings(pid)
            if not antinu:
                nc_weight += (gL * gL + gR * gR * omy2) * q
                nc_weight += (gR * gR + gL * gL * omy2) * qb
            else:
                nc_weight += (gR * gR + gL * gL * omy2) * q
                nc_weight += (gL * gL + gR * gR * omy2) * qb
        nc = pref * prop_z * nc_weight
        return cc * GEV2_TO_CM2, nc * GEV2_TO_CM2

    def _d2_np(self, E, x, y):
        s = 2.0 * M_N * E
        q2 = x * y * s
        if q2 < self.cfg.qmin_gev**2:
            return 0.0
        p = self._pdf_components(x, q2)
        sum_xf = sum(p.values())
        d2_gev_m2 = (
            (4.0 / math.pi)
            * self.cfg.c4_gev_m3**2
            * x
            * s**2
            * y
            * (1.0 - y)
            * sum_xf
        )
        return d2_gev_m2 * GEV2_TO_CM2

    def dsdy_sm(self, E, y, antinu=False):
        bounds = self._x_bounds(E, y)
        if bounds is None:
            return 0.0, 0.0
        xmin, xmax = bounds
        cc = self._integrate_logx(lambda x: self._d2_sm_components(E, x, y, antinu)[0], xmin, xmax)
        nc = self._integrate_logx(lambda x: self._d2_sm_components(E, x, y, antinu)[1], xmin, xmax)
        return cc, nc

    def dsdy_np(self, E, y, lambda_tev, cut="shat"):
        bounds = self._x_bounds(E, y, lambda_tev=lambda_tev, cut=cut)
        if bounds is None:
            return 0.0
        xmin, xmax = bounds
        return self._integrate_logx(lambda x: self._d2_np(E, x, y), xmin, xmax)

    def spectra(self, E, lambdas_tev=(1.5, 2.0, 2.44, 3.0), antinu=False, cut="shat"):
        edges, y, dy = self.y_grid()
        sm_cc = np.zeros_like(y)
        sm_nc = np.zeros_like(y)
        np_by_lambda = {float(lam): np.zeros_like(y) for lam in lambdas_tev}
        for i, yy in enumerate(y):
            sm_cc[i], sm_nc[i] = self.dsdy_sm(E, yy, antinu=antinu)
            for lam in np_by_lambda:
                np_by_lambda[lam][i] = self.dsdy_np(E, yy, lam, cut=cut)
        return dict(edges=edges, y=y, dy=dy, sm_cc=sm_cc, sm_nc=sm_nc, np=np_by_lambda)


def probability_from_dsdy(dsdy, dy):
    weights = np.asarray(dsdy, dtype=float) * np.asarray(dy, dtype=float)
    sigma = float(np.sum(weights))
    if sigma <= 0:
        return sigma, np.zeros_like(weights)
    return sigma, weights / sigma


def mean_y(dsdy, y, dy):
    dsdy = np.asarray(dsdy, dtype=float)
    y = np.asarray(y, dtype=float)
    dy = np.asarray(dy, dtype=float)
    norm = float(np.sum(dsdy * dy))
    if norm <= 0:
        return float("nan")
    return float(np.sum(y * dsdy * dy) / norm)
