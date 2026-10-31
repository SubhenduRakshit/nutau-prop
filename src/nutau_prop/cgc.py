"""CGC/saturation utilities for UHE neutrino Earth transport.

The primary benchmark is a literature-calibrated rcBK+DLA 2015 suppression
surrogate.  It rescales SM CC and NC total cross sections by an energy-dependent
factor while retaining the conditional SM y distributions.  This is deliberate:
the cited CGC papers publish total/structure-function predictions but not a
machine-readable full d sigma / dy table suitable for exact transport.

References
----------
J. L. Albacete, J. I. Illana, A. Soto-Ontoso,
Phys. Rev. D 92, 014027 (2015), arXiv:1505.06583.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import expm_multiply

from .constants import NA
from .prem import chord_profile
from .chirality_transport_expm import _SpeciesOp


# Anchors encode only published qualitative/quantitative information used in
# the surrogate: BK ~ DGLAP through ~1e8 GeV; at 1e9 GeV the DLA-BK prediction
# is about 1/1.5 of DGLAP; by ~1e12 GeV the suppression is at least ~1/2.5;
# at 1e14 GeV the difference reaches ~1/4.5.  Interpolation is in log E-log R.
_RCBK_DLA_E = np.asarray([1.0e4, 1.0e8, 1.0e9, 1.0e12, 1.0e14], float)
_RCBK_DLA_R = np.asarray([1.0,    1.0,    1.0/1.5, 1.0/2.5, 1.0/4.5], float)


def rcbk_dla_2015_ratio(energy_gev):
    """Return sigma_CGC/sigma_DGLAP for the v0.88 rcBK+DLA surrogate.

    Values outside the anchor range are held fixed at the endpoint value.
    The object returned has the same shape as ``energy_gev``.
    """
    e = np.asarray(energy_gev, float)
    if np.any(e <= 0):
        raise ValueError("energies must be positive")
    lr = np.interp(
        np.log(e), np.log(_RCBK_DLA_E), np.log(_RCBK_DLA_R),
        left=np.log(_RCBK_DLA_R[0]), right=np.log(_RCBK_DLA_R[-1]),
    )
    return np.exp(lr)


def rescale_sm_table(table, ratio):
    """Rescale CC and NC totals, retaining the conditional SM y kernels."""
    out = {k: np.array(v, copy=True) if isinstance(v, np.ndarray) else v
           for k, v in table.items()}
    r = np.asarray(ratio, float)
    e = np.asarray(out["energies"], float)
    if r.ndim == 0:
        r = np.full_like(e, float(r))
    if r.shape != e.shape:
        raise ValueError("ratio shape must match table energies")
    if np.any((r <= 0) | (r > 1.5)):
        raise ValueError("unphysical CGC rescaling ratio")
    out["sigma_cc"] = np.asarray(out["sigma_cc"], float) * r
    out["sigma_nc_sm"] = np.asarray(out["sigma_nc_sm"], float) * r
    out["cgc_ratio"] = r
    return out


class LightFlavorTransport:
    """PREM CC absorption + NC energy migration for e/mu neutrinos.

    Oscillations and charged-lepton regeneration are neglected, appropriate for
    the UHE regime treated here.  Neutrino and antineutrino are evolved
    independently with the same validated discrete interaction operator used by
    the tau transport package.
    """

    def __init__(self, table_nu, table_antinu):
        self.nu = _SpeciesOp(table_nu)
        self.anu = _SpeciesOp(table_antinu)
        if not np.allclose(self.nu.E, self.anu.E):
            raise ValueError("nu/antinu energy grids differ")
        self.E = self.nu.E.copy()
        self.Hn = self.nu.sm_block().tocsr()
        self.Ha = self.anu.sm_block().tocsr()

    @staticmethod
    def _advance(Hrow, x, rho, h_cm):
        alpha = float(NA * rho * h_cm)
        if alpha == 0.0:
            return np.asarray(x, float)
        return np.asarray(expm_multiply(Hrow.T * alpha, np.asarray(x, float)), float)

    def propagate(self, cosz_abs, phi_nu0, phi_antinu0, max_step_km=48.0,
                  profile=None):
        n = np.asarray(phi_nu0, float).copy()
        a = np.asarray(phi_antinu0, float).copy()
        prof = chord_profile(float(cosz_abs), max_step_km=float(max_step_km)) \
            if profile is None else profile
        for ds_km, rho in zip(prof["ds_km"], prof["rho"]):
            h = float(ds_km) * 1.0e5
            n = self._advance(self.Hn, n, rho, h)
            a = self._advance(self.Ha, a, rho, h)
        return dict(nu=n, antinu=a, profile=prof)
