"""Exact constant-density tau energy-loss + decay evolution.

This module fixes the finite-beta_tau operator-splitting weakness identified
in the v0.15.1 audit.  For continuous loss dE/ds = - beta rho E and decay
length L(E)=c tau E/m_tau, both the surviving-tau energy and the decay
probability can be evolved analytically over a constant-density segment.

For a tau with initial energy E0,
  E(s) = E0 exp(-b s),  b = beta rho,
and with a = m_tau/(c tau E0),
  U(s) = integral ds/L(E(s)) = a/b [exp(b s)-1]
for b>0 (U=a s when b=0).  Hence P_surv=exp(-U).

The decay-energy distribution is integrated in the cumulative variable
  t = 1-exp(-u),
where u is decay optical depth.  Since exp(-u) du = dt, Gauss-Legendre
quadrature in t directly integrates the physical decay probability without
requiring metre-scale spatial steps.
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss

from .constants import M_TAU, C_TAU_CM
from .kernels import deposit_logE


class CoupledTauEvolution:
    """Linear tau propagator for one constant-density path segment."""

    def __init__(self, energies, beta_tau=0.85e-6, nquad=24):
        self.E = np.asarray(energies, dtype=float)
        self.beta = float(beta_tau)
        self.nodes, self.weights = leggauss(int(nquad))

    def analytic_parent(self, E0, rho, ds_cm):
        """Return (survival, E_end, decay_optical_depth) for one tau."""
        E0 = float(E0); rho = float(rho); ds_cm = float(ds_cm)
        if ds_cm <= 0:
            return 1.0, E0, 0.0
        a = M_TAU / (C_TAU_CM * E0)
        b = self.beta * rho
        if b == 0.0:
            U = a * ds_cm
            E_end = E0
        else:
            x = b * ds_cm
            U = (a / b) * np.expm1(x)
            E_end = E0 * np.exp(-x)
        # Avoid harmless underflow warnings in extremely opaque segments.
        survival = float(np.exp(-min(float(U), 745.0)))
        return survival, float(E_end), float(U)

    def advance(self, tau, rho, ds_cm):
        """Advance a tau spectrum through one constant-density segment.

        Returns
        -------
        tau_out : ndarray
            Surviving taus after continuous energy loss and decay.
        decay_parent : ndarray
            Distribution of tau energies *at the instant of decay*.  Applying
            the existing tau-decay kernels to this vector produces daughter
            neutrino spectra with the correct pre-decay tau energy.
        """
        tau = np.asarray(tau, dtype=float)
        tau_out = np.zeros_like(tau)
        decay_parent = np.zeros_like(tau)
        if ds_cm <= 0 or not np.any(tau > 0):
            return tau.copy(), decay_parent

        rho = float(rho); ds_cm = float(ds_cm)
        b = self.beta * rho

        for j, (E0, population) in enumerate(zip(self.E, tau)):
            if population <= 0:
                continue
            survival, E_end, U = self.analytic_parent(E0, rho, ds_cm)
            if survival > 0:
                deposit_logE(tau_out, self.E, E_end, population * survival)

            pdec = -np.expm1(-min(U, 745.0))
            if pdec <= 0:
                continue

            if b == 0.0:
                # All decays occur at the same tau energy when there is no
                # continuous loss; this branch is exact and preserves the
                # historical beta_tau=0 limit.
                deposit_logE(decay_parent, self.E, E0, population * pdec)
                continue

            a = M_TAU / (C_TAU_CM * E0)
            # t = 1-exp(-u) is uniform under the physical decay measure:
            # exp(-u) du = dt.  Integrate t from 0 to pdec.
            t = 0.5 * (self.nodes + 1.0) * pdec
            w = 0.5 * pdec * self.weights
            u = -np.log1p(-t)
            Edec = E0 / (1.0 + (b / a) * u)
            for ee, ww in zip(Edec, w):
                if ww > 0:
                    deposit_logE(decay_parent, self.E, float(ee),
                                 population * float(ww))

        return tau_out, decay_parent

    def advance_with_kernels(self, tau, rho, ds_cm, Kdec, Ksecondary=None):
        """Convenience wrapper returning surviving tau and decay daughters."""
        tau_out, decay_parent = self.advance(tau, rho, ds_cm)
        primary = decay_parent @ Kdec
        secondary = None if Ksecondary is None else decay_parent @ Ksecondary
        return tau_out, primary, secondary, decay_parent
