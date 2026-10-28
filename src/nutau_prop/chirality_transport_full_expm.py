"""v0.20 candidate: one sparse generator for neutrinos + taus.

This removes the A(h/2) exp(Bh) A(h/2) spatial splitting used in v0.19.
On the discrete energy grid, each constant-density PREM slab is evolved as

    x(s+h) = x(s) exp[G(rho) h]

with one linear generator containing
  * SM/NP neutrino interaction and redistribution,
  * CC tau production,
  * continuous tau energy loss discretised as conservative down-logE drift,
  * tau decay,
  * regenerated nu_tau and secondary nu_mu daughters.

The tau drift generator is a first-order finite-volume/upwind discretisation in
log E.  It is a true semigroup on the chosen energy grid, unlike repeatedly
projecting an exact continuous shift back to the grid after every spatial step.
Energy-grid convergence must therefore be validated separately before use in
production.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import bmat, csr_matrix, diags, lil_matrix
from scipy.sparse.linalg import expm_multiply

from .constants import NA, M_TAU, C_TAU_CM
from .prem import chord_profile
from .kernels import tau_decay_kernel, tau_secondary_numu_kernel
from .chirality_transport_expm import _SpeciesOp, _scaled_kernel


class TauHelicityFlipTransportFullExpm:
    """Combined-generator candidate for Case 6/7/resolved tau-flavor models."""

    def __init__(self, table_nu, table_antinu, beta_tau=0.85e-6,
                 secondary_branching=0.18):
        self.nu = _SpeciesOp(table_nu)
        self.anu = _SpeciesOp(table_antinu)
        if not np.allclose(self.nu.E, self.anu.E):
            raise ValueError("nu/antinu energy grids differ")
        self.E = self.nu.E.copy()
        self.N = len(self.E)
        self.beta = float(beta_tau)
        self.Kdec_minus = tau_decay_kernel(self.E, -1.0)
        self.Kdec_plus = tau_decay_kernel(self.E, +1.0)
        self.Ksec_minus = tau_secondary_numu_kernel(
            self.E, -1.0, secondary_branching)
        self.Ksec_plus = tau_secondary_numu_kernel(
            self.E, +1.0, secondary_branching)
        self._cache = {}

    def _interaction_sigma(self, mode, use_np):
        if mode not in ("majorana_flip", "dirac_sink"):
            raise ValueError("mode must be majorana_flip or dirac_sink")
        N = self.N
        Z = csr_matrix((N, N))
        Gn = self.nu.same_block(use_np)
        Ga = self.anu.same_block(use_np)
        Gmn = self.nu.sm_block()
        Gma = self.anu.sm_block()
        NPA = self.anu.np_transfer() if use_np and mode == "majorana_flip" else Z
        NPN = self.nu.np_transfer() if use_np and mode == "majorana_flip" else Z
        Tn = self.nu.cc_to_tau()
        Ta = self.anu.cc_to_tau()
        # row-state order: n, a, mn, ma, tm, tp
        return bmat([
            [Gn,  NPN, Z,   Z,   Tn, Z],
            [NPA, Ga,  Z,   Z,   Z,  Ta],
            [Z,   Z,   Gmn, Z,   Z,  Z],
            [Z,   Z,   Z,   Gma, Z,  Z],
            [Z,   Z,   Z,   Z,   Z,  Z],
            [Z,   Z,   Z,   Z,   Z,  Z],
        ], format="csr")

    def _tau_loss_per_rho(self):
        """Tau log-energy drift generator divided by density rho.

        For d ln(E)/ds = -beta*rho, population in log-E cell j jumps to
        j-1 at rate beta*rho/Delta lnE.  The lowest cell leaks below the
        represented energy range, matching the historical out-of-support loss.
        """
        lE = np.log(self.E)
        dl = np.diff(lE)
        if np.any(dl <= 0):
            raise ValueError("energy grid must be strictly increasing")
        rate = np.empty(self.N, float)
        rate[0] = self.beta / dl[0]
        rate[1:] = self.beta / dl
        Q = lil_matrix((self.N, self.N), dtype=float)
        Q.setdiag(-rate)
        for j in range(1, self.N):
            Q[j, j-1] = rate[j]
        return Q.tocsr()

    def _tau_decay_blocks(self, secondary=True):
        d = M_TAU / (C_TAU_CM * self.E)  # decay rate per cm
        Q = -diags(d, 0, format="csr")
        Dm = _scaled_kernel(d, self.Kdec_minus)
        Dp = _scaled_kernel(d, self.Kdec_plus)
        Sm = _scaled_kernel(d, self.Ksec_minus) if secondary else csr_matrix((self.N,self.N))
        Sp = _scaled_kernel(d, self.Ksec_plus) if secondary else csr_matrix((self.N,self.N))
        return Q, Dm, Dp, Sm, Sp

    def _generators(self, mode, use_np=True, secondary=True):
        """Return (G_rho, G_0) such that G(rho)=rho*G_rho+G_0 [cm^-1]."""
        key = (mode, bool(use_np), bool(secondary))
        if key in self._cache:
            return self._cache[key]
        N = self.N
        Z = csr_matrix((N, N))
        Hsig = self._interaction_sigma(mode, use_np)

        L = self._tau_loss_per_rho()
        Tloss = bmat([
            [Z,Z,Z,Z,Z,Z], [Z,Z,Z,Z,Z,Z],
            [Z,Z,Z,Z,Z,Z], [Z,Z,Z,Z,Z,Z],
            [Z,Z,Z,Z,L,Z], [Z,Z,Z,Z,Z,L],
        ], format="csr")
        G_rho = (NA * Hsig + Tloss).tocsr()

        Q, Dm, Dp, Sm, Sp = self._tau_decay_blocks(secondary)
        Tdec = bmat([
            [Z,Z,Z,Z,Z,Z], [Z,Z,Z,Z,Z,Z],
            [Z,Z,Z,Z,Z,Z], [Z,Z,Z,Z,Z,Z],
            [Dm,Z,Z,Sm,Q,Z],
            [Z,Dp,Sp,Z,Z,Q],
        ], format="csr")
        G_0 = Tdec.tocsr()
        self._cache[key] = (G_rho, G_0)
        return G_rho, G_0

    def generator(self, rho, mode="majorana_flip", use_np=True, secondary=True):
        Gr, G0 = self._generators(mode, use_np, secondary)
        return (float(rho) * Gr + G0).tocsr()

    @staticmethod
    def _advance(Grow, x, h_cm):
        if h_cm == 0:
            return np.asarray(x, float)
        return np.asarray(expm_multiply(Grow.T * float(h_cm), np.asarray(x, float)), float)

    def propagate(self, cosz_abs, phi_nu0, phi_antinu0,
                  mode="majorana_flip", use_np=True, max_step_km=12.0,
                  track_secondary_numu=True, profile=None):
        n = np.array(phi_nu0, float, copy=True)
        a = np.array(phi_antinu0, float, copy=True)
        mn = np.zeros_like(n); ma = np.zeros_like(n)
        tm = np.zeros_like(n); tp = np.zeros_like(n)
        prof = chord_profile(cosz_abs, max_step_km=max_step_km) if profile is None else profile
        N = self.N
        x = np.concatenate([n,a,mn,ma,tm,tp])
        for ds_km, rho in zip(prof["ds_km"], prof["rho"]):
            G = self.generator(rho, mode, use_np, track_secondary_numu)
            x = self._advance(G, x, float(ds_km)*1e5)
        n,a,mn,ma,tm,tp = [x[i*N:(i+1)*N] for i in range(6)]
        return dict(nu_tau=n, antinu_tau=a,
                    nu_mu_secondary=mn, antinu_mu_secondary=ma,
                    tau_minus=tm, tau_plus=tp, profile=prof)
