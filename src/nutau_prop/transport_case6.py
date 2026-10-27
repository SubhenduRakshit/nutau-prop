"""Flavor-aware transport for tau-specific Case-6 NSI.

v0.17 replaces the old full-step tau energy-loss-then-decay ordering by exact
constant-density coupled tau loss+decay evolution inside a Strang-style
splitting. Secondary nu_mu produced in tau decay is propagated with an
explicitly supplied SM-only table and never inherits the tau-specific NSI.
"""
from __future__ import annotations

import numpy as np

from .constants import NA
from .prem import chord_profile
from .kernels import interaction_kernel, tau_decay_kernel, tau_secondary_numu_kernel
from .tau_evolution import CoupledTauEvolution


class TauTransportFlavorAware:
    def __init__(self, tau_table, secondary_numu_table, beta_tau=0.85e-6,
                 tau_polarization=-1.0, secondary_branching=0.18):
        self.E = np.asarray(tau_table["energies"])
        self.scc = np.asarray(tau_table["sigma_cc"])
        self.snc = np.asarray(tau_table["sigma_nc"])
        self.Kcc = interaction_kernel(self.E, tau_table["y"], tau_table["p_cc"])
        self.Knc = interaction_kernel(self.E, tau_table["y"], tau_table["p_nc"])
        self.Kdec = tau_decay_kernel(self.E, tau_polarization)
        self.Kdec_numu = tau_secondary_numu_kernel(
            self.E, tau_polarization, secondary_branching
        )
        self.beta = float(beta_tau)
        self.tau_evolution = CoupledTauEvolution(self.E, beta_tau=self.beta)

        e2 = np.asarray(secondary_numu_table["energies"])
        if not np.allclose(self.E, e2, rtol=0, atol=0):
            raise ValueError("tau and secondary-numu energy grids must match")
        self.mucc = np.asarray(secondary_numu_table["sigma_cc"])
        self.munc = np.asarray(secondary_numu_table["sigma_nc"])
        self.muKnc = interaction_kernel(
            self.E, secondary_numu_table["y"], secondary_numu_table["p_nc"]
        )

    @staticmethod
    def _step_flux(nu, rho, h, scc, snc, Knc):
        stot = scc + snc
        pint = 1.0 - np.exp(-NA * rho * stot * h)
        interacted = nu * pint
        survive = nu - interacted
        fcc = np.divide(scc, stot, out=np.zeros_like(stot), where=stot > 0)
        fnc = np.divide(snc, stot, out=np.zeros_like(stot), where=stot > 0)
        return survive + (interacted * fnc) @ Knc, interacted * fcc

    def _tau_half_step(self, tau, nu, numu, rho, half_h, track_secondary_numu):
        tau, dnutau, dnumu, _ = self.tau_evolution.advance_with_kernels(
            tau, rho, half_h, self.Kdec,
            self.Kdec_numu if track_secondary_numu else None,
        )
        nu = nu + dnutau
        if track_secondary_numu:
            numu = numu + dnumu
        return tau, nu, numu

    def propagate(self, cosz_abs, phi_nutau0, phi_tau0=None, max_step_km=10.0,
                  max_p=0.025, track_secondary_numu=True):
        nu = np.array(phi_nutau0, float, copy=True)
        tau = np.zeros_like(nu) if phi_tau0 is None else np.array(phi_tau0, float, copy=True)
        numu = np.zeros_like(nu)
        prof = chord_profile(cosz_abs, max_step_km=max_step_km)

        for ds_km, rho in zip(prof["ds_km"], prof["rho"]):
            ds_cm = ds_km * 1e5
            rate_tau = np.max(NA * rho * (self.scc + self.snc))
            rate_mu = np.max(NA * rho * (self.mucc + self.munc)) if track_secondary_numu else 0.0
            nsub = max(1, int(np.ceil(max(rate_tau, rate_mu) * ds_cm / max_p)))
            h = ds_cm / nsub
            half_h = 0.5 * h

            for _ in range(nsub):
                # A(h/2): exact coupled tau energy loss + decay.
                tau, nu, numu = self._tau_half_step(
                    tau, nu, numu, rho, half_h, track_secondary_numu
                )

                # B(h): propagate neutrino species. Secondary muons are SM-only.
                nu, cc_parent = self._step_flux(
                    nu, rho, h, self.scc, self.snc, self.Knc
                )
                if track_secondary_numu:
                    numu, _ = self._step_flux(
                        numu, rho, h, self.mucc, self.munc, self.muKnc
                    )
                tau = tau + cc_parent @ self.Kcc

                # A(h/2): evolve newly produced tau as well.
                tau, nu, numu = self._tau_half_step(
                    tau, nu, numu, rho, half_h, track_secondary_numu
                )

        return dict(nu=nu, tau=tau, numu=numu, profile=prof)
