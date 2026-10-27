import numpy as np

from .constants import NA
from .prem import chord_profile
from .kernels import interaction_kernel, tau_decay_kernel, tau_secondary_numu_kernel
from .tau_evolution import CoupledTauEvolution


class TauTransport:
    """nu_tau/tau transport with exact finite-beta tau evolution.

    v0.17 uses Strang-style splitting inside each matter substep:
      tau evolution for h/2 -> neutrino interaction for h -> tau evolution for h/2.
    The tau half-steps solve continuous energy loss and decay analytically in
    constant density, so tau decay no longer has to be resolved by metre-scale
    geometric stepping.
    """

    def __init__(self, table, beta_tau=0.85e-6, tau_polarization=-1.0,
                 secondary_branching=0.18):
        self.E = np.asarray(table["energies"])
        self.scc = np.asarray(table["sigma_cc"])
        self.snc = np.asarray(table["sigma_nc"])
        self.Kcc = interaction_kernel(self.E, table["y"], table["p_cc"])
        self.Knc = interaction_kernel(self.E, table["y"], table["p_nc"])
        self.Kdec = tau_decay_kernel(self.E, tau_polarization)
        self.Kdec_numu = tau_secondary_numu_kernel(
            self.E, tau_polarization, secondary_branching
        )
        self.beta = float(beta_tau)
        self.tau_evolution = CoupledTauEvolution(self.E, beta_tau=self.beta)

    def _nu_step(self, nu, rho, h):
        stot = self.scc + self.snc
        pint = 1.0 - np.exp(-NA * rho * stot * h)
        interacted = nu * pint
        survive = nu - interacted
        fcc = np.divide(self.scc, stot, out=np.zeros_like(stot), where=stot > 0)
        fnc = np.divide(self.snc, stot, out=np.zeros_like(stot), where=stot > 0)
        return survive + (interacted * fnc) @ self.Knc, interacted * fcc

    def _tau_half_step(self, tau, nu, numu, rho, half_h, track_secondary_numu):
        tau, dnutau, dnumu, _ = self.tau_evolution.advance_with_kernels(
            tau, rho, half_h, self.Kdec,
            self.Kdec_numu if track_secondary_numu else None,
        )
        nu = nu + dnutau
        if track_secondary_numu:
            numu = numu + dnumu
        return tau, nu, numu

    def propagate(self, cosz_abs, phi_nu0, phi_tau0=None,
                  max_step_km=10.0, max_p=0.025,
                  phi_numu0=None, track_secondary_numu=False):
        nu = np.array(phi_nu0, float, copy=True)
        tau = np.zeros_like(nu) if phi_tau0 is None else np.array(phi_tau0, float, copy=True)
        numu = np.zeros_like(nu) if phi_numu0 is None else np.array(phi_numu0, float, copy=True)
        prof = chord_profile(cosz_abs, max_step_km=max_step_km)

        for ds_km, rho in zip(prof["ds_km"], prof["rho"]):
            ds_cm = ds_km * 1e5
            rate_max = float(np.max(NA * rho * (self.scc + self.snc)))
            nsub = max(1, int(np.ceil(rate_max * ds_cm / max_p)))
            h = ds_cm / nsub
            half_h = 0.5 * h

            for _ in range(nsub):
                # A(h/2): exact coupled tau loss+decay.
                tau, nu, numu = self._tau_half_step(
                    tau, nu, numu, rho, half_h, track_secondary_numu
                )

                # B(h): neutrino interaction/redistribution.
                nu, cc_parent = self._nu_step(nu, rho, h)
                if phi_numu0 is not None or track_secondary_numu:
                    numu, _ = self._nu_step(numu, rho, h)
                tau = tau + cc_parent @ self.Kcc

                # A(h/2): evolve both pre-existing and newly produced taus.
                tau, nu, numu = self._tau_half_step(
                    tau, nu, numu, rho, half_h, track_secondary_numu
                )

        return dict(nu=nu, tau=tau, numu=numu, profile=prof)


class NuOnlyTransport:
    """CC attenuation + NC regeneration for nu_mu/e."""

    def __init__(self, table):
        self.E = np.asarray(table["energies"])
        self.scc = np.asarray(table["sigma_cc"])
        self.snc = np.asarray(table["sigma_nc"])
        self.Knc = interaction_kernel(self.E, table["y"], table["p_nc"])

    def propagate(self, cosz_abs, phi0, max_step_km=10.0, max_p=0.025):
        nu = np.array(phi0, float, copy=True)
        prof = chord_profile(cosz_abs, max_step_km=max_step_km)
        for ds_km, rho in zip(prof["ds_km"], prof["rho"]):
            ds_cm = ds_km * 1e5
            stot = self.scc + self.snc
            rate_max = float(np.max(NA * rho * stot))
            nsub = max(1, int(np.ceil(rate_max * ds_cm / max_p)))
            h = ds_cm / nsub
            for _ in range(nsub):
                pint = 1.0 - np.exp(-NA * rho * stot * h)
                inter = nu * pint
                survive = nu - inter
                fnc = np.divide(self.snc, stot, out=np.zeros_like(stot), where=stot > 0)
                nu = survive + (inter * fnc) @ self.Knc
        return dict(nu=nu, profile=prof)
