"""Chirality-aware transport screens for dimension-7 momentum-growing NSI.

This module is deliberately a physics bracketing tool.  The low-energy
dimension-7 currents are chirality flipping for Dirac neutrinos.  For Majorana
neutrinos the corresponding right-helicity state is active (antineutrino-like),
and the tensor Case-5 operator is necessarily flavor off diagonal.

Implemented brackets:
  * Dirac sink: NP scattering removes an active neutrino into a RH state.
  * Majorana diagonal flip: nu_tau <-> antinu_tau with NP energy degradation.
  * Majorana tau-mu transition: nu_tau <-> antinu_mu and
    antinu_tau <-> nu_mu, including the reverse transitions required by the
    Hermitian flavor-off-diagonal coupling.

v0.15.1 numerical-audit corrections:
  * SM and BSM comparison runs use the same substep schedule.  The schedule is
    set by the full SM+NP interaction rate even when use_np=False.  This removes
    the spurious possibility that a pure Dirac sink appears above the matched
    SM solution solely because the two runs used different numerical
    resolutions.
  * secondary nu_mu/antinu_mu produced by tau decay in the helicity-flip solver
    are subsequently propagated with SM CC/NC interactions instead of merely
    being accumulated.

SM CC/NC transport and tau decay/regeneration use the already validated
PREM/DRS machinery conventions of the project.
"""
from __future__ import annotations

import numpy as np

from .constants import NA, M_TAU, C_TAU_CM
from .prem import chord_profile
from .kernels import (
    interaction_kernel,
    tau_decay_kernel,
    tau_secondary_numu_kernel,
    deposit_logE,
)


class _Species:
    def __init__(self, table):
        self.E = np.asarray(table["energies"], float)
        self.scc = np.asarray(table["sigma_cc"], float)
        self.snc = np.asarray(table["sigma_nc_sm"], float)
        self.snp = np.asarray(table.get("sigma_np", np.zeros_like(self.scc)), float)
        self.Kcc = interaction_kernel(self.E, table["y"], table["p_cc"])
        self.Knc = interaction_kernel(self.E, table["y"], table["p_nc_sm"])
        self.Knp = interaction_kernel(self.E, table["y"], table["p_np"])

    def step(self, phi, rho, h, use_np=True):
        snp = self.snp if use_np else np.zeros_like(self.snp)
        stot = self.scc + self.snc + snp
        pint = 1.0 - np.exp(-NA * rho * stot * h)
        inter = phi * pint
        survive = phi - inter
        fcc = np.divide(self.scc, stot, out=np.zeros_like(stot), where=stot > 0)
        fnc = np.divide(self.snc, stot, out=np.zeros_like(stot), where=stot > 0)
        fnp = np.divide(snp, stot, out=np.zeros_like(stot), where=stot > 0)
        same = survive + (inter * fnc) @ self.Knc
        cc_parent = inter * fcc
        np_daughter = (inter * fnp) @ self.Knp if use_np else np.zeros_like(phi)
        return same, cc_parent, np_daughter

    def matched_step_rate(self, rho):
        """Rate used to define a common SM/BSM numerical substep schedule."""
        return float(np.max(NA * rho * (self.scc + self.snc + self.snp)))


class _TauDecayMixin:
    def _init_tau_decay(self, energies, beta_tau=0.85e-6, secondary_branching=0.18):
        self.E = np.asarray(energies, float)
        self.beta = float(beta_tau)
        self.Kdec_minus = tau_decay_kernel(self.E, -1.0)
        self.Kdec_plus = tau_decay_kernel(self.E, +1.0)
        self.Ksec_minus = tau_secondary_numu_kernel(self.E, -1.0, secondary_branching)
        self.Ksec_plus = tau_secondary_numu_kernel(self.E, +1.0, secondary_branching)

    def _tau_loss(self, tau, rho, ds_cm):
        if tau.sum() == 0 or self.beta == 0:
            return tau
        factor = np.exp(-self.beta * rho * ds_cm)
        out = np.zeros_like(tau)
        for j, E in enumerate(self.E):
            deposit_logE(out, self.E, E * factor, tau[j])
        return out

    def _decay(self, tau_minus, tau_plus, h):
        ldec = C_TAU_CM * self.E / M_TAU
        pdec = 1.0 - np.exp(-h / ldec)
        dm = tau_minus * pdec
        dp = tau_plus * pdec
        tau_minus = tau_minus - dm
        tau_plus = tau_plus - dp
        return tau_minus, tau_plus, dm, dp


class TauHelicityFlipTransport(_TauDecayMixin):
    """Tau-flavor operator: Dirac sink or Majorana nu<->antinu flip."""

    def __init__(self, table_nu, table_antinu, beta_tau=0.85e-6):
        self.nu = _Species(table_nu)
        self.anu = _Species(table_antinu)
        if not np.allclose(self.nu.E, self.anu.E):
            raise ValueError("nu/antinu energy grids differ")
        self._init_tau_decay(self.nu.E, beta_tau)

    def propagate(self, cosz_abs, phi_nu0, phi_antinu0, mode="majorana_flip",
                  use_np=True, max_step_km=12.0, max_p=0.02,
                  track_secondary_numu=True):
        if mode not in ("majorana_flip", "dirac_sink"):
            raise ValueError("mode must be majorana_flip or dirac_sink")
        n = np.array(phi_nu0, float, copy=True)
        a = np.array(phi_antinu0, float, copy=True)
        tm = np.zeros_like(n); tp = np.zeros_like(n)
        mn = np.zeros_like(n); ma = np.zeros_like(n)
        prof = chord_profile(cosz_abs, max_step_km=max_step_km)
        for ds_km, rho in zip(prof["ds_km"], prof["rho"]):
            ds_cm = ds_km * 1e5
            # IMPORTANT: use the full SM+NP rate for the substep schedule even
            # when use_np=False, so SM and BSM runs are compared at identical
            # numerical resolution.
            maxrate = max(self.nu.matched_step_rate(rho), self.anu.matched_step_rate(rho))
            nsub = max(1, int(np.ceil(maxrate * ds_cm / max_p)))
            h = ds_cm / nsub
            for _ in range(nsub):
                n0, ccn, npn = self.nu.step(n, rho, h, use_np=use_np)
                a0, cca, npa = self.anu.step(a, rho, h, use_np=use_np)

                # Existing secondary muon-neutrino flux propagates through SM
                # CC/NC interactions.  The NP operator in this class is tau
                # flavor only, so secondary muon neutrinos use use_np=False.
                if track_secondary_numu:
                    mn, _ccmn, _ = self.nu.step(mn, rho, h, use_np=False)
                    ma, _ccma, _ = self.anu.step(ma, rho, h, use_np=False)

                if use_np and mode == "majorana_flip":
                    n = n0 + npa
                    a = a0 + npn
                else:
                    n, a = n0, a0

                tm += ccn @ self.nu.Kcc
                tp += cca @ self.anu.Kcc
                tm = self._tau_loss(tm, rho, h)
                tp = self._tau_loss(tp, rho, h)
                tm, tp, dm, dp = self._decay(tm, tp, h)
                n += dm @ self.Kdec_minus
                a += dp @ self.Kdec_plus
                if track_secondary_numu:
                    ma += dm @ self.Ksec_minus
                    mn += dp @ self.Ksec_plus
        return dict(nu_tau=n, antinu_tau=a, nu_mu_secondary=mn,
                    antinu_mu_secondary=ma, tau_minus=tm, tau_plus=tp, profile=prof)


class TauMuonTransitionTransport(_TauDecayMixin):
    """Case-5 tau-mu transition bracket.

    Majorana mode includes both directions of the off-diagonal transition:
      nu_tau -> antinu_mu, antinu_tau -> nu_mu,
      nu_mu  -> antinu_tau, antinu_mu  -> nu_tau.
    Dirac-sink mode removes the NP-interacted active state instead.
    """

    def __init__(self, table_nu, table_antinu, beta_tau=0.85e-6):
        self.nu = _Species(table_nu)
        self.anu = _Species(table_antinu)
        if not np.allclose(self.nu.E, self.anu.E):
            raise ValueError("nu/antinu energy grids differ")
        self._init_tau_decay(self.nu.E, beta_tau)

    def propagate(self, cosz_abs, tau_nu0, tau_antinu0, mu_nu0, mu_antinu0,
                  mode="majorana_taumu", use_np=True,
                  max_step_km=12.0, max_p=0.02):
        if mode not in ("majorana_taumu", "dirac_sink"):
            raise ValueError("mode must be majorana_taumu or dirac_sink")
        tn = np.array(tau_nu0, float, copy=True)
        ta = np.array(tau_antinu0, float, copy=True)
        mn = np.array(mu_nu0, float, copy=True)
        ma = np.array(mu_antinu0, float, copy=True)
        tm = np.zeros_like(tn); tp = np.zeros_like(tn)
        prof = chord_profile(cosz_abs, max_step_km=max_step_km)
        for ds_km, rho in zip(prof["ds_km"], prof["rho"]):
            ds_cm = ds_km * 1e5
            # Matched SM/BSM numerical resolution, as above.
            maxrate = max(self.nu.matched_step_rate(rho), self.anu.matched_step_rate(rho))
            nsub = max(1, int(np.ceil(maxrate * ds_cm / max_p)))
            h = ds_cm / nsub
            for _ in range(nsub):
                tn0, cctn, nptn = self.nu.step(tn, rho, h, use_np=use_np)
                ta0, ccta, npta = self.anu.step(ta, rho, h, use_np=use_np)
                mn0, _ccmn, npmn = self.nu.step(mn, rho, h, use_np=use_np)
                ma0, _ccma, npma = self.anu.step(ma, rho, h, use_np=use_np)
                if use_np and mode == "majorana_taumu":
                    tn = tn0 + npma
                    ta = ta0 + npmn
                    mn = mn0 + npta
                    ma = ma0 + nptn
                else:
                    tn, ta, mn, ma = tn0, ta0, mn0, ma0
                tm += cctn @ self.nu.Kcc
                tp += ccta @ self.anu.Kcc
                tm = self._tau_loss(tm, rho, h)
                tp = self._tau_loss(tp, rho, h)
                tm, tp, dm, dp = self._decay(tm, tp, h)
                tn += dm @ self.Kdec_minus
                ta += dp @ self.Kdec_plus
                # tau- -> anti-nu_mu; tau+ -> nu_mu.  These daughters enter the
                # propagated muon-neutrino pools and are evolved on later steps.
                ma += dm @ self.Ksec_minus
                mn += dp @ self.Ksec_plus
        return dict(nu_tau=tn, antinu_tau=ta, nu_mu=mn, antinu_mu=ma,
                    tau_minus=tm, tau_plus=tp, profile=prof)
