"""v0.20.3 candidate: mixed-grid combined-generator transport.

Neutrinos remain on the existing DIS table energy grid while taus use a finer
nested log-energy grid.  A single sparse generator per constant-density slab
contains neutrino interactions, CC tau production, tau energy loss, tau decay,
and regenerated daughter neutrinos.  This preserves the exact discrete
semigroup property while reducing tau-advection diffusion without rebuilding
DIS tables.

State order for TauHelicityFlipTransportMixedExpm:
    nu_tau[Nn], antinu_tau[Nn], nu_mu[Nn], antinu_mu[Nn],
    tau_minus[Nt], tau_plus[Nt]
where Nt = 1 + tau_refine*(Nn-1).
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import bmat, csr_matrix, diags, lil_matrix
from scipy.sparse.linalg import expm_multiply

from .constants import NA, M_TAU, C_TAU_CM
from .prem import chord_profile
from .kernels import tau_nutau_pdf, tau_secondary_numu_pdf
from .chirality_transport_expm import _SpeciesOp, _scaled_kernel


def _refined_log_grid(E, factor):
    E = np.asarray(E, float)
    if int(factor) < 1:
        raise ValueError("tau_refine must be >=1")
    l = np.log(E)
    out = [l[0]]
    for a, b in zip(l[:-1], l[1:]):
        out.extend(np.linspace(a, b, int(factor)+1)[1:].tolist())
    return np.exp(np.asarray(out, float))


def _deposit_row(row, log_grid, eout, weight):
    if weight == 0.0:
        return
    u = np.log(float(eout))
    if u < log_grid[0] or u > log_grid[-1]:
        return
    j = int(np.searchsorted(log_grid, u))
    if j == 0:
        row[0] += weight
    elif j >= len(log_grid):
        row[-1] += weight
    else:
        f = (u-log_grid[j-1])/(log_grid[j]-log_grid[j-1])
        row[j-1] += weight*(1.0-f)
        row[j] += weight*f


def _rect_interaction_kernel(E_parent, y, p_y, E_daughter):
    """Conditional parent-grid -> daughter-grid redistribution kernel."""
    Ep = np.asarray(E_parent, float)
    Ed = np.asarray(E_daughter, float)
    ld = np.log(Ed)
    K = np.zeros((len(Ep), len(Ed)), float)
    for j, E in enumerate(Ep):
        row = K[j]
        for yy, w in zip(y, p_y[j]):
            _deposit_row(row, ld, E*(1.0-float(yy)), float(w))
    return K


def _decay_kernel_rect(E_tau, E_nu, polarization, secondary=False,
                       branching=0.18, nz=800):
    z = np.linspace(1e-6, 1.0-1e-6, int(nz))
    if secondary:
        f = tau_secondary_numu_pdf(z, polarization, branching)
        norm = np.trapezoid(f, z)
        if norm > 0:
            f *= branching/norm
    else:
        f = tau_nutau_pdf(z, polarization)
        f /= np.trapezoid(f, z)
    w = f*(z[1]-z[0])
    ln = np.log(np.asarray(E_nu, float))
    K = np.zeros((len(E_tau), len(E_nu)), float)
    for j, E in enumerate(E_tau):
        row = K[j]
        for zz, ww in zip(z, w):
            _deposit_row(row, ln, E*float(zz), float(ww))
    return csr_matrix(K)


class TauHelicityFlipTransportMixedExpm:
    """Mixed-grid combined-generator transport for Case 6/7/resolved models."""

    def __init__(self, table_nu, table_antinu, beta_tau=0.85e-6,
                 secondary_branching=0.18, tau_refine=16):
        self.nu = _SpeciesOp(table_nu)
        self.anu = _SpeciesOp(table_antinu)
        if not np.allclose(self.nu.E, self.anu.E):
            raise ValueError("nu/antinu energy grids differ")
        self.E = self.nu.E.copy()
        self.Nn = len(self.E)
        self.tau_refine = int(tau_refine)
        self.Etau = _refined_log_grid(self.E, self.tau_refine)
        self.Nt = len(self.Etau)
        self.beta = float(beta_tau)
        self.secondary_branching = float(secondary_branching)

        # Rectangular CC production: neutrino grid -> tau grid.
        Kcc_n = _rect_interaction_kernel(self.E, table_nu["y"],
                                          table_nu["p_cc"], self.Etau)
        Kcc_a = _rect_interaction_kernel(self.E, table_antinu["y"],
                                          table_antinu["p_cc"], self.Etau)
        self.Tn = csr_matrix(self.nu.scc[:, None] * Kcc_n)
        self.Ta = csr_matrix(self.anu.scc[:, None] * Kcc_a)

        # Rectangular tau-decay kernels: tau grid -> neutrino grid.
        self.Kdec_minus = _decay_kernel_rect(self.Etau, self.E, -1.0, False)
        self.Kdec_plus  = _decay_kernel_rect(self.Etau, self.E, +1.0, False)
        self.Ksec_minus = _decay_kernel_rect(self.Etau, self.E, -1.0, True,
                                              self.secondary_branching)
        self.Ksec_plus  = _decay_kernel_rect(self.Etau, self.E, +1.0, True,
                                              self.secondary_branching)
        self._cache = {}

    def _interaction_sigma(self, mode, use_np):
        if mode not in ("majorana_flip", "dirac_sink"):
            raise ValueError("mode must be majorana_flip or dirac_sink")
        Nn, Nt = self.Nn, self.Nt
        Znn = csr_matrix((Nn, Nn))
        Znt = csr_matrix((Nn, Nt))
        Ztn = csr_matrix((Nt, Nn))
        Ztt = csr_matrix((Nt, Nt))
        Gn = self.nu.same_block(use_np)
        Ga = self.anu.same_block(use_np)
        Gmn = self.nu.sm_block()
        Gma = self.anu.sm_block()
        NPA = self.anu.np_transfer() if use_np and mode == "majorana_flip" else Znn
        NPN = self.nu.np_transfer() if use_np and mode == "majorana_flip" else Znn
        return bmat([
            [Gn,  NPN, Znn, Znn, self.Tn, Znt],
            [NPA, Ga,  Znn, Znn, Znt, self.Ta],
            [Znn, Znn, Gmn, Znn, Znt, Znt],
            [Znn, Znn, Znn, Gma, Znt, Znt],
            [Ztn, Ztn, Ztn, Ztn, Ztt, Ztt],
            [Ztn, Ztn, Ztn, Ztn, Ztt, Ztt],
        ], format="csr")

    def _tau_loss_per_rho(self):
        l = np.log(self.Etau)
        dl = np.diff(l)
        if np.any(dl <= 0):
            raise ValueError("tau energy grid must be strictly increasing")
        rate = np.empty(self.Nt, float)
        rate[0] = self.beta/dl[0]
        rate[1:] = self.beta/dl
        Q = lil_matrix((self.Nt, self.Nt), dtype=float)
        Q.setdiag(-rate)
        for j in range(1, self.Nt):
            Q[j, j-1] = rate[j]
        return Q.tocsr()

    def _tau_decay_blocks(self, secondary=True):
        d = M_TAU/(C_TAU_CM*self.Etau)
        Q = -diags(d, 0, format="csr")
        Dm = diags(d, 0, format="csr") @ self.Kdec_minus
        Dp = diags(d, 0, format="csr") @ self.Kdec_plus
        if secondary:
            Sm = diags(d, 0, format="csr") @ self.Ksec_minus
            Sp = diags(d, 0, format="csr") @ self.Ksec_plus
        else:
            Sm = csr_matrix((self.Nt, self.Nn))
            Sp = csr_matrix((self.Nt, self.Nn))
        return Q, Dm, Dp, Sm, Sp

    def _generators(self, mode, use_np=True, secondary=True):
        key = (mode, bool(use_np), bool(secondary))
        if key in self._cache:
            return self._cache[key]
        Nn, Nt = self.Nn, self.Nt
        Znn = csr_matrix((Nn, Nn)); Znt = csr_matrix((Nn, Nt))
        Ztn = csr_matrix((Nt, Nn)); Ztt = csr_matrix((Nt, Nt))
        Hsig = self._interaction_sigma(mode, use_np)
        L = self._tau_loss_per_rho()
        Tloss = bmat([
            [Znn,Znn,Znn,Znn,Znt,Znt], [Znn,Znn,Znn,Znn,Znt,Znt],
            [Znn,Znn,Znn,Znn,Znt,Znt], [Znn,Znn,Znn,Znn,Znt,Znt],
            [Ztn,Ztn,Ztn,Ztn,L,Ztt],    [Ztn,Ztn,Ztn,Ztn,Ztt,L],
        ], format="csr")
        G_rho = (NA*Hsig + Tloss).tocsr()

        Q,Dm,Dp,Sm,Sp = self._tau_decay_blocks(secondary)
        Tdec = bmat([
            [Znn,Znn,Znn,Znn,Znt,Znt], [Znn,Znn,Znn,Znn,Znt,Znt],
            [Znn,Znn,Znn,Znn,Znt,Znt], [Znn,Znn,Znn,Znn,Znt,Znt],
            [Dm,Ztn,Ztn,Sm,Q,Ztt],      [Ztn,Dp,Sp,Ztn,Ztt,Q],
        ], format="csr")
        G0 = Tdec.tocsr()
        self._cache[key] = (G_rho, G0)
        return G_rho, G0

    def generator(self, rho, mode="majorana_flip", use_np=True, secondary=True):
        Gr, G0 = self._generators(mode, use_np, secondary)
        return (float(rho)*Gr + G0).tocsr()

    @staticmethod
    def _advance(Grow, x, h_cm):
        if h_cm == 0:
            return np.asarray(x, float)
        return np.asarray(expm_multiply(Grow.T*float(h_cm), np.asarray(x, float)), float)

    def propagate(self, cosz_abs, phi_nu0, phi_antinu0,
                  mode="majorana_flip", use_np=True, max_step_km=12.0,
                  track_secondary_numu=True, profile=None):
        n = np.array(phi_nu0, float, copy=True)
        a = np.array(phi_antinu0, float, copy=True)
        mn = np.zeros(self.Nn); ma = np.zeros(self.Nn)
        tm = np.zeros(self.Nt); tp = np.zeros(self.Nt)
        prof = chord_profile(cosz_abs, max_step_km=max_step_km) if profile is None else profile
        x = np.concatenate([n,a,mn,ma,tm,tp])
        for ds_km, rho in zip(prof["ds_km"], prof["rho"]):
            G = self.generator(rho, mode, use_np, track_secondary_numu)
            x = self._advance(G, x, float(ds_km)*1e5)
        i0=0; i1=self.Nn; n=x[i0:i1]
        i0=i1; i1+=self.Nn; a=x[i0:i1]
        i0=i1; i1+=self.Nn; mn=x[i0:i1]
        i0=i1; i1+=self.Nn; ma=x[i0:i1]
        i0=i1; i1+=self.Nt; tm=x[i0:i1]
        tp=x[i1:i1+self.Nt]
        return dict(nu_tau=n, antinu_tau=a,
                    nu_mu_secondary=mn, antinu_mu_secondary=ma,
                    tau_minus=tm, tau_plus=tp,
                    tau_energies=self.Etau, profile=prof)
