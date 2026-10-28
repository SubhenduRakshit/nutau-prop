"""Sparse-exponential chirality transport.

For each constant-density PREM slice, the neutrino interaction/redistribution
operator is linear.  We therefore evolve it exactly on the discretised energy
grid with scipy.sparse.linalg.expm_multiply instead of enforcing an
'at most one interaction per max_p substep' rule.

Tau loss+decay retains the validated CoupledTauEvolution exact half-step, so
one matter slice is composed as A(h/2) exp(B h) A(h/2).  Numerical convergence
is then controlled by the PREM/matter slice size, not by an interaction
probability substep parameter.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import bmat, csr_matrix, diags
from scipy.sparse.linalg import expm_multiply

from .constants import NA
from .prem import chord_profile
from .kernels import interaction_kernel, tau_decay_kernel, tau_secondary_numu_kernel
from .tau_evolution import CoupledTauEvolution


def _scaled_kernel(sigma, K):
    """Row-vector block diag(sigma) @ K as CSR."""
    return csr_matrix(np.asarray(sigma, float)[:, None] * np.asarray(K, float))


class _SpeciesOp:
    def __init__(self, table):
        self.E = np.asarray(table["energies"], float)
        self.scc = np.asarray(table["sigma_cc"], float)
        self.snc = np.asarray(table["sigma_nc_sm"], float)
        self.snp = np.asarray(table.get("sigma_np", np.zeros_like(self.scc)), float)
        self.Kcc = interaction_kernel(self.E, table["y"], table["p_cc"])
        self.Knc = interaction_kernel(self.E, table["y"], table["p_nc_sm"])
        self.Knp = interaction_kernel(self.E, table["y"], table["p_np"])

    def same_block(self, use_np):
        loss = self.scc + self.snc + (self.snp if use_np else 0.0)
        return -diags(loss, 0, format="csr") + _scaled_kernel(self.snc, self.Knc)

    def sm_block(self):
        return -diags(self.scc + self.snc, 0, format="csr") + _scaled_kernel(self.snc, self.Knc)

    def np_transfer(self):
        return _scaled_kernel(self.snp, self.Knp)

    def cc_to_tau(self):
        return _scaled_kernel(self.scc, self.Kcc)


class _TauMixin:
    def _init_tau(self, E, beta_tau=0.85e-6, secondary_branching=0.18):
        self.E = np.asarray(E, float)
        self.Kdec_minus = tau_decay_kernel(self.E, -1.0)
        self.Kdec_plus = tau_decay_kernel(self.E, +1.0)
        self.Ksec_minus = tau_secondary_numu_kernel(self.E, -1.0, secondary_branching)
        self.Ksec_plus = tau_secondary_numu_kernel(self.E, +1.0, secondary_branching)
        self.tau_evolution = CoupledTauEvolution(self.E, beta_tau=float(beta_tau))

    def _tau_half(self, tm, tp, n, a, mn, ma, rho, h, secondary=True):
        tm, dn, dma, _ = self.tau_evolution.advance_with_kernels(
            tm, rho, h, self.Kdec_minus, self.Ksec_minus if secondary else None)
        tp, da, dmn, _ = self.tau_evolution.advance_with_kernels(
            tp, rho, h, self.Kdec_plus, self.Ksec_plus if secondary else None)
        n = n + dn; a = a + da
        if secondary:
            ma = ma + dma; mn = mn + dmn
        return tm, tp, n, a, mn, ma

    @staticmethod
    def _advance_linear(Hrow, x, rho, h_cm):
        # Hrow acts on row vectors: dx_row/ds = NA*rho*x_row*Hrow.
        # expm_multiply acts on columns, hence Hrow.T.
        alpha = float(NA * rho * h_cm)
        if alpha == 0.0:
            return x
        return np.asarray(expm_multiply(Hrow.T * alpha, np.asarray(x, float)), float)


class TauHelicityFlipTransportExpm(_TauMixin):
    """Exact-discrete interaction evolution for Case 6/7/resolved tau-flavor models."""

    def __init__(self, table_nu, table_antinu, beta_tau=0.85e-6):
        self.nu = _SpeciesOp(table_nu); self.anu = _SpeciesOp(table_antinu)
        if not np.allclose(self.nu.E, self.anu.E):
            raise ValueError("nu/antinu energy grids differ")
        self.N = len(self.nu.E); self._init_tau(self.nu.E, beta_tau)
        self._cache = {}

    def _operator(self, mode, use_np):
        key=(mode,bool(use_np))
        if key in self._cache: return self._cache[key]
        if mode not in ("majorana_flip","dirac_sink"):
            raise ValueError("mode must be majorana_flip or dirac_sink")
        Z=csr_matrix((self.N,self.N))
        Gn=self.nu.same_block(use_np); Ga=self.anu.same_block(use_np)
        Gmn=self.nu.sm_block(); Gma=self.anu.sm_block()
        NPA=self.anu.np_transfer() if use_np and mode=="majorana_flip" else Z
        NPN=self.nu.np_transfer() if use_np and mode=="majorana_flip" else Z
        Tn=self.nu.cc_to_tau(); Ta=self.anu.cc_to_tau()
        # row-state order: n, a, mn, ma, tm, tp
        B=[[Gn, NPN, Z, Z, Tn, Z],
           [NPA, Ga, Z, Z, Z, Ta],
           [Z, Z, Gmn, Z, Z, Z],
           [Z, Z, Z, Gma, Z, Z],
           [Z, Z, Z, Z, Z, Z],
           [Z, Z, Z, Z, Z, Z]]
        H=bmat(B,format="csr")
        self._cache[key]=H
        return H

    def propagate(self, cosz_abs, phi_nu0, phi_antinu0, mode="majorana_flip",
                  use_np=True, max_step_km=12.0, track_secondary_numu=True):
        n=np.array(phi_nu0,float,copy=True); a=np.array(phi_antinu0,float,copy=True)
        mn=np.zeros_like(n); ma=np.zeros_like(n); tm=np.zeros_like(n); tp=np.zeros_like(n)
        H=self._operator(mode,use_np)
        prof=chord_profile(cosz_abs,max_step_km=max_step_km)
        N=self.N
        for ds_km,rho in zip(prof["ds_km"],prof["rho"]):
            h=ds_km*1e5; hh=0.5*h
            tm,tp,n,a,mn,ma=self._tau_half(tm,tp,n,a,mn,ma,rho,hh,track_secondary_numu)
            x=np.concatenate([n,a,mn,ma,tm,tp])
            x=self._advance_linear(H,x,rho,h)
            n,a,mn,ma,tm,tp=[x[i*N:(i+1)*N] for i in range(6)]
            tm,tp,n,a,mn,ma=self._tau_half(tm,tp,n,a,mn,ma,rho,hh,track_secondary_numu)
        return dict(nu_tau=n,antinu_tau=a,nu_mu_secondary=mn,antinu_mu_secondary=ma,
                    tau_minus=tm,tau_plus=tp,profile=prof)


class TauMuonTransitionTransportExpm(_TauMixin):
    """Exact-discrete interaction evolution for Case-5 tau<->mu brackets."""

    def __init__(self, table_nu, table_antinu, beta_tau=0.85e-6):
        self.nu=_SpeciesOp(table_nu); self.anu=_SpeciesOp(table_antinu)
        if not np.allclose(self.nu.E,self.anu.E):
            raise ValueError("nu/antinu energy grids differ")
        self.N=len(self.nu.E); self._init_tau(self.nu.E,beta_tau); self._cache={}

    def _operator(self,mode,use_np):
        key=(mode,bool(use_np))
        if key in self._cache: return self._cache[key]
        if mode not in ("majorana_taumu","dirac_sink"):
            raise ValueError("mode must be majorana_taumu or dirac_sink")
        Z=csr_matrix((self.N,self.N))
        Gn=self.nu.same_block(use_np); Ga=self.anu.same_block(use_np)
        Nn=self.nu.np_transfer() if use_np and mode=="majorana_taumu" else Z
        Na=self.anu.np_transfer() if use_np and mode=="majorana_taumu" else Z
        Tn=self.nu.cc_to_tau(); Ta=self.anu.cc_to_tau()
        # row-state order: tn, ta, mn, ma, tm, tp
        # tn->ma via Nn; ta->mn via Na; mn->ta via Nn; ma->tn via Na.
        B=[[Gn, Z, Z, Nn, Tn, Z],
           [Z, Ga, Na, Z, Z, Ta],
           [Z, Nn, Gn, Z, Z, Z],
           [Na, Z, Z, Ga, Z, Z],
           [Z, Z, Z, Z, Z, Z],
           [Z, Z, Z, Z, Z, Z]]
        H=bmat(B,format="csr"); self._cache[key]=H; return H

    def propagate(self,cosz_abs,tau_nu0,tau_antinu0,mu_nu0,mu_antinu0,
                  mode="majorana_taumu",use_np=True,max_step_km=12.0):
        tn=np.array(tau_nu0,float,copy=True); ta=np.array(tau_antinu0,float,copy=True)
        mn=np.array(mu_nu0,float,copy=True); ma=np.array(mu_antinu0,float,copy=True)
        tm=np.zeros_like(tn); tp=np.zeros_like(tn); H=self._operator(mode,use_np)
        prof=chord_profile(cosz_abs,max_step_km=max_step_km); N=self.N
        for ds_km,rho in zip(prof["ds_km"],prof["rho"]):
            h=ds_km*1e5; hh=0.5*h
            tm,tp,tn,ta,mn,ma=self._tau_half(tm,tp,tn,ta,mn,ma,rho,hh,True)
            x=np.concatenate([tn,ta,mn,ma,tm,tp])
            x=self._advance_linear(H,x,rho,h)
            tn,ta,mn,ma,tm,tp=[x[i*N:(i+1)*N] for i in range(6)]
            tm,tp,tn,ta,mn,ma=self._tau_half(tm,tp,tn,ta,mn,ma,rho,hh,True)
        return dict(nu_tau=tn,antinu_tau=ta,nu_mu=mn,antinu_mu=ma,
                    tau_minus=tm,tau_plus=tp,profile=prof)
