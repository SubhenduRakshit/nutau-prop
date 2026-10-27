"""Flavor-aware transport for tau-specific Case-6 NSI.

Unlike the legacy TauTransport, secondary nu_mu produced in tau decay is
propagated with an explicitly supplied SM-only nu_mu table.  This matters for
a tau-specific Case-6 operator: the new NC interaction must not be inherited
by the secondary muon-neutrino channel.
"""
from __future__ import annotations

import numpy as np

from .constants import NA, M_TAU, C_TAU_CM
from .prem import chord_profile
from .kernels import interaction_kernel, tau_decay_kernel, tau_secondary_numu_kernel, deposit_logE


class TauTransportFlavorAware:
    def __init__(self, tau_table, secondary_numu_table, beta_tau=0.85e-6,
                 tau_polarization=-1.0, secondary_branching=0.18):
        self.E=np.asarray(tau_table["energies"])
        self.scc=np.asarray(tau_table["sigma_cc"])
        self.snc=np.asarray(tau_table["sigma_nc"])
        self.Kcc=interaction_kernel(self.E,tau_table["y"],tau_table["p_cc"])
        self.Knc=interaction_kernel(self.E,tau_table["y"],tau_table["p_nc"])
        self.Kdec=tau_decay_kernel(self.E,tau_polarization)
        self.Kdec_numu=tau_secondary_numu_kernel(self.E,tau_polarization,secondary_branching)
        self.beta=float(beta_tau)

        e2=np.asarray(secondary_numu_table["energies"])
        if not np.allclose(self.E,e2,rtol=0,atol=0):
            raise ValueError("tau and secondary-numu energy grids must match")
        self.mucc=np.asarray(secondary_numu_table["sigma_cc"])
        self.munc=np.asarray(secondary_numu_table["sigma_nc"])
        self.muKnc=interaction_kernel(self.E,secondary_numu_table["y"],secondary_numu_table["p_nc"])

    def _tau_energy_loss(self,tau,rho,ds_cm):
        if tau.sum()==0 or self.beta==0:
            return tau
        factor=np.exp(-self.beta*rho*ds_cm)
        out=np.zeros_like(tau)
        for j,E in enumerate(self.E):
            deposit_logE(out,self.E,E*factor,tau[j])
        return out

    @staticmethod
    def _step_flux(nu,rho,h,scc,snc,Knc):
        stot=scc+snc
        pint=1.0-np.exp(-NA*rho*stot*h)
        interacted=nu*pint
        survive=nu-interacted
        fcc=np.divide(scc,stot,out=np.zeros_like(stot),where=stot>0)
        fnc=np.divide(snc,stot,out=np.zeros_like(stot),where=stot>0)
        return survive+(interacted*fnc)@Knc, interacted*fcc

    def propagate(self,cosz_abs,phi_nutau0,phi_tau0=None,max_step_km=10.0,
                  max_p=0.025,track_secondary_numu=True):
        nu=np.array(phi_nutau0,float,copy=True)
        tau=np.zeros_like(nu) if phi_tau0 is None else np.array(phi_tau0,float,copy=True)
        numu=np.zeros_like(nu)
        prof=chord_profile(cosz_abs,max_step_km=max_step_km)
        for ds_km,rho in zip(prof["ds_km"],prof["rho"]):
            ds_cm=ds_km*1e5
            rate_tau=np.max(NA*rho*(self.scc+self.snc))
            rate_mu=np.max(NA*rho*(self.mucc+self.munc)) if track_secondary_numu else 0.0
            nsub=max(1,int(np.ceil(max(rate_tau,rate_mu)*ds_cm/max_p)))
            h=ds_cm/nsub
            for _ in range(nsub):
                nu,cc_parent=self._step_flux(nu,rho,h,self.scc,self.snc,self.Knc)
                tau += cc_parent @ self.Kcc
                if track_secondary_numu:
                    numu,_=self._step_flux(numu,rho,h,self.mucc,self.munc,self.muKnc)
                tau=self._tau_energy_loss(tau,rho,h)
                ldec=C_TAU_CM*self.E/M_TAU
                pdec=1.0-np.exp(-h/ldec)
                dec=tau*pdec
                tau-=dec
                nu += dec @ self.Kdec
                if track_secondary_numu:
                    numu += dec @ self.Kdec_numu
        return dict(nu=nu,tau=tau,numu=numu,profile=prof)
