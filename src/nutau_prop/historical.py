"""
Historical regression machinery.

The RRZIterator follows the structure of Rakshit-Reya's Eqs. (16)-(24):
the flux is represented by an effective Z(E,X), the initial Z^(0) already
contains tau-decay regeneration, and later iterations recompute NC + tau-decay
source terms from the previous flux.

It is implemented directly on the discretized transport kernels rather than
retyping the nested analytic y-integrals.  This preserves the same integral
equation while providing an implementation independent of the stepwise
production solver.
"""
import numpy as np
from .constants import NA, M_TAU, C_TAU_CM, R_EARTH_KM
from .prem import chord_profile
from .kernels import interaction_kernel, tau_decay_kernel, tau_secondary_numu_kernel

def historical_flux_per_logE(E,gamma,cutoff_gev=1e8):
    """Return E*dPhi/dE up to normalization."""
    E=np.asarray(E,float)
    if float(gamma)==1.0:
        return 1.0/(1.0+E/float(cutoff_gev))**2
    return (E/1e3)**(1.0-float(gamma))

def nadir_to_mu(theta_deg):
    return float(np.cos(np.deg2rad(theta_deg)))

def interp_logE(E,y,target):
    return float(np.interp(np.log10(target),np.log10(E),y))

class RRZIterator:
    """
    Decay-aware Z iteration for one particle species (nu or antinu).

    Z^(0): tau-decay regeneration generated from a purely attenuated neutrino
           beam, with Z_nutau^(0)=0, matching the physical input prescription
           of RR Eq. (21).

    Z^(n+1): lambda/X * integral dX [S_NC^n + S_tau^n]/F_n^n,
             which is the discretized source form of RR Eqs. (18),(19),(22).

    Tau interactions and radiative loss are omitted in this historical mode,
    as in RR/DRS below 1e8 GeV.
    """
    def __init__(self,table,polarization):
        self.E=np.asarray(table["energies"])
        self.scc=np.asarray(table["sigma_cc"])
        self.snc=np.asarray(table["sigma_nc"])
        self.stot=self.scc+self.snc
        self.lam=1.0/(NA*self.stot)
        self.Kcc=interaction_kernel(self.E,table["y"],table["p_cc"])
        self.Knc=interaction_kernel(self.E,table["y"],table["p_nc"])
        self.Kdec=tau_decay_kernel(self.E,polarization)

    def _geometry(self,mu,depth_steps):
        L=2*R_EARTH_KM*float(mu)
        max_step=max(0.5,L/max(1,int(depth_steps)))
        p=chord_profile(mu,max_step_km=max_step)
        dX=p["rho"]*p["ds_km"]*1e5
        X=np.r_[0.0,np.cumsum(dX)]
        return X,dX,p["rho"]

    def _tau_profile(self,F,X,dX,rho):
        """Solve dT/dX = S_CC - T/lambda_dec on each column-depth segment."""
        T=np.zeros_like(F)
        rcc=NA*self.scc
        for k in range(len(dX)):
            Fmid=0.5*(F[k]+F[k+1])
            Scc=(Fmid*rcc)@self.Kcc
            ldec=rho[k]*(C_TAU_CM*self.E/M_TAU)  # g/cm^2
            a=np.exp(-dX[k]/ldec)
            T[k+1]=T[k]*a + Scc*ldec*(1-a)
        return T

    def _z_from_sources(self,F,T,X,dX,rho,include_nc=True):
        z=np.zeros_like(F)
        integ=np.zeros(len(self.E))
        rnc=NA*self.snc
        tiny=max(float(np.max(F[0]))*1e-280,1e-300)
        for k in range(len(dX)):
            Fmid=0.5*(F[k]+F[k+1])
            if include_nc:
                Snc=(Fmid*rnc)@self.Knc
            else:
                Snc=np.zeros_like(Fmid)
            Tmid=0.5*(T[k]+T[k+1])
            ldec=rho[k]*(C_TAU_CM*self.E/M_TAU)
            Sdec=(Tmid/ldec)@self.Kdec
            ratio=(Snc+Sdec)/np.maximum(Fmid,tiny)
            integ += ratio*dX[k]
            z[k+1]=self.lam*integ/X[k+1]
        return z

    def _flux_from_z(self,phi0,Z,X):
        expo=-(X[:,None]/self.lam[None,:])*(1-Z)
        return phi0[None,:]*np.exp(np.clip(expo,-700,700))

    def solve(self,mu,phi0,depth_steps=900,n_iter=3):
        X,dX,rho=self._geometry(mu,depth_steps)

        # Pure absorption used only to build decay-aware Z^(0).
        Fabs=phi0[None,:]*np.exp(np.clip(-X[:,None]/self.lam[None,:],-700,0))
        T0=self._tau_profile(Fabs,X,dX,rho)
        Z0=self._z_from_sources(Fabs,T0,X,dX,rho,include_nc=False)

        Zs=[Z0]
        Fs=[self._flux_from_z(phi0,Z0,X)]
        Ts=[self._tau_profile(Fs[0],X,dX,rho)]

        for _ in range(int(n_iter)):
            Znew=self._z_from_sources(Fs[-1],Ts[-1],X,dX,rho,include_nc=True)
            Fnew=self._flux_from_z(phi0,Znew,X)
            Tnew=self._tau_profile(Fnew,X,dX,rho)
            Zs.append(Znew); Fs.append(Fnew); Ts.append(Tnew)

        return dict(X=X,Z=Zs,F=Fs,T=Ts)

def relative_max(a,b,mask):
    r=np.divide(a,b,out=np.ones_like(a),where=mask & (b>0))
    return float(np.max(np.abs(r[mask]-1))) if np.any(mask) else np.inf
