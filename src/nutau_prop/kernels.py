import numpy as np

MTAU = 1.77686
MPI = 0.13957
MRHO = 0.77526
MA1 = 1.230

def deposit_logE(row, energies, eout, weight):
    if weight == 0 or eout < energies[0] or eout > energies[-1]:
        return
    lE=np.log(energies); u=np.log(eout); j=np.searchsorted(lE,u)
    if j==0:
        row[0]+=weight
    elif j>=len(energies):
        row[-1]+=weight
    else:
        f=(u-lE[j-1])/(lE[j]-lE[j-1])
        row[j-1]+=weight*(1-f); row[j]+=weight*f

def interaction_kernel(energies, y, p_y):
    """Conditional redistribution kernel, parent log-E bin -> daughter log-E bin."""
    K=np.zeros((len(energies),len(energies)))
    for j,E in enumerate(energies):
        for yy,w in zip(y,p_y[j]):
            deposit_logE(K[j],energies,E*(1-yy),float(w))
    return K

def _two_body_hadronic(z, mass, polarization, vector=False):
    r=(mass/MTAU)**2
    support=(z <= (1-r))
    g0=np.where(support,1.0/(1-r),0.0)
    if vector:
        analyzing=(1-2*r)/(1+2*r)
        # Table-I formula is:
        # g1 = -[(2z-1+r)/(1-r)] [(1-2r)/(1+2r)] theta(...)
        g1=np.where(support,
                    -((2*z-1+r)/(1-r))*analyzing,0.0)
    else:
        # pion:
        # g1 = -(2z-1+r)/(1-r)^2 theta(...)
        g1=np.where(support,-(2*z-1+r)/(1-r)**2,0.0)
    return g0 + polarization*g1

def tau_nutau_pdf(z, polarization=-1.0):
    r"""
    Inclusive dn/dz for tau -> nu_tau + X from Dutta, Reno & Sarcevic,
    PRD 62, 123001 (2000), Appendix A, Table I.

    Branching fractions used in that historical calculation:
      e nu nu  0.18
      mu nu nu 0.18
      pi        0.12
      rho       0.26
      a1        0.13
      other X   0.13
    They sum to one.
    """
    z=np.asarray(z,dtype=float)
    P=float(polarization)

    # Each leptonic mode has the same nu_tau distribution.
    g0_l=5/3 - 3*z*z + (4/3)*z**3
    g1_l=1/3 - 3*z*z + (8/3)*z**3
    lept=0.36*(g0_l + P*g1_l)

    pion=0.12*_two_body_hadronic(z,MPI,P,vector=False)
    rho =0.26*_two_body_hadronic(z,MRHO,P,vector=True)
    a1  =0.13*_two_body_hadronic(z,MA1,P,vector=True)

    # Multiprong remainder, Table I.
    other=0.13*np.where(z<=0.3,1/0.3,0.0)

    return np.clip(lept+pion+rho+a1+other,0.0,None)

def tau_decay_pdf(z, polarization=-1.0):
    """Backward-compatible name; now returns the full inclusive tau->nu_tau spectrum."""
    return tau_nutau_pdf(z, polarization)

def tau_decay_kernel(energies, polarization=-1.0, nz=800):
    z=np.linspace(1e-6,1-1e-6,nz)
    f=tau_nutau_pdf(z,polarization)
    f/=np.trapezoid(f,z)
    dz=z[1]-z[0]
    K=np.zeros((len(energies),len(energies)))
    for j,E in enumerate(energies):
        for zz,ff in zip(z,f):
            deposit_logE(K[j],energies,E*zz,float(ff*dz))
    return K

def tau_secondary_numu_pdf(z, polarization=-1.0, branching=0.18):
    r"""
    Rakshit-Reya Eq. (31), identical to the DRS leptonic secondary spectrum:
      dn/dz = B [2-6z^2+4z^3 + P(-2+12z-18z^2+8z^3)].
    tau- -> anti-nu_mu: P=-1
    tau+ -> nu_mu:      P=+1
    """
    z=np.asarray(z)
    base=2 - 6*z*z + 4*z**3
    pol=-2 + 12*z - 18*z*z + 8*z**3
    return branching*np.clip(base + polarization*pol,0,None)

def tau_secondary_numu_kernel(energies, polarization=-1.0, branching=0.18, nz=800):
    z=np.linspace(1e-6,1-1e-6,nz)
    f=tau_secondary_numu_pdf(z,polarization,branching)
    norm=np.trapezoid(f,z)
    if norm>0:
        f*=branching/norm
    dz=z[1]-z[0]
    K=np.zeros((len(energies),len(energies)))
    for j,E in enumerate(energies):
        for zz,ff in zip(z,f):
            deposit_logE(K[j],energies,E*zz,float(ff*dz))
    return K
