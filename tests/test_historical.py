import numpy as np
from nutau_prop.kernels import tau_nutau_pdf,tau_secondary_numu_pdf
from nutau_prop.historical import historical_flux_per_logE,nadir_to_mu

def test_inclusive_nutau_decay_normalization():
    z=np.linspace(1e-6,1-1e-6,200000)
    for P in (-1.0,1.0):
        val=np.trapezoid(tau_nutau_pdf(z,P),z)
        assert abs(val-1.0)<4e-3

def test_secondary_numu_branching():
    z=np.linspace(1e-6,1-1e-6,30000)
    for P in (-1.0,1.0):
        assert abs(np.trapezoid(tau_secondary_numu_pdf(z,P,0.18),z)-0.18)<3e-4

def test_historical_fluxes():
    E=np.array([1e3,1e4,1e5])
    f1=historical_flux_per_logE(E,1,1e8)
    f2=historical_flux_per_logE(E,2,1e8)
    assert f1[0]>f1[-1]
    assert np.allclose(f2,[1,.1,.01])

def test_nadir_mapping():
    assert abs(nadir_to_mu(0)-1)<1e-12
    assert abs(nadir_to_mu(60)-.5)<1e-12
