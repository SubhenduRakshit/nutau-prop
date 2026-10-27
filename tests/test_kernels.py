import numpy as np
from nutau_prop.kernels import tau_decay_pdf,tau_decay_kernel

def test_decay_pdf_positive():
    z=np.linspace(0,1,1001)
    assert np.all(tau_decay_pdf(z,-1)>=0)

def test_decay_kernel_probability():
    E=np.geomspace(1e3,1e7,80)
    K=tau_decay_kernel(E,-1)
    assert np.all(K.sum(axis=1)<=1.001)
    assert np.all(K.sum(axis=1)>=0)
