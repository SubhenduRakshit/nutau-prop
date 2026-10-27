import numpy as np
from pathlib import Path
from nutau_prop.icecube import (
    parse_total_effective_area,powerlaw_integral,solid_angle_from_cosz,
    poisson_asimov_q
)

def test_effective_area_logE_cosz(tmp_path):
    p=tmp_path/"effective_area.nu_tau.txt"
    p.write_text(
        "# log10Elo log10Ehi coszlo coszhi Aeff_m2\n"
        "4 5 -1 -0.5 2.0\n"
        "5 6 -0.5 0 3.0\n"
    )
    d=parse_total_effective_area(p)
    assert np.allclose(d.E_lo_GeV,[1e4,1e5])
    assert np.allclose(d.E_hi_GeV,[1e5,1e6])
    assert np.allclose(d.Aeff_m2,[2,3])

def test_effective_area_degree_zenith(tmp_path):
    p=tmp_path/"effective_area.nu_tau.txt"
    p.write_text("1e4 1e5 90 180 2.0\n")
    d=parse_total_effective_area(p)
    assert d.cosz_lo.iloc[0] < -0.999
    assert abs(d.cosz_hi.iloc[0]) < 1e-12

def test_powerlaw_integral_positive():
    x=powerlaw_integral(2e-18,1e5,2.5,1e4,1e5)
    assert x>0

def test_solid_angle_half_sky():
    assert abs(solid_angle_from_cosz(-1,0)-2*np.pi)<1e-12

def test_poisson_q_identity():
    x=np.array([1.,2.,3.])
    assert abs(poisson_asimov_q(x,x))<1e-12
