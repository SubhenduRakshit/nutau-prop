"""Checks of normalization invariants, not historical numerical answers."""
import math
import pytest
from scipy.integrate import quad
from nutau_prop.pkr_resolved import (partonic_dsdy,invisible_width,
    allowed_gq_squared,quark_width_coefficient)

def test_massless_crossing_and_common_cutoff_equivalence():
    for y in (.001,.1,.5,.9,.999):
        d=partonic_dsdy(1e5,y,50,.002,.03,'4.4')
        assert partonic_dsdy(1e5,y,50,.002,.03,'4.3') == pytest.approx(4*d)
        assert partonic_dsdy(1e5,y,50,.001,.03,'4.3') == pytest.approx(d)
    assert invisible_width(50,.001,'4.3') == pytest.approx(invisible_width(50,.002,'4.4'))

def test_contact_limit_and_endpoints():
    s,M,a,g=100.,1e5,.001,.03
    result=quad(lambda y:partonic_dsdy(s,y,M,a,g),0,1,epsabs=1e-40)[0]
    expected=(a*g)**2*s*s/(24*math.pi*M**4)
    assert result == pytest.approx(expected,rel=2*s/M**2,abs=0)
    assert partonic_dsdy(s,0,M,a,g)==0
    assert partonic_dsdy(s,1,M,a,g)==0
    assert partonic_dsdy(s,.5,M,a,g)>0

@pytest.mark.parametrize('op',['4.4','4.3'])
@pytest.mark.parametrize('a',[.0001,.001,.01])
@pytest.mark.parametrize('dijet',[None,.03,.08])
def test_exact_coupling_bound_is_feasible_and_saturated(op,a,dijet):
    M,pmax=50.,.0003
    z=allowed_gq_squared(M,a,op,pmax,dijet)
    I=invisible_width(M,a,op); B=quark_width_coefficient(M)
    P=z*I/(B*z+I); D=z-P
    assert P<=pmax*(1+1e-12)
    assert B*z+I<=.1*M*(1+1e-12)
    if dijet is not None: assert D<=dijet**2*(1+1e-12)
    violations=[]
    zz=z*(1+1e-5); PP=zz*I/(B*zz+I)
    violations.extend([PP>pmax,B*zz+I>.1*M,zz>4*math.pi])
    if dijet is not None: violations.append(zz-PP>dijet**2)
    assert any(violations)
