import math
from nutau_prop.icecube import poisson_asimov_q

def test_significance_scales_as_sqrt_exposure():
    # Poisson Asimov q scales linearly if all counts scale by k.
    true=[11.0,22.0]
    test=[10.0,20.0]
    q1=poisson_asimov_q(true,test)
    k=5.0
    q2=poisson_asimov_q([k*x for x in true],[k*x for x in test])
    assert abs(q2/k-q1)<1e-12
    assert abs(math.sqrt(q2)/math.sqrt(q1)-math.sqrt(k))<1e-12
