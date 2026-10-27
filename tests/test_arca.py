import math
import numpy as np
from nutau_prop.arca import fig7_response, calibrate_fig7_response

def _cfg():
    return {
      "arca":{
        "reference_flux_phi0_at_1GeV":1.2e-8,
        "reference_gamma":2.0,
        "reference_track_rate_per_year":220.0,
        "reference_shower_rate_per_year":60.0,
      }
    }

def test_fig7_positive():
    E=np.logspace(4,6,20)
    R=fig7_response(E)
    for d in ("shower","track"):
        for k,v in R[d].items():
            if k.endswith("cos_bins"): continue
            if hasattr(v,"__len__"):
                assert np.all(np.asarray(v)>0)

def test_calibration_finite_positive():
    a,b=calibrate_fig7_response(_cfg())
    assert np.isfinite(a) and a>0
    assert np.isfinite(b) and b>0
