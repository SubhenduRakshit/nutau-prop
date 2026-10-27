import math

import numpy as np
import pandas as pd

from scripts import run_public_ehe_audit_v0_9 as audit


def test_icecube_area_units_and_flavor_sum():
    _, areas = audit.response_areas()
    assert abs(areas[0, 0] - 12700.0) < 1.0
    assert np.allclose(areas[:, 0], areas[:, 1:].sum(axis=1), rtol=0.01)


def test_published_reference_flux_normalization():
    center = 3.162e8
    limit = 7.535e-9
    reference_flux = limit / (3.0 * center * audit.ENERGY)
    center_index = np.argmin(abs(audit.ENERGY - center))
    assert math.isclose(3.0 * audit.ENERGY[center_index] ** 2 * reference_flux[center_index], limit, rel_tol=0.001)


def test_benchmark_and_direct_track_topology():
    summary = pd.read_csv(audit.RESULTS / "case6_public_audit_v0_9.csv")
    assert np.allclose(summary.KM3_brighttrack_direct_Case6, summary.KM3_brighttrack_SM)
    assert np.all(summary.KM3_NC_equivalent_proxy > summary.KM3_brighttrack_SM)
    assert np.all(np.isfinite(summary.public_window_allflavor_scale_proxy))
    assert np.all(summary.public_window_allflavor_scale_proxy > 0.0)
