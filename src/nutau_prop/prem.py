import numpy as np
from .constants import R_EARTH_KM

def density_g_cm3(r_km):
    """PREM density parametrization, r = radius from Earth's center."""
    r = np.asarray(r_km, dtype=float)
    x = r / R_EARTH_KM
    rho = np.zeros_like(x)

    masks = [
        r < 1221.5,
        (r >= 1221.5) & (r < 3480.0),
        (r >= 3480.0) & (r < 5701.0),
        (r >= 5701.0) & (r < 5771.0),
        (r >= 5771.0) & (r < 5971.0),
        (r >= 5971.0) & (r < 6151.0),
        (r >= 6151.0) & (r < 6346.6),
        (r >= 6346.6) & (r < 6356.0),
        (r >= 6356.0) & (r < 6368.0),
        (r >= 6368.0) & (r <= R_EARTH_KM),
    ]
    formulas = [
        13.0885 - 8.8381*x**2,
        12.5815 - 1.2638*x - 3.6426*x**2 - 5.5281*x**3,
        7.9565 - 6.4761*x + 5.5283*x**2 - 3.0807*x**3,
        5.3197 - 1.4836*x,
        11.2494 - 8.0298*x,
        7.1089 - 3.8045*x,
        2.6910 + 0.6924*x,
        np.full_like(x, 2.9),
        np.full_like(x, 2.6),
        np.full_like(x, 1.02),
    ]
    for m, f in zip(masks, formulas):
        rho[m] = f[m] if np.ndim(f) else f
    return rho.item() if rho.ndim == 0 else rho

def chord_profile(cos_zenith_abs, max_step_km=10.0):
    """
    Earth chord for |cos(theta_z)| in [0,1]; 1 = vertical through diameter.
    Returns segment midpoints with ds[km], rho[g/cm^3], radius[km].
    """
    mu = float(cos_zenith_abs)
    if not 0.0 <= mu <= 1.0:
        raise ValueError("|cos(theta_z)| must lie in [0,1].")
    if mu == 0:
        return dict(ds_km=np.array([]), rho=np.array([]), r_km=np.array([]))
    L = 2.0 * R_EARTH_KM * mu
    b = R_EARTH_KM * np.sqrt(max(0.0, 1.0-mu*mu))
    n = max(1, int(np.ceil(L/max_step_km)))
    edges = np.linspace(-L/2.0, L/2.0, n+1)
    mids = 0.5*(edges[:-1]+edges[1:])
    ds = np.diff(edges)
    r = np.sqrt(b*b + mids*mids)
    rho = density_g_cm3(r)
    return dict(ds_km=ds, rho=rho, r_km=r)

def column_depth_g_cm2(cos_zenith_abs, max_step_km=2.0):
    p = chord_profile(cos_zenith_abs, max_step_km=max_step_km)
    return float(np.sum(p["rho"] * p["ds_km"] * 1.0e5))
