from dataclasses import dataclass
from .constants import M_W, SIN2W_DEFAULT, CM2_TO_GEV2_INV

@dataclass(frozen=True)
class BSMModel:
    kind: str = "sm"  # sm | charge_radius | derivative_em
    r2_cm2: float = 0.0
    eps_em_inf: float = 0.0
    mediator_mass_gev: float = 0.0
    antinu_sign: float = 1.0

def charge_radius_eps_em(r2_cm2, sin2w=SIN2W_DEFAULT):
    """
    Defines delta g_V^q = Q_q * eps_em.
    Literature convention used here:
      delta g_V^q = -(2/3) Q_q M_W^2 sin^2(theta_W) <r^2>.
    <r^2> is converted from cm^2 to GeV^-2.
    """
    r2_gev2_inv = r2_cm2 * CM2_TO_GEV2_INV
    return -(2.0/3.0) * M_W**2 * sin2w * r2_gev2_inv

def delta_gv(q_charge, q2, model: BSMModel, sin2w=SIN2W_DEFAULT, antinu=False):
    if model.kind == "sm":
        return 0.0
    sign = model.antinu_sign if antinu else 1.0
    if model.kind == "charge_radius":
        eps = charge_radius_eps_em(model.r2_cm2, sin2w)
        return sign * q_charge * eps
    if model.kind == "derivative_em":
        m2 = model.mediator_mass_gev**2
        ff = q2/(q2+m2) if (q2+m2) > 0 else 0.0
        return sign * q_charge * model.eps_em_inf * ff
    raise ValueError(f"Unknown BSM kind: {model.kind}")
