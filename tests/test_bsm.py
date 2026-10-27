from nutau_prop.bsm import BSMModel,delta_gv,charge_radius_eps_em

def test_derivative_forward_blind():
    m=BSMModel(kind="derivative_em",eps_em_inf=0.05,mediator_mass_gev=0.01)
    assert delta_gv(2/3,0.0,m)==0.0

def test_charge_radius_zero():
    assert charge_radius_eps_em(0.0)==0.0
