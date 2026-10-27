from nutau_prop.prem import density_g_cm3,column_depth_g_cm2

def test_prem_center():
    assert 13.0 < density_g_cm3(0.0) < 13.2

def test_vertical_column():
    # Rakshit-Reya quote ~1.1e10 g/cm^2.
    x=column_depth_g_cm2(1.0,max_step_km=1.0)
    assert 1.0e10 < x < 1.2e10
