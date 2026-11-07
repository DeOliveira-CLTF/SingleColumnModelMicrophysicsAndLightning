import numpy as np
from scm_lightning.evaluation.mynn_audit_products import supported_reach,saturation_specific


def test_reach_requires_both_interfaces_and_never_bridges_inactive_gap():
    mask=np.array([[1,1,0],[1,0,1],[1,1,1],[1,1,1]],bool)
    assert supported_reach(mask,[100,300,600],np.array([200,200,700,np.nan])).tolist()==[True,False,False,False]


def test_surface_saturation_varies_with_skin_temperature_and_pressure():
    value=saturation_specific(np.array([280.,300.]),np.array([1e5,1e5]))
    np.testing.assert_allclose(value,[.006187804812187,.02228140411649],rtol=.001)
    assert value[1]>value[0]
