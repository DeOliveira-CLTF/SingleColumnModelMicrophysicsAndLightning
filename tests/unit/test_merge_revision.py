import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from scm_lightning.evaluation.merge import regional_field,metrics


def test_missing_precipitation_cell_never_becomes_zero_or_a_renormalized_box(tmp_path):
    path=tmp_path/'fixture.tif';a=np.ones((20,20),dtype='float64');a[10,10]=-9999.9
    with rasterio.open(path,'w',driver='GTiff',height=20,width=20,count=2,dtype='float64',transform=from_origin(-60, -31,.1,.1)) as d:
        d.write(a,1);d.write(np.full_like(a,99),2)
        d.update_tags(1,GRIB_ELEMENT='PREC',GRIB_UNIT='[kg/(m^2)]',GRIB_VALID_TIME='1544662800')
        d.update_tags(2,GRIB_ELEMENT='PRMSL',GRIB_UNIT='[Pa]',GRIB_VALID_TIME='1544662800')
    amount,field,meta=regional_field(path,[-59.5,-32.5,-58.5,-31.5],'2018-12-13T01:00Z')
    assert np.isnan(amount) and meta['covered_area_fraction']<1
    assert field.max()==1  # station count band never used
    with pytest.raises(ValueError,match='valid-time'):regional_field(path,[-59.5,-32.5,-58.5,-31.5],'2018-12-13T02:00Z')


def test_precipitation_metrics_share_actual_hourly_support_and_no_fake_timing():
    f=pd.DataFrame({'time':['2018-12-13T00:00Z','2018-12-13T01:00Z','2018-12-13T02:00Z'],
        'MERGE_mm':[1.,np.nan,3.],'ERA5_mm':[4.,5.,1.],'TH_C0_mm':[4.,5.,1.],'WS_C0_mm':[4.,5.,1.]})
    table,mask=metrics(f)
    np.testing.assert_array_equal(mask,[True,False,True])
    assert table.n_shared_hourly_intervals.eq(2).all()
    assert table.peak_time_error_hours.eq(-2).all()
    assert not table.full_period_accumulation_valid.any()


def test_constant_series_has_no_defined_peak_timing():
    f=pd.DataFrame({'time':['2018-12-13T00:00Z','2018-12-13T01:00Z'],
        'MERGE_mm':[1.,1.],'ERA5_mm':[1.,1.],'TH_C0_mm':[1.,1.],'WS_C0_mm':[1.,1.]})
    table,_=metrics(f);assert table.peak_time_error_hours.isna().all()
