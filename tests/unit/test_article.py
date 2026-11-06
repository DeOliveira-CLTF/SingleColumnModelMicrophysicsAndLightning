import numpy as np
import pandas as pd
import pytest
from scm_lightning.evaluation.article import statistics,performance


def test_article_spread_excludes_control_and_uses_sample_std():
    frame=pd.DataFrame({'time':['2018-12-14T01:00:00Z']*41,'method':['IxG']*41,
                        'configuration':['TH_C0']*41,'member':range(41),
                        'prediction_A':[1.e9]+list(range(1,41))})
    row=statistics(frame,'prediction_A').iloc[0]
    assert row.n==40 and row['mean']==20.5 and row['median']==20.5
    assert row['std']==pytest.approx(np.sqrt(40*41/12))
    assert row.minimum==1 and row.maximum==40
    assert row.p05==pytest.approx(2.95) and row.p95==pytest.approx(38.05)


def test_article_count_metrics_and_undefined_correlations():
    value=performance([0.,2.],[0.,4.])
    assert value['rmse']==pytest.approx(np.sqrt(2.))
    assert value['bias']==1. and value['mae']==1.
    assert value['accumulated_error']==2. and value['peak_time_error_minutes']==0
    flat=performance([1.,2.,3.],[0.,0.,0.])
    assert np.isnan(flat['pearson']) and np.isnan(flat['spearman'])
    missing=performance([1.,2.],[np.nan,np.nan])
    assert missing['n']==0 and np.isnan(missing['rmse'])


def test_article_peak_ties_and_rmse_reverse_triangle():
    a=np.array([0.,2.,2.,0.]);b=np.array([3.,0.,0.,3.]);obs=np.array([0.,0.,4.,0.])
    aa=performance(obs,a);bb=performance(obs,b)
    assert aa['peak_time_error_minutes']==-15
    assert abs(aa['rmse']-bb['rmse'])<=np.sqrt(np.mean((a-b)**2))
