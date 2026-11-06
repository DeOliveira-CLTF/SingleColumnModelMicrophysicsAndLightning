"""Validate local article bindings independently of the rendering functions."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scm_lightning.evaluation.article import OUT,METHODS
from scm_lightning.provenance import ROOT,sha256

pytestmark=pytest.mark.skipif(not (OUT/'audit/current_statistics.json').exists(),reason='Local article derivatives are not distributed')


def test_article_members_scales_splits_and_csv_metrics():
    p=pd.read_csv(OUT/'data/calibrated_member_predictions.csv.gz');p.time=pd.to_datetime(p.time,utc=True)
    observations=pd.read_csv(OUT/'data/observations_15min.csv');observations.time=pd.to_datetime(observations.time,utc=True)
    table=pd.read_csv(OUT/'tables/Table_03_observational_performance.csv')
    assert len(table)==18 and set(table.method)==set(METHODS)
    assert (p.groupby(['method','configuration','time']).member.nunique()==41).all()
    for row in table.itertuples():
        f=p[(p.method==row.method)&(p.configuration==row.configuration)&(p.member!=0)&(p.split=='evaluation')]
        members=f.pivot(index='time',columns='member',values='prediction_A')
        assert members.shape==(92,40) and members.notna().all().all()
        y=observations.set_index('time').loc[members.index,'count'].values
        mean=np.sum(members.to_numpy(),axis=1)/40
        assert np.sqrt(np.sum((mean-y)**2)/92)==pytest.approx(row.rmse,rel=1e-12)
        assert (mean-y).sum()==pytest.approx(row.accumulated_error,abs=1e-8)
        np.testing.assert_allclose(f.prediction_A,f.predictor*row.coefficient,rtol=1e-12)
        dev=p[(p.method==row.method)&(p.configuration==row.configuration)&(p.member!=0)&(p.split=='development')].groupby('time').predictor.mean()
        yd=observations.set_index('time').loc[dev.index,'count'].values
        assert np.dot(dev,yd)/np.dot(dev,dev)==pytest.approx(row.coefficient,rel=1e-12)
        assert members.index.min()==pd.Timestamp('2018-12-14T01:00:00Z')


def test_article_statistics_and_signed_contrasts():
    stats=pd.read_csv(OUT/'data/ensemble_statistics_A.csv')
    assert (stats.n==40).all() and len(stats)==192*3*6
    assert (stats.minimum<=stats['mean']).all() and (stats['mean']<=stats.maximum).all()
    assert (stats.p05<=stats.p95).all()
    contrasts=pd.read_csv(OUT/'data/rmse_contrasts.csv')
    np.testing.assert_allclose(contrasts.delta_rmse,contrasts.rmse_a-contrasts.rmse_b)
    assert (abs(contrasts.delta_rmse)<=contrasts.forecast_RMS_distance+1e-8).all()
    assert (abs(contrasts.delta_rmse)<=contrasts.rmse_a+contrasts.rmse_b+1e-8).all()


def test_article_figure_csv_bindings_and_observation_accumulation():
    stats=pd.read_csv(OUT/'data/ensemble_statistics_A.csv')
    bound=pd.read_csv(OUT/'data/Figure_03_ensemble_predictions_statistics.csv')
    expected=stats[stats.configuration.isin(['TH_C0','WS_C0'])]
    expected=expected[pd.to_datetime(expected.time,utc=True)>=pd.Timestamp('2018-12-14T01:00Z')]
    pd.testing.assert_frame_equal(bound.reset_index(drop=True),expected.reset_index(drop=True))
    o=pd.read_csv(OUT/'data/Figure_02_observations.csv')
    np.testing.assert_array_equal(o.IC+o.CG,o['count'])
    np.testing.assert_array_equal(o.cumulative_total,o['count'].cumsum())
    assert o.covered.all() and len(o)==192
    table=pd.read_csv(OUT/'tables/Table_03_observational_performance_ranked.csv')
    assert table.rmse.is_monotonic_increasing


def test_article_exports_are_vector_pdf_and_300dpi_png():
    from PIL import Image
    from pypdf import PdfReader
    figures=list((OUT/'figures').glob('Figure_*.png'))
    assert len(figures)==5
    for path in figures:
        with Image.open(path) as img:assert min(img.info['dpi'])>=299.9
        pdf=PdfReader(str(path.with_suffix('.pdf')))
        assert len(pdf.pages)==1
        objects=pdf.pages[0]['/Resources'].get('/XObject',{})
        if hasattr(objects,'get_object'):objects=objects.get_object()
        for obj in objects.values():assert obj.get_object().get('/Subtype')!='/Image'


def test_original_graphics_and_observations_are_preserved():
    original=ROOT/'results/argentina_20181213/figures'
    identity=json.loads((original/'manifest.json').read_text())
    for figure in identity['figures']:assert sha256(original/figure['path'])==figure['sha256']
    source=json.loads((OUT/'audit/source_identity.json').read_text())
    assert sha256(ROOT/'results/argentina_20181213/observations/counts_15min.csv')==source['observations_sha256']
