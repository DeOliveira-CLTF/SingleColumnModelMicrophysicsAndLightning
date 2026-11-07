import numpy as np
import pandas as pd
import pytest
from scm_lightning.observations.profiles import read_campaign,read_csv,forcing_eligibility
from scm_lightning.calibration.independent import fit
from scm_lightning.evaluation.independent import metric_values
from test_independent_events import inputs


def test_one_minute_fit_uses_real_one_minute_support_not_fifteen():
    c,p,o,p5,o5=inputs();c['interval_minutes']=1
    clock=pd.date_range(c['start'],c['end'],freq='min',inclusive='left').astype(str)
    p=pd.DataFrame([dict(time=t,method=m,configuration='TH_C0',member=0,predictor=float(i+1))
        for m in ['CAPExP','IxG','McCaulF3'] for i,t in enumerate(clock)])
    o=pd.DataFrame({'time':clock,'count':np.arange(1,len(clock)+1)*2.,'covered':True})
    coef,baseline=fit(c,p,o,p5,o5)
    assert set(coef.loc[coef.method!='McCaulF3_LS_1D_adapted','n_intervals'])=={2880}
    assert set(coef.loc[coef.method!='McCaulF3_LS_1D_adapted','fit_interval_minutes'])=={1}
    assert coef.set_index('method').loc['McCaulF3_LS_1D_adapted','fit_interval_minutes']==5


def test_one_minute_peak_timing_is_one_minute_not_fifteen():
    r=metric_values([0.,5.,0.],[0.,0.,5.],1.,minutes=1)
    assert r['peak_time_error_minutes']==1


def test_campaign_utc_rollover_si_units_and_balloon_speed_not_w(tmp_path):
    p=tmp_path/'edt_20181214_0002.txt'
    p.write_text('Station name\tCSU_atmos1\nBalloon release date and time\t2018-12-14T00:02:11\n'
        'Elapsed time HeightMSL P Temp Dewp RH Speed Dir AscRate Lat Lon GpsHeightMSL\n'
        's m hPa C C percent m/s deg m/s deg deg m\n'
        '0 100 1000 20 15 70 10 90 0 -31 -64 100\n'
        '60 400 965 18 13 72 12 180 5 -31.01 -64.01 400\n',encoding='utf-8')
    f,_=read_campaign(p)
    assert f.time.iloc[0]==pd.Timestamp('2018-12-14T00:02:11Z')
    assert f.pressure_Pa.iloc[0]==100000 and f.temperature_K.iloc[0]==293.15
    assert 0<f.rv_kg_kg_dry.iloc[0]<.02
    assert f.u_m_s.iloc[0]==pytest.approx(-10) and f.v_m_s.iloc[1]==pytest.approx(12)
    assert 'w_LS' not in f and 'w_EDMF' not in f
    assert f.height_geopotential_m.iloc[1]<f.height_msl_m.iloc[1]


def test_distant_future_sonde_cannot_become_case_initial_condition(tmp_path):
    p=tmp_path/'2019102812-83928.csv'
    pd.DataFrame({'time':['2019-10-28 11:34:00']*3,'latitude':[-29.78]*3,'longitude':[-57.03]*3,
        'pressure_hPa':[1000,500,50],'temperature_C':[20,-10,-55],
        'geopotential height_m':[74,5500,21000],'mixing ratio_g/kg':[10,2,.001]}).to_csv(p,index=False)
    f,_=read_csv(p)
    c={'area':[-56.75,-32.25,-55.75,-31.25],'start':'2019-10-28T00:00Z'}
    result=forcing_eligibility(f,c)
    assert result['distance_to_column_center_km']>200 and result['age_at_case_start_hours']<0
    assert not result['initial_condition_eligible'] and not result['boundary_condition_eligible']


def test_missing_time_is_not_invented_from_nominal_filename(tmp_path):
    p=tmp_path/'2018121412-87344.csv'
    pd.DataFrame({'pressure_hPa':[1000,900],'temperature_C':[20,15],'geopotential height_m':[0,1000],
        'mixing ratio_g/kg':[10,8]}).to_csv(p,index=False)
    f,_=read_csv(p);assert f.time.isna().all() and f.latitude.isna().all()


def test_exact_cape_cache_does_not_round_profiles_or_fabricate_values(tmp_path):
    from scm_lightning.lightning.cape_cache import CapeCache
    memo=CapeCache(tmp_path/'cape.sqlite3');calls=[]
    p=np.array([100000.,90000.]);t=np.array([300.,290.]);r=np.array([.01,.005])
    assert memo.key(p,t,r)!=memo.key(p.astype('float32'),t.astype('float32'),r.astype('float32'))
    def provider(p,t,r,layer_selection):calls.append(t.copy());return float(t[0])
    assert memo.get(p,t,r,provider)==300.
    assert memo.get(p,t.copy(),r,provider)==300. and len(calls)==1
    changed=t.copy();changed[0]=np.nextafter(changed[0],np.inf)
    assert memo.get(p,changed,r,provider)==changed[0] and len(calls)==2
    memo.seed(p,t+1,r,np.nan);assert np.isnan(memo.get(p,t+1,r,provider))
    memo.close()
    second=CapeCache(tmp_path/'cape.sqlite3')
    assert second.get(p,t,r,provider)==300. and len(calls)==2
    second.close()


def test_pulse_cases_cannot_be_mistaken_for_initial_condition_members():
    from scm_lightning.evaluation.minute import member
    with pytest.raises(ValueError,match='pulse scenarios'):
        member({'folder':'runs/argentina_20181213/updraft_sensitivity/weak/TH_C0','configuration':'TH_C0','member':1,'sha256':'irrelevant'})
