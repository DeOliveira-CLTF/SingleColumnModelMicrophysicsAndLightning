from scm_lightning.config import load_config
"""Cross-file scientific data contracts; no ranking or interpretation report."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
from scm_lightning.provenance import ROOT,write_json

def main():
    out=ROOT/'results/argentina_20181213';c=load_config(ROOT/'configs/cases/argentina_20181213.yaml')
    ensemble=ROOT/c['ensemble'].get('active_directory','runs/argentina_20181213/_execution/ensemble')
    plan=json.loads((ensemble/'manifest.json').read_text());runs=json.loads((ensemble/'run_checks.json').read_text())
    n=c['ensemble']['perturbed_members']
    assert len(runs)==6*(n+1) and all(r['physical_preflight_passed'] for r in runs)
    o=pd.read_csv(out/'observations/counts_15min.csv');types=pd.read_csv(out/'observations/counts_by_type_15min.csv')
    totals=types[types.type!='total'].groupby('time')['count'].sum()
    np.testing.assert_array_equal(totals.values,o['count'].values)
    assert len(o)==192 and o.covered.all() and not o['count'].isna().any()
    h=pd.read_csv(out/'evaluation/heldout_predictions.csv');h['time']=pd.to_datetime(h.time,utc=True)
    assert len(h)==6*11*(n+1)*92 and (h.split=='evaluation').all()
    assert h.time.min()>=pd.Timestamp(c['split']['evaluation_start']) and h.time.max()<pd.Timestamp(c['end'])
    assert (h.groupby(['time','method','configuration']).member.nunique()==n+1).all()
    coeff=pd.read_csv(out/'evaluation/coefficients.csv')
    assert not coeff.duplicated(['method','configuration']).any()
    assert len(coeff)==66 and (coeff.loc[coeff.status=='fitted','scale']>=0).all()
    assert (coeff.loc[coeff.status=='fitted','n_development']==96).all()
    assert h[h.method.isin(['PW10','IFluxT','McCaul2009','PR92_wmax'])].prediction_A.isna().all()
    ref=coeff[coeff.configuration==c['reference_configuration']].set_index('method').scale
    for method,f in h.groupby('method'):
        np.testing.assert_allclose(f.reference_scale,ref[method],equal_nan=True)
        for configuration,g in f.groupby('configuration'):
            expected=coeff[(coeff.method==method)&(coeff.configuration==configuration)].scale.iloc[0]
            np.testing.assert_allclose(g.scale,expected,equal_nan=True)
    stats=pd.read_csv(out/'evaluation/ensemble_statistics.csv')
    assert (stats.loc[stats['mean'].notna(),'n']==n).all()
    assert stats.loc[stats['mean'].notna(),'std'].ge(0).all()
    assert (stats.loc[stats.n>0,'q05']<=stats.loc[stats.n>0,'q95']).all()
    inheritance={}
    if n>20:
        oldplan=json.loads((ROOT/'runs/argentina_20181213/_execution/ensemble/manifest.json').read_text())
        oldruns=json.loads((ROOT/'runs/argentina_20181213/_execution/ensemble/run_checks.json').read_text())
        old_draws=np.load(ROOT/'runs/argentina_20181213/_execution/ensemble/realized_perturbations.npz')
        new_draws=np.load(ensemble/'realized_perturbations.npz')
        np.testing.assert_array_equal(old_draws['delta'],new_draws['delta'][:20])
        np.testing.assert_array_equal(old_draws['weights'],new_draws['weights'][:20])
        old_hash={(r['configuration'],r['member']):r['sha256'] for r in oldruns}
        for run in runs:
            if run['member']<=20:assert run['sha256']==old_hash[(run['configuration'],run['member'])]
        old_inputs={(r['configuration'],r['member']):r['forcing_sha256'] for r in oldplan['runs']}
        for run in plan['runs']:
            if run['member']<=20:assert run['forcing_sha256']==old_inputs[(run['configuration'],run['member'])]
        inheritance={'first_20_draws_bit_identical':True,'first_20_forcings_hash_identical':True,'first_20_outputs_hash_identical':True}
    write_json(out/'tests/delivery_contracts.json',{'status':'passed','ensemble_members':n,'configurations':6,
        'heldout_intervals':92,'development_intervals':96,'gap_intervals':4,'observation_bins':192,
        'coverage_confirmed_by_user':True,'unknown_methods_remain_nan':True,'coefficient_application_checked':True,
        'input_output_inheritance':inheritance})
    print('Delivery contracts passed')

if __name__=='__main__':main()
