from scm_lightning.config import load_config
"""Machine-readable delivery/provenance, quality tables and unresolved inputs."""
from pathlib import Path
import json
import re
import sys
import pandas as pd
from scm_lightning.evaluation import ensemble_statistics
from scm_lightning.provenance import ROOT,git,sha256,write_json

def main():
    out=ROOT/'results/argentina_20181213';tables=out/'tables';tables.mkdir(exist_ok=True)
    c=load_config(ROOT/'configs/cases/argentina_20181213.yaml')
    p=pd.read_csv(out/'predictors_15min.csv.gz')
    ensemble_statistics(p,value='predictor').to_csv(tables/'ensemble_raw_statistics_full_period.csv',index=False)
    availability=p.groupby(['method','configuration']).agg(rows=('predictor','size'),valid_rows=('predictor','count'),members=('member','nunique')).reset_index()
    registry=pd.DataFrame(load_config(ROOT/'configs/lightning/schemes.yaml')['methods'])
    availability=availability.merge(registry[['id','version','units','scm','requirements']],left_on='method',right_on='id')
    availability.to_csv(tables/'parameterization_availability.csv',index=False)
    data=[]
    for record in json.loads((ROOT/'data/intermediate/input_inventory.json').read_text()):
        data.append({'type':'provided_observation','source':record['source'],'status':record.get('status','verified_copy'),**record})
    for record in json.loads((ROOT/'data/raw/era5/inventory.json').read_text()):
        data.append({'type':'ERA5_EDA','source':record['path'],'status':record['status'],'sha256':record.get('sha256')})
    pd.DataFrame(data).to_csv(tables/'data_availability_cache.csv',index=False)
    ensemble=ROOT/c['ensemble'].get('active_directory','runs/argentina_20181213/_execution/ensemble')
    models=json.loads((ensemble/'run_checks.json').read_text())
    pd.DataFrame(models).to_csv(tables/'member_physical_numerical_checks.csv',index=False)
    baseline='5c1655c6d23514d9a993b61b208dc9454f0e36ef'
    protected=['src/module_mp_thompson.F','src/wsm6','src/mynn']
    core_diff=git('diff',baseline,'--',*protected)
    old=git('show',baseline+':src/mp_offline_driver.F90')
    new=(ROOT/'src/mp_offline_driver.F90').read_text()
    pattern=r'  subroutine apply_forcing_and_advance\(.*?  end subroutine apply_forcing_and_advance'
    integrator_same=re.search(pattern,old,re.S).group()==re.search(pattern,new,re.S).group()
    if core_diff or not integrator_same:
        raise ValueError('Protected physics or integrator differs from audited baseline')
    pending={
        'CAPExP':'Original isentropic definition remains conditional. CAPExP_SBCAPE15_1D_adaptado now uses MetPy pseudoadiabatic CAPE from simulated profiles, diagnosed every 15min; separate control-only exploratory outputs.',
        'convective_velocity_methods':'Original CRM methods remain conditional. Separate LS/EDMF 1D adaptations are computed in vertical_diagnostics. EDMF w=mf_w/rhoz, inactive plume masked; grid-hydrometeor fluxes explicitly assume plume representativeness. No LS+EDMF sum.',
        'storm_geometry_methods':'Radar reflectivity, storm volumes, terminal fall speeds/divergence required by remaining Barthe formulas are unavailable; no silent substitutions.',
        'IGRA':'Alternative station inventory documented; alternative profile acquisition/parser not implemented. Primary Wyoming and three provided Córdoba files are cached.',
        'LEONA_manuscript':'No manuscript file path supplied; case metadata follows user dates/image and meteorological support manifest.',
        'budget_scope':'Entry-density microphysics water residual is retained; full variable-density atmospheric water/energy closure is not certified.'}
    write_json(out/'pending_dependencies.json',pending)
    write_json(out/'execution_manifest.json',{'baseline':baseline,'head':git('rev-parse','HEAD'),
        'commits':git('log','--oneline',baseline+'..HEAD').splitlines(),'config':c,
        'protected_physics_diff_empty':not bool(core_diff),'integrator_subroutine_unchanged':integrator_same,
        'ensemble_runs':len(models),'all_member_preflights_passed':all(r.get('physical_preflight_passed',False) for r in models),
        'case_controls_sha256':sha256(out/'tests/case_controls.json'),
        'instrumentation_invariance_sha256':sha256(out/'tests/driver_checks.json'),
        'tests_xml_sha256':sha256(out/'tests/all.xml'),'ensemble_stability_sha256':sha256(out/'ensemble_stability.json'),
        'native_product_index_sha256':sha256(out/'native_product_index.json'),
        'native_utc_clock_check_sha256':sha256(out/'tests/native_utc_clock.json'),
        'preferred_native_outputs':'results/argentina_20181213/simulations/<configuration>/<member>/output_utc.nc',
        'ensemble_realizations':str(ensemble/'realized_perturbations.npz'),
        'ensemble_realizations_sha256':sha256(ensemble/'realized_perturbations.npz'),
        'source_config_sha256':sha256(ROOT/'configs/cases/argentina_20181213.yaml'),'pending_dependencies':pending,
        'redistribution':'private observations and derivatives remain local; no push or external publication'})

if __name__=='__main__':main()
