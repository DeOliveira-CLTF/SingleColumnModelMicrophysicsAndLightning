"""Machine-readable delivery/provenance, quality tables and unresolved inputs."""
from pathlib import Path
import json
import re
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import pandas as pd
from leona.evaluate import ensemble_statistics
from leona.provenance import ROOT,git,sha256,write_json

def main():
    out=ROOT/'products/leona';tables=out/'tables';tables.mkdir(exist_ok=True)
    c=json.loads((ROOT/'configs/leona.json').read_text())
    p=pd.read_csv(out/'predictors_15min.csv.gz')
    ensemble_statistics(p,value='predictor').to_csv(tables/'ensemble_raw_statistics_full_period.csv',index=False)
    availability=p.groupby(['method','configuration']).agg(rows=('predictor','size'),valid_rows=('predictor','count'),members=('member','nunique')).reset_index()
    registry=pd.DataFrame(json.loads((ROOT/'configs/parameterizations.json').read_text())['methods'])
    availability=availability.merge(registry[['id','version','units','scm','requirements']],left_on='method',right_on='id')
    availability.to_csv(tables/'parameterization_availability.csv',index=False)
    data=[]
    for record in json.loads((ROOT/'data_private/leona/manifest.json').read_text()):
        data.append({'type':'provided_observation','source':record['source'],'status':record.get('status','verified_copy'),**record})
    for record in json.loads((ROOT/'data_private/leona/era5/inventory.json').read_text()):
        data.append({'type':'ERA5_EDA','source':record['path'],'status':record['status'],'sha256':record.get('sha256')})
    pd.DataFrame(data).to_csv(tables/'data_availability_cache.csv',index=False)
    ensemble=ROOT/c['ensemble'].get('active_directory','work/leona/ensemble')
    models=json.loads((ensemble/'run_checks.json').read_text())
    pd.DataFrame(models).to_csv(tables/'member_physical_numerical_checks.csv',index=False)
    baseline='5c1655c6d23514d9a993b61b208dc9454f0e36ef'
    protected=['src/module_mp_thompson.F','src/wsm6','src/mynn','src/mynn_offline_coupling.F90']
    core_diff=git('diff',baseline,'--',*protected)
    old=git('show',baseline+':src/mp_offline_driver.F90')
    new=(ROOT/'src/mp_offline_driver.F90').read_text()
    pattern=r'  subroutine apply_forcing_and_advance\(.*?  end subroutine apply_forcing_and_advance'
    integrator_same=re.search(pattern,old,re.S).group()==re.search(pattern,new,re.S).group()
    if core_diff or not integrator_same:
        raise ValueError('Protected physics or integrator differs from audited baseline')
    pending={
        'CAPExP':'Published functional product implemented; native near-surface isentropic CAPE diagnostic/provider not implemented. Unavailable in native SCM postprocessing.',
        'convective_velocity_methods':'PW10, IFluxT, McCaul2009 and PR92_wmax require resolved convective velocity and collocated hydrometeors; wa/mf_w are not substitutes.',
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
        'ensemble_realizations':str(ensemble/'realized_perturbations.npz'),
        'ensemble_realizations_sha256':sha256(ensemble/'realized_perturbations.npz'),
        'source_config_sha256':sha256(ROOT/'configs/leona.json'),'pending_dependencies':pending,
        'redistribution':'private observations and derivatives remain local; no push or external publication'})

if __name__=='__main__':main()
