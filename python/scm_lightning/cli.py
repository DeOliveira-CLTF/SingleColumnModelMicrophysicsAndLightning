from scm_lightning.config import load_config
"""python -m scm_lightning.cli --help; no interpretation reports are generated."""
import argparse
import json
from pathlib import Path
import sys
import pandas as pd
from scm_lightning.provenance import ROOT,write_json,sha256,audit

def process_lightning(config, sample=False):
    from scm_lightning.observations import read_events,group_flashes,aggregate
    folder=ROOT/config.get('paths',{}).get('entln','data/raw/entln')
    last=pd.Timestamp(config['end']).normalize()
    if config.get('total_anchor_counts_causal_end',False):last-=pd.Timedelta(days=1)
    days=pd.date_range(pd.Timestamp(config['start']).normalize()-pd.Timedelta(days=1),last,freq='D')
    paths=[folder/'2018-11-25.csv'] if sample else [folder/f'{day:%Y-%m-%d}.csv' for day in days]
    out=ROOT/config.get('paths',{}).get('results','results/argentina_20181213')/('reader_sample' if sample else 'observations')
    out.mkdir(parents=True,exist_ok=True)
    e,d=read_events(paths,config['flag_mapping'])
    f,a=group_flashes(e,config['grouping'],progress=lambda status:write_json(out/'progress.json',status))
    f.to_csv(out/'flashes.csv.gz',index=False)
    a.to_csv(out/'event_membership.csv.gz',index=False)
    d.to_csv(out/'duplicate_rows.csv.gz',index=False)
    manifest={'inputs':[{'path':str(p),'sha256':sha256(p)} for p in paths],
        'events':len(e),'duplicates':len(d),'flashes':len(f),'grouping':config['grouping'],
        'flag_mapping':config['flag_mapping'],'flag_mapping_confirmed':config['flag_mapping_confirmed'],
        'classification_status':'confirmed' if config['flag_mapping_confirmed'] else 'provisional_export_convention_unconfirmed',
        'location':'first detected event; not known physical flash origin',
        'tie_break':'earliest eligible anchor; simultaneous anchors ordered by event_id',
        'duration_limit_ms':700,'radius_km':10,'inclusive_boundaries':True,
        'grouping_domain':'all events in cached files; adjacent days included; no final area selection',
        'outer_file_boundary_policy':('one-day preceding margin; first-event total-flash counts are causal; future events cannot change earlier anchors; end multiplicity/type may be censored'
            if config.get('total_anchor_counts_causal_end',False) else 'first/last 1.4s outside target case; target has full adjacent-day grouping margin'),
        'coverage_status':'unknown until explicit coverage intervals supplied',
        'restricted_data':True}
    write_json(out/'manifest.json',manifest)
    if not sample and config['area'] is not None and config['observation_coverage'] is not None:
        if config['area_selection']['status'] != 'fixed_before_lightning_evaluation':
            raise ValueError('Area must be selected independently before evaluation')
        coverage=config['observation_coverage']
        if isinstance(coverage,str) and coverage=='full_case':
            coverage=[True]*len(pd.date_range(config['start'],config['end'],freq=f"{config['interval_minutes']}min",inclusive='left'))
        counts=aggregate(f,config['start'],config['end'],config['area'],coverage,config['interval_minutes'])
        counts.to_csv(out/'counts_15min.csv',index=False)
        five=aggregate(f,config['start'],config['end'],config['area'],
                       [bool(v) for v in coverage for _ in range(3)],5)
        five.to_csv(out/'counts_5min.csv',index=False)
        manifest.update(area_wsen=config['area'],coverage_status='user_confirmed' if all(coverage) else 'explicit_mask',
                        coverage_confirmation=config.get('coverage_confirmation'),counts_sha256=sha256(out/'counts_15min.csv'))
        write_json(out/'manifest.json',manifest)
    print(json.dumps({'events':len(e),'flashes':len(f),'products':str(out)}))

def query_soundings(config):
    from scm_lightning.download import sounding_inventory,wyoming_profile,sounding_quality
    cache=ROOT/'data/raw/radiosondes/wyoming'
    rows=[]
    for date in pd.date_range(config['start'],config['end'],freq='3h',inclusive='left'):
        try:
            inv=sounding_inventory(date,cache/'inventory')
        except Exception as exc:
            rows.extend({'time':str(date),'station':s,'status':'server_or_schema_failure','reason':str(exc)} for s in config['stations'])
            continue
        for station in config['stations']:
            available=[s for s in inv['stations'] if str(s['stationid'])==station]
            if not available:
                rows.append({'time':str(date),'station':station,'status':'absent_from_inventory','fallback':'IGRA/campaign archive not yet queried'})
            for item in available:
                try:
                    wyoming_profile(date,station,item['src'],cache/'profiles')
                    path=cache/'profiles'/f'{date:%Y%m%d%H}-{station}-{item["src"]}.csv'
                    _,quality=sounding_quality(path)
                    rows.append({'time':str(date),'station':station,'source':item['src'],**quality})
                except Exception as exc:
                    rows.append({'time':str(date),'station':station,'source':item['src'],'status':'profile_fetch_or_schema_failure','reason':str(exc)})
    target=ROOT/'results/argentina_20181213/soundings'
    target.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(target/'availability.csv',index=False)
    write_json(target/'endpoint_verification.json',{'interface':'https://weather.uwyo.edu/upperair/sounding.shtml',
        'inventory':'https://weather.uwyo.edu/cgi-bin/sounding_json','profile':'https://weather.uwyo.edu/wsgi/sounding',
        'verified_from':'upperair/radiosonde.js; live inventory and profile requests','rows':len(rows)})

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--config',default=str(ROOT/'configs/cases/argentina_20181213.yaml'))
    sub=ap.add_subparsers(dest='command',required=True)
    for name in ['audit','cache-inputs','soundings','prepare-era5','download-era5','lightning','sample-lightning','status']:
        sub.add_parser(name)
    ev=sub.add_parser('evaluate'); ev.add_argument('predictors');ev.add_argument('observations');ev.add_argument('output')
    forc=sub.add_parser('forcing');forc.add_argument('pl');forc.add_argument('sl');forc.add_argument('output')
    control=sub.add_parser('controls');control.add_argument('executable');control.add_argument('forcing');control.add_argument('output')
    ens=sub.add_parser('prepare-ensemble');ens.add_argument('forcing');ens.add_argument('eda_profiles');ens.add_argument('control_checks');ens.add_argument('output');ens.add_argument('--inherit')
    eda=sub.add_parser('collocate-eda');eda.add_argument('forcing');eda.add_argument('pl');eda.add_argument('sl');eda.add_argument('output')
    run=sub.add_parser('run-ensemble');run.add_argument('executable');run.add_argument('manifest')
    post=sub.add_parser('postprocess');post.add_argument('simulation');post.add_argument('configuration');post.add_argument('member',type=int);post.add_argument('output')
    args=ap.parse_args(); config=load_config(args.config)
    if args.command=='audit': audit(ROOT/'results/argentina_20181213/audit_current.json')
    elif args.command=='cache-inputs':
        from scm_lightning.download import local_inputs
        local_inputs(config,ROOT/'data')
    elif args.command=='soundings': query_soundings(config)
    elif args.command in ['prepare-era5','download-era5']:
        from scm_lightning.download import acquire_era5,acquire_era5_parallel
        cache=ROOT/config.get('paths',{}).get('era5','data/raw/era5')
        if args.command=='download-era5':acquire_era5_parallel(config,cache)
        else:acquire_era5(config,cache)
    elif args.command in ['lightning','sample-lightning']: process_lightning(config,args.command=='sample-lightning')
    elif args.command=='evaluate':
        if config.get('role')=='evaluation_only':
            raise ValueError('2018 is evaluation-only; use the independent frozen-2019 evaluation entry point')
        from scm_lightning.evaluation import evaluate,ensemble_statistics
        if config['area'] is None or config['area_selection']['status']!='fixed_before_lightning_evaluation':
            raise ValueError('Area must be fixed from meteorology before evaluation')
        if config['observation_coverage'] is None:
            raise ValueError('Observation coverage required; file presence is not complete coverage')
        p=pd.read_csv(args.predictors);o=pd.read_csv(args.observations)
        results=evaluate(p,o,config['reference_configuration'],config['split']['development_end'],config['split']['evaluation_start'])
        out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
        for name,value in zip(['coefficients','heldout_predictions','metrics','fixed_scale_skill_contrasts','physical_trajectories'],results):
            value.to_csv(out/(name+'.csv'),index=False)
        ensemble_statistics(results[1]).to_csv(out/'ensemble_statistics.csv',index=False)
        write_json(out/'manifest.json',{'configuration':config,'predictors_sha256':sha256(args.predictors),'observations_sha256':sha256(args.observations),'spread_scope':'initial_condition_spread_only'})
    elif args.command=='forcing':
        import numpy as np
        from scm_lightning.forcing import build
        build(config,args.pl,args.sl,args.output,np.arange(0,20001,10))
    elif args.command=='controls':
        from scm_lightning.ensemble.experiments import run_controls
        run_controls(args.executable,args.forcing,args.output)
    elif args.command=='prepare-ensemble':
        from scm_lightning.ensemble import prepare
        prepare(config,args.forcing,args.eda_profiles,args.control_checks,args.output,args.inherit)
    elif args.command=='collocate-eda':
        from scm_lightning.ensemble import collocate
        collocate(config,args.forcing,args.pl,args.sl,args.output)
    elif args.command=='run-ensemble':
        from scm_lightning.ensemble import run_prepared
        run_prepared(args.executable,args.manifest)
    elif args.command=='postprocess':
        from scm_lightning.lightning.postprocess import native_diagnostics,integrate_bins
        native=native_diagnostics(args.simulation,config['start'])
        out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
        native.to_csv(out/'native_diagnostics.csv.gz',index=False)
        bins=integrate_bins(native,config['start'],config['end'],config['interval_minutes'])
        bins['configuration']=args.configuration;bins['member']=args.member
        bins.to_csv(out/'predictors.csv',index=False)
    elif args.command=='status':
        active=ROOT/config['ensemble'].get('active_directory','runs/argentina_20181213/_execution/ensemble')
        checks=json.loads((active/'run_checks.json').read_text()) if (active/'run_checks.json').exists() else []
        stability=ROOT/'results/argentina_20181213/ensemble_stability.json'
        status={'area':config['area_selection']['status'],'flag_convention_confirmed':config['flag_mapping_confirmed'],
            'coverage_provided':config['observation_coverage'] is not None,
            'era5_files':{p.name:p.stat().st_size for p in (ROOT/'data/raw/era5').glob('*.nc')},
            'private_input_manifest':(ROOT/'data/intermediate/input_inventory.json').exists(),
            'ensemble_members':config['ensemble']['perturbed_members'],'completed_runs':len(checks),
            'expected_runs':6*(config['ensemble']['perturbed_members']+1),
            'all_finished_preflights_passed':all(r.get('physical_preflight_passed',False) for r in checks) if checks else False,
            'ensemble_stability':json.loads(stability.read_text()).get('all_available_stable') if stability.exists() else None,
            'heldout_evaluation_available':(ROOT/'results/argentina_20181213/evaluation/manifest.json').exists(),
            'pending_dependencies_file':'results/argentina_20181213/pending_dependencies.json'}
        if config.get('role')=='evaluation_only':
            status.update(role='evaluation_only',independent_calibration_event='argentina_20191028',
                independent_coefficients_available=(ROOT/config['independent_calibration']/'frozen.json').exists(),
                independent_evaluation_available=(ROOT/'results/argentina_20181213/article/independent_manifest.json').exists(),
                pending_dependencies_file='results/argentina_20181213/article/pending_dependencies.json')
        write_json(ROOT/'results/argentina_20181213/status.json',status);print(json.dumps(status,indent=2))

if __name__=='__main__': main()
