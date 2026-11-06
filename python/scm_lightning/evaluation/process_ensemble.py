from scm_lightning.config import load_config
from scm_lightning.provenance import resolve_path
"""Postprocess finished members; evaluate only after the complete paired set exists."""
from pathlib import Path
import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor,as_completed
import pandas as pd
from scm_lightning.provenance import ROOT,write_json,sha256
from scm_lightning.lightning.postprocess import native_diagnostics,integrate_bins
from scm_lightning.evaluation import evaluate,ensemble_statistics

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--evaluate',action='store_true');ap.add_argument('--evaluation-output');args=ap.parse_args()
    c=load_config(ROOT/'configs/cases/argentina_20181213.yaml')
    if args.evaluate and c.get('role')=='evaluation_only':
        raise ValueError('2018 observations cannot fit coefficients; use independent evaluate')
    ensemble=ROOT/c['ensemble'].get('active_directory','runs/argentina_20181213/_execution/ensemble')
    plan=json.loads((ensemble/'manifest.json').read_text())
    for attempt in range(5):
        try:
            checks=json.loads((ensemble/'run_checks.json').read_text());break
        except json.JSONDecodeError:
            if attempt==4:raise
            time.sleep(.2)
    out=ROOT/'results/argentina_20181213';members=out/'members';members.mkdir(parents=True,exist_ok=True)
    all_predictions=[]
    def one(run):
        if not run.get('physical_preflight_passed',False):
            raise ValueError(f'Member fails physical preflight: {run["configuration"]}/{run["member"]}')
        folder=members/run['configuration']/f'{run["member"]:02}';folder.mkdir(parents=True,exist_ok=True)
        rawfolder=resolve_path(run['folder']);rawfolder=rawfolder if rawfolder.is_absolute() else ROOT/rawfolder
        source=rawfolder/'output.nc';target=folder/'predictors.csv';tag=folder/'manifest.json'
        if target.exists():
            if not tag.exists() or json.loads(tag.read_text())['simulation_sha256']!=run['sha256']:
                raise ValueError('Existing member products use a different simulation')
        else:
            if sha256(source)!=run['sha256']:raise ValueError('Changed finished model output')
            native=native_diagnostics(source,c['start'])
            native.to_csv(folder/'native_diagnostics.csv.gz',index=False)
            pred=integrate_bins(native,c['start'],c['end'],c['interval_minutes'])
            pred['configuration']=run['configuration'];pred['member']=run['member']
            pred.to_csv(target,index=False)
            write_json(tag,{'simulation_sha256':run['sha256'],'configuration':run['configuration'],
                'member':run['member'],'origin_utc':c['start'],'predictor_units':'published raw proxy units integrated over seconds; scale fits interval counts',
                'density_basis':'dry-air density; MPAS Rd=287, Rv=461.6, kappa=2/7',
                'diagnostics_mutate_simulation':False})
        return pd.read_csv(target)
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(one,run) for run in checks]
        for future in as_completed(futures):
            all_predictions.append(future.result())
            write_json(out/'postprocess_progress.json',{'finished_members':len(checks),'expected_members':len(plan['runs']),
                'postprocessed_members':len(all_predictions)})
    write_json(out/'postprocess_progress.json',{'finished_members':len(checks),'expected_members':len(plan['runs']),
        'postprocessed_members':len(all_predictions)})
    if args.evaluate:
        if len(checks)!=len(plan['runs']):raise ValueError('Wait for complete paired ensemble')
        p=pd.concat(all_predictions,ignore_index=True)
        o=pd.read_csv(out/'observations/counts_15min.csv')
        p.to_csv(out/'predictors_15min.csv.gz',index=False)
        result=evaluate(p,o,c['reference_configuration'],c['split']['development_end'],c['split']['evaluation_start'])
        dest=ROOT/args.evaluation_output if args.evaluation_output else out/'evaluation';dest.mkdir(parents=True,exist_ok=False)
        for name,value in zip(['coefficients','heldout_predictions','metrics','fixed_scale_skill_contrasts','physical_trajectories'],result):
            value.to_csv(dest/(name+'.csv'),index=False)
        ensemble_statistics(result[1]).to_csv(dest/'ensemble_statistics.csv',index=False)
        write_json(dest/'manifest.json',{'config':c,'ensemble_manifest_sha256':sha256(ensemble/'manifest.json'),
            'observations_sha256':sha256(out/'observations/counts_15min.csv'),
            'split_semantics':'all operations fitted on development only; same contiguous split for every method/configuration/member',
            'spread_scope':'initial condition only; not total predictive uncertainty'})

if __name__=='__main__':main()
