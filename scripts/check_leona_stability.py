"""Prefix ensemble stability, independent of lightning observations and errors."""
from pathlib import Path
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from leona.provenance import ROOT,write_json,sha256

def main():
    root=ROOT/'products/leona';criterion=json.loads((ROOT/'configs/ensemble_stability.json').read_text())
    config=json.loads((ROOT/'configs/leona.json').read_text())
    final_size=config['ensemble']['perturbed_members']
    sizes=[final_size//4,final_size//2,3*final_size//4,final_size]
    small_size=sizes[-2]
    files=sorted((root/'members').glob('*/*/predictors.csv'))
    total=[]
    for p in files:
        f=pd.read_csv(p);f=f[f.member!=0]
        total.append(f.groupby(['method','configuration','member']).predictor.sum(min_count=len(f.time.unique())).reset_index())
    all_values=pd.concat(total)
    rows=[];checks=[]
    for (method,configuration),f in all_values.groupby(['method','configuration']):
        f=f.sort_values('member');a=f.predictor.to_numpy()
        if not np.isfinite(a).any():
            checks.append({'method':method,'configuration':configuration,'status':'unavailable'});continue
        if len(a)!=final_size or not np.isfinite(a).all():raise ValueError('Complete paired ensemble required')
        summary={}
        for size in sizes:
            values=a[:size];q=np.quantile(values,[.05,.5,.95])
            stats={'mean':float(values.mean()),'std':float(values.std(ddof=1)),'q05':float(q[0]),'median':float(q[1]),'q95':float(q[2])}
            summary[size]=stats;rows.append({'method':method,'configuration':configuration,'members':size,**stats})
        small,large=summary[small_size],summary[final_size]
        scale=large['std']
        mean_ok=abs(small['mean']-large['mean'])<=2*scale*np.sqrt(1/small_size-1/final_size)+1e-30
        std_ok=abs(small['std']-scale)<=criterion['std_change_relative_limit']*scale+1e-30
        width=large['q95']-large['q05']
        quantile_ok=all(abs(small[k]-large[k])<=criterion['quantile_change_fraction_of_final_05_95_range']*width+1e-30 for k in ['q05','median','q95'])
        checks.append({'method':method,'configuration':configuration,'status':'degenerate_zero' if np.all(a==0) else ('stable' if mean_ok and std_ok and quantile_ok else 'increase_members'),
            'mean_ok':bool(mean_ok),'std_ok':bool(std_ok),'quantiles_ok':bool(quantile_ok)})
    pd.DataFrame(rows).to_csv(root/'tables/ensemble_convergence.csv',index=False)
    write_json(root/'ensemble_stability.json',{'criterion':criterion,'actual_prefix_sizes':sizes,'criterion_sha256':sha256(ROOT/'configs/ensemble_stability.json'),
        'all_available_stable':all(r['status']!='increase_members' for r in checks),'checks':checks})

if __name__=='__main__':main()
