"""Archive native products without moving or replacing original run files."""
from pathlib import Path
import json
import os
import shutil
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from leona.experiments import CONFIGURATIONS
from leona.provenance import ROOT,write_json,sha256

def main():
    out=ROOT/'products/leona';c=json.loads((ROOT/'configs/leona.json').read_text())
    ensemble=ROOT/c['ensemble'].get('active_directory','work/leona/ensemble')
    plan=json.loads((ensemble/'manifest.json').read_text());checks=json.loads((ensemble/'run_checks.json').read_text())
    by_run={(r['configuration'],r['member']):r for r in checks};records=[]
    def archive(source,target,digest=None):
        source=Path(source).resolve();target=Path(target).resolve()
        assert source.is_relative_to(ROOT) and target.is_relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():
            if target.stat().st_ino!=source.stat().st_ino and sha256(target)!=sha256(source):raise FileExistsError(target)
            method='existing_verified'
        else:
            try:os.link(source,target);method='hardlink'
            except OSError:shutil.copy2(source,target);method='copy'
        expected=digest or sha256(source)
        if method=='copy' and sha256(target)!=expected:raise IOError('Native archive checksum mismatch')
        records.append({'source':str(source.relative_to(ROOT)),'product':str(target.relative_to(ROOT)),
            'bytes':target.stat().st_size,'method':method,'sha256':expected})
    for run in plan['runs']:
        folder=Path(run['folder']);folder=folder if folder.is_absolute() else ROOT/folder
        target=out/'simulations'/run['configuration']/f'{run["member"]:02}'
        for name in ['output.nc','output.nc.water_budget.csv','namelist.mp_offline','execution.log','run_identity.json']:
            archive(folder/name,target/name,by_run[(run['configuration'],run['member'])]['sha256'] if name=='output.nc' else None)
    for configuration in CONFIGURATIONS:
        source=ROOT/'work/leona/controls'/configuration
        for p in source.iterdir():
            if p.is_file():archive(p,out/'controls'/configuration/p.name)
    unique_inputs={r['member']:Path(r['forcing']) for r in plan['runs']}
    for member,p in unique_inputs.items():archive(p,out/'ensemble_inputs'/f'member_{member:02}.nc')
    for p in (ROOT/'work/leona/case_forcing_v2').iterdir():
        if p.is_file():archive(p,out/'forcing'/p.name)
    archive(ensemble/'realized_perturbations.npz',out/f'ensemble_realizations_{c["ensemble"]["perturbed_members"]}.npz')
    archive(ROOT/'work/leona/ensemble/realized_perturbations.npz',out/'ensemble_realizations_20.npz')
    archive(ROOT/'work/leona/eda_profiles.npz',out/'eda_profiles.npz')
    write_json(out/'native_product_index.json',{'immutability':'Treat native files as immutable. Hardlinks share bytes with preserved run files; never edit either path in place.',
        'restricted':'Local only; no observations or derivatives published','files':records})

if __name__=='__main__':main()
