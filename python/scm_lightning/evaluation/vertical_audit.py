"""Run read-only vertical instrumentation and compare every saved native variable."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import numpy as np
import xarray as xr
from scm_lightning.download import NETCDF_LOCK
from scm_lightning.ensemble.experiments import CONFIGURATIONS, execute, namelist
from scm_lightning.provenance import ROOT, sha256, write_json


def configuration_path(configuration):
    return configuration.replace('TH_', 'thompson_').replace('WS_', 'wsm6_')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', default='build/vertical_audit/mp_offline_driver')
    parser.add_argument('--output', default='runs/argentina_20181213/_diagnostic_audit_v2')
    args = parser.parse_args()
    executable = (ROOT/args.executable).resolve()
    root = ROOT/args.output
    forcing = ROOT/'data/processed/forcing/LEONA_SCM.nc'

    def one(configuration):
        folder = root/configuration
        folder.mkdir(parents=True, exist_ok=True)
        output = folder/'output.nc'
        text = namelist(configuration, forcing, output).replace(
            'write_water_budget=.true.,',
            'write_water_budget=.false.,write_vertical_diagnostics=.true.,')
        nml = folder/'namelist.mp_offline'
        if output.exists():
            identity = json.loads((folder/'audit_identity.json').read_text())
            if identity != {'executable': sha256(executable), 'forcing': sha256(forcing), 'namelist': text}:
                raise ValueError('Changed diagnostic audit inputs; use a fresh output folder')
        else:
            nml.write_text(text, encoding='utf-8')
            write_json(folder/'audit_identity.json', {
                'executable': sha256(executable), 'forcing': sha256(forcing), 'namelist': text})
            execute(executable, folder)
        baseline = ROOT/'runs/argentina_20181213'/configuration_path(configuration)/'control/output.nc'
        with NETCDF_LOCK:
            with xr.open_dataset(baseline, decode_times=False) as original, xr.open_dataset(output, decode_times=False) as actual:
                if set(original.variables) != set(actual.variables):
                    raise AssertionError('Native output schema changed')
                for name in original.variables:
                    np.testing.assert_array_equal(actual[name].values, original[name].values, err_msg=configuration+'/'+name)
                names = sorted(original.variables)
            with xr.open_dataset(str(output)+'.vertical.nc', decode_times=False) as diagnostics:
                assert diagnostics.attrs['mynn_edmf_type'] == 0
                assert diagnostics.attrs['mynn_uedmf_mp'] == 0
                assert np.all(diagnostics.graupel_plume_flux_available.values == 0)
                assert diagnostics.sizes['time'] == 5760
        return {'configuration': configuration, 'native_variables_bit_identical': True,
                'variables': names, 'baseline_sha256': sha256(baseline), 'instrumented_sha256': sha256(output),
                'vertical_file': str(output.relative_to(ROOT))+'.vertical.nc',
                'vertical_sha256': sha256(str(output)+'.vertical.nc'),
                'mynn_edmf_type': 0, 'mynn_uedmf_mp': 0,
                'graupel_plume_flux': 'unavailable: unified explicit plume microphysics is disabled'}

    with ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(one, CONFIGURATIONS))
    write_json(ROOT/'results/argentina_20181213/tests/vertical_instrumentation.json', {
        'scope': 'six full 48-hour controls compared against saved pre-instrumentation runs; all native variables and timesteps',
        'all_physical_outputs_bit_identical': all(r['native_variables_bit_identical'] for r in results),
        'runs': results})
    print('Six full controls: every native variable bit-identical to the saved baseline.')


if __name__ == '__main__':
    main()
