from pathlib import Path
import json,re
R=Path.cwd();P=R/'python/scm_lightning';C='argentina_20181213'
replacements={'products/leona':'results/'+C,'data_private/leona/era5':'data/raw/era5','data_private/leona/raw/entln':'data/raw/entln','data_private/leona/raw/soundings':'data/raw/radiosondes','data_private/leona/processed/soundings':'data/processed/profiles','data_private/leona/wyoming':'data/raw/radiosondes/wyoming','data_private/leona/forcing_source':'data/intermediate/era5','data_private/leona':'data/intermediate/input_cache','work/leona/case_forcing_v2':'data/processed/forcing','work/leona/build':'build/leona','work/leona':'runs/'+C+'/_execution','era5_forcing/DEPHY-SCM':'external/DEPHY-SCM','configs/leona.json':'configs/cases/'+C+'.yaml','configs/parameterizations.json':'configs/lightning/schemes.yaml','configs/ensemble_stability.json':'configs/ensemble/stability.yaml','configs/budget_contract.json':'configs/evaluation/budget_contract.yaml','configs/sounding_sources.json':'configs/cases/sounding_sources.yaml'}
for p in list(P.rglob('*.py'))+list((R/'tests').rglob('*.py')):
 s=p.read_text(encoding='utf-8')
 for a,b in replacements.items():s=s.replace(a,b)
 p.write_text(s,encoding='utf-8')
# JSON is a strict subset of YAML 1.2: retain the existing validated parser.
for old,new in [('leona','cases/'+C),('parameterizations','lightning/schemes'),('ensemble_stability','ensemble/stability'),('budget_contract','evaluation/budget_contract'),('sounding_sources','cases/sounding_sources')]:
 src=R/'configs'/f'{old}.json';data=json.loads(src.read_text())
 if old=='leona':
  data['ensemble']['active_directory']='runs/'+C+'/_execution/ensemble40'
  data['station']={'id':'87344','icao':'SACO','latitude':-31.32,'longitude':-64.22,'source':'https://www.data.jma.go.jp/tcc/tcc/products/climate/extlist_weekly/extlist_20261014.html','launch_coordinate_policy':'CSV launch locations plotted separately from station metadata'}
  data['plot_extent_requested']=data.pop('plot_extent_case');data['plot_extent_policy']='intersection_with_finite_ERA5_coverage_no_extrapolation'
 dest=R/'configs'/f'{new}.yaml';dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf-8');src.unlink()
# ROOT is now under python/scm_lightning.
p=P/'provenance.py';s=p.read_text().replace('parents[1]','parents[2]');p.write_text(s)
# Historical manifest references are resolved without rewriting their scientific provenance.
with (P/'provenance.py').open('a',encoding='utf-8') as f:f.write('''\n\ndef resolve_path(value):
    """Resolve immutable pre-reorganization manifest paths into the current layout."""
    value = str(value).replace('\\\\', '/')
    root = ROOT.as_posix()
    if value.startswith(root + '/'):
        value = value[len(root)+1:]
    elif value.startswith('/mnt/c/git_repos/SingleColumnModelMicrophysicsAndLightning/'):
        value = value.split('SingleColumnModelMicrophysicsAndLightning/', 1)[1]
    for short, long in [('TH', 'thompson'), ('WS', 'wsm6')]:
        for mode in ['C0', 'CL', 'CF']:
            cfg = short+'_'+mode
            dest = 'runs/argentina_20181213/'+long+'_'+mode
            for old in ['work/leona/ensemble/', 'work/leona/ensemble40/']:
                prefix = old+cfg+'/'
                if value.startswith(prefix):
                    return ROOT/dest/'members'/value[len(prefix):]
            prefix = 'work/leona/controls/'+cfg
            if value == prefix or value.startswith(prefix+'/'):
                return ROOT/(dest+'/control'+value[len(prefix):])
    mapping = '''+repr(replacements)+'''
    for old, new in sorted(mapping.items(), key=lambda item: -len(item[0])):
        if value == old or value.startswith(old+'/'):
            return ROOT/(new+value[len(old):])
    return Path(value) if Path(value).is_absolute() else ROOT/value
''')
# Resolve path fields where old immutable manifests are consumed.
for p in P.rglob('*.py'):
 if p.name=='provenance.py':continue
 s=p.read_text(encoding='utf-8')
 if "Path(run['folder'])" in s or "Path(run['forcing'])" in s or "Path(r['forcing'])" in s:
  s='from scm_lightning.provenance import resolve_path\n'+s
  s=s.replace("Path(run['folder'])","resolve_path(run['folder'])").replace("Path(run['forcing'])","resolve_path(run['forcing'])").replace("Path(r['forcing'])","resolve_path(r['forcing'])")
 p.write_text(s,encoding='utf-8')
# Make builds just the active SCM. Prototype plume core remains a real MYNN dependency.
p=R/'Makefile';s=p.read_text();s=s.replace('all: $(BUILD)/mp_offline_driver $(BUILD)/ensemble_driver','all: $(BUILD)/mp_offline_driver')
s=re.sub(r'ENSEMBLE_DRIVER_OBJS :=.*?(?=\.PHONY)', '',s,flags=re.S)
s=re.sub(r'\$\(BUILD\)/ensemble_driver:.*?(?=# ---- compile)', '',s,flags=re.S)
s=s.replace('$(BUILD)/%.o: $(SRC)/%.F90 | $(BUILD)','$(BUILD)/%.o: $(SRC)/diagnostics/%.F90 | $(BUILD)\n\t$(FC) $(FFLAGS) -c $< -o $@\n$(BUILD)/%.o: $(SRC)/%.F90 | $(BUILD)')
s=re.sub(r'^\$\(BUILD\)/(?:ensemble_driver|era5_bridge|saturation_perturbation|updraft_forcing|output_ensemble_netcdf)\.o:.*\n','',s,flags=re.M)
s='# Active Argentina SCM; historical ensemble recipes: archive/goamazon.\n'+s[s.index('FC      :='):]
s+='\n.PHONY: mp_offline_driver test reproduce\nmp_offline_driver: $(BUILD)/mp_offline_driver\ntest:\n\tpython -m pytest tests -q\nreproduce:\n\tpython scripts/reproduce.py --help\n'
p.write_text(s)
(R/'.gitmodules').write_text('[submodule "DEPHY-SCM"]\n\tpath = external/DEPHY-SCM\n\turl = https://github.com/GdR-DEPHY/DEPHY-SCM.git\n')
(R/'pyproject.toml').write_text('''[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "scm-lightning"
version = "0.2.0"
description = "Argentina 13-14 December 2018 single-column microphysics and lightning"
requires-python = ">=3.10"
dependencies = ["numpy>=1.24", "pandas>=2.2", "xarray>=2024.1", "netCDF4>=1.6", "scipy>=1.10", "matplotlib>=3.8", "metpy>=1.5", "requests>=2.31", "cdsapi>=0.7"]

[project.optional-dependencies]
test = ["pytest>=8"]

[project.scripts]
scm-lightning = "scm_lightning.cli:main"

[tool.setuptools.packages.find]
where = ["python"]

[tool.pytest.ini_options]
pythonpath = ["python"]
testpaths = ["tests"]
''')
(R/'environment.yml').write_text('name: scm-lightning\nchannels: [conda-forge]\ndependencies:\n  - python=3.12\n  - numpy\n  - pandas>=2.2\n  - xarray\n  - netcdf4\n  - scipy\n  - matplotlib\n  - metpy>=1.5\n  - requests\n  - cdsapi\n  - pytest\n  - pip\n  - pip:\n      - -e .\n# Build Fortran in Linux/WSL with gfortran and libnetcdff-dev.\n')
(R/'.gitignore').write_text('''# Local/private inputs, model runs and generated scientific products.
/data/raw/
/data/intermediate/
/data/processed/
/runs/
/results/
/build/
*.egg-info/
__pycache__/
.pytest_cache/
*.pyc
''')
(R/'python/scm_lightning/calibration/__init__.py').write_text('from scm_lightning.evaluation import fit_scale\n__all__ = ["fit_scale"]\n')
for kind in ['regression','fixtures']:(R/'tests'/kind).mkdir(exist_ok=True)
# Keep the migration recipe with technical history rather than as a ninth entry point.
(R/'archive/migrations').mkdir(parents=True,exist_ok=True)
(R/'.restructure.py').rename(R/'archive/migrations/restructure_20181213.py')
print('Configured package and manifest path resolver')
