from pathlib import Path
import ast,json,re,shutil,subprocess
R=Path.cwd().resolve(); log=[]
def move(a,b):
 s=(R/a).resolve();d=(R/b).resolve()
 assert s.is_relative_to(R) and d.is_relative_to(R) and s!=R and d!=R
 if not s.exists():
  if d.exists():log.append({"old":a,"new":b})
  return
 if d.exists():
  if s.is_file() and d.is_file() and s.read_bytes()==d.read_bytes():
   log.append({'old':a,'new':b,'source_retained':'open in another application'});return
  raise FileExistsError(d)
 d.parent.mkdir(parents=True,exist_ok=True);shutil.move(str(s),str(d));log.append({'old':a,'new':b})
# Extract only the process category definitions; no legacy runtime imports.
tree=ast.parse((R/'scripts/executar_processos_figura7.py').read_text(encoding='utf-8'))
category='\n\n'.join(ast.get_source_segment((R/'scripts/executar_processos_figura7.py').read_text(encoding='utf-8'),n) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in ['GROUPS','CATEGORY'] for t in n.targets))
mods={'acquire':'download','forcing':'forcing','contracts':'forcing/contracts','spatial':'forcing/spatial','lightning':'observations','proxies':'lightning','postprocess':'lightning/postprocess','ensemble':'ensemble','experiments':'ensemble/experiments','evaluate':'evaluation','figures':'plotting/basic','io':'io','provenance':'provenance','cli':'cli'}
packages={'download','forcing','observations','lightning','ensemble','evaluation'}
for old,new in mods.items():
 target='python/scm_lightning/'+new+('/__init__.py' if new in packages else '.py')
 move('leona/'+old+'.py',target)
move('leona/__init__.py','python/scm_lightning/__init__.py')
flows={'assemble_leona_era5':'forcing/assemble','select_leona_support':'forcing/select_support','finalize_leona_observations':'observations/finalize','process_leona_controls':'evaluation/process_controls','process_leona_ensemble':'evaluation/process_ensemble','figures_leona':'plotting/figures','check_leona_controls':'evaluation/check_controls','check_leona_delivery':'evaluation/check_delivery','check_leona_stability':'evaluation/check_stability','check_leona_driver':'evaluation/check_driver','archive_leona_products':'evaluation/archive_products','normalize_leona_native_time':'evaluation/normalize_time','finalize_leona_execution':'evaluation/finalize_execution','prepare_leona_products':'evaluation/prepare_products'}
for old,new in flows.items():move('scripts/'+old+'.py','python/scm_lightning/'+new+'.py')
for p in (R/'python/scm_lightning').rglob('*.py'):
 text=p.read_text(encoding='utf-8')
 # Relative imports in the old flat package become explicit absolute imports.
 text=re.sub(r'from \.(\w+) import',lambda m:'from scm_lightning.'+mods.get(m[1],m[1]).replace('/','.')+' import',text)
 text=re.sub(r'from leona\.(\w+) import',lambda m:'from scm_lightning.'+mods.get(m[1],m[1]).replace('/','.')+' import',text)
 text=text.replace("__import__('leona.forcing'","__import__('scm_lightning.forcing'")
 text=text.replace('from executar_processos_figura7 import CATEGORY','from scm_lightning.evaluation.process_categories import CATEGORY')
 text=re.sub(r"sys.path.insert\(0,\s*str\(Path\(__file__\).resolve\(\).parents\[1\]\)\)\n",'',text)
 p.write_text(text,encoding='utf-8')
(R/'python/scm_lightning/evaluation/process_categories.py').write_text(category+'\n',encoding='utf-8')
for folder in ['plotting','calibration']:
 (R/'python/scm_lightning'/folder).mkdir(exist_ok=True);(R/'python/scm_lightning'/folder/'__init__.py').touch()
for p in (R/'tests').glob('test*.py'):
 text=p.read_text(encoding='utf-8');text=re.sub(r'from leona\.(\w+) import',lambda m:'from scm_lightning.'+mods[m[1]].replace('/','.')+' import',text);p.write_text(text,encoding='utf-8')
move('tests/test_leona.py','tests/unit/test_contracts.py');move('tests/test_leona_forcing.py','tests/integration/test_forcing.py')
# Exact physical run folders move once. Immutable historical manifests retain original paths.
for short,long in [('TH','thompson'),('WS','wsm6')]:
 for mode in ['C0','CL','CF']:
  key=short+'_'+mode;dest='runs/argentina_20181213/'+long+'_'+mode
  move('work/leona/controls/'+key,dest+'/control')
  for ens in ['ensemble','ensemble40']:
   parent=R/'work/leona'/ens/key
   if parent.exists():
    for member in list(parent.iterdir()):move(member.relative_to(R).as_posix(),dest+'/members/'+member.name)
move('work/leona/build','build/leona')
move('work/leona/case_forcing_v2','data/processed/forcing')
move('work/leona','runs/argentina_20181213/_execution')
move('products/leona','results/argentina_20181213')
for a,b in [('era5','data/raw/era5'),('raw/entln','data/raw/entln'),('raw/soundings','data/raw/radiosondes'),('processed/soundings','data/processed/profiles'),('wyoming','data/raw/radiosondes/wyoming'),('forcing_source','data/intermediate/era5')]:move('data_private/leona/'+a,b)
move('data_private/leona/manifest.json','data/intermediate/input_inventory.json')
move('era5_forcing/DEPHY-SCM','external/DEPHY-SCM')
# Historical materials, including prior outputs and recipes, stay outside the active workflow.
for name in ['data_era5','data_goamazon','era5_forcing','figuras_artigo','figuras_mp_offline','resultados']:
 move(name,'archive/goamazon/'+name)
for p in list((R/'data').glob('*.nc')):move(p.relative_to(R).as_posix(),'archive/goamazon/data/'+p.name)
for name in ['andes_portillo_thompson.nml','andes_portillo_wsm6.nml','ensemble.goamazon.nml.sample','ensemble.nml.sample','Makefile_old','namelist.mp_offline','output_mp_offline.nc','plot_mp_offline.py','Projeto_artigo_microfisica_raios_Amazonia.docx','README_ANALISE_ARTIGO.md','requirements-analysis.txt','requirements-leona.txt']:
 move(name,'archive/goamazon/'+name)
for p in list((R/'scripts').glob('*.py')):move(p.relative_to(R).as_posix(),'archive/goamazon/scripts/'+p.name)
move('src/Makefile','archive/goamazon/src/Makefile')
for name in ['ensemble_driver.F90','era5_bridge.F90','output_ensemble_netcdf.F90','saturation_perturbation.F90','updraft_forcing.F90','mpas_atmphys_driver_microphysics_REFERENCE.F']:
 move('src/'+name,'archive/goamazon/src/'+name)
move('src/process_diagnostics.F90','src/diagnostics/process_diagnostics.F90')
(R/'docs').mkdir(exist_ok=True)
(R/'docs/path_migration.json').write_text(json.dumps({'baseline':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),'moves':log,'historical_manifests':'preserved byte-for-byte; use resolve_path for old references'},indent=2)+'\n')
print('Moved',len(log),'paths; no native model outputs modified')
