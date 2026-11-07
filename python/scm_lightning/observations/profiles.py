"""Station-aware sounding harmonization; campaign files take priority.

Observations are not silently transplanted into the fixed SCM column. All
raw coordinates, UTC sample times and incomplete profiles remain explicit.
"""
import argparse
from io import StringIO
from pathlib import Path
import re
import hashlib
import numpy as np
import pandas as pd
from scm_lightning.config import load_config
from scm_lightning.observations import haversine_km
from scm_lightning.provenance import ROOT,sha256,write_json

OUT=ROOT/'data/processed/radiosondes'
EPSILON=287.04/461.5


def read_text(path):
    raw=Path(path).read_bytes()
    try:return raw.decode('utf-8-sig')
    except UnicodeDecodeError:return raw.decode('cp1252')


def read_campaign(path):
    lines=read_text(path).splitlines();meta={}
    for line in lines:
        if '\t' in line:
            key,value=line.split('\t',1);meta[key.strip()]=value.strip()
    i=next(i for i,line in enumerate(lines) if line.startswith('Elapsed time'))
    columns=['elapsed_s','height_msl_m','pressure_hPa','temperature_C','dewpoint_C','RH_percent',
        'wind_speed_m_s','wind_direction_deg','balloon_ascent_m_s','latitude','longitude','gps_height_msl_m']
    f=pd.read_csv(StringIO('\n'.join(lines[i+2:])),sep=r'\s+',names=columns,na_values=['/////','////','///'],engine='python')
    for name in columns:f[name]=pd.to_numeric(f[name],errors='coerce')
    release=pd.Timestamp(meta['Balloon release date and time'],tz='UTC')
    f['time']=release+pd.to_timedelta(f.elapsed_s,unit='s')
    # GPS/geometric height -> geopotential height, without altering raw height.
    f['height_geopotential_m']=6371000*f.height_msl_m/(6371000+f.height_msl_m)
    station=meta.get('Station name','CSU_mobile')
    return finish(f,station,'CSU_RELAMPAGO_in_projetoCHUVA','GPS geometric height converted to geopotential'),meta


def read_csv(path):
    f=pd.read_csv(path)
    f=f.rename(columns={'geopotential height_m':'height_geopotential_m','dew point temperature_C':'dewpoint_C',
        'mixing ratio_g/kg':'mixing_ratio_g_kg','wind speed_m/s':'wind_speed_m_s','wind direction_degree':'wind_direction_deg','relative humidity_%':'RH_percent'})
    required=['pressure_hPa','temperature_C','height_geopotential_m']
    if not set(required)<=set(f):raise ValueError('Unsupported sounding CSV schema')
    for name in f.columns.difference(['time']):f[name]=pd.to_numeric(f[name],errors='coerce')
    if 'time' in f:f.time=pd.to_datetime(f.time,utc=True,errors='coerce')
    else:f['time']=pd.NaT
    match=re.search(r'-(\d{5})',Path(path).name)
    station=match.group(1) if match else 'unidentified'
    return finish(f,station,'Wyoming_download' if 'wyoming' in Path(path).parts else 'user_supplied_CSV','provided geopotential height'),{}


def finish(f,station,source,height_definition):
    f=f.copy()
    f['pressure_Pa']=100*f.pressure_hPa;f['temperature_K']=f.temperature_C+273.15
    if 'mixing_ratio_g_kg' in f:f['rv_kg_kg_dry']=f.mixing_ratio_g_kg/1000
    else:
        td=f.dewpoint_C.where(f.dewpoint_C.between(-150,70))
        vapor=611.2*np.exp(17.67*td/(td+243.5))
        f['rv_kg_kg_dry']=EPSILON*vapor/(f.pressure_Pa-vapor)
    f['theta_K']=f.temperature_K*(1e5/f.pressure_Pa)**(2/7)
    for name in ['wind_speed_m_s','wind_direction_deg']:
        if name not in f:f[name]=np.nan
    f['u_m_s']=-f.wind_speed_m_s*np.sin(np.deg2rad(f.wind_direction_deg))
    f['v_m_s']=-f.wind_speed_m_s*np.cos(np.deg2rad(f.wind_direction_deg))
    for name in ['latitude','longitude']:
        if name not in f:f[name]=np.nan
    valid=f.pressure_Pa.between(100,110000)&f.temperature_K.between(150,350)&f.rv_kg_kg_dry.ge(0)&f.height_geopotential_m.notna()
    f['thermodynamic_valid']=valid
    if 'height_msl_m' not in f:f['height_msl_m']=f.height_geopotential_m
    # Preserve descent in the normalized file but flag the ascent used in plots.
    maximum=f.height_geopotential_m.idxmax()
    f['ascending_segment']=f.index<=maximum
    f['station']=station;f['source']=source;f['height_definition']=height_definition
    return f


def forcing_eligibility(f,case):
    valid=f[f.thermodynamic_valid & f.ascending_segment]
    first=f.iloc[0];lat,lon=float(first.latitude),float(first.longitude)
    west,south,east,north=case['area'];inside=west<=lon<=east and south<=lat<=north
    launch=pd.Timestamp(f.time.dropna().iloc[0]) if f.time.notna().any() else None
    age=(pd.Timestamp(case['start'])-launch).total_seconds()/3600 if launch is not None else None
    coverage=bool(len(valid) and valid.height_geopotential_m.max()-first.height_geopotential_m>=20000)
    initial=bool(inside and age is not None and 0<=age<=6 and coverage)
    return {'inside_fixed_column':inside,'distance_to_column_center_km':haversine_km(lat,lon,(south+north)/2,(west+east)/2),
        'distance_to_box_km':haversine_km(lat,lon,np.clip(lat,south,north),np.clip(lon,west,east)),
        'age_at_case_start_hours':age,'covers_20km_AGL':coverage,'initial_condition_eligible':initial,
        'boundary_condition_eligible':False,'reason':('candidate IC; full pressure/terrain consistency still required' if initial else
            'outside fixed column and/or not a pre-initialization full-depth profile'),
        'boundary_reason':'a single flight is not complete time-varying forcing or horizontal-advection information'}


def prepare():
    OUT.mkdir(parents=True,exist_ok=True)
    paths=sorted((ROOT/'data/raw/projetoCHUVA').glob('edt_*.txt'))+sorted((ROOT/'data/raw/radiosondes').glob('*.csv'))+sorted((ROOT/'data/raw/radiosondes/wyoming/profiles').glob('*.csv'))
    rows=[];seen={}
    for path in paths:
        try:
            f,rawmeta=read_campaign(path) if path.suffix=='.txt' else read_csv(path)
            canonical=f[['time','pressure_Pa','temperature_K','rv_kg_kg_dry','latitude','longitude']].to_csv(index=False)
            identity=hashlib.sha256(canonical.encode()).hexdigest()
            if identity in seen:
                rows.append({'file':str(path),'status':'duplicate_profile','duplicate_of':seen[identity],'sha256':sha256(path)})
                continue
            seen[identity]=str(path)
            first=f.iloc[0];launch=pd.Timestamp(f.time.dropna().iloc[0]) if f.time.notna().any() else None
            stem=path.stem;target=OUT/(stem+'.csv.gz');f.to_csv(target,index=False)
            row={'file':str(path),'normalized':str(target),'sha256':sha256(path),'normalized_sha256':sha256(target),
                'status':'available','source':first.source,'station':str(first.station),'launch_UTC':str(launch) if launch else None,
                'latitude':float(first.latitude) if pd.notna(first.latitude) else None,
                'longitude':float(first.longitude) if pd.notna(first.longitude) else None,'rows':len(f),
                'valid_ascent_rows':int((f.thermodynamic_valid&f.ascending_segment).sum()),
                'top_geopotential_m':float(f.height_geopotential_m.max()),'min_pressure_hPa':float(f.pressure_hPa.min()),
                'height_definition':first.height_definition,'actual_per_level_times':bool(f.time.nunique()>1),
                'balloon_ascent_is_atmospheric_w':False,'raw_metadata':rawmeta,'cases':{}}
            for event in ['argentina_20181213','argentina_20191028']:
                c=load_config(ROOT/'configs/cases'/(event+'.yaml'))
                if launch is not None and pd.Timestamp(c['start'])-pd.Timedelta(hours=6)<=launch<pd.Timestamp(c['end']):
                    row['cases'][event]=forcing_eligibility(f,c)
            write_json(target.with_suffix('.json'),row);rows.append(row)
        except Exception as exc:rows.append({'file':str(path),'status':'failed','reason':str(exc)})
    write_json(OUT/'inventory.json',{'source_priority':['CSU_RELAMPAGO_in_projetoCHUVA','user_supplied_CSV','Wyoming_download'],
        'station_profiles_never_merged':True,'source_time_policy':'actual sample UTC; never infer launch from nominal filename',
        'profiles':rows,'existing_SCM_initial_and_boundary_conditions':'ERA5 only; observations have not been assimilated',
        'initial_condition_policy':'inside fixed box, launch no later than initial time, age <=6h, >=20km AGL coverage',
        'boundary_policy':'same-location observations with complete temporal coverage and explicitly specified nudging; no automatic transplantation'})
    table=[]
    for row in rows:
        for event,decision in row.get('cases',{}).items():
            table.append({k:row[k] for k in ['file','normalized','source','station','launch_UTC','latitude','longitude','top_geopotential_m','min_pressure_hPa']}|{'event':event}|decision)
    pd.DataFrame(table).to_csv(OUT/'case_availability_and_forcing_eligibility.csv',index=False)


def selected(event):
    import json
    inventory=json.loads((OUT/'inventory.json').read_text())
    eligible=[r for r in inventory['profiles'] if r.get('status')=='available' and event in r.get('cases',{})]
    if event=='argentina_20181213' and any(r['source']=='CSU_RELAMPAGO_in_projetoCHUVA' for r in eligible):
        return [r for r in eligible if r['source']=='CSU_RELAMPAGO_in_projetoCHUVA']
    if event=='argentina_20191028':return [r for r in eligible if r['source']=='user_supplied_CSV']
    return eligible


def compare_era5():
    """Compare preferred observed profiles with ERA5 at their own locations/times.

    This is not a comparison of distant soundings against the fixed SCM column.
    Statistics use a common log-pressure grid, avoiding high-resolution sampling
    receiving more weight than mandatory/significant-level profiles.
    """
    import xarray as xr
    from scm_lightning.forcing import normalize
    rows=[]
    for event in ['argentina_20181213','argentina_20191028']:
        c=load_config(ROOT/'configs/cases'/(event+'.yaml'))
        cache=ROOT/c.get('paths',{}).get('era5','data/raw/era5')
        for record in selected(event):
            f=pd.read_csv(record['normalized']);f=f[f.thermodynamic_valid&f.ascending_segment].copy()
            f.time=pd.to_datetime(f.time,utc=True)
            # Remove duplicate pressure levels only in the common-grid derivative.
            source=f.groupby('pressure_Pa',as_index=False).agg(temperature_K=('temperature_K','mean'),
                rv_kg_kg_dry=('rv_kg_kg_dry','mean'),u_m_s=('u_m_s','mean'),v_m_s=('v_m_s','mean'),
                latitude=('latitude','mean'),longitude=('longitude','mean'),height_geopotential_m=('height_geopotential_m','mean'),
                elapsed=('time',lambda t:float((t.mean()-pd.Timestamp(record['launch_UTC'])).total_seconds()))).sort_values('pressure_Pa')
            pgrid=100*np.array([50,70,100,125,150,175,200,225,250,300,350,400,450,500,550,600,650,700,750,775,800,825,850,875,900,925,950,975,1000])
            pgrid=pgrid[(pgrid>=source.pressure_Pa.min())&(pgrid<=source.pressure_Pa.max())]
            if len(pgrid)<2:continue
            q=pd.DataFrame({'pressure_Pa':pgrid})
            for name in source.columns.difference(['pressure_Pa']):
                q[name]=np.interp(np.log(pgrid),np.log(source.pressure_Pa),source[name])
            clock=pd.Timestamp(record['launch_UTC'])+pd.to_timedelta(q.elapsed,unit='s')
            days=pd.date_range(clock.min().floor('D'),clock.max().ceil('D'),freq='D')
            parts=[];surface=[];inputs=[]
            for day in days:
                for kind,holder in [('pl',parts),('sl',surface)]:
                    path=cache/f'reanalysis_{kind}_{day:%Y%m%d}.nc'
                    if not path.exists():continue
                    with xr.open_dataset(path) as original:
                        d=normalize(original).sel(latitude=slice(q.latitude.min()-.5,q.latitude.max()+.5),
                            longitude=slice(q.longitude.min()-.5,q.longitude.max()+.5))
                        holder.append(d[['t','q','u','v'] if kind=='pl' else ['sp']].load())
                    inputs.append({'path':str(path),'sha256':sha256(path)})
            if not parts or not surface:continue
            coords={'time':xr.DataArray(clock.dt.tz_localize(None).to_numpy(),dims='sample'),
                'latitude':xr.DataArray(q.latitude.values,dims='sample'),'longitude':xr.DataArray(q.longitude.values,dims='sample')}
            pl=xr.concat(parts,dim='time').sortby('time').interp(**coords,level=xr.DataArray(pgrid/100,dims='sample'))
            sl=xr.concat(surface,dim='time').sortby('time').interp(**coords)
            above=pgrid<sl.sp.values
            for obs,model in [('temperature_K',pl.t.values),('rv_kg_kg_dry',(pl.q/(1-pl.q)).values),('u_m_s',pl.u.values),('v_m_s',pl.v.values)]:
                q['ERA5_'+obs]=np.where(above,model,np.nan)
                q['obs_minus_ERA5_'+obs]=q[obs]-q['ERA5_'+obs]
            q['time_UTC']=clock.astype(str);q['station']=record['station'];q['event']=event
            target=OUT/(Path(record['file']).stem+'_ERA5_common_pressure.csv');q.to_csv(target,index=False)
            summary={'event':event,'station':record['station'],'launch_UTC':record['launch_UTC'],
                'source':record['source'],'profile_file':record['normalized'],'comparison_file':str(target),
                'levels':len(q),'support':'own station and observed ascent time/trajectory; common log-pressure grid',
                'comparison_is_at_fixed_SCM_column':False}
            for name in ['temperature_K','rv_kg_kg_dry','u_m_s','v_m_s']:
                e=q['obs_minus_ERA5_'+name].dropna()
                summary['n_'+name]=len(e);summary['RMSE_'+name]=float(np.sqrt(np.mean(e**2))) if len(e) else None
            write_json(target.with_suffix('.json'),{'profile_source_sha256':record['sha256'],'ERA5_inputs':inputs,
                'pressure_grid_Pa':pgrid.tolist(),'horizontal_support':'each sounding own location, not SCM box',
                'time_support':'observed flight time at each pressure; repeated reported time retained where flight timing is unavailable',
                'extrapolation':False,'terrain':'pressure below interpolated ERA5 surface is unavailable',
                'mixing_ratio_basis':'kg/kg dry air; ERA5 q converted using q/(1-q)'})
            rows.append(summary)
    pd.DataFrame(rows).to_csv(OUT/'ERA5_profile_comparison_summary.csv',index=False)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--compare-era5',action='store_true');a=ap.parse_args();prepare()
    if a.compare_era5:compare_era5()


if __name__=='__main__':main()
