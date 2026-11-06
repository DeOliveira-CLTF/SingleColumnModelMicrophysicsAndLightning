"""Group all event types jointly before final spatial and temporal selection."""
from collections import deque,defaultdict
from itertools import product
from pathlib import Path
import hashlib
import math
import numpy as np
import pandas as pd
EVENT_COLUMNS = ['time', 'latitude', 'longitude', 'current_kA', 'kind']

def read_events(paths, flag_mapping):
    if set(flag_mapping.values()) != {'IC','CG'}:
        raise ValueError('Explicit export flag convention required')
    frames = []
    for path in paths:
        f = pd.read_csv(path)
        t = pd.to_datetime(f[['year','month','day','hour','minute','second']], utc=True)
        if not (f.nanosec.between(0,999999999) & (f.nanosec%1==0)).all():
            raise ValueError('Invalid nanosecond')
        f['time'] = t + pd.to_timedelta(f.nanosec.astype('int64'), unit='ns')
        f['kind'] = f.ic_cg_flag.astype(str).map(flag_mapping)
        if f.kind.isna().any():
            raise ValueError('Unknown IC/CG flag')
        f = f.rename(columns={'corrente_pico_kA':'current_kA'})
        f['source'], f['source_row'] = Path(path).name, np.arange(len(f)) + 2
        if not (f.latitude.between(-90,90) & f.longitude.between(-180,180)).all():
            raise ValueError('Invalid coordinates')
        frames.append(f[EVENT_COLUMNS + ['source','source_row']])
    events = pd.concat(frames, ignore_index=True).sort_values(EVENT_COLUMNS + ['source','source_row'], kind='stable')
    duplicates = events[events.duplicated(EVENT_COLUMNS, keep='first')].copy()
    events = events.drop_duplicates(EVENT_COLUMNS).reset_index(drop=True)
    events['event_id'] = [hashlib.sha256('|'.join(map(str,r)).encode()).hexdigest()[:24]
                         for r in events[EVENT_COLUMNS].itertuples(index=False, name=None)]
    return events, duplicates

def haversine_km(lat1, lon1, lat2, lon2):
    a,b = math.radians(lat1),math.radians(lat2)
    h = math.sin((b-a)/2)**2 + math.cos(a)*math.cos(b)*math.sin(math.radians(lon2-lon1)/2)**2
    return 12742 * math.asin(math.sqrt(min(1.,max(0.,h))))

def group_flashes(events, version='chronis2015_first_v1', progress=None):
    if version not in ('chronis2015_first_v1','azambuja2017_v1'):
        raise ValueError('Unknown grouping version')
    f = events.sort_values(['time','event_id'], kind='stable').copy()
    active, flashes, assigned = deque(), [], []
    buckets=defaultdict(set)
    offsets=list(product([-1,0,1],repeat=3))
    def cell(lat,lon):
        a,b=math.radians(lat),math.radians(lon)
        # A 10-km geocentric Cartesian cell, valid even near the dateline/poles.
        return tuple(math.floor(v) for v in (637.1*math.cos(a)*math.cos(b),637.1*math.cos(a)*math.sin(b),637.1*math.sin(a)))
    limit = 700_000_000 if version.startswith('chronis') else 1_400_000_000
    for number,row in enumerate(f.itertuples(),1):
        ns = row.time.value
        while active and ns - flashes[active[0]]['start_ns'] > limit:
            old=active.popleft();g=flashes[old]
            key=cell(g['latitude'],g['longitude']);buckets[key].discard(old)
            if not buckets[key]:del buckets[key]
        selected = None
        key=cell(row.latitude,row.longitude)
        candidates=set()
        for delta in offsets:
            candidates.update(buckets.get(tuple(a+b for a,b in zip(key,delta)),()))
        # Earliest eligible anchor wins. No polarity separation or neighbor chaining.
        for i in sorted(candidates):
            g = flashes[i]
            if version.startswith('azambuja') and ns - g['last_ns'] > 700_000_000:
                continue
            if haversine_km(g['latitude'],g['longitude'],row.latitude,row.longitude) <= 10 + 1e-10:
                selected = i
                break
        if selected is None:
            selected = len(flashes)
            flashes.append({'flash_id':version+':'+row.event_id, 'start_ns':ns,'last_ns':ns,
                'latitude':row.latitude,'longitude':row.longitude,'multiplicity':0,'kind':'IC','version':version})
            active.append(selected)
            buckets[key].add(selected)
        g = flashes[selected]
        g['last_ns'], g['multiplicity'] = ns, g['multiplicity']+1
        if row.kind == 'CG':
            g['kind'] = 'CG'
        assigned.append(g['flash_id'])
        if progress is not None and number%100000==0:
            progress({'processed_events':number,'total_events':len(f),'flashes_so_far':len(flashes)})
    f['flash_id'] = assigned
    out = pd.DataFrame(flashes, columns=['flash_id','start_ns','last_ns','latitude','longitude','multiplicity','kind','version'])
    out['start'] = pd.to_datetime(out.pop('start_ns'), utc=True, unit='ns')
    out['end'] = pd.to_datetime(out.pop('last_ns'), utc=True, unit='ns')
    return out, f

def aggregate(flashes, start, end, box, covered, minutes=15):
    from scm_lightning.forcing.contracts import area_km2
    edges = pd.date_range(start,end,freq=f'{minutes}min',inclusive='left')
    if len(covered) != len(edges):
        raise ValueError('Coverage must match bins')
    w,s,e,n = box
    f = flashes[flashes.longitude.between(w,e) & flashes.latitude.between(s,n)]
    t = f.start[(f.start >= pd.Timestamp(start)) & (f.start < pd.Timestamp(end))]
    counts = t.dt.floor(f'{minutes}min').value_counts().reindex(edges,fill_value=0).astype(float)
    counts.loc[~np.asarray(covered,dtype=bool)] = np.nan
    result = pd.DataFrame({'time':edges,'covered':covered,'count':counts.values})
    result['rate_per_min'] = result['count']/minutes
    result['density_km2'] = result['count']/area_km2(box)
    result['cumulative'] = result['count'].cumsum(skipna=False)
    return result
