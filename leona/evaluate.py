"""Contiguous temporal holdout, nonnegative scale and fixed-reference contrasts."""
import numpy as np
import pandas as pd

def split_times(times, development_end, evaluation_start):
    t = pd.DatetimeIndex(times)
    a,b = pd.Timestamp(development_end),pd.Timestamp(evaluation_start)
    if b <= a:
        raise ValueError('Positive temporal separation required')
    return np.where(t < a,'development',np.where(t >= b,'evaluation','gap'))

def fit_scale(x,y):
    x,y = np.asarray(x,float),np.asarray(y,float)
    ok = np.isfinite(x)&np.isfinite(y)
    if not ok.any():
        return np.nan,'no_valid_development_data'
    if np.any(x[ok] < 0) or np.any(y[ok] < 0):
        raise ValueError('Counts/predictors must be nonnegative')
    norm = x[ok]@x[ok]
    if norm == 0:
        return np.nan,'identically_zero_unrecoverable'
    return max(0.,float(x[ok]@y[ok]/norm)),'fitted'

def metrics(y,p):
    y,p = np.asarray(y,float),np.asarray(p,float)
    ok = np.isfinite(y)&np.isfinite(p)
    y,p = y[ok],p[ok]
    if not len(y):
        return {'n':0,'rmse':np.nan,'mae':np.nan,'correlation':np.nan,'peak_count_error':np.nan,'peak_time_error_bins':np.nan,'accumulated_error':np.nan}
    return {'n':len(y),'rmse':float(np.sqrt(np.mean((p-y)**2))),'mae':float(np.mean(abs(p-y))),
        'correlation':float(np.corrcoef(y,p)[0,1]) if len(y)>1 and y.std()>0 and p.std()>0 else np.nan,
        'peak_count_error':float(p.max()-y.max()),
        'peak_time_error_bins':int(np.flatnonzero(ok)[p.argmax()]-np.flatnonzero(ok)[y.argmax()]),
        'accumulated_error':float(p.sum()-y.sum())}

def evaluate(predictors, observations, reference, development_end, evaluation_start):
    """Input predictor is COUNT-support per bin, member; fit ensemble mean only.

    columns: time,method,configuration,member,predictor. Observation: time,count.
    All members at a timestamp share a split. Same reference coefficient in B.
    """
    keys = ['time','method','configuration','member']
    if predictors.duplicated(keys).any() or observations.duplicated('time').any():
        raise ValueError('Duplicate prediction/observation keys')
    members=predictors.groupby(['time','method','configuration']).member.agg(lambda x:tuple(sorted(x)))
    if members.nunique()>1:
        raise ValueError('Member support differs between times/configurations/methods')
    p = predictors.copy()
    p['time'] = pd.to_datetime(p.time,utc=True)
    o = observations.copy()
    o['time'] = pd.to_datetime(o.time,utc=True)
    p = p.merge(o[['time','count']],on='time',validate='many_to_one')
    p['split'] = split_times(p.time,development_end,evaluation_start)
    # Member 0 is the unperturbed control. Fit on mean of synthetic members if present.
    ensemble = p[p.member != 0] if (p.member != 0).any() else p
    mean = ensemble.groupby(['time','method','configuration','split'],as_index=False).agg(predictor=('predictor','mean'),count=('count','first'))
    coefficients = []
    for (method,configuration),f in mean.groupby(['method','configuration']):
        d = f[f.split == 'development']
        scale,status = fit_scale(d.predictor,d['count'])
        coefficients.append({'method':method,'configuration':configuration,'scale':scale,'status':status,
            'development_end':development_end,'evaluation_start':evaluation_start,'fit_member_statistic':'perturbed_mean_or_control',
            'n_development':int((d.predictor.notna()&d['count'].notna()).sum())})
    c = pd.DataFrame(coefficients)
    a = p.merge(c[['method','configuration','scale','status']].rename(columns={'status':'fit_status'}),on=['method','configuration'],validate='many_to_one')
    ref = c[c.configuration == reference][['method','scale','status']].rename(columns={'scale':'reference_scale','status':'reference_fit_status'})
    if len(ref) != p.method.nunique():
        raise ValueError('Reference configuration missing for a method')
    a = a.merge(ref,on='method',validate='many_to_one')
    a['prediction_A'] = a.predictor*a.scale
    a['prediction_B'] = a.predictor*a.reference_scale
    # A known zero predictor is zero for every finite nonnegative scale, even
    # when the development scale is unidentifiable. Nonzero held-out predictors
    # remain unavailable with an unidentified scale; missing input never becomes 0.
    a.loc[(a.predictor==0)&(a.fit_status=='identically_zero_unrecoverable'),'prediction_A']=0.
    a.loc[(a.predictor==0)&(a.reference_fit_status=='identically_zero_unrecoverable'),'prediction_B']=0.
    held = a[a.split == 'evaluation'].copy()
    rows=[]
    for (method,configuration),f in held.groupby(['method','configuration']):
        f = f[f.member != 0] if (f.member != 0).any() else f
        g = f.groupby('time').agg(count=('count','first'),prediction_A=('prediction_A','mean'),prediction_B=('prediction_B','mean'))
        for comparison in ['A','B']:
            rows.append({'method':method,'configuration':configuration,'comparison':comparison,**metrics(g['count'],g['prediction_'+comparison])})
    m = pd.DataFrame(rows)
    contrasts = m[m.comparison == 'B'].merge(m[(m.configuration == reference)&(m.comparison == 'B')][['method','rmse']].rename(columns={'rmse':'reference_rmse'}),on='method')
    contrasts['rmse_change'] = contrasts.rmse-contrasts.reference_rmse
    contrasts['rmse_improvement'] = -contrasts.rmse_change
    physical_members=held[held.member!=0] if (held.member!=0).any() else held
    fixed = physical_members.groupby(['time','method','configuration'],as_index=False).prediction_B.mean()
    baseline = fixed[fixed.configuration == reference][['time','method','prediction_B']].rename(columns={'prediction_B':'reference_prediction'})
    fixed = fixed.merge(baseline,on=['time','method'])
    fixed['physical_change'] = fixed.prediction_B-fixed.reference_prediction
    return c,held,m,contrasts,fixed

def ensemble_statistics(predictions, value='prediction_A'):
    f = predictions[predictions.member != 0]
    return f.groupby(['time','method','configuration'])[value].agg(
        n='count',mean='mean',std='std',median='median',q05=lambda x:x.quantile(.05),q25=lambda x:x.quantile(.25),
        q75=lambda x:x.quantile(.75),q95=lambda x:x.quantile(.95)).reset_index()
