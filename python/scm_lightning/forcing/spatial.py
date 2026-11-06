"""Spherical cell/box overlap weights, with explicit fixed-area support."""
import numpy as np
import xarray as xr

def box_weights(latitude,longitude,box):
    w,s,e,n=box
    def edges(a):
        a=np.asarray(a,float)
        if len(a)<2 or np.any(np.diff(a)<=0):raise ValueError('Increasing regional coordinate required')
        return np.r_[a[0]-(a[1]-a[0])/2,(a[:-1]+a[1:])/2,a[-1]+(a[-1]-a[-2])/2]
    la,lo=edges(latitude),edges(longitude)
    if lo[0]>w or lo[-1]<e or la[0]>s or la[-1]<n:
        raise ValueError('Fixed box outside source grid-cell support')
    dy=np.maximum(0,np.sin(np.deg2rad(np.minimum(la[1:],n)))-np.sin(np.deg2rad(np.maximum(la[:-1],s))))
    dx=np.deg2rad(np.maximum(0,np.minimum(lo[1:],e)-np.maximum(lo[:-1],w)))
    return xr.DataArray(6371.**2*dy[:,None]*dx[None,:],dims=['latitude','longitude'],coords={'latitude':np.asarray(latitude),'longitude':np.asarray(longitude)})

def subset(data,box):
    weight=box_weights(data.latitude.values,data.longitude.values,box)
    iy=np.flatnonzero(weight.sum('longitude').values>0)
    ix=np.flatnonzero(weight.sum('latitude').values>0)
    return data.isel(latitude=iy,longitude=ix),weight.isel(latitude=iy,longitude=ix)
