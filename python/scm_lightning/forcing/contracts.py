"""Explicit SI contracts. No model state is mutated by diagnostics."""
import numpy as np
# DEPHY Cpd=3.5*Rd and MPAS cp=3.5*rgas: same theta/T exponent.
KAPPA = 2. / 7.

def exner(p):
    p = np.asarray(p)
    if np.any(p <= 0):
        raise ValueError('Pressure must be positive Pa')
    return (p / 1.e5) ** KAPPA

def theta_from_temperature(t, p):
    return np.asarray(t) / exner(p)

def q_to_r(q):
    q = np.asarray(q)
    if np.any((q < 0) | (q >= 1)):
        raise ValueError('Specific humidity outside [0,1)')
    return q / (1 - q)

def q_tendency_to_r(dq, q):
    return np.asarray(dq) / (1 - np.asarray(q)) ** 2

def horizontal_advection(field, u, v):
    lat = np.deg2rad(field.latitude)
    dx = field.differentiate('longitude') * 180 / np.pi / (6371000 * np.cos(lat))
    dy = field.differentiate('latitude') * 180 / np.pi / 6371000
    return -(u * dx + v * dy)

def layer_integral(concentration, dz):
    if np.any(np.asarray(dz) <= 0):
        raise ValueError('Layer thickness must be positive m')
    return np.sum(np.asarray(concentration) * dz, axis=-1)

def isotherm(t, z, target, values):
    """Lowest descending crossing; missing crossing remains NaN."""
    t, z, values = map(np.asarray, (t, z, values))
    if np.any(np.diff(z) <= 0):
        raise ValueError('Height must increase')
    hits = np.flatnonzero((t[:-1] >= target) & (t[1:] < target))
    if not len(hits):
        return np.nan
    k = hits[0]
    a = (target - t[k]) / (t[k+1] - t[k])
    return (1-a)*values[k] + a*values[k+1]

def precipitation_rate(rainnc, seconds):
    """rainnc includes total solid/liquid water equivalent; never add snownc."""
    a, dt = np.diff(rainnc), np.diff(seconds)
    if np.any(a < -1.e-7) or np.any(dt <= 0):
        raise ValueError('Reset or nonmonotonic precipitation/time')
    return np.maximum(a, 0) / dt

def area_km2(box):
    w,s,e,n = box
    if not (-180 <= w < e <= 180 and -90 <= s < n <= 90):
        raise ValueError('Expected west,south,east,north')
    return 6371.**2 * np.deg2rad(e-w) * (np.sin(np.deg2rad(n))-np.sin(np.deg2rad(s)))
