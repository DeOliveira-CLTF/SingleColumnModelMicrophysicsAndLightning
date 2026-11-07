"""Exact-input CAPE memoization; no rounding, interpolation or physical changes."""
import hashlib
import importlib.metadata
from pathlib import Path
import sqlite3
import numpy as np
from scm_lightning.provenance import sha256


class CapeCache:
    def __init__(self,path):
        path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
        self.connection=sqlite3.connect(path,timeout=120)
        self.connection.execute('PRAGMA journal_mode=WAL')
        self.connection.execute('PRAGMA synchronous=NORMAL')
        self.connection.execute('CREATE TABLE IF NOT EXISTS cape (namespace TEXT,key BLOB,value REAL,PRIMARY KEY(namespace,key)) WITHOUT ROWID')
        self.connection.commit()
        self.namespace='exact_dtype_v2:most_cape:'+importlib.metadata.version('MetPy')+':'+sha256(Path(__file__).with_name('adaptations.py'))
        self.pending={};self.hits=0;self.misses=0

    def key(self,p,t,r):
        # Preserve both values and input precision; one changed bit is a new key.
        h=hashlib.sha256()
        for values in [p,t,r]:
            a=np.ascontiguousarray(values);h.update(a.dtype.str.encode());h.update(np.asarray(a.shape,dtype='<i8').tobytes());h.update(a.tobytes())
        return h.digest()

    def seed(self,p,t,r,value):
        self.pending[self.key(p,t,r)]=float(value) if np.isfinite(value) else None
        if len(self.pending)>=256:self.flush()

    def get(self,p,t,r,provider):
        key=self.key(p,t,r)
        if key in self.pending:
            self.hits+=1;value=self.pending[key]
            return np.nan if value is None else value
        row=self.connection.execute('SELECT value FROM cape WHERE namespace=? AND key=?',(self.namespace,key)).fetchone()
        if row is not None:
            self.hits+=1;return np.nan if row[0] is None else float(row[0])
        self.misses+=1;value=provider(p,t,r,layer_selection='most_cape')
        self.seed(p,t,r,value)
        return value

    def flush(self):
        if self.pending:
            self.connection.executemany('INSERT OR IGNORE INTO cape VALUES (?,?,?)',[(self.namespace,k,v) for k,v in self.pending.items()])
            self.connection.commit();self.pending.clear()

    def close(self):
        self.flush();self.connection.close()
