"""Numeric Fortran CSV reader, including omitted E in three-digit exponents."""
import pandas as pd

def numeric_csv(path,columns=None):
    return numeric_frame(pd.read_csv(path),columns)

def numeric_frame(f,columns=None):
    for name in columns or list(f.columns):
        if not pd.api.types.is_numeric_dtype(f[name]):
            text=f[name].astype(str).str.strip().str.replace('[dD]','E',regex=True)
            text=text.str.replace(r'(?<=\d)([+-]\d{3})$',r'E\1',regex=True)
            f[name]=pd.to_numeric(text,errors='raise')
    return f
