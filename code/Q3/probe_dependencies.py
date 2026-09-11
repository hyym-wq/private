import importlib.util
for p in ['numpy','scipy','numba','openpyxl']:
    s=importlib.util.find_spec(p)
    print(p,s.origin if s else None,flush=True)
