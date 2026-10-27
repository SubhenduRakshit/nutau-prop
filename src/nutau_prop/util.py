import os, hashlib, json
from pathlib import Path
import numpy as np

def nproc_default():
    n=os.cpu_count() or 2
    return max(1,min(8,n//2 if n>2 else 1))

def powerlaw_per_logE(E,gamma=2.5,E0=1e5):
    # dPhi/dlnE = E dPhi/dE; arbitrary normalization.
    return (E/E0)**(1-gamma)

def cache_key(obj):
    s=json.dumps(obj,sort_keys=True,default=str).encode()
    return hashlib.sha256(s).hexdigest()[:16]
