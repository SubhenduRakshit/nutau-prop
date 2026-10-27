from pathlib import Path
import yaml

def load_config(path):
    with open(path,"r",encoding="utf-8") as f:
        return yaml.safe_load(f)

def energy_grid(cfg):
    import numpy as np
    e=cfg["energy_grid"]
    return np.geomspace(float(e["emin_gev"]),float(e["emax_gev"]),int(e["n"]))

def y_kwargs(cfg):
    x=cfg["integration"]
    return dict(
        y_min=float(x["y_min"]),
        y_split=float(x["y_split"]),
        nlog=int(x["y_log_bins"]),
        nlin=int(x["y_linear_bins"]),
    )
