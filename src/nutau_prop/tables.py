import json
from pathlib import Path
import numpy as np
from .pdfs import LHAPDFNucleon
from .cross_sections import DISCalculator, DISConfig
from .bsm import BSMModel

def build_table(
    energies,
    pdf_set="CT18NLO",
    member=0,
    dis_config=None,
    model=None,
    y_kwargs=None,
    antinu=False,
):
    pdf = LHAPDFNucleon(pdf_set,member)
    calc = DISCalculator(pdf,dis_config or DISConfig())
    model = model or BSMModel()
    y_kwargs = y_kwargs or {}
    cc, nc, pcc, pnc = [], [], [], []
    y_ref = None
    dy_ref = None
    for E in energies:
        a = calc.spectrum(E,"CC",antinu=antinu,model=BSMModel(),**y_kwargs)
        b = calc.spectrum(E,"NC",antinu=antinu,model=model,**y_kwargs)
        if y_ref is None:
            y_ref,dy_ref = a["y"],a["dy"]
        cc.append(a["sigma"]); nc.append(b["sigma"])
        pcc.append(a["p"]); pnc.append(b["p"])
    return dict(
        energies=np.asarray(energies),
        y=np.asarray(y_ref),
        dy=np.asarray(dy_ref),
        sigma_cc=np.asarray(cc),
        sigma_nc=np.asarray(nc),
        p_cc=np.asarray(pcc),
        p_nc=np.asarray(pnc),
    )

def save_table(path, table, metadata=None):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path, **table)
    if metadata is not None:
        path.with_suffix(".json").write_text(json.dumps(metadata,indent=2))

def load_table(path):
    z=np.load(path)
    return {k:z[k] for k in z.files}
