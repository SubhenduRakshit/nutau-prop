#!/usr/bin/env python3
import argparse
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from nutau_prop.config import load_config, energy_grid, y_kwargs
from nutau_prop.cross_sections import DISConfig
from nutau_prop.bsm import BSMModel
from nutau_prop.tables import build_table, save_table
from nutau_prop.util import nproc_default

def _chunk_job(payload):
    E,pdf_set,member,dc_dict,model_dict,ykw,antinu=payload
    return build_table(E,pdf_set,member,DISConfig(**dc_dict),BSMModel(**model_dict),ykw,antinu)

def _combine(parts):
    return dict(
        energies=np.concatenate([p["energies"] for p in parts]),
        y=parts[0]["y"], dy=parts[0]["dy"],
        sigma_cc=np.concatenate([p["sigma_cc"] for p in parts]),
        sigma_nc=np.concatenate([p["sigma_nc"] for p in parts]),
        p_cc=np.concatenate([p["p_cc"] for p in parts],axis=0),
        p_nc=np.concatenate([p["p_nc"] for p in parts],axis=0))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--config",default="config/default.yaml")
    ap.add_argument("--out",required=True)
    ap.add_argument("--kind",choices=["sm","charge_radius","derivative_em"],default="sm")
    ap.add_argument("--r2",type=float,default=0.0)
    ap.add_argument("--eps-em-inf",type=float,default=0.0)
    ap.add_argument("--mx",type=float,default=0.01)
    ap.add_argument("--antinu",action="store_true")
    ap.add_argument("--workers",type=int,default=nproc_default())
    args=ap.parse_args()
    c=load_config(args.config); E=energy_grid(c); ic=c["integration"]
    dc=dict(q2_min_gev2=float(ic["q2_min_gev2"]),w2_min_gev2=float(ic["w2_min_gev2"]),
            nx=int(ic["nx"]),sin2w=float(c["physics"]["sin2thetaW"]))
    model=dict(kind=args.kind,r2_cm2=args.r2,eps_em_inf=args.eps_em_inf,
               mediator_mass_gev=args.mx,
               antinu_sign=float(c["physics"]["antinu_charge_radius_sign"]))
    workers=max(1,min(int(args.workers),len(E)))
    if workers==1:
        t=_chunk_job((E,c["pdf"]["set"],c["pdf"]["member"],dc,model,y_kwargs(c),args.antinu))
    else:
        chunks=[x for x in np.array_split(E,workers) if len(x)]
        payloads=[(ch,c["pdf"]["set"],c["pdf"]["member"],dc,model,y_kwargs(c),args.antinu)
                  for ch in chunks]
        with ProcessPoolExecutor(max_workers=workers) as ex:
            t=_combine(list(ex.map(_chunk_job,payloads)))
    save_table(args.out,t,dict(kind=args.kind,r2_cm2=args.r2,eps_em_inf=args.eps_em_inf,
                              mediator_mass_gev=args.mx,pdf=c["pdf"]["set"],
                              member=c["pdf"]["member"],workers=workers,antinu=args.antinu))
    print(args.out)
if __name__=="__main__": main()
