#!/usr/bin/env python3
"""Refine all a intervals and directly verify maxima; no detector events."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar
from pkr_normalization_fold_v0_86 import OUT, ROOT, fold, collider_refs
from nutau_prop.pkr_resolved import allowed_gq_squared,operator_factor,invisible_width,quark_width_coefficient

def main():
    frame=pd.read_csv(OUT/'base_folds.csv'); refs=collider_refs(); rows=[]
    for M in (50.,100.,220.):
        q=frame[frame.M_GeV==M].sort_values('a_GeV_inv')
        la=np.log(q.a_GeV_inv.to_numpy())
        logf=PchipInterpolator(la,np.log(q.ratio44_gq1.to_numpy()))
        for op in ('4.4','4.3'):
            for label,dijet in [('global_floor_sufficient',.03),('global_weakest_relaxed',.13),
                               ('approx_mass_value',{50.:.038,100.:.083,220.:.116}[M])]:
                for rcap in (.8,1.,1.2):
                    def negative(loga):
                        a=math.exp(loga)
                        z=allowed_gq_squared(M,a,op,rcap/refs[M]['K'],dijet)
                        return -z*operator_factor(op)*math.exp(float(logf(loga)))
                    # Every interval, not only the first crossing or one local extremum.
                    candidates=[(negative(x),x) for x in la]
                    for lo,hi in zip(la[:-1],la[1:]):
                        result=minimize_scalar(negative,bounds=(lo,hi),method='bounded',options={'xatol':1e-10})
                        candidates.append((result.fun,result.x))
                    val,loga=min(candidates); a=math.exp(loga)
                    z=allowed_gq_squared(M,a,op,rcap/refs[M]['K'],dijet)
                    I=invisible_width(M,a,op); B=quark_width_coefficient(M); P=z*I/(B*z+I)
                    rows.append(dict(operator=op,M_GeV=M,dijet_screen=label,r_cap=rcap,
                        a_GeV_inv=a,gq=math.sqrt(z),interpolated_ratio100=-val,
                        monojet_r=refs[M]['K']*P,Gamma_over_M=(B*z+I)/M,
                        gq_equiv_dijet=math.sqrt(max(0,z-P))))
    tasks=list(dict.fromkeys((1e8,r['M_GeV'],r['a_GeV_inv'],1/r['a_GeV_inv'],n)
                            for r in rows for n in (72,112)))
    with ProcessPoolExecutor(max_workers=14) as pool:
        verified=list(pool.map(fold,tasks))
    lookup={(r['M_GeV'],r['a_GeV_inv'],r['nquad']):r for r in verified}
    for row in rows:
        low=lookup[(row['M_GeV'],row['a_GeV_inv'],72)]
        high=lookup[(row['M_GeV'],row['a_GeV_inv'],112)]
        ratio=high['ratio44_gq1']*row['gq']**2*operator_factor(row['operator'])
        row.update(direct_ratio100=ratio,interpolation_relative_error=abs(ratio/row['interpolated_ratio100']-1),
                   quadrature_relative_change=abs(high['ratio44_gq1']/low['ratio44_gq1']-1))
    out=pd.DataFrame(rows); out.to_csv(OUT/'refined_screens.csv',index=False)
    pd.DataFrame(verified).to_csv(OUT/'refined_verification_folds.csv',index=False)
    assert out.quadrature_relative_change.max()<.01
    assert out.interpolation_relative_error.max()<.01
    print(out[(out.r_cap==1)&(out.dijet_screen!='global_floor_sufficient')].to_string(index=False))
    print('max relative quadrature change',out.quadrature_relative_change.max())
    print('max interpolation error',out.interpolation_relative_error.max())
    (OUT/'refinement_validation.json').write_text(json.dumps(dict(version='v0.86',status='PASS',
        worker_processes=14,all_intervals_optimized=True,
        max_relative_quadrature_change=float(out.quadrature_relative_change.max()),
        max_relative_interpolation_error=float(out.interpolation_relative_error.max()),
        note='Numerical maxima under specified screens; not a general exclusion theorem.'),indent=2)+'\n')

if __name__=='__main__': main()
