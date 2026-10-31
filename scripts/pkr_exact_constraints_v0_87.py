#!/usr/bin/env python3
"""Exact CMS table, analytic coupling envelope, and public response inventory.

No automatic Git actions and no substitute detector likelihood claimed.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import hashlib
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize_scalar,brentq
from pkr_normalization_fold_v0_86 import ROOT,fold,collider_refs
from nutau_prop.pkr_resolved import allowed_gq_squared,invisible_width,quark_width_coefficient,operator_factor

OUT=ROOT/'results/v0.87_constraints'
TABLE_URL='https://www.hepdata.net/download/table/ins3133152/Figure%202a/json'
EHE_URL='https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:10.7910/DVN/JHK49D'

def fetch(url,path):
    if not path.exists():
        response=requests.get(url,timeout=45); response.raise_for_status()
        data=response.json(); path.write_text(json.dumps(data,indent=2)+'\n')
    return json.loads(path.read_text())

def inventory():
    directory=ROOT/'cache/HESE-7-year-data-release/HESE-7-year-data-release/resources/data'
    records=[]
    for name in ('HESE_mc_truth.json','HESE_mc_observable.json','HESE_mc_flux.json'):
        p=directory/name; raw=p.read_bytes(); data=json.loads(raw)
        records.append(dict(file=name,sha256=hashlib.sha256(raw).hexdigest(),
                            fields={k:len(v) for k,v in data.items()}))
    union={k for record in records for k in record['fields']}
    result=dict(files=records,has_interaction_type='interactionType' in union,
        has_true_inelasticity=bool(union&{'y','bjorkenY','inelasticity','BjorkenY'}),
        has_true_hadronic_energy=bool(union&{'hadronicEnergy','trueHadronicEnergy','trueDepositedEnergy'}),
        explanation='Reconstructed deposited energy is not true yE. The old HESE loader also discards interactionType.')
    # Constructive non-identifiability proof, no detector assumptions hidden.
    sm=np.array([.5,.5]); np_shape=np.array([.75,.25])
    A=np.array([.2,.8]); B=np.array([.8,.2])
    assert np.isclose(A@sm,B@sm) and not np.isclose(A@np_shape,B@np_shape)
    result['response_counterexample']=dict(SM_y_probabilities=sm.tolist(),BSM_y_probabilities=np_shape.tolist(),
        response_A=A.tolist(),response_B=B.tolist(),SM_acceptance_A=float(A@sm),SM_acceptance_B=float(B@sm),
        BSM_acceptance_A=float(A@np_shape),BSM_acceptance_B=float(B@np_shape))
    try:
        ehe=fetch(EHE_URL,OUT/'ehe_release_metadata.json')
        version=ehe['data']['latestVersion']
        result['EHE']=dict(url=EHE_URL,version=str(version.get('versionNumber')),files=[
            dict(name=f['dataFile']['filename'],description=f.get('description',''),
                 id=f['dataFile']['id']) for f in version['files']])
    except Exception as exc:
        result['EHE']=dict(url=EHE_URL,error=str(exc))
    (OUT/'icecube_input_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    table=fetch(TABLE_URL,OUT/'cms_dijet_figure2a.json')
    bounds={float(row['x'][0]['value']):float(row['y'][0]['value']) for row in table['values']}
    pd.DataFrame([dict(M_GeV=m,observed_gq_limit=g) for m,g in bounds.items()]).to_csv(OUT/'exact_dijet_limits.csv',index=False)
    refs=collider_refs(); base=pd.read_csv(ROOT/'results/v0.86_normalization/base_folds.csv')
    rows=[]; targets=[]
    for M in (50.,100.,220.):
        assert M in bounds, 'Do not silently interpolate mass-specific bounds'
        q=base[base.M_GeV==M].sort_values('a_GeV_inv')
        la=np.log(q.a_GeV_inv.to_numpy()); logf=PchipInterpolator(la,np.log(q.ratio44_gq1.to_numpy()))
        dijet=bounds[M]
        for op in ('4.4','4.3'):
            fac=operator_factor(op)
            def F(x): return fac*math.exp(float(logf(x)))
            for cap in (.8,1.,1.2,1e100):
                def objective(x):
                    a=math.exp(x)
                    return -F(x)*allowed_gq_squared(M,a,op,cap/refs[M]['K'],dijet)
                candidates=[(objective(x),x) for x in la]
                for lo,hi in zip(la[:-1],la[1:]):
                    r=minimize_scalar(objective,bounds=(lo,hi),method='bounded',options={'xatol':1e-11})
                    candidates.append((r.fun,r.x))
                val,x=min(candidates); a=math.exp(x)
                z=allowed_gq_squared(M,a,op,cap/refs[M]['K'],dijet)
                I=invisible_width(M,a,op); B=quark_width_coefficient(M)
                rows.append(dict(operator=op,M_GeV=M,dijet_gq_observed=dijet,r_cap=cap,
                    screen='dijet_width_only' if cap==1e100 else 'dijet_and_monojet',
                    a_GeV_inv=a,gq=math.sqrt(z),interpolated_ratio100=-val,
                    monojet_r=refs[M]['K']*z*I/(B*z+I),Gamma_over_M=(B*z+I)/M))
            # At fixed a, minimal gq for 5% gives minimal collider signals.
            def gap(x):
                return allowed_gq_squared(M,math.exp(x),op,1e100,dijet)-.05/F(x)
            def required_r(x):
                a=math.exp(x); z=.05/F(x)
                I=invisible_width(M,a,op); B=quark_width_coefficient(M)
                return refs[M]['K']*z*I/(B*z+I)
            grid=np.linspace(la[0],la[-1],6001); values=np.array([gap(x) for x in grid])
            crossings=[brentq(gap,x,y) for x,y,g,h in zip(grid[:-1],grid[1:],values[:-1],values[1:]) if g*h<0]
            boundaries=[la[0],*crossings,la[-1]]; possible=[]
            for lo,hi in zip(boundaries[:-1],boundaries[1:]):
                if gap((lo+hi)/2)<0: continue
                # Optimize every base-grid segment inside each feasible interval.
                knots=np.r_[lo,la[(la>lo)&(la<hi)],hi]
                possible.extend((required_r(x),x) for x in knots)
                for x,y in zip(knots[:-1],knots[1:]):
                    opt=minimize_scalar(required_r,bounds=(x,y),method='bounded')
                    possible.append((opt.fun,opt.x))
            if possible:
                r,x=min(possible); a=math.exp(x); z=.05/F(x)
                targets.append(dict(operator=op,M_GeV=M,status='TARGET_REQUIRES_MONOJET_R',
                    min_r_for_5pct=r,a_GeV_inv=a,gq=math.sqrt(z),dijet_gq_observed=dijet))
            else:
                targets.append(dict(operator=op,M_GeV=M,status='NO_DIJET_WIDTH_VALID_5PCT_POINT',
                    min_r_for_5pct=None,a_GeV_inv=None,gq=None,dijet_gq_observed=dijet))
    tasks=list(dict.fromkeys((E,row['M_GeV'],row['a_GeV_inv'],1/row['a_GeV_inv'],n)
                            for row in rows for E in (1e8,2.2e8) for n in (72,112)))
    tasks+=list(dict.fromkeys((1e8,row['M_GeV'],row['a_GeV_inv'],1/row['a_GeV_inv'],n)
                             for row in targets if row['a_GeV_inv'] is not None for n in (72,112)))
    print(f'{len(tasks)} direct folds on 14 processes; exact CMS mass entries, no detector rerun',flush=True)
    with ProcessPoolExecutor(max_workers=14) as pool: verifications=list(pool.map(fold,tasks))
    lookup={(r['E_GeV'],r['M_GeV'],r['a_GeV_inv'],r['nquad']):r for r in verifications}
    for row in rows:
        for E,label in ((1e8,'100'),(2.2e8,'220')):
            low=lookup[(E,row['M_GeV'],row['a_GeV_inv'],72)]
            high=lookup[(E,row['M_GeV'],row['a_GeV_inv'],112)]
            row['ratio'+label]=row['gq']**2*operator_factor(row['operator'])*high['ratio44_gq1']
            row['quadrature_change'+label]=abs(high['ratio44_gq1']/low['ratio44_gq1']-1)
        row['interpolation_error']=abs(row['ratio100']/row['interpolated_ratio100']-1)
    output=pd.DataFrame(rows); output.to_csv(OUT/'exact_dijet_envelopes.csv',index=False)
    for row in targets:
        if row['a_GeV_inv'] is not None:
            high=lookup[(1e8,row['M_GeV'],row['a_GeV_inv'],112)]
            row['verified_target_ratio']=row['gq']**2*operator_factor(row['operator'])*high['ratio44_gq1']
            assert abs(row['verified_target_ratio']/.05-1)<.01
    pd.DataFrame(targets).to_csv(OUT/'minimum_r_for_5pct.csv',index=False)
    pd.DataFrame(verifications).to_csv(OUT/'direct_folds.csv',index=False)
    assert output.interpolation_error.max()<.01
    assert max(output.quadrature_change100.max(),output.quadrature_change220.max())<.01
    inv=inventory()
    print(output[output.r_cap==1].to_string(index=False));print(pd.DataFrame(targets).to_string(index=False))
    print('HESE true inelasticity available:',inv['has_true_inelasticity'])
    summary=dict(version='v0.87',status='EXACT_DIJET_INPUT_AND_CONDITIONAL_COLLIDER_ENVELOPE_COMPLETE',
        table_doi=table['doi'],source_url=TABLE_URL,workers=14,
        max_interpolation_error=float(output.interpolation_error.max()),
        max_quadrature_change=float(max(output.quadrature_change100.max(),output.quadrature_change220.max())),
        limits='NWA with calibrated monojet shapes; not a combined confidence level or full IceCube fit')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

if __name__=='__main__':main()
