"""Derive listed-union endpoints from accepted CRC allocations and fixed graph controls."""
import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
import analyse_fraction as shared

def load(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))

def normalise_calls(data):
    result={}
    for r in data['calls']:
        key=(r['sample_id'],r['gene'])
        value=('MEASURED' if r['state'].startswith('MEASURED') else r['state'],r['value'])
        shared.require(key not in result or result[key]==value,'Unique consistent source call')
        result[key]=value
    return result

def summarise_pair(old,new,calls):
    totals=[]
    for values in (old,new):
        t=dict.fromkeys(('R','H','Z','U'),F())
        for key,w in values.items():
            shared.require(key in calls,'Selected identity source call supplied')
            state,y=calls[key]
            t['R']+=w
            t[('H' if y else 'Z') if state=='MEASURED' else 'U']+=w
        shared.require(t['R']==t['H']+t['Z']+t['U'],'Pair mass conservation')
        shared.require(t['R']>0,'Pooled listed denominator positive')
        t['p']=t['H']/t['R']
        totals.append(t)
    known,unknown=F(),{}
    mass_known,mass_unknown=F(),{}
    for key in old.keys()|new.keys():
        c=new.get(key,F())/totals[1]['R']-old.get(key,F())/totals[0]['R']
        d=new.get(key,F())-old.get(key,F())
        state,y=calls[key]
        if state=='MEASURED': known+=c*y;mass_known+=d*y
        else: unknown[key]=c;mass_unknown[key]=d
    lo,hi=shared.optimise(known,unknown)
    mlo,mhi=shared.optimise(mass_known,mass_unknown)
    shared.require(known==totals[1]['p']-totals[0]['p'],'Fraction known-part reconciliation')
    shared.require(mass_known==totals[1]['H']-totals[0]['H'],'Capture known-part reconciliation')
    return {**{k+str(i):v for i,t in enumerate(totals) for k,v in t.items()},
            'delta_R':totals[1]['R']-totals[0]['R'],'delta_H':mass_known,'observed_delta':known,
            'sharp_lower':lo,'sharp_upper':hi,'capture_lower':mlo,'capture_upper':mhi}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--archive',type=Path,required=True)
    ap.add_argument('--controls',type=Path,required=True)
    ap.add_argument('--union',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    shared.require(not a.output.exists(),'Fresh output directory')
    files=[a.archive/'results/crc_precision/allocations.json',a.archive/'results/crc_prodigy_allocations.json',
           a.archive/'data/crc/categories.json',a.archive/'data/crc/selected_labels.json',
           a.controls/'ALLOCATION_FREEZE.json',a.controls/'CONTRAST_FREEZE.json',a.controls/'ALL_SELECTED_LABELS.json',
           a.controls/'population_state_levels.csv',a.union/'population_levels.csv']
    dawn,pro,cat,labels,freeze,plan,control_labels=map(load,files[:7])
    roles=set(cat['roles'])
    calls=normalise_calls(labels)
    donor_counts=Counter(str(m['donor_id']) for m in dawn['models'])
    models={m['model_id']:m for m in dawn['models']}
    shared.require(len(models)==85 and len(donor_counts)==83,'Complete CRC population')
    weights={m['sample_id']:F(1,83*donor_counts[str(m['donor_id'])]) for m in models.values()}
    previous=shared.read(a.union/'population_levels.csv')
    budget_rows=[]
    for method,source in (('PRODIGY',pro),('DawnRank',dawn)):
        for k in (1,5,10,20):
            state_values=[defaultdict(F),defaultdict(F)]
            selected={}
            for r in source['allocations']:
                cell=r['cell'].replace('_','')
                if r['budget']!=k or cell not in ('G0H0','G1H1'): continue
                t=0 if cell=='G0H0' else 1
                model=models[r['model_id']]
                shared.require(r['sample_id']==model['sample_id'],'Budget model/sample identity')
                record_key=(r['model_id'],t,str(r['seed']))
                shared.require(record_key not in selected,'No duplicate budget state')
                selected[record_key]=True
                w={g:F(x) for g,x in r['weights'].items()}
                shared.require(all(F()<=x<=1 for x in w.values()),'Allocation range')
                query=set(model['query'])
                shared.require(set(w)<=query,'Allocation within fixed query')
                capacity=min(k,len(query)) if method=='DawnRank' else min(k,len(set(r['menu'])&query))
                shared.require(sum(w.values(),F())==capacity,'Actual offered-budget capacity')
                multiplier=weights[r['sample_id']]/(5 if method=='PRODIGY' else 1)
                for g,v in w.items():
                    if g in roles: state_values[t][r['sample_id'],g]+=multiplier*v
            shared.require(len(selected)==85*2*(5 if method=='PRODIGY' else 1),'All original model-seed-states')
            if method=='PRODIGY':
                expected={str(s) for s in (104729,130363,155921,196613,228017)}
                shared.require(all({seed for mid,t,seed in selected if (mid,t)==(m,j)}==expected for m in models for j in (0,1)),'Five original seeds in each model/state')
            result=summarise_pair(*state_values,calls)
            if k==10:
                for t in (0,1):
                    accepted=next(r for r in previous if (r['context'],r['policy'],r['state'])==('CRC',method,str(t)))
                    shared.require(all(result[f+str(t)]==F(accepted[f]) for f in ('R','H','Z','U')),'k10 accepted union reconciliation')
            budget_rows.append(dict(context='CRC',policy=method,budget=k,**result))
    # Fixed graph-control population and all frozen contrasts; no fitting/reranking.
    gcalls=normalise_calls(control_labels)
    gmodels=freeze['models']
    shared.require(len(gmodels)==85 and len({m['donor_id'] for m in gmodels})==83,'Complete graph-control population')
    gw={m['sample_id']:F(m['alpha']) for m in gmodels}
    shared.require(gw==weights,'Identical graph-control donor weighting')
    groles={g for g,x in freeze['annotations'].items() if any(c in ('ONCOGENE','TSG','DUAL_ONCOGENE_TSG','OTHER_OR_UNSPECIFIED') for c in x)}
    expected_roles=roles&set(freeze['annotations'])
    shared.require(groles==expected_roles,'Identical role union on common domain')
    states={}
    level_rows=[]
    accepted_levels=shared.read(a.controls/'population_state_levels.csv')
    for name,allocs in freeze['states'].items():
        shared.require(set(allocs)==set(weights),'Complete state roster')
        values={}
        for sample,alloc in allocs.items():
            shared.require(sum((F(x) for x in alloc.values()),F())==min(10,len(freeze['queries'][sample])),'Graph state capacity')
            for g,w in alloc.items():
                if g in groles: values[sample,g]=gw[sample]*F(w)
        states[name]=values
        result=summarise_pair(values,values,gcalls)
        for f in ('R','H','Z','U'):
            accepted=sum((F(r[f]) for r in accepted_levels if r['state']==name and r['unit']=='population' and r['category'] in ('ONCOGENE','TSG','DUAL_ONCOGENE_TSG','OTHER_OR_UNSPECIFIED')),F())
            shared.require(result[f+'0']==accepted,'Accepted control-state role sum')
        level_rows.append(dict(state=name,**{f:result[f+'0'] for f in ('R','H','Z','U','p')}))
    control_rows=[]
    for name,p in plan['contrasts'].items():
        result=summarise_pair(states[p['old_state']],states[p['new_state']],gcalls)
        control_rows.append(dict(contrast=name,pair=p['pair'],mode=p['mode'],old_state=p['old_state'],new_state=p['new_state'],**result))
    shared.require(len(control_rows)==27 and len(states)==36,'Complete frozen graph-control roster')
    baseline=summarise_pair(states['G0_T700_H0'],states['G1_T700_H1'],gcalls)
    b10=next(r for r in budget_rows if r['policy']=='DawnRank' and r['budget']==10)
    shared.require(all(v==b10[k] for k,v in baseline.items()),'Independent graph input matches original joint DawnRank union endpoints')
    a.output.mkdir(parents=True)
    shared.save(a.output/'budget_endpoints.csv',budget_rows)
    shared.save(a.output/'graph_control_endpoints.csv',control_rows)
    shared.save(a.output/'graph_state_levels.csv',level_rows)
    notes={'scope':'CRC PRODIGY and DawnRank budgets 1,5,10,20; all 27 existing frozen graph contrasts (18 DawnRank,9 degree), k10. Fixed-cohort descriptive endpoints only.',
           'not_extended':'No ESCA secondary budgets or degree secondary budgets; current available allocation archive has k10 only for these policies. No missing calls imputed and no new donor-envelope family.',
           'status':'PASS','checks':shared.CHECKS,'native_fits':0,'rerankings':0,'inputs':[shared.bind(p) for p in files],
           'code':[shared.bind(Path(__file__)),shared.bind(Path(__file__).with_name('analyse_fraction.py'))],
           'protocol':shared.bind(Path(__file__).with_name('PROTOCOL.md')),
           'outputs':[shared.bind(p) for p in sorted(a.output.iterdir())]}
    with (a.output/'manifest.json').open('x',encoding='utf-8') as h: json.dump(notes,h,indent=2)
    for r in budget_rows: print(r['policy'],r['budget'],'deltaR/H/p',*[round(float(r[k]),6) for k in ('delta_R','delta_H','observed_delta')],'sharp pp',round(float(100*r['sharp_lower']),6),round(float(100*r['sharp_upper']),6))
    for mode in ('H0','H1','DEGREE'):
        rs=[r for r in control_rows if r['mode']==mode]
        print('graph',mode,'R range',min(float(r['delta_R']) for r in rs),max(float(r['delta_R']) for r in rs),'observed p pp',min(float(100*r['observed_delta']) for r in rs),max(float(100*r['observed_delta']) for r in rs),'negative sharp',sum(r['sharp_upper']<0 for r in rs),'/',len(rs))
    print('PASS',shared.CHECKS)

if __name__=='__main__': main()
