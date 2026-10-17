"""Independent arithmetic checks; no imports from the analysis implementations."""
import argparse
import csv
import json
import math
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as h: return list(csv.DictReader(h))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--union',type=Path,required=True)
    ap.add_argument('--category',type=Path,required=True)
    ap.add_argument('--analysis',type=Path,required=True)
    ap.add_argument('--controls',type=Path,required=True)
    ap.add_argument('--archive',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    assert not a.output.exists()
    ledger=rows(a.union/'model_gene_ledger.csv')
    donor_levels=rows(a.union/'donor_levels.csv')
    category_levels=rows(a.category/'donor_category_levels.csv')
    bounds=rows(a.analysis/'results/fraction_bounds.csv')
    dec=rows(a.analysis/'results/decomposition_summary.csv')
    exhaustive=[];decomp=[]
    for b in bounds:
        context,policy=b['context'],b['policy']
        selected=[r for r in ledger if (r['context'],r['policy'])==(context,policy)]
        R=[sum((F(r[k])*F(r['population_weight']) for r in selected),F()) for k in ('old','new')]
        coeff=defaultdict(F);labels={}
        for r in selected:
            key=(r['sample_id'],r['gene'])
            coeff[key]+=F(r['population_weight'])*(F(r['new'])*R[0]-F(r['old'])*R[1])/(R[0]*R[1])
            labels[key]=None if r['label_state']!='MEASURED' else int(r['label_value'])
        observed=sum((v*labels[k] for k,v in coeff.items() if labels[k] is not None),F())
        unknown=[v for k,v in coeff.items() if labels[k] is None]
        den=math.lcm(*(v.denominator for v in unknown))
        ints=[v.numerator*(den//v.denominator) for v in unknown]
        # Exhaust every actual completion by a Gray-code walk; one integer changes at a time.
        low=high=running=previous=0
        for n in range(1,1<<len(ints)):
            gray=n^(n>>1);changed=gray^previous;bit=changed.bit_length()-1
            running+=ints[bit] if gray&changed else -ints[bit]
            low=min(low,running);high=max(high,running);previous=gray
        assert observed+F(low,den)==F(b['sharp_lower'])
        assert observed+F(high,den)==F(b['sharp_upper'])
        exhaustive.append(dict(context=context,policy=policy,unknown_count=len(ints),enumerated_completions=1<<len(ints),status='PASS'))
        # Reconstruct partitions from accepted category/donor tables, not the new ledger totals.
        for partition in ('donor','role'):
            values=defaultdict(lambda: [[F(),F()],[F(),F()]])
            source=donor_levels if partition=='donor' else category_levels
            for r in source:
                if r['context']!=context or r['policy'].replace('DawnRank_strict','DawnRank')!=policy: continue
                if partition=='role' and r['category'] not in ('ONCOGENE','TSG','DUAL_ONCOGENE_TSG','OTHER_OR_UNSPECIFIED'): continue
                unit=r['donor_id'] if partition=='donor' else r['category']
                t=int(r['state']);values[unit][t][0]+=F(r['R']);values[unit][t][1]+=F(r['H'])
            denom=[sum((v[t][0] for v in values.values()),F()) for t in (0,1)]
            p=[sum((v[t][1] for v in values.values()),F())/denom[t] for t in (0,1)]
            centre=(p[0]+p[1])/2
            within=weight=entry=exit=F()
            for v in values.values():
                r0,h0=v[0];r1,h1=v[1]
                a0,a1=r0/denom[0],r1/denom[1]
                if r0 and r1:
                    f0,f1=h0/r0,h1/r1
                    within+=(a1+a0)*(f1-f0)/2
                    weight+=(a1-a0)*(f1+f0-p[1]-p[0])/2
                elif r1: entry+=h1/denom[1]-centre*a1
                elif r0: exit+=centre*a0-h0/denom[0]
            expected=next(r for r in dec if (r['context'],r['policy'],r['partition'])==(context,policy,partition))
            assert (within,weight,entry,exit)==tuple(F(expected[k]) for k in ('within','reweight','entry','exit'))
            assert sum((within,weight,entry,exit),F())==p[1]-p[0]
            decomp.append(dict(context=context,policy=policy,partition=partition,units=len(values),status='PASS'))
    # Direct old/new completed numerators on every frozen control; no optimiser helper is used.
    freeze=json.loads((a.controls/'ALLOCATION_FREEZE.json').read_text(encoding='utf-8-sig'))
    callrows=json.loads((a.controls/'ALL_SELECTED_LABELS.json').read_text(encoding='utf-8-sig'))['calls']
    calls={(r['sample_id'],r['gene']):r for r in callrows}
    alpha={r['sample_id']:F(r['alpha']) for r in freeze['models']}
    cancer={g for g,c in freeze['annotations'].items() if 'UNLISTED' not in c}
    controls=rows(a.analysis/'sensitivity_results/graph_control_endpoints.csv')
    witnesses=[]
    for r in controls:
        arms=[]
        for state in (r['old_state'],r['new_state']):
            arms.append({(s,g):F(w)*alpha[s] for s,v in freeze['states'][state].items() for g,w in v.items() if g in cancer})
        den=[sum(arm.values(),F()) for arm in arms]
        got=[]
        for upper in (False,True):
            h=[F(),F()]
            for key in arms[0].keys()|arms[1].keys():
                call=calls[key]
                w0,w1=arms[0].get(key,F()),arms[1].get(key,F())
                y=call['value'] if call['state']=='MEASURED' else int((w1*den[0]>w0*den[1]) if upper else (w1*den[0]<w0*den[1]))
                h[0]+=w0*y;h[1]+=w1*y
            got.append(h[1]/den[1]-h[0]/den[0])
        assert got==[F(r['sharp_lower']),F(r['sharp_upper'])]
        witnesses.append(dict(contrast=r['contrast'],status='PASS'))
    groups=[]
    for group in ('original','low_count','high_count','threshold'):
        selected=[r for r in controls if r['pair'].startswith(group) and r['mode']!='DEGREE']
        groups.append(dict(group=group,n=len(selected),**{field+'_min':min(float(F(r[field])) for r in selected) for field in ('delta_R','delta_H','observed_delta','sharp_lower','sharp_upper')},
                           **{field+'_max':max(float(F(r[field])) for r in selected) for field in ('delta_R','delta_H','observed_delta','sharp_lower','sharp_upper')},
                           positive_returned=sum(F(r['delta_R'])>0 for r in selected),negative_completed_fraction=sum(F(r['sharp_upper'])<0 for r in selected)))
    # Independently reconstruct all budget numerators from each seed's raw source allocations.
    labels=json.loads((a.archive/'data/crc/selected_labels.json').read_text(encoding='utf-8-sig'))['calls']
    bcalls={(r['sample_id'],r['gene']):r for r in labels}
    categories=json.loads((a.archive/'data/crc/categories.json').read_text(encoding='utf-8-sig'))
    source_paths={'PRODIGY':a.archive/'results/crc_prodigy_allocations.json','DawnRank':a.archive/'results/crc_precision/allocations.json'}
    sources={p:json.loads(f.read_text(encoding='utf-8-sig')) for p,f in source_paths.items()}
    budgets=rows(a.analysis/'sensitivity_results/budget_endpoints.csv');budget_witnesses=[]
    for b in budgets:
        p,k=b['policy'],int(b['budget']);arms=[defaultdict(F),defaultdict(F)]
        nseed=5 if p=='PRODIGY' else 1
        for r in sources[p]['allocations']:
            cell=r['cell'].replace('_','')
            if r['budget']!=k or cell not in ('G0H0','G1H1'): continue
            t=cell=='G1H1'
            for g,w in r['weights'].items():
                if g in categories['roles']: arms[t][r['sample_id'],g]+=F(w)*alpha[r['sample_id']]/nseed
        den=[sum(arm.values(),F()) for arm in arms]
        got=[]
        for upper in (False,True):
            h=[F(),F()]
            for key in arms[0].keys()|arms[1].keys():
                call=bcalls[key];w0,w1=arms[0].get(key,F()),arms[1].get(key,F())
                y=call['value'] if call['state'].startswith('MEASURED') else int((w1*den[0]>w0*den[1]) if upper else (w1*den[0]<w0*den[1]))
                h[0]+=w0*y;h[1]+=w1*y
            got.append(h[1]/den[1]-h[0]/den[0])
        assert den==[F(b['R0']),F(b['R1'])]
        assert got==[F(b['sharp_lower']),F(b['sharp_upper'])]
        budget_witnesses.append(dict(policy=p,budget=k,status='PASS'))
    result=dict(status='PASS',actual_completion_enumeration=exhaustive,accepted_source_decompositions=decomp,control_endpoint_witnesses=witnesses,budget_endpoint_witnesses=budget_witnesses,DawnRank_control_groups=groups)
    with a.output.open('x',encoding='utf-8') as h: json.dump(result,h,indent=2)
    print(json.dumps({'status':'PASS','exhaustive':exhaustive,'graph_groups':groups},indent=2))

if __name__=='__main__': main()
