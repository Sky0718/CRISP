"""Finite existing-data completed-fraction bounds and midpoint decomposition."""
import argparse
import csv
import hashlib
import itertools
import json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

CHECKS = 0

def require(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise AssertionError(message)

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))

def save(path, rows):
    rows = [{**row, **{key+'_decimal':float(value) for key,value in row.items() if isinstance(value,F)}} for row in rows]
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)

def bind(path):
    return dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest(), bytes=path.stat().st_size)

def optimise(known, unknown):
    return known+sum((min(F(),c) for c in unknown.values()), F()), known+sum((max(F(),c) for c in unknown.values()), F())

def decompose(groups, context, policy, partition):
    sums = [sum((v[t][0] for v in groups.values()), F()) for t in (0,1)]
    p = [sum((v[t][1] for v in groups.values()), F())/sums[t] for t in (0,1)]
    centre = sum(p)/2
    detail = []
    support_counts = defaultdict(int)
    changes = defaultdict(int)
    total = dict.fromkeys(('within','reweight','entry','exit'), F())
    for name, values in sorted(groups.items()):
        r = [values[t][0] for t in (0,1)]
        h = [values[t][1] for t in (0,1)]
        a = [r[t]/sums[t] for t in (0,1)]
        f = [h[t]/r[t] if r[t] else None for t in (0,1)]
        terms = dict.fromkeys(total, F())
        if r[0] and r[1]:
            support = 'common'
            terms['within'] = (a[0]+a[1])*(f[1]-f[0])/2
            terms['reweight'] = ((f[0]+f[1])/2-centre)*(a[1]-a[0])
            changes['increase' if f[1]>f[0] else 'decrease' if f[1]<f[0] else 'tie'] += 1
        elif r[1]:
            support = 'entrant'; terms['entry'] = a[1]*(f[1]-centre)
        elif r[0]:
            support = 'exit'; terms['exit'] = -a[0]*(f[0]-centre)
        else:
            support = 'neither'
        require(sum(terms.values(),F()) == (h[1]/sums[1]-h[0]/sums[0])-centre*(a[1]-a[0]), 'Unit midpoint identity')
        support_counts[support] += 1
        for k in total: total[k] += terms[k]
        detail.append(dict(context=context,policy=policy,partition=partition,unit=name,support=support,R0=r[0],R1=r[1],H0=h[0],H1=h[1],a0=a[0],a1=a[1],f0=f[0],f1=f[1],**terms))
    require(sum(total.values(),F()) == p[1]-p[0], 'Exact midpoint summation')
    summary = dict(context=context,policy=policy,partition=partition,p0=p[0],p1=p[1],delta=p[1]-p[0],**total,
                   **{k+'_units':support_counts[k] for k in ('common','entrant','exit','neither')},
                   **{k+'_common_units':changes[k] for k in ('increase','decrease','tie')})
    return summary, detail

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--union',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a = ap.parse_args()
    require(not a.output.exists(), 'Use fresh output directory')
    ledger_path = a.union/'model_gene_ledger.csv'
    levels_path = a.union/'population_levels.csv'
    donors_path = a.union/'donor_levels.csv'
    ledger, existing, donor_existing = map(read, (ledger_path,levels_path,donors_path))
    fixture = {'a':F(-2,3), 'b':F(1,4), 'c':F()}
    exact = [F(2,7)+sum((c*y for c,y in zip(fixture.values(),bits)),F()) for bits in itertools.product((0,1),repeat=3)]
    require(optimise(F(2,7),fixture)==(min(exact),max(exact)), 'Enumerated sharp-bound fixture')
    bounds, coefficients, summaries, details, role_levels = [], [], [], [], []
    for context, policy in dict.fromkeys((r['context'],r['policy']) for r in existing):
        rows = [r for r in ledger if (r['context'],r['policy'])==(context,policy)]
        totals = [{k:F() for k in ('R','H','Z','U')} for t in (0,1)]
        donor_groups, role_groups = {}, {}
        for r in donor_existing:
            if (r['context'],r['policy'])==(context,policy):
                donor_groups.setdefault(r['donor_id'],[[F(),F()],[F(),F()]])
        for r in rows:
            role_groups.setdefault(r['role'],[[F(),F()],[F(),F()]])
            require(r['role'] in ('ONCOGENE','TSG','DUAL_ONCOGENE_TSG','OTHER_OR_UNSPECIFIED'), 'Exclusive cancer roles')
            for t,k in enumerate(('old','new')):
                w = F(r[k])*F(r['population_weight'])
                totals[t]['R'] += w
                field = ('H' if int(r['label_value']) else 'Z') if r['label_state']=='MEASURED' else 'U'
                totals[t][field] += w
                donor_groups[r['donor_id']][t][0] += w
                role_groups[r['role']][t][0] += w
                if field=='H':
                    donor_groups[r['donor_id']][t][1] += w
                    role_groups[r['role']][t][1] += w
        for t in (0,1):
            accepted = next(r for r in existing if (r['context'],r['policy'],r['state'])==(context,policy,str(t)))
            require(all(totals[t][k]==F(accepted[k]) for k in totals[t]), 'Accepted population mass reconciliation')
            require(totals[t]['R']==sum((totals[t][k] for k in ('H','Z','U')),F()), 'Population conservation')
        merged, calls = defaultdict(F), {}
        for r in rows:
            key = (r['sample_id'],r['gene'])
            call = (r['label_state'],r['label_value'])
            require(key not in calls or calls[key]==call, 'Shared identity has one source call')
            calls[key] = call
            merged[key] += F(r['population_weight'])*(F(r['new'])/totals[1]['R']-F(r['old'])/totals[0]['R'])
        known, unknown = F(), {}
        for key,c in merged.items():
            state,y = calls[key]
            if state=='MEASURED': known += c*int(y)
            else: unknown[key] = c
            coefficients.append(dict(context=context,policy=policy,sample_id=key[0],gene=key[1],label_state=state,label_value=y,coefficient=c,
                                     minimising_label=(int(c<0) if state!='MEASURED' else int(y)),
                                     maximising_label=(int(c>0) if state!='MEASURED' else int(y))))
        observed = totals[1]['H']/totals[1]['R']-totals[0]['H']/totals[0]['R']
        require(known==observed, 'Known coefficient sum equals observed fraction change')
        lo,hi = optimise(known,unknown)
        attained=[]
        for sign in (-1,1):
            complete = [F(),F()]
            for r in rows:
                key=(r['sample_id'],r['gene']); state,y=calls[key]
                label=int(y) if state=='MEASURED' else int(merged[key]*sign>0)
                for t,k in enumerate(('old','new')):
                    complete[t] += F(r['population_weight'])*F(r[k])*label
            attained.append(complete[1]/totals[1]['R']-complete[0]/totals[0]['R'])
        require(attained==[lo,hi], 'Explicit shared-label assignments attain both bounds')
        independent_lo=totals[1]['H']/totals[1]['R']-(totals[0]['H']+totals[0]['U'])/totals[0]['R']
        independent_hi=(totals[1]['H']+totals[1]['U'])/totals[1]['R']-totals[0]['H']/totals[0]['R']
        require(independent_lo<=lo<=observed<=hi<=independent_hi, 'Independent-arm containment')
        bounds.append(dict(context=context,policy=policy,**{k+str(t):v for t in (0,1) for k,v in totals[t].items()},
                           observed_delta=observed,sharp_lower=lo,sharp_upper=hi,unknown_identities=len(unknown),
                           unknown_nonzero_coefficients=sum(c!=0 for c in unknown.values()),independent_lower=independent_lo,independent_upper=independent_hi))
        for partition,groups in (('donor',donor_groups),('role',role_groups)):
            summary,detail=decompose(groups,context,policy,partition)
            require(summary['delta']==observed,'Decomposition agrees with observed contrast')
            summaries.append(summary);details.extend(detail)
        for role,vals in role_groups.items():
            for t in (0,1):
                role_levels.append(dict(context=context,policy=policy,role=role,state=t,R=vals[t][0],H=vals[t][1],p=vals[t][1]/vals[t][0] if vals[t][0] else None,
                                        union_weight=vals[t][0]/totals[t]['R']))
    require(len(bounds)==5 and len(summaries)==10, 'All five context-policy rows retained')
    a.output.mkdir(parents=True)
    for name,rows in [('fraction_bounds.csv',bounds),('fraction_coefficients.csv',coefficients),('decomposition_summary.csv',summaries),('decomposition_units.csv',details),('role_levels.csv',role_levels)]:
        save(a.output/name,rows)
    manifest=dict(status='PASS',checks=CHECKS,design='Post-exposure, fixed-cohort descriptive analysis; no new resampling family',native_fits=0,
                  inputs=[bind(p) for p in (ledger_path,levels_path,donors_path)],code=bind(Path(__file__)),protocol=bind(Path(__file__).with_name('PROTOCOL.md')),
                  outputs=[bind(p) for p in sorted(a.output.iterdir())])
    with (a.output/'manifest.json').open('x',encoding='utf-8') as h: json.dump(manifest,h,indent=2)
    for r in bounds: print(r['context'],r['policy'],'observed pp',float(100*r['observed_delta']),'sharp pp',float(100*r['sharp_lower']),float(100*r['sharp_upper']))
    for r in summaries: print(r['context'],r['policy'],r['partition'],{k:round(float(100*r[k]),6) for k in ('within','reweight','entry','exit')})
    print('PASS checks',CHECKS)

if __name__=='__main__': main()
