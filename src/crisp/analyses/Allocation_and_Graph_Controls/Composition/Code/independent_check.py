"""Independent reconciliation from literal identities and derived summaries."""
import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--source',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))

F=Fraction
ledger=read(a.source/'functional_content/ledger.csv')
donors=sorted({r['donor_id'] for r in ledger})
expected=read(a.source/'Category_Sensitivity/donor_category_levels.csv')
mapping={r['model_id']:r['donor_id'] for r in ledger}
model_counts=Counter(mapping.values())
assert len(mapping)==85 and len(donors)==83
totals=defaultdict(lambda:[F(0),F(0)])
for r in ledger:
    method=r['method'].replace('_strict','')
    for c in ('BOTH_CORE','ONCOGENE'):
        if (c=='BOTH_CORE' and r['core']!=c) or (c=='ONCOGENE' and r['role']!=c):
            continue
        assert r['state']=='MEASURED'
        label=F(r['value'])
        for t,w in enumerate(('old_weight','new_weight')):
            weight=F(r[w])/model_counts[r['donor_id']]
            z=totals[method,c,r['donor_id'],t]
            z[0]+=weight
            z[1]+=weight*label
count=0
for r in expected:
    if r['context']!='CRC' or r['policy']=='degree' or r['category'] not in ('BOTH_CORE','ONCOGENE'):
        continue
    key=r['policy'].replace('_strict',''),r['category'],r['donor_id'],int(r['state'])
    assert totals[key]==[F(r['R']),F(r['H'])]
    count+=1
assert count==664
detail=read(a.output/'ratio_decomposition_donors.csv')
summary=read(a.output/'ratio_decomposition_summary.csv')
common=[]
checks=[]
for method in ('PRODIGY','DawnRank'):
    for c in ('BOTH_CORE','ONCOGENE'):
        support=Counter()
        changes=[]
        mass=[sum((totals[method,c,d,t][0] for d in donors),F(0)) for t in (0,1)]
        pos=[sum((totals[method,c,d,t][1] for d in donors),F(0)) for t in (0,1)]
        p0,p1=pos[0]/mass[0],pos[1]/mass[1]
        for d in donors:
            r0,h0=totals[method,c,d,0]
            r1,h1=totals[method,c,d,1]
            group='both' if r0 and r1 else 'exits' if r0 else 'entrants' if r1 else 'neither'
            support[group]+=1
            if group=='both':
                changes.append(h1/r1-h0/r0)
                if method=='PRODIGY' and c=='BOTH_CORE':
                    common.append({'donor_id':d,'R0':str(r0),'R1':str(r1),'p0':str(h0/r0),'p1':str(h1/r1)})
        got=next(r for r in summary if r['policy']==method and r['category']==c and r['measure']=='binary_fraction')
        assert F(got['ratio0'])==p0 and F(got['ratio1'])==p1 and F(got['difference'])==p1-p0
        assert all(int(got[k])==support[k] for k in ('both','exits','entrants','neither'))
        assert sum((F(got[k]) for k in ('within_common','reweight_common','entrant','exit')),F(0))==p1-p0
        checks.append({'method':method,'category':c,'support':dict(support),'within_common_positive':sum(x>0 for x in changes),
                       'within_common_negative':sum(x<0 for x in changes),'within_common_exact_zero':sum(x==0 for x in changes)})
assert len(common)==14 and all(r['p0']==r['p1'] for r in common)

# Reconstruct the whole continuous family from exported donor sufficient
# statistics using direct donor indexing, sample-variance scales and sorting.
cont=read(a.output/'donor_continuous_levels.csv')
cv={(r['policy'],r['category'],r['donor_id'],int(r['state'])):(float(F(r['R'])),float(F(r['S']))) for r in cont}
idx=read(a.source/'Category_Sensitivity/donor_index.csv')
order=[r['donor_id'] for r in idx if r['context']=='CRC']
dr=read(a.source/'Category_Sensitivity/colorectal_donor_draws.csv')
draws=np.array([[int(r['position_'+str(i+1)]) for i in range(83)] for r in dr])
endpoints=read(a.output/'continuous_composition_endpoints.csv')
replicates=[]
for r in endpoints:
    states=[np.array([cv[r['policy'],r['category'],d,t] for d in order]) for t in (0,1)]
    r0,s0=states[0].mean(axis=0)
    r1,s1=states[1].mean(axis=0)
    if r['endpoint']=='Delta_S':
        est=s1-s0
        influences=states[1][:,1]-states[0][:,1]
        boot=(states[1][:,1]-states[0][:,1])[draws].mean(axis=1)
    else:
        est=s1/r1-s0/r0
        influences=(states[1][:,1]-s1/r1*states[1][:,0])/r1-(states[0][:,1]-s0/r0*states[0][:,0])/r0
        bootmeans=[st[draws].mean(axis=1) for st in states]
        boot=bootmeans[1][:,1]/bootmeans[1][:,0]-bootmeans[0][:,1]/bootmeans[0][:,0]
    scale=np.std(influences,ddof=1)/np.sqrt(83)
    assert abs(est-float(F(r['estimate'])))<2e-12
    assert abs(scale-float(r['scale']))<2e-12
    replicates.append(boot)
replicates=np.array(replicates).T
saved=np.load(a.output/'continuous_resampled_endpoints.npy')
assert np.max(abs(replicates-saved))<2e-12
est=np.array([float(F(r['estimate'])) for r in endpoints])
scale=np.array([float(r['scale']) for r in endpoints])
maximum=np.max(abs(replicates-est)/scale,axis=1)
sorted_max=np.sort(maximum)
# Equivalent independent explicit linear interpolation for n=10,000, p=.95.
q=.95*sorted_max[9499]+.05*sorted_max[9500]
assert abs(q-float(endpoints[0]['critical_value']))<2e-12
q_lower,q_upper=sorted_max[9456],sorted_max[9542]
coverage=sum(math.exp(math.lgamma(10001)-math.lgamma(k+1)-math.lgamma(10001-k)+k*math.log(.95)+(10000-k)*math.log(.05)) for k in range(9457,9543))
mc_endpoints=[]
for r in endpoints:
    theta=float(F(r['estimate']))
    scale=float(r['scale'])
    assert abs(float(r['lower'])-(theta-float(q)*scale))<2e-12
    assert abs(float(r['upper'])-(theta+float(q)*scale))<2e-12
    def classification(critical):
        return 'negative' if theta+critical*scale<0 else 'positive' if theta-critical*scale>0 else 'crosses_zero'
    assert classification(q_lower)==classification(q_upper)==r['direction']
    mc_endpoints.append({'policy':r['policy'],'category':r['category'],'endpoint':r['endpoint'],
                        'outer_lower_using_q_upper':theta-float(q_upper)*scale,
                        'outer_upper_using_q_upper':theta+float(q_upper)*scale})
report={'status':'PASS','literal_donor_category_states_checked':count,
        'binary_support_checks':checks,'prodigy_shared_core_common_donors':common,
        'continuous_endpoints_checked':12,'continuous_resampled_values_checked':120000,
        'maximum_resample_difference':float(np.max(abs(replicates-saved))),
        'critical_value':float(q),'new_fits':0,
        'monte_carlo_audit':{'order_statistics_one_based':[9457,9543],
                            'coverage':coverage,'quantile_lower':float(q_lower),'quantile_upper':float(q_upper),
                            'outer_endpoints':mc_endpoints}}
(a.output/'INDEPENDENT_CHECK.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
