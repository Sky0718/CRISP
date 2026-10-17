"""Independent checks from previously accepted category tables and coefficient ledgers."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import numpy as np

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--reference',type=Path,required=True)
ap.add_argument('--results',type=Path,required=True)
a=ap.parse_args()
checks=0

def req(v,msg):
    global checks
    checks+=1
    if not v: raise AssertionError(msg)

def csvread(p):
    with p.open(encoding='utf-8-sig',newline='') as h:return list(csv.DictReader(h))

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def near(x,y,msg):req(abs(float(x)-float(y))<2e-12,msg)

old=csvread(a.reference/'Category_Sensitivity/donor_category_levels.csv')
roles={'ONCOGENE','TSG','DUAL_ONCOGENE_TSG','OTHER_OR_UNSPECIFIED'}
levels={}
for r in old:
    if r['category'] not in roles:continue
    key=(r['context'],r['policy'].replace('DawnRank_strict','DawnRank'),r['donor_id'],r['state'])
    levels.setdefault(key,{f:F() for f in ['R','H','Z','U']})
    for f in ['R','H','Z','U']:levels[key][f]+=F(r[f])
new=csvread(a.results/'donor_levels.csv')
for r in new:
    expected=levels[r['context'],r['policy'],r['donor_id'],r['state']]
    req(all(F(r[f])==expected[f] for f in expected),'Independent role-sum mass')

contrast_rows=csvread(a.results/'donor_contrasts.csv')
contrasts={(r['context'],r['contrast'],r['donor_id']):r for r in contrast_rows}
# A different pre-existing ESCA generator supplies complete selected model-gene contrast records.
esca=csvread(a.reference/'Oesophageal_Identity_Analysis/model_gene_coefficients.csv')
old_esca=defaultdict(lambda:dict.fromkeys(('known','lower','upper'),F()))
aliases={'Dawn_joint':'DawnRank','degree_update':'degree','Theta':'DawnRank_minus_degree'}
D=58
for r in esca:
    if r['role'] not in roles:continue
    key=('ESCA',aliases[r['contrast']],r['donor_id'])
    w=F(r['population_coefficient'])*D
    v=old_esca[key]
    if r['label_state'] in ('MEASURED_POSITIVE','MEASURED_NEGATIVE'):
        k=w*int(r['label_value'])
        for f in v:v[f]+=k
    else:
        v['lower']+=min(w,0);v['upper']+=max(w,0)
for key,r in contrasts.items():
    if key[0]!='ESCA':continue
    req(all(F(r[f])==old_esca[key][f] for f in ('known','lower','upper')),
        'Independent accepted ESCA model-gene coefficients')
# CRC independent accepted identity ledger is already seed-averaged and records exact donor alpha.
old_crc=defaultdict(lambda:dict.fromkeys(('known','lower','upper'),F()))
for r in csvread(a.reference/'functional_content/ledger.csv'):
    if r['role'] not in roles:continue
    p=r['method'].replace('DawnRank_strict','DawnRank')
    v=old_crc['CRC',p,r['donor_id']]
    w=F(r['coefficient'])*83
    if r['state']=='MEASURED':
        k=w*int(r['value'])
        for f in v:v[f]+=k
    else:
        v['lower']+=min(w,0);v['upper']+=max(w,0)
for key,r in contrasts.items():
    if key[0]=='CRC' and key[1] in ('PRODIGY','DawnRank'):
        req(all(F(r[f])==old_crc[key][f] for f in ('known','lower','upper')),
            'Independent accepted CRC identity coefficients')

endpoints=csvread(a.results/'endpoints.csv')
details=json.loads((a.results/'family_details.json').read_text())
leaveouts=csvread(a.results/'leave_one_donor_out.csv')
loo={(r['context'],r['contrast'],r['endpoint'],r['omitted_donor']):r for r in leaveouts}
comparisons=[]
for context,D,name in [('CRC',83,'colorectal_donor_draws.csv'),('ESCA',58,'oesophageal_donor_draws.csv')]:
    donors=details[context]['donor_order']
    dr=csvread(a.reference/'Category_Sensitivity'/name)
    draws=np.array([[int(r['position_'+str(i+1)]) for i in range(D)] for r in dr])
    recorded=csvread(a.results/(context+'_resamples.csv'))
    records=[r for r in endpoints if r['context']==context]
    sampled=[];scales=[];estimates=[]
    for r in records:
        p,e=r['contrast'],r['endpoint']
        if e=='delta_H_over_R':
            ar=np.array([[float(levels[context,p,d,t][f]) for d in donors]
                         for t,f in [('0','R'),('0','H'),('1','R'),('1','H')]])
            means=ar.mean(axis=1)
            p0,p1=means[1]/means[0],means[3]/means[2]
            est=p1-p0
            infl=(ar[3]-p1*ar[2])/means[2]-(ar[1]-p0*ar[0])/means[0]
            scale=math.sqrt(np.sum(infl**2)/(D*(D-1)))
            # Direct index-and-average evaluation, rather than multiplication by count matrix.
            sums=np.array([v[draws].mean(axis=1) for v in ar])
            values=sums[3]/sums[2]-sums[1]/sums[0]
            for i,d in enumerate(donors):
                left=(ar.sum(axis=1)-ar[:,i])/(D-1)
                actual=left[3]/left[2]-left[1]/left[0]
                near(actual,F(loo[context,p,e,d]['estimate']),'Independent leave-one-donor ratio')
        else:
            if e in ('delta_R','delta_H'):
                field=e[-1]
                x=np.array([float(levels[context,p,d,'1'][field]-levels[context,p,d,'0'][field]) for d in donors])
            else:
                x=np.array([float(F(contrasts[context,p,d][e])) for d in donors])
            est=x.mean();scale=x.std(ddof=1)/math.sqrt(D)
            values=x[draws].mean(axis=1)
            for i,d in enumerate(donors):
                near((x.sum()-x[i])/(D-1),F(loo[context,p,e,d]['estimate']),'Independent leave-one-donor linear')
        near(est,r['estimate_decimal'],'Independent endpoint')
        near(scale,r['scale'],'Independent fixed donor scale')
        truth=np.array([float(row[p+'__'+e]) for row in recorded])
        req(np.max(np.abs(values-truth))<2e-12,'All stored resamples independent recomputation')
        sampled.append(values);scales.append(scale);estimates.append(est)
    maxima=np.max(abs(np.array(sampled).T-np.array(estimates))/np.array(scales),axis=1)
    q=float(np.quantile(maxima,.95,method='linear'))
    near(q,details[context]['q'],'Independent family calibration')
    for r in records:
        near(float(r['estimate_decimal'])-q*float(r['scale']),r['envelope_lower'],'Independent lower envelope')
        near(float(r['estimate_decimal'])+q*float(r['scale']),r['envelope_upper'],'Independent upper envelope')
    comparisons.append({'context':context,'coordinates':len(records),'draws':len(draws),'q':q})

manifest=json.loads((a.results/'manifest.json').read_text())
for r in manifest['source_files']:
    req(sha(Path(r['path']))==r['sha256'],'Unchanged input '+r['path'])
for name,r in manifest['outputs'].items():
    req(sha(a.results/name)==r['sha256'],'Unchanged recorded output '+name)
report={'status':'PASS','independent_checks':checks,'methods':['Accepted disjoint-role donor sums',
        'Independent accepted CRC and ESCA model-gene signed coefficients',
        'Direct indexing of 10,000 donor draws, not count-matrix implementation',
        'Independent standard-deviation/delta-method scales, q and leave-one-donor values',
        'All bound and input/output hashes'], 'families':comparisons,
        'verification_script_sha256':sha(Path(__file__))}
with (a.results/'independent_verification.json').open('x',encoding='utf-8') as h:json.dump(report,h,indent=2)
print(json.dumps(report,indent=2))
