"""Cancer-listed union allocation analysis with read-only inputs."""
import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np

FIELDS = ('R', 'H', 'Z', 'U', 'absent', 'explicit_null', 'identity_conflict', 'call_conflict')
ROLES = ('ONCOGENE', 'TSG', 'DUAL_ONCOGENE_TSG', 'OTHER_OR_UNSPECIFIED')
CONFIG = {
    'CRC': {'degree_allocations': 'results/crc_precision/degree_allocations.json',
            'allocations': [{'method': 'PRODIGY', 'path': 'results/crc_prodigy_allocations.json'},
                            {'method': 'DawnRank', 'path': 'results/crc_precision/allocations.json'}],
            'labels': ['data/crc/selected_labels.json'], 'categories': 'data/crc/categories.json'},
    'ESCA': {'allocations': 'results/esca_allocations/allocations.json',
             'labels': ['results/esca_evaluation/selected_labels.json'],
             'categories': 'results/esca_scope/scope.json'}}
checks = 0
inputs = set()


def require(ok, message):
    global checks
    checks += 1
    if not ok:
        raise AssertionError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    inputs.add(Path(path).resolve())
    with Path(path).open(encoding='utf-8-sig', newline='') as h:
        return list(csv.DictReader(h))


def write_csv(path, rows):
    encoded = []
    for row in rows:
        encoded.append({k: str(v) if isinstance(v, F) else v for k, v in row.items()})
    keys = list(dict.fromkeys(k for r in encoded for k in r))
    with path.open('x', encoding='utf-8', newline='') as h:
        w = csv.DictWriter(h, fieldnames=keys)
        w.writeheader()
        w.writerows(encoded)


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    return value


def json_save(path, value):
    with path.open('x', encoding='utf-8') as h:
        json.dump(serial(value), h, indent=2, allow_nan=False)


def zero_mass():
    return dict.fromkeys(FIELDS, F())


def quantile_mc_interval(maxima):
    """Exact binomial order-statistic interval for the unknown 0.95 quantile."""
    n, p = len(maxima), .95
    logp = np.array([math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
                    + k*math.log(p)+(n-k)*math.log1p(-p) for k in range(n+1)])
    probs = np.exp(logp-np.max(logp))
    probs /= probs.sum()
    cdf = probs.cumsum()
    lo = int(np.searchsorted(cdf, .025))
    hi = int(np.searchsorted(cdf, .975)) + 1
    ordered = np.sort(maxima)
    return {'lower_order_one_based': lo, 'upper_order_one_based': hi,
            'lower': float(ordered[max(0, lo-1)]), 'upper': float(ordered[min(n-1, hi-1)]),
            'binomial_coverage': float(cdf[hi-1] - (cdf[lo-1] if lo else 0))}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--archive', type=Path, required=True)
    ap.add_argument('--crisp-source', type=Path, required=True)
    ap.add_argument('--reference', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    out, archive, ref = a.output.resolve(), a.archive.resolve(), a.reference.resolve()
    require(not out.exists(), 'Output directory must be new')
    require(not out.is_relative_to(archive) and not out.is_relative_to(ref), 'Output outside input trees')
    sys.path.insert(0, str(a.crisp_source.resolve()))
    from crisp.statistics.inputs import load_crc, load_esca, load_categories
    from crisp.statistics.category_sensitivity import endpoint as ratio_endpoint, count_draws
    from crisp.statistics.paired import combine, bound, SEEDS
    for module in ('inputs.py', 'category_sensitivity.py', 'paired.py', 'endpoint.py'):
        inputs.add(a.crisp_source.resolve()/'crisp/statistics'/module)
    references = json.loads((ref/'Colorectal_Functional_Analysis/Code/INPUT_REFERENCES.json').read_text())
    inputs.add(ref/'Colorectal_Functional_Analysis/Code/INPUT_REFERENCES.json')
    for r in references['inputs']:
        if r['root'] == 'archive':
            require(sha(archive/r['path']) == r['sha256'], 'Accepted archive hash: '+r['path'])
    previous = read_csv(ref/'Category_Sensitivity/model_category_levels.csv')
    previous_donors = read_csv(ref/'Category_Sensitivity/donor_category_levels.csv')
    donor_order = read_csv(ref/'Category_Sensitivity/donor_index.csv')
    overview = read_csv(ref/'Cancer_Allocation_Overview/Figure_Overview_Source.csv')
    old_ledger = read_csv(ref/'functional_content/ledger.csv')
    previous_models = {}
    for row in previous:
        if row['category'] not in ROLES:
            continue
        key = (row['context'], row['policy'].replace('DawnRank_strict', 'DawnRank'), row['model_id'], row['state'])
        previous_models.setdefault(key, zero_mass())
        for f in FIELDS:
            previous_models[key][f] += F(row[f])
    previous_donor_union = {}
    for row in previous_donors:
        if row['category'] not in ROLES:
            continue
        key = (row['context'], row['policy'].replace('DawnRank_strict','DawnRank'), row['donor_id'], row['state'])
        previous_donor_union.setdefault(key, zero_mass())
        for f in FIELDS:
            previous_donor_union[key][f] += F(row[f])
    all_models, all_donors, all_contrasts, all_population = [], [], [], []
    all_endpoints, all_envelopes, all_loo, all_flows, all_concentration = [], [], [], [], []
    matrices, families, roster_rows, gene_rows = {}, {}, [], []
    for context in ('CRC', 'ESCA'):
        config = CONFIG[context]
        for k in ('degree_allocations', 'categories'):
            if k in config:
                inputs.add(archive/config[k])
        for s in config['labels']:
            inputs.add(archive/s)
        if context == 'CRC':
            for s in config['allocations']:
                inputs.add(archive/s['path'])
            models, allocations, means, degrees, calls, policies = load_crc(archive, config)
            policy_names = ('PRODIGY', 'DawnRank', 'degree')
            require(len(allocations) == 85*4*6, 'CRC complete four-cell five-seed plus deterministic population')
        else:
            inputs.add(archive/config['allocations'])
            models, policies, calls = load_esca(archive, config)
            policy_names = ('DawnRank', 'degree')
        cats = load_categories(archive/config['categories'], models, context)
        model_counts = Counter(m['donor_id'] for m in models)
        donors = [r['donor_id'] for r in sorted((r for r in donor_order if r['context']==context),
                                               key=lambda r:int(r['donor_index']))]
        require(donors == sorted(model_counts), 'Exact recorded donor order '+context)
        D = len(donors)
        expected = (85,83) if context=='CRC' else (59,58)
        require((len(models),D)==expected, 'Complete population '+context)
        draw_name = 'colorectal_donor_draws.csv' if context=='CRC' else 'oesophageal_donor_draws.csv'
        draw_rows = read_csv(ref/'Category_Sensitivity'/draw_name)
        draws = np.array([[int(r['position_'+str(i+1)]) for i in range(D)] for r in draw_rows])
        require([int(r['draw_index']) for r in draw_rows]==list(range(10000)), 'All ordered draws')
        counts = count_draws(draws, D)
        require(counts.shape==(10000,D) and np.all(counts.sum(axis=1)==D), 'Valid draw counts')
        donor_levels = {(p,d,t): zero_mass() for p in policy_names for d in donors for t in ('0','1')}
        vectors = {(p,d): defaultdict(F) for p in policy_names for d in donors}
        by_model = {}
        for model in models:
            mid, sample, donor = (model[k] for k in ('model_id','sample_id','donor_id'))
            roster_rows.append({'context':context, 'model_id':mid,'sample_id':sample,'donor_id':donor,
                                'models_in_donor':model_counts[donor]})
            wm = F(1,model_counts[donor])
            for p in policy_names:
                cells = ('G0','G1') if p=='degree' else ('G0H0','G1H1')
                pair = []
                for t,cell in zip(('0','1'),cells):
                    v = {g:w for g,w in policies[mid][p][cell].items() if cats['role'][g] in ROLES}
                    pair.append(v)
                    mass = zero_mass()
                    for g,w in v.items():
                        state,y = calls[sample,g]
                        mass['R'] += w
                        if state=='MEASURED':
                            mass['H' if y else 'Z'] += w
                        else:
                            mass['U'] += w
                            mass[{'ABSENT':'absent','EXPLICIT_NULL':'explicit_null',
                                  'IDENTITY_CONFLICT':'identity_conflict','CONFLICTING_CALLS':'call_conflict'}[state]] += w
                    require(mass['R']==mass['H']+mass['Z']+mass['U'], 'Model mass conservation')
                    require(mass==previous_models[context,p,mid,t], 'Independent accepted model-total reconciliation')
                    all_models.append({'context':context,'policy':p,'model_id':mid,'sample_id':sample,
                                       'donor_id':donor,'state':t,**mass})
                    for f in FIELDS:
                        donor_levels[p,donor,t][f] += wm*mass[f]
                vec = combine([(1,pair[1]),(-1,pair[0])])
                by_model[p,mid] = {(sample,g):c for g,c in vec.items()}
                for g in sorted(set(pair[0])|set(pair[1])):
                    old,new=pair[0].get(g,F()),pair[1].get(g,F())
                    state,y=calls[sample,g]
                    gene_rows.append({'context':context,'policy':p,'model_id':mid,'sample_id':sample,
                                      'donor_id':donor,'gene':g,'role':cats['role'][g],
                                      'label_state':state,'label_value':y,'old':old,'new':new,
                                      'difference':new-old,'donor_model_weight':wm,
                                      'population_weight':wm/D})
                    if context=='CRC' and p!='degree':
                        matches=[r for r in old_ledger if r['method'].replace('DawnRank_strict','DawnRank')==p
                                 and r['model_id']==mid and r['gene']==g]
                        require(len(matches)==1,'Complete independent CRC identity ledger')
                        r=matches[0]
                        require(F(r['old_weight'])==old and F(r['new_weight'])==new,
                                'Independent CRC identity allocation match')
                for g,c in vec.items():
                    vectors[p,donor][sample,g] += wm*c
        donor_bound = {}
        for p in policy_names:
            for d in donors:
                for t in ('0','1'):
                    mass=donor_levels[p,d,t]
                    require(mass==previous_donor_union[context,p,d,t], 'Independent donor-total reconciliation')
                    all_donors.append({'context':context,'policy':p,'donor_id':d,'state':t,
                                       'models_in_donor':model_counts[d],**mass})
                b=bound(vectors[p,d],calls)
                donor_bound[p,d]=b
                require(b['known']==donor_levels[p,d,'1']['H']-donor_levels[p,d,'0']['H'],
                        'Known donor contrast agrees with levels')
                all_contrasts.append({'context':context,'contrast':p,'donor_id':d,**b})
            for t in ('0','1'):
                mass={f:sum((donor_levels[p,d,t][f] for d in donors),F())/D for f in FIELDS}
                all_population.append({'context':context,'policy':p,'state':t,**mass,
                                       'pooled_positive_fraction':mass['H']/mass['R']})
                if context=='CRC' and p!='degree':
                    for field,endpoint in [('R','returned'),('H','known_positive')]:
                        match=[r for r in overview if r['method'].replace('DawnRank_strict','DawnRank')==p and r['state']==('G0H0' if t=='0' else 'G1H1')
                               and r['category']=='Cancer-listed' and r['endpoint']==endpoint]
                        require(len(match)==1 and F(match[0]['exact'])==mass[field], 'CRC overview exact match')
        for p in policy_names[:-1]:
            label=p+'_minus_degree'
            for d in donors:
                vectors[label,d]=combine([(1,vectors[p,d]),(-1,vectors['degree',d])])
                b=bound(vectors[label,d],calls)
                donor_bound[label,d]=b
                all_contrasts.append({'context':context,'contrast':label,'donor_id':d,**b})
        # Recompute joint labels before bounding for selected resamples (including duplicates).
        for p in tuple(policy_names)+tuple(p+'_minus_degree' for p in policy_names[:-1]):
            for idx in (0,1,4999,9999):
                aggregate=combine([(F(int(counts[idx,i]),D),vectors[p,d]) for i,d in enumerate(donors)])
                exact=bound(aggregate,calls)
                for f in ('known','lower','upper'):
                    linear=sum((F(int(counts[idx,i]),D)*donor_bound[p,d][f] for i,d in enumerate(donors)),F())
                    require(exact[f]==linear,'Sharp bound remains exact under donor resampling')
        records, sampled, loo_values = [],[],[]
        def add_linear(p,name,x):
            value=sum(x,F())/D
            centred=[v-value for v in x]
            variance=sum((v*v for v in centred),F())/(D*(D-1))
            scale=math.sqrt(float(variance))
            vals=counts@np.array([float(v) for v in x])/D
            rec={'context':context,'contrast':p,'endpoint':name,'estimate':value,
                 'estimate_decimal':float(value),'scale_squared':variance,'scale':scale,
                 'constant':variance==0,'influences':[str(v) for v in centred]}
            records.append(rec);sampled.append(vals)
            loo_values.append([(sum(x,F())-v)/(D-1) for v in x])
        for p in policy_names:
            old=[donor_levels[p,d,'0'] for d in donors]
            new=[donor_levels[p,d,'1'] for d in donors]
            add_linear(p,'delta_R',[b['R']-a['R'] for a,b in zip(old,new)])
            add_linear(p,'delta_H',[b['H']-a['H'] for a,b in zip(old,new)])
            arr=[[r[f] for r in rows] for rows,f in ((old,'R'),(old,'H'),(new,'R'),(new,'H'))]
            e=ratio_endpoint(*arr,'delta_p')
            require(e['status'] in ('AVAILABLE','CONSTANT'),'Nondegenerate ratio linearisation')
            totals=counts@np.array([[float(v) for v in x] for x in arr]).T
            require(np.all(totals[:,[0,2]]>0),'No undefined pooled denominator in stored draws')
            vals=totals[:,3]/totals[:,2]-totals[:,1]/totals[:,0]
            records.append({'context':context,'contrast':p,'endpoint':'delta_H_over_R',
                            'estimate':F(e['value']),'estimate_decimal':e['value_float'],
                            'scale_squared':F(e['scale_squared']),'scale':e['scale'],
                            'constant':e['constant_certificate'],'influences':e['influences']})
            sampled.append(vals)
            sums=[sum(x,F()) for x in arr]
            ratio_loo=[]
            for i in range(D):
                v=[s-x[i] for s,x in zip(sums,arr)]
                require(v[0]>0 and v[2]>0,'Positive leave-one-out pooled denominators')
                ratio_loo.append(v[3]/v[2]-v[1]/v[0])
            loo_values.append(ratio_loo)
            for f in ('lower','upper'):
                add_linear(p,f,[donor_bound[p,d][f] for d in donors])
        for p in policy_names[:-1]:
            label=p+'_minus_degree'
            for f in ('lower','upper'):
                add_linear(label,f,[donor_bound[label,d][f] for d in donors])
        M=np.array(sampled).T
        require(M.shape==(10000,19 if context=='CRC' else 12),'Complete fixed analysis family')
        varying=[]
        for j,r in enumerate(records):
            if r['constant']:
                require(np.max(np.abs(M[:,j]-r['estimate_decimal']))<1e-12,'Constant coordinate certificate')
            else:
                require(r['scale']>0,'Positive varying scale')
                varying.append(j)
        maxes=np.max(np.abs(M[:,varying]-np.array([records[j]['estimate_decimal'] for j in varying]))/
                     np.array([records[j]['scale'] for j in varying]),axis=1)
        q=float(np.quantile(maxes,.95,method='linear'))
        family={'context':context,'models':len(models),'donors':D,'draws':10000,'coordinates':len(records),
                'constant_coordinates':len(records)-len(varying),'q':q,
                'q_monte_carlo_interval':quantile_mc_interval(maxes),
                'donor_order':donors,'endpoints':records}
        for j,r in enumerate(records):
            r['q']=q;r['envelope_lower']=r['estimate_decimal']-q*r['scale']
            r['envelope_upper']=r['estimate_decimal']+q*r['scale']
            r['direction']='positive' if r['envelope_lower']>0 else 'negative' if r['envelope_upper']<0 else 'crosses_zero'
            all_endpoints.append({k:v for k,v in r.items() if k!='influences'})
            for d,value in zip(donors,loo_values[j]):
                all_loo.append({'context':context,'contrast':r['contrast'],'endpoint':r['endpoint'],
                                'omitted_donor':d,'estimate':value,'estimate_decimal':float(value),
                                'change_from_full':value-r['estimate']})
        families[context]=family
        matrices[context]=[{'draw_index':i,**{r['contrast']+'__'+r['endpoint']:float(M[i,j])
                          for j,r in enumerate(records)},'maximum_standardised_deviation':float(maxes[i])}
                           for i in range(10000)]
        for p in tuple(policy_names)+tuple(p+'_minus_degree' for p in policy_names[:-1]):
            lo=next(r for r in records if r['contrast']==p and r['endpoint']=='lower')
            hi=next(r for r in records if r['contrast']==p and r['endpoint']=='upper')
            all_envelopes.append({'context':context,'contrast':p,'sharp_lower':lo['estimate'],
                                  'sharp_upper':hi['estimate'],'outer_lower':lo['envelope_lower'],
                                  'outer_upper':hi['envelope_upper'],'q':q})
        for p in policy_names:
            b=[donor_bound[p,d] for d in donors]
            positive=sorted([(d,r['known']) for d,r in zip(donors,b) if r['known']>0],key=lambda x:(-x[1],x[0]))
            negative=sorted([(d,-r['known']) for d,r in zip(donors,b) if r['known']<0],key=lambda x:(-x[1],x[0]))
            posmass=sum((v for d,v in positive),F());negmass=sum((v for d,v in negative),F())
            summary={'context':context,'policy':p,'positive_donors':len(positive),'negative_donors':len(negative),
                     'zero_donors':sum(r['known']==0 for r in b),'gross_positive_mass':posmass/D,
                     'gross_negative_mass':negmass/D,'net_mass':(posmass-negmass)/D,
                     'largest_positive_donor':positive[0][0] if positive else None,
                     'largest_positive_donor_contribution':positive[0][1]/D if positive else F(),
                     'largest_negative_donor':negative[0][0] if negative else None,
                     'largest_negative_donor_contribution':-negative[0][1]/D if negative else F(),
                     'top_three_positive_fraction':sum((v for d,v in positive[:3]),F())/posmass if posmass else None}
            all_concentration.append(summary)
    out.mkdir(parents=True)
    tables={'model_levels.csv':all_models,'donor_levels.csv':all_donors,'donor_contrasts.csv':all_contrasts,
            'population_levels.csv':all_population,'endpoints.csv':all_endpoints,
            'bound_envelopes.csv':all_envelopes,'leave_one_donor_out.csv':all_loo,
            'donor_concentration.csv':all_concentration,'model_roster.csv':roster_rows,
            'model_gene_ledger.csv':gene_rows,**{c+'_resamples.csv':v for c,v in matrices.items()}}
    for name,rows in tables.items():
        write_csv(out/name,rows)
    json_save(out/'family_details.json',families)
    manifest={'status':'PASS','checks':checks,'design':'Post-exposure finite exploratory union analysis',
              'native_fits':0,'rerankings':0,'prodigy_seeds':SEEDS,
              'python':platform.python_version(),'numpy':np.__version__,
              'protocol_sha256':sha(Path(__file__).with_name('PROTOCOL.md')),
              'code_sha256':sha(__file__),
              'source_files':[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(inputs)],
              'outputs':{p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(out.iterdir())},
              'population_levels':all_population,'bound_envelopes':all_envelopes,
              'donor_concentration':all_concentration,
              'family_summaries':{c:{k:v for k,v in f.items() if k not in ('donor_order','endpoints')} for c,f in families.items()}}
    json_save(out/'manifest.json',manifest)
    print(json.dumps(serial({'checks':checks,'families':manifest['family_summaries'],
                             'levels':all_population,'endpoints':all_endpoints,
                             'bounds':all_envelopes,'donors':all_concentration}),indent=2))


if __name__=='__main__':
    main()
