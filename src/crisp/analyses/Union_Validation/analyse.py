"""Finite analysis extensions from existing immutable study records.

Python 3.10+, NumPy. See PROTOCOL.md. No native model fitting or reranking.
"""
import argparse
import csv
import hashlib
import itertools
import json
import shutil
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

import numpy as np

POLICIES = ('PRODIGY', 'DawnRank', 'degree')
ROLES = ('ONCOGENE', 'TSG', 'DUAL_ONCOGENE_TSG', 'OTHER_OR_UNSPECIFIED')
CHECKS = 0


def require(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise AssertionError(message)


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as h:
        return list(csv.DictReader(h))


def write(path, rows):
    encoded = []
    for row in rows:
        out = {}
        for k, v in row.items():
            out[k] = str(v) if isinstance(v, F) else v
            if isinstance(v, F):
                out[k+'_decimal'] = float(v)
        encoded.append(out)
    fields = list(dict.fromkeys(k for row in encoded for k in row))
    with path.open('x', encoding='utf-8', newline='') as h:
        w = csv.DictWriter(h, fieldnames=fields)
        w.writeheader()
        w.writerows(encoded)


def digest(path):
    return dict(bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    source, output = a.source.resolve(), a.output.resolve()
    require(not output.exists(), 'Use fresh output directory')
    require(not output.is_relative_to(source), 'Never write into immutable source')
    paths = {
        'union_ledger.csv': 'Union_Sensitivity_Analysis/inputs/union/model_gene_ledger.csv',
        'union_levels.csv': 'Union_Sensitivity_Analysis/inputs/union/population_levels.csv',
        'union_donor_levels.csv': 'Union_Sensitivity_Analysis/inputs/union/donor_levels.csv',
        'continuous_identities.csv': 'Colorectal_Functional_Analysis/current_continuous_selected_identities.csv',
        'donor_index.csv': 'Category_Sensitivity/donor_index.csv',
        'accepted_fraction_bounds.csv': 'Union_Sensitivity_Analysis/expected/results/fraction_bounds.csv',
        'accepted_role_decomposition.csv': 'Union_Sensitivity_Analysis/expected/results/decomposition_units.csv',
    }
    inputs = {k: read(source/v) for k, v in paths.items()}
    ledger, levels = inputs['union_ledger.csv'], inputs['union_levels.csv']
    continuous = {(r['sample_id'], r['gene']): r for r in inputs['continuous_identities.csv']}
    require(len(continuous) == len(inputs['continuous_identities.csv']) == 2287, '2287 unique continuous source identities')
    lfc = {}
    for key, r in continuous.items():
        value = F(r['supplied_LFC']) if r['supplied_LFC'] else None
        if value is not None:
            raw = json.loads(r['raw_representations'])
            require(bool(raw), 'LFC has retained source rows')
            for v in raw:
                require(v[3] == key[0] and v[4] == r['retained_library'], 'LFC sample-library binding')
                require(F(v[5]) == value, 'LFC agrees with literal raw decimal')
        lfc[key] = value
    donors = tuple(r['donor_id'] for r in inputs['donor_index.csv'] if r['context'] == 'CRC')
    require(len(donors) == len(set(donors)) == 83, '83 unique CRC donors')
    crc_models = {r['model_id']: r['donor_id'] for r in continuous.values()}
    donor_counts = Counter(crc_models.values())
    require(len(crc_models) == 85 and set(donor_counts) == set(donors), 'Complete CRC 85/83 identity roster')
    for r in ledger:
        require(r['role'] in ROLES, 'Only cancer-listed union roles')
        require(F(r['difference']) == F(r['new'])-F(r['old']), 'Immutable allocation difference')
        if r['context'] == 'CRC':
            require(crc_models[r['model_id']] == r['donor_id'], 'CRC model/donor ledger link')
            require(F(r['population_weight']) == F(1, 83*donor_counts[r['donor_id']]), 'CRC equal donor/model weights')
            sr = continuous[r['sample_id'], r['gene']]
            require(sr['model_id'] == r['model_id'] and sr['donor_id'] == r['donor_id'], 'Continuous source exact identity')
            require(sr['role'] == r['role'], 'Continuous cancer role matches union')

    # A: continuous sums and ratios, keeping every donor including zero allocation.
    dvals = {(pol, d, t): dict.fromkeys(('R', 'A', 'S'), F())
             for pol in POLICIES for d in donors for t in (0, 1)}
    mvals = {(pol, m, t): dict.fromkeys(('R', 'A', 'S'), F())
             for pol in POLICIES for m in crc_models for t in (0, 1)}
    identities = []
    for r in ledger:
        if r['context'] != 'CRC':
            continue
        pol, mid, d = r['policy'], r['model_id'], r['donor_id']
        value = lfc[r['sample_id'], r['gene']]
        within = F(r['donor_model_weight'])
        require(within == F(1, donor_counts[d]), 'Within-donor model weight')
        for t, field in enumerate(('old', 'new')):
            w = F(r[field])
            for group, factor in ((mvals[pol, mid, t], F(1)), (dvals[pol, d, t], within)):
                group['R'] += factor*w
                if value is not None:
                    group['A'] += factor*w
                    group['S'] += factor*w*(-value)
        identities.append({k: r[k] for k in ('policy','model_id','sample_id','donor_id','gene','role','old','new','population_weight')} | {
            'source_LFC': value, 'depletion': -value if value is not None else None,
            'continuous_state': continuous[r['sample_id'],r['gene']]['continuous_state'],
            'retained_library': continuous[r['sample_id'],r['gene']]['retained_library'],
            'source_rows': continuous[r['sample_id'],r['gene']]['source_rows']})
    def expand(v):
        return v | {'missing':v['R']-v['A'], 'coverage':v['A']/v['R'] if v['R'] else None,
                    'mean_depletion':v['S']/v['A'] if v['A'] else None,
                    'complete_return_sum': v['S'] if v['R']==v['A'] else None}
    population = []
    for pol in POLICIES:
        for t in (0,1):
            v = {k:sum((dvals[pol,d,t][k] for d in donors),F())/83 for k in ('R','A','S')}
            expected = next(r for r in levels if (r['context'],r['policy'],r['state']) == ('CRC',pol,str(t)))
            require(v['R'] == F(expected['R']), 'Continuous union R matches accepted binary denominator')
            population.append(dict(policy=pol,state=t,**expand(v)))
    family = [(pol,metric) for pol in POLICIES for metric in ('R','A','coverage','S','mean_depletion')]
    states = np.asarray([[[[float(dvals[pol,d,t][k]) for k in ('R','A','S')]
                          for t in (0,1)] for d in donors] for pol in POLICIES])
    rng = np.random.Generator(np.random.PCG64(470016))
    draw_ids = rng.integers(0,83,size=(10000,83))
    boot_states = states[:,draw_ids,:,:].mean(axis=2)
    def endpoints(values):
        r, assessed, depletion = (values[...,i] for i in range(3))
        require(bool(np.all(assessed > 0) and np.all(r > 0)), 'Positive resampled union denominators')
        five = np.stack((r,assessed,assessed/r,depletion,depletion/assessed),axis=-1)
        return (five[...,1,:]-five[...,0,:])
    boot = endpoints(boot_states).transpose(1,0,2).reshape(10000,15)
    point = endpoints(states.mean(axis=1)).reshape(15)
    sd = boot.std(axis=0,ddof=1)
    active = sd > 0
    standard = np.zeros_like(boot)
    standard[:,active] = (boot[:,active]-point[active])/sd[active]
    maxima = np.max(np.abs(standard),axis=1)
    multiplier = float(np.quantile(maxima,.95,method='linear'))
    endpoints_rows = [dict(policy=pol,metric=metric,estimate=point[i],bootstrap_sd=sd[i],
                           lower=point[i]-multiplier*sd[i],upper=point[i]+multiplier*sd[i],
                           family_size=15,draws=10000,seed=470016,multiplier=multiplier)
                      for i,(pol,metric) in enumerate(family)]

    # B: maximum gain after k flips with exact coefficients, independent of their order in input.
    bounds, gain_rows, witnesses, curves, breakdowns, role_levels, role_terms = [], [], [], [], [], [], []
    accepted = {(r['context'],r['policy']):r for r in inputs['accepted_fraction_bounds.csv']}
    for context, pol in dict.fromkeys((r['context'],r['policy']) for r in levels):
        records = [r for r in ledger if (r['context'],r['policy']) == (context,pol)]
        totals = [{k:F() for k in ('R','H','Z','U')} for _ in (0,1)]
        roles = {(role,t):dict.fromkeys(('R','H','Z','U'),F()) for role in ROLES for t in (0,1)}
        calls, coeff, metadata = {}, defaultdict(F), {}
        for r in records:
            for t,k in enumerate(('old','new')):
                w = F(r['population_weight'])*F(r[k])
                field = ('H' if int(r['label_value']) else 'Z') if r['label_state']=='MEASURED' else 'U'
                totals[t]['R'] += w; totals[t][field] += w
                roles[r['role'],t]['R'] += w; roles[r['role'],t][field] += w
        for r in records:
            identity = r['sample_id'],r['gene']
            call = r['label_state'],r['label_value']
            require(identity not in calls or calls[identity] == call, 'One binary call per immutable identity')
            calls[identity] = call
            metadata[identity] = r
            coeff[identity] += F(r['population_weight'])*(F(r['new'])/totals[1]['R']-F(r['old'])/totals[0]['R'])
        upper = F()
        positive = []
        measured = 0
        for key,c in coeff.items():
            state,y = calls[key]
            if state == 'MEASURED':
                measured += 1
                upper += c*int(y)
                gain = c*(1-2*int(y))
                r = metadata[key]
                gain_rows.append(dict(context=context,policy=pol,sample_id=key[0],gene=key[1],
                                      model_id=r['model_id'],donor_id=r['donor_id'],role=r['role'],
                                      original_label=y,coefficient=c,gain=gain))
                if gain > 0:
                    positive.append((gain,key))
            else:
                upper += max(c,F())
        require(upper == F(accepted[context,pol]['sharp_upper']), 'Accepted exact upper bound reconciled')
        positive.sort(key=lambda z:(-z[0],z[1]))
        cumulative = F()
        minimal = 0 if upper >= 0 else None
        threshold_before, threshold_after = upper, upper
        chosen = []
        curves.append(dict(context=context,policy=pol,flips=0,gain=F(),upper=upper))
        for k,(gain,key) in enumerate(positive,1):
            before = upper+cumulative
            cumulative += gain
            curves.append(dict(context=context,policy=pol,flips=k,gain=cumulative,upper=upper+cumulative))
            if minimal is None and upper+cumulative >= 0:
                minimal = k; threshold_before = before; threshold_after = upper+cumulative
        require(minimal is not None, 'Finite reachability retained; cannot silently drop unreachable case')
        chosen = positive[:minimal]
        if minimal:
            require(threshold_before < 0 <= threshold_after, 'Exact minimality threshold')
        counts = Counter()
        for rank,(gain,key) in enumerate(chosen,1):
            r = metadata[key]; state,y = calls[key]
            witnesses.append(dict(context=context,policy=pol,rank=rank,sample_id=key[0],gene=key[1],
                                  model_id=r['model_id'],donor_id=r['donor_id'],role=r['role'],
                                  original_label=y,flipped_label=1-int(y),coefficient=coeff[key],gain=gain))
            for dimension,value in [('role',r['role']),('donor',r['donor_id']),('direction',y+'->'+str(1-int(y)))]:
                counts[dimension,value] += 1
        for (dimension,value),count in counts.items():
            breakdowns.append(dict(context=context,policy=pol,dimension=dimension,value=value,count=count))
        # Independent two-state reconstruction of attained bound after witness flips.
        chosen_ids = {key for gain,key in chosen}
        complete = [F(),F()]
        for r in records:
            key=(r['sample_id'],r['gene']); state,y=calls[key]
            label = (1-int(y) if key in chosen_ids else int(y)) if state=='MEASURED' else int(coeff[key]>0)
            for t,k in enumerate(('old','new')):
                complete[t] += F(r['population_weight'])*F(r[k])*label
        attained = complete[1]/totals[1]['R']-complete[0]/totals[0]['R']
        require(attained == threshold_after, 'Direct state-level witness reaches exact bound')
        max_alpha = max(F(r['population_weight']) for r in records)
        loose_cap=max_alpha*max(1/totals[0]['R'],1/totals[1]['R'])
        conservative = 0 if upper >= 0 else (-upper/loose_cap).__ceil__()
        require(minimal >= conservative, 'Exact count respects universal conservative lower bound')
        bounds.append(dict(context=context,policy=pol,R0=totals[0]['R'],R1=totals[1]['R'],upper_before=upper,
                           measured_identities=measured,positive_gain_identities=len(positive),minimal_flips=minimal,
                           fraction_measured=F(minimal,measured),conservative_count_lower_bound=conservative,
                           maximal_upper_at_k_minus_1=threshold_before,upper_at_k=threshold_after,
                           witness_donors=len({metadata[key]['donor_id'] for gain,key in chosen}),
                           witness_models=len({metadata[key]['model_id'] for gain,key in chosen})))
        # C: exhaustive four-role context and exact within-role/reweighting contributions.
        p0,p1 = (totals[t]['H']/totals[t]['R'] for t in (0,1))
        center=(p0+p1)/2
        terms=[]
        for role in ROLES:
            rr=[roles[role,t]['R'] for t in (0,1)]
            hh=[roles[role,t]['H'] for t in (0,1)]
            aa=[rr[t]/totals[t]['R'] for t in (0,1)]
            ff=[hh[t]/rr[t] if rr[t] else None for t in (0,1)]
            within=reweight=entry=exit=F()
            if all(rr):
                within=(aa[0]+aa[1])*(ff[1]-ff[0])/2
                reweight=((ff[0]+ff[1])/2-center)*(aa[1]-aa[0])
            elif rr[1]: entry=aa[1]*(ff[1]-center)
            elif rr[0]: exit=-aa[0]*(ff[0]-center)
            total=within+reweight+entry+exit; terms.append(total)
            role_terms.append(dict(context=context,policy=pol,role=role,within=within,reweight=reweight,entry=entry,exit=exit,total=total))
            for t in (0,1):
                role_levels.append(dict(context=context,policy=pol,role=role,state=t,**roles[role,t],
                                        observed_fraction=ff[t],union_allocation_share=aa[t]))
        require(sum(terms,F()) == p1-p0, 'Four-role terms reconcile union fraction change')
    # Enumerated finite fixture: greedy maximum for each flip count equals all subset enumeration.
    fixture=[F(2,3),F(-1,2),F(1,5),F(1,5)]
    for k in range(5):
        exhaustive=max(sum((fixture[i] for i in combo),F()) for n in range(k+1) for combo in itertools.combinations(range(4),n))
        require(exhaustive==sum(sorted((x for x in fixture if x>0),reverse=True)[:k],F()),'Enumerated minimal-flip fixture')

    output.mkdir(parents=True)
    (output/'inputs').mkdir()
    provenance=[]
    for name,rel in paths.items():
        shutil.copy2(source/rel,output/'inputs'/name)
        provenance.append(dict(file=name,source_relative=rel,**digest(source/rel)))
    for name, rows in (
        ('continuous_identity_records.csv',identities),
        ('continuous_model_levels.csv',[dict(policy=pol,model_id=mid,donor_id=crc_models[mid],state=t,**expand(v)) for (pol,mid,t),v in mvals.items()]),
        ('continuous_donor_levels.csv',[dict(policy=pol,donor_id=d,state=t,**expand(v)) for (pol,d,t),v in dvals.items()]),
        ('continuous_population_levels.csv',population),('continuous_endpoints.csv',endpoints_rows),
        ('flip_summary.csv',bounds),('flip_all_measured_gains.csv',gain_rows),('flip_minimal_witness.csv',witnesses),
        ('flip_curves.csv',curves),('flip_witness_breakdown.csv',breakdowns),
        ('four_role_levels.csv',role_levels),('four_role_decomposition.csv',role_terms)):
        write(output/name,rows)
    np.save(output/'donor_draw_indices.npy',draw_ids,allow_pickle=False)
    np.save(output/'continuous_resampled_endpoints.npy',boot,allow_pickle=False)
    np.save(output/'continuous_standardised_maxima.npy',maxima,allow_pickle=False)
    report=dict(status='PASS',checks=CHECKS,source_inputs=provenance,
                python=sys.version,numpy=np.__version__,random_generator='PCG64',random_seed=470016,
                draws=10000,donor_ids=donors,family=[dict(policy=pol,metric=metric) for pol,metric in family],
                multiplier=multiplier,new_native_fits=0,protocol=digest(Path(__file__).with_name('PROTOCOL.md')),
                analysis_code=digest(Path(__file__)))
    (output/'ANALYSIS_MANIFEST.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    for row in population:
        print('CONTINUOUS',row['policy'],row['state'],{k:float(v) if isinstance(v,F) else v for k,v in row.items() if k not in ('policy','state')})
    for row in endpoints_rows:
        print('ENVELOPE',row)
    for row in bounds:
        print('FLIPS',row['context'],row['policy'],row['minimal_flips'],'measured',row['measured_identities'],'donors',row['witness_donors'],'upper',float(row['upper_before']*100),float(row['upper_at_k']*100))
    print('PASS',CHECKS,'checks')


if __name__ == '__main__':
    main()
