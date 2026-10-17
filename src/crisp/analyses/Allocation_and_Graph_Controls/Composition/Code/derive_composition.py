"""Donor-composition and continuous-score summaries from allocation ledgers.

Requires NumPy. Definitions and resampling families are in ANALYSIS_DESIGN.txt.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

import numpy as np

POLICIES = ('PRODIGY', 'DawnRank', 'degree')
CATS = ('BOTH_CORE', 'ONCOGENE')
ROLES = ('ONCOGENE', 'TSG', 'DUAL_ONCOGENE_TSG', 'OTHER_OR_UNSPECIFIED', 'UNLISTED')
MAP = {'DawnRank_strict': 'DawnRank'}
checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(message)


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, records):
    encoded = []
    for row in records:
        e = {}
        for k, v in row.items():
            if isinstance(v, F):
                e[k] = str(v)
                e[k + '_decimal'] = float(v)
            else:
                e[k] = v
        encoded.append(e)
    fields = list(dict.fromkeys(k for row in encoded for k in row))
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(encoded)


def decompose(policy, category, values, numerator, measure):
    """Exact symmetric centered decomposition. Undefined donor ratios stay None."""
    total = {t: sum((values[d, t]['R'] for d in DONORS), F()) for t in (0, 1)}
    require(all(v > 0 for v in total.values()), 'Positive pooled denominators')
    p = {t: sum((values[d, t][numerator] for d in DONORS), F()) / total[t] for t in (0, 1)}
    center = (p[0] + p[1]) / 2
    details = []
    counts = Counter()
    terms = defaultdict(F)
    for donor in DONORS:
        r0, r1 = (values[donor, t]['R'] for t in (0, 1))
        a0, a1 = r0 / total[0], r1 / total[1]
        f0 = values[donor, 0][numerator] / r0 if r0 else None
        f1 = values[donor, 1][numerator] / r1 if r1 else None
        support = 'both' if r0 and r1 else 'exit' if r0 else 'entrant' if r1 else 'neither'
        counts[support] += 1
        x = dict.fromkeys(('within_common', 'reweight_common', 'entrant', 'exit'), F())
        if support == 'both':
            x['within_common'] = (a0 + a1) / 2 * (f1 - f0)
            x['reweight_common'] = ((f0 + f1) / 2 - center) * (a1 - a0)
        elif support == 'entrant':
            x['entrant'] = a1 * (f1 - center)
        elif support == 'exit':
            x['exit'] = -a0 * (f0 - center)
        for k, value in x.items():
            terms[k] += value
        # Raw uncentered contribution supplies an independent per-donor identity.
        direct = values[donor, 1][numerator] / total[1] - values[donor, 0][numerator] / total[0]
        require(sum(x.values(), F()) == direct - center * (a1 - a0), 'Donor decomposition identity')
        details.append({'policy': policy, 'category': category, 'measure': measure,
                        'donor_id': donor, 'support': support, 'R0': r0, 'R1': r1,
                        'numerator0': values[donor, 0][numerator],
                        'numerator1': values[donor, 1][numerator],
                        'donor_ratio0': f0, 'donor_ratio1': f1, **x})
    require(sum(terms.values(), F()) == p[1] - p[0], 'Exact full-cohort decomposition identity')
    require(sum(counts.values()) == 83, 'All 83 support-status donors retained')
    summary = {'policy': policy, 'category': category, 'measure': measure,
               'ratio0': p[0], 'ratio1': p[1], 'difference': p[1]-p[0],
               'both': counts['both'], 'entrants': counts['entrant'],
               'exits': counts['exit'], 'neither': counts['neither'], **terms}
    return summary, details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    require(not output.is_relative_to(source), 'Never write into source')
    paths = {name: source / name for name in (
        'functional_content/ledger.csv', 'functional_content/population_levels.csv',
        'Category_Sensitivity/donor_category_levels.csv', 'Category_Sensitivity/donor_index.csv',
        'Category_Sensitivity/colorectal_donor_draws.csv',
        'Colorectal_Functional_Analysis/current_continuous_model_levels.csv',
        'Colorectal_Functional_Analysis/current_continuous_selected_identities.csv',
        'Colorectal_Functional_Analysis/current_continuous_summary.csv')}
    source_hashes = {str(p.relative_to(source)): sha(p) for p in paths.values()}
    index = [r for r in rows(paths['Category_Sensitivity/donor_index.csv']) if r['context'] == 'CRC']
    require([int(r['donor_index']) for r in index] == list(range(83)), 'Stored donor index 0..82')
    global DONORS
    DONORS = tuple(r['donor_id'] for r in index)
    require(len(set(DONORS)) == 83, '83 unique donors')
    binary = {}
    for r in rows(paths['Category_Sensitivity/donor_category_levels.csv']):
        if r['context'] != 'CRC':
            continue
        key = (MAP.get(r['policy'], r['policy']), r['category'], r['donor_id'], int(r['state']))
        require(key not in binary, 'Unique donor-category key')
        binary[key] = {k: F(r[k]) for k in ('R', 'H', 'Z', 'U')}
        require(binary[key]['R'] == sum((binary[key][k] for k in ('H', 'Z', 'U')), F()), 'Binary mass conservation')
    # Complete role table includes degree from accepted donor levels; methods
    # are additionally reconstructed from literal model-gene allocation ledger.
    ledger = rows(paths['functional_content/ledger.csv'])
    model_map = {}
    literal_roles = defaultdict(lambda: dict.fromkeys(('R', 'H', 'Z', 'U'), F()))
    literal_score = defaultdict(lambda: dict.fromkeys(('R', 'S', 'A'), F()))
    source_identity = rows(paths['Colorectal_Functional_Analysis/current_continuous_selected_identities.csv'])
    lfc = {(r['sample_id'], r['gene']): F(r['supplied_LFC']) if r['supplied_LFC'] else None for r in source_identity}
    require(len(lfc) == len(source_identity) == 2287, 'Complete selected LFC identities')
    for r in ledger:
        policy = MAP.get(r['method'], r['method'])
        mid, donor = r['model_id'], r['donor_id']
        require(mid not in model_map or model_map[mid] == donor, 'Model-donor consistency')
        model_map[mid] = donor
        require(r['role'] in ROLES, 'Exhaustive original five roles')
        alpha = F(r['alpha'])
        for t, field in enumerate(('old_weight', 'new_weight')):
            weight = F(r[field])
            v = literal_roles[policy, r['role'], t]
            v['R'] += alpha*weight
            group = 'H' if r['state'] == 'MEASURED' and F(r['value']) == 1 else 'Z' if r['state'] == 'MEASURED' else 'U'
            v[group] += alpha*weight
            for cat in CATS:
                if (cat == 'BOTH_CORE' and r['core'] != cat) or (cat == 'ONCOGENE' and r['role'] != cat):
                    continue
                z = literal_score[policy, mid, cat, t]
                z['R'] += weight
                value = lfc[r['sample_id'], r['gene']]
                if value is not None:
                    z['A'] += weight
                    z['S'] += weight*(-value)
    counts = Counter(model_map.values())
    require(len(model_map) == 85 and set(counts) == set(DONORS), 'Full 85/83 population')
    for r in ledger:
        require(F(r['alpha']) == F(1, 83*counts[r['donor_id']]), 'Equal donor/model weights')
    role_rows = []
    for policy in POLICIES:
        for t in (0, 1):
            pooled = {role: {field: sum((binary[policy, role, d, t][field] for d in DONORS), F())/83
                             for field in ('R', 'H', 'Z', 'U')} for role in ROLES}
            total = sum((v['R'] for v in pooled.values()), F())
            expected_total = sum((binary[policy, 'ALL', d, t]['R'] for d in DONORS), F())/83
            require(total == expected_total, 'Five roles exhaust returned mass')
            listed = total - pooled['UNLISTED']['R']
            for role, vals in pooled.items():
                if policy != 'degree':
                    require(vals == literal_roles[policy, role, t], 'Independent literal-ledger role reconstruction')
                role_rows.append({'policy': policy, 'state': t, 'role': role, **vals,
                                  'all_returned': total, 'cancer_listed_returned': listed,
                                  'share_all_returned': vals['R']/total,
                                  'share_cancer_listed': vals['R']/listed if role != 'UNLISTED' else None,
                                  'donors_with_allocation': sum(binary[policy, role, d, t]['R'] > 0 for d in DONORS)})
    reference = rows(paths['functional_content/population_levels.csv'])
    lookup_roles = {(r['policy'], r['role'], r['state']): r for r in role_rows}
    for r in reference:
        if r['partition'] != 'role':
            continue
        got = lookup_roles[MAP.get(r['method'], r['method']), r['category'], {'G0H0': 0, 'G1H1': 1}[r['cell']]]
        for field, old in (('R','returned'), ('H','known_positive'), ('Z','known_zero'), ('U','unknown')):
            require(got[field] == F(r[old]), 'Accepted population role table exact match')
    # Aggregate complete continuous model rows to donors, no donor-wise ratio averaging.
    continuous = defaultdict(lambda: dict.fromkeys(('R', 'A', 'S'), F()))
    seen = Counter()
    model_scores = rows(paths['Colorectal_Functional_Analysis/current_continuous_model_levels.csv'])
    for r in model_scores:
        if r['category'] not in CATS or r['cell'] not in ('G0H0','G1H1','G0','G1'):
            continue
        t = int(r['cell'] in ('G1H1', 'G1'))
        policy, cat, mid, donor = r['method'], r['category'], r['model_id'], r['donor_id']
        require(model_map[mid] == donor, 'Continuous model donor matches')
        require(F(r['continuous_missing']) == 0, 'Focal category completely continuously measured')
        vals = {'R': F(r['returned']), 'A': F(r['continuous_assessed']), 'S': F(r['observed_negative_LFC_sum'])}
        if policy != 'degree':
            independent = literal_score[policy, mid, cat, t]
            require(vals['R'] == independent['R'] and vals['A'] == independent['A'], 'Continuous literal allocation match')
            require(abs(float(vals['S'] - independent['S'])) < 2e-12, 'Continuous score independently reconstructed')
            vals = independent
        require(vals['R'] == binary[policy, cat, donor, t]['R'] or counts[donor] > 1, 'Single-model category R matches')
        for f, value in vals.items():
            continuous[policy, cat, donor, t][f] += value/counts[donor]
        seen[policy, cat, t] += 1
    require(set(seen.values()) == {85} and len(seen) == 12, 'Every continuous policy/class/state has 85 models')
    cont_rows = []
    for policy in POLICIES:
        for cat in CATS:
            for d in DONORS:
                for t in (0, 1):
                    v = continuous[policy, cat, d, t]
                    require(v['R'] == v['A'] == binary[policy,cat,d,t]['R'], 'Donor focal continuous denominator matches')
                    cont_rows.append({'policy': policy,'category':cat,'donor_id':d,'state':t, **v})
    summaries, details = [], []
    for policy in POLICIES:
        for cat in CATS:
            for measure, sourcevals, numerator in (('binary_fraction', binary, 'H'), ('negative_LFC_mean', continuous, 'S')):
                vals = {(d,t): sourcevals[policy,cat,d,t] for d in DONORS for t in (0,1)}
                s, z = decompose(policy,cat,vals,numerator,measure)
                summaries.append(s)
                details.extend(z)
    drawrows = rows(paths['Category_Sensitivity/colorectal_donor_draws.csv'])
    require(len(drawrows) == 10000, 'All original 10000 draws')
    require([int(r['draw_index']) for r in drawrows] == list(range(10000)), 'Original draw order')
    draws = np.array([[int(r['position_'+str(j+1)]) for j in range(83)] for r in drawrows], dtype=np.int64)
    require(draws.min() >= 0 and draws.max() < 83, 'Valid draw indices')
    multiplicities = np.zeros((10000,83), dtype=np.int64)
    for b, d in enumerate(draws):
        multiplicities[b] = np.bincount(d, minlength=83)
    endpoint_rows, draw_endpoints = [], []
    for policy in POLICIES:
        for cat in CATS:
            vals = {(d,t): continuous[policy,cat,d,t] for d in DONORS for t in (0,1)}
            RR = [sum((vals[d,t]['R'] for d in DONORS), F())/83 for t in (0,1)]
            SS = [sum((vals[d,t]['S'] for d in DONORS), F())/83 for t in (0,1)]
            QQ = [SS[t]/RR[t] for t in (0,1)]
            array = {f: np.array([[float(vals[d,t][f]) for d in DONORS] for t in (0,1)]) for f in ('R','S')}
            sample = {f: multiplicities @ array[f].T/83 for f in ('R','S')}
            require(np.all(sample['R'] > 0), 'No undefined continuous resample denominator')
            direct = {f: array[f][:,draws].mean(axis=2).T for f in ('R','S')}
            require(all(np.max(np.abs(sample[f]-direct[f])) < 2e-12 for f in ('R','S')), 'Independent every-draw direct/multiplicity mean agreement')
            for endpoint in ('Delta_S','Delta_Q'):
                theta = SS[1]-SS[0] if endpoint == 'Delta_S' else QQ[1]-QQ[0]
                influence = [(vals[d,1]['S']-vals[d,0]['S']-theta) if endpoint == 'Delta_S' else
                             (vals[d,1]['S']-QQ[1]*vals[d,1]['R'])/RR[1] - (vals[d,0]['S']-QQ[0]*vals[d,0]['R'])/RR[0] for d in DONORS]
                require(sum(influence, F()) == 0, 'Exact centered continuous influence')
                var = sum((z*z for z in influence), F())/(83*82)
                require(var > 0, 'All planned continuous scales nondegenerate')
                scale = float(var)**0.5
                if endpoint == 'Delta_S':
                    rep = sample['S'][:,1]-sample['S'][:,0]
                    alt = direct['S'][:,1]-direct['S'][:,0]
                else:
                    rep = sample['S'][:,1]/sample['R'][:,1]-sample['S'][:,0]/sample['R'][:,0]
                    alt = direct['S'][:,1]/direct['R'][:,1]-direct['S'][:,0]/direct['R'][:,0]
                require(np.max(np.abs(rep-alt)) < 2e-12, 'Every continuous resampled endpoint independently checked')
                endpoint_rows.append({'policy':policy,'category':cat,'endpoint':endpoint,
                                      'estimate':theta,'scale':scale,'R0':RR[0],'R1':RR[1],
                                      'S0':SS[0],'S1':SS[1],'Q0':QQ[0],'Q1':QQ[1]})
                draw_endpoints.append(rep)
    reps = np.array(draw_endpoints).T
    theta = np.array([float(r['estimate']) for r in endpoint_rows])
    scales = np.array([r['scale'] for r in endpoint_rows])
    maxima = np.max(np.abs(reps-theta)/scales, axis=1)
    critical = float(np.quantile(maxima, .95, method='linear'))
    for r in endpoint_rows:
        r['critical_value'] = critical
        r['lower'] = float(r['estimate'])-critical*r['scale']
        r['upper'] = float(r['estimate'])+critical*r['scale']
        r['direction'] = 'negative' if r['upper'] < 0 else 'positive' if r['lower'] > 0 else 'crosses_zero'
    # All donor leave-one-out ratios supplement support accounting without resampling selection.
    loo = []
    for policy in POLICIES:
        for cat in CATS:
            for numerator, sourcevals, measure in (('H',binary,'binary_fraction'),('S',continuous,'negative_LFC_mean')):
                for omitted in DONORS:
                    totals = [{f:sum((sourcevals[policy,cat,d,t][f] for d in DONORS if d != omitted),F()) for f in ('R',numerator)} for t in (0,1)]
                    require(all(z['R'] > 0 for z in totals), 'Leave-one-out denominator remains positive')
                    delta = totals[1][numerator]/totals[1]['R']-totals[0][numerator]/totals[0]['R']
                    loo.append({'policy':policy,'category':cat,'measure':measure,'omitted_donor':omitted,'difference':delta})
    output.mkdir(parents=True, exist_ok=True)
    outputs = {'role_overview.csv':role_rows, 'donor_continuous_levels.csv':cont_rows,
               'ratio_decomposition_summary.csv':summaries,'ratio_decomposition_donors.csv':details,
               'continuous_composition_endpoints.csv':endpoint_rows, 'leave_one_donor_out_ratios.csv':loo}
    for name, records in outputs.items():
        write_csv(output/name, records)
    np.save(output/'continuous_resampled_endpoints.npy', reps)
    np.save(output/'continuous_standardised_maxima.npy', maxima)
    require(all(sha(source / k) == h for k,h in source_hashes.items()), 'All immutable inputs unchanged')
    report = {'status':'PASS','checks':checks,'models':85,'donors':83,'prodigy_seed_count':5,
              'seed_semantics':'Original five-seed means, strict deterministic DawnRank, full-query degree unchanged',
              'post_exposure':True,'native_fits':0,'rerankings':0,
              'family_coordinates':12,'draws':10000,'undefined_draws':0,'critical_value':critical,
              'family':'3 policies x 2 categories x Delta_S / Delta_Q; distinct from earlier binary families',
              'uncertainty':'Exploratory conditional donor-composition sensitivity; not population or assay-error intervals',
              'design_sha256':sha(Path(__file__).with_name('ANALYSIS_DESIGN.txt')),
              'sources':source_hashes,'outputs':{name:{'rows':len(records),'sha256':sha(output/name)} for name,records in outputs.items()}}
    (output/'COMPOSITION_MANIFEST.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({'checks':checks,'critical_value':critical,'endpoints':[{k:float(v) if isinstance(v,F) else v for k,v in r.items()} for r in endpoint_rows],
                      'decomposition':[{k:float(v) if isinstance(v,F) else v for k,v in r.items()} for r in summaries]},indent=2))


if __name__ == '__main__':
    main()
