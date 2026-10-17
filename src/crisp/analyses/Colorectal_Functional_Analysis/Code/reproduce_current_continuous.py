"""Reproduce stored continuous tables from supplied source records and archived allocations.

Standard-library-only Python 3.10+. No fit, reranking, network or random draw.
Usage:
  python reproduce_current_continuous.py --archive-root ../research-data --source-data ../Supplementary_Data/Colorectal_Functional_Analysis

Default source-data directory is the Code directory's parent. It contains
current_continuous_selected_identities.csv, which retains unmodified source Table6
values, sample/library identities and source row representations, plus the expected
summary/contrast outputs. --source-data can point elsewhere. The archive root is
the published dataset root containing data/crc and results/crc_precision.

Every summary is independently recalculated using exact fractions of the supplied
decimal LFC representations. Missing values have no imputation or finite bounds.
The source assay preparation and its scaling convention are upstream inputs; this
script verifies the recorded values against their retained raw representations.
"""
import argparse
import csv
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

def rows(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))
def jread(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))
def require(ok, message):
    if not ok:
        raise AssertionError(message)
def weightdict(v):
    return {g: F(w) for g, w in (v.items() if isinstance(v, dict) else v)}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive-root', type=Path, required=True)
    p.add_argument('--source-data', type=Path, default=Path(__file__).resolve().parent.parent)
    a = p.parse_args()
    expected = rows(a.source_data / 'current_continuous_summary.csv')
    identity_rows = rows(a.source_data / 'current_continuous_selected_identities.csv')
    identities = {(r['sample_id'], r['gene']): r for r in identity_rows}
    require(len(identities) == len(identity_rows) == 2287, 'Complete unique selected identities')
    values = {}
    for key, r in identities.items():
        if not r['supplied_LFC']:
            values[key] = None
            continue
        value = F(r['supplied_LFC'])
        raw = json.loads(r['raw_representations'])
        require(bool(raw), 'Raw continuous source evidence')
        for representation in raw:
            require(representation[3] == key[0] and representation[4] == r['retained_library'], 'Source sample/library')
            require(F(representation[5]) == value, 'Supplied LFC unchanged')
        values[key] = value
    degree = jread(a.archive_root / 'results/crc_precision/degree_allocations.json')
    models = degree['models']
    counts = Counter(m['donor_id'] for m in models)
    require(len(models) == 85 and len(counts) == 83, 'Full85/83 population')
    by_mid = {m['model_id']: m for m in models}
    vectors = {}
    for method, path in (
        ('PRODIGY', 'results/crc_prodigy_allocations.json'),
        ('DawnRank', 'results/crc_precision/allocations.json')):
        doc = jread(a.archive_root / path)
        accum = defaultdict(dict)
        for r in doc['allocations']:
            if r['budget'] != 10 or r.get('method', method) != method:
                continue
            require(r.get('numerically_certified', True), 'Certified allocation')
            cell = r.get('cell', r.get('state')).replace('_', '')
            seed = str(r['seed']) if method == 'PRODIGY' else 'deterministic'
            key = (method, r['model_id'], cell)
            v = weightdict(r['weights'])
            require(seed not in accum[key] or accum[key][seed] == v, 'Seed agreement')
            accum[key][seed] = v
        seed_set = {'104729', '130363', '155921', '196613', '228017'} if method == 'PRODIGY' else {'deterministic'}
        for m in models:
            for cell in ('G0H0', 'G0H1', 'G1H0', 'G1H1'):
                key = (method, m['model_id'], cell)
                sv = accum[key]
                require(set(sv) == seed_set, 'Complete original seeds')
                genes = set().union(*(v.keys() for v in sv.values()))
                vectors[key] = {g: sum((v.get(g, F()) for v in sv.values()), F()) / len(sv) for g in genes}
    degree_rows = {(r['model_id'], r.get('state', r.get('cell'))): weightdict(r['weights'])
                   for r in degree.get('allocations', []) if r['budget'] == 10}
    for m in models:
        for cell in ('G0', 'G1'):
            vectors[('degree', m['model_id'], cell)] = weightdict(m['full_degree'][cell]) if 'full_degree' in m else degree_rows[(m['model_id'], cell)]
    calculated = {}
    for r in expected:
        R = A = S = F()
        for mid, m in by_mid.items():
            alpha = F(1, 83 * counts[m['donor_id']])
            v = vectors[(r['method'], mid, r['cell'])]
            require(sum(v.values(), F()) <= 10, 'Capacity')
            for g, w in v.items():
                if not w:
                    continue
                key = (m['sample_id'], g)
                source = identities[key]
                c = r['category']
                keep = c == 'ALL' or (c == 'BOTH_CORE' and source['core'] == c) or (c == 'ONCOGENE' and source['role'] == c)
                if not keep:
                    continue
                R += alpha*w
                if values[key] is not None:
                    A += alpha*w
                    S += alpha*w*(-values[key])
        result = {'returned': R, 'continuous_assessed': A, 'continuous_missing': R-A,
                  'observed_negative_LFC_sum': S,
                  'pooled_measured_negative_LFC_mean': S/A if A else None}
        for field, value in result.items():
            require(value is None and not r[field] or value is not None and abs(float(value)-float(F(r[field]))) < 2e-12,
                    'Summary disagreement: ' + str((r['method'], r['cell'], r['category'], field)))
        require(bool(r['complete_return_negative_LFC_sum']) == (R == A), 'Complete-return missingness')
        calculated[(r['method'], r['cell'], r['category'])] = result
    contrasts = {
        'graph_given_H0': {'G1H0': 1, 'G0H0': -1},
        'graph_given_H1': {'G1H1': 1, 'G0H1': -1},
        'input_given_G0': {'G0H1': 1, 'G0H0': -1},
        'input_given_G1': {'G1H1': 1, 'G1H0': -1},
        'joint': {'G1H1': 1, 'G0H0': -1},
        'interaction': {'G1H1': 1, 'G1H0': -1, 'G0H1': -1, 'G0H0': 1},
        'graph_update': {'G1': 1, 'G0': -1}}
    expected_contrasts = rows(a.source_data / 'current_continuous_contrasts.csv')
    for r in expected_contrasts:
        for field in ('returned', 'continuous_assessed', 'continuous_missing', 'observed_negative_LFC_sum'):
            value = sum((sign*calculated[(r['method'], cell, r['category'])][field]
                         for cell, sign in contrasts[r['contrast']].items()), F())
            require(abs(float(value)-float(F(r[field]))) < 2e-12, 'Contrast disagreement')
    print(json.dumps({'status': 'PASS', 'models': 85, 'donors': 83,
                      'independently_reconstructed_vectors': len(vectors),
                      'verified_summary_rows': len(expected),
                      'verified_contrast_rows': len(expected_contrasts),
                      'maximum_absolute_tolerance': 2e-12, 'new_fits': 0}))

if __name__ == '__main__':
    main()
