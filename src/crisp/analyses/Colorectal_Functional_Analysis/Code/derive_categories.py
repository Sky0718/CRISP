"""Finite descriptive derivations from accepted allocations; no native fits or reranking.

Preserves four cells, 85 models / 83 donors, five PRODIGY seeds and strict
DawnRank allocations. Exact rational arithmetic is primary; decimal columns
are provided for plotting only. Existing sources are never modified.
"""
import csv
import argparse
import hashlib
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--archive', type=Path, required=True,
                    help='Extracted base research-data archive containing data/ and results/.')
parser.add_argument('--crisp-source', type=Path,
                    help='Source directory containing the crisp package.')
parser.add_argument('--reference-tables', type=Path,
                    default=Path(__file__).resolve().parents[2] / 'functional_content',
                    help='Accepted functional_content tables supplied beside Colorectal_Functional_Analysis.')
parser.add_argument('--output', type=Path, required=True,
                    help='New output directory; do not point to original archived sources.')
parser.add_argument('--input-references', type=Path, required=True,
                    help='Input-reference JSON supplied with the supplementary data.')
args = parser.parse_args()
Z = args.archive.resolve()
OUT = args.output.resolve()
if args.crisp_source:
    sys.path.insert(0, str(args.crisp_source.resolve()))
from crisp.statistics.inputs import load_crc, load_categories

CONFIG = {
    'degree_allocations': 'results/crc_precision/degree_allocations.json',
    'allocations': [
        {'method': 'PRODIGY', 'path': 'results/crc_prodigy_allocations.json'},
        {'method': 'DawnRank', 'path': 'results/crc_precision/allocations.json'}],
    'labels': ['data/crc/selected_labels.json']}
FILES = [CONFIG['degree_allocations'], *(x['path'] for x in CONFIG['allocations']),
         *CONFIG['labels'], 'data/crc/categories.json']
METHODS = ('PRODIGY', 'DawnRank')
CELLS = ('G0H0', 'G0H1', 'G1H0', 'G1H1')
CATS = ('ALL', 'BOTH_CORE', 'ONCOGENE')
FIELDS = ('returned', 'positive', 'zero', 'unknown', 'absent', 'explicit_null',
          'identity_conflict', 'call_conflict')
FLOW = ('old', 'retained', 'gain', 'loss', 'new', 'net')
GROUPS = ('positive', 'zero', 'unknown')
CONTRASTS = {
    'graph_given_H0': {'G1H0': 1, 'G0H0': -1},
    'graph_given_H1': {'G1H1': 1, 'G0H1': -1},
    'input_given_G0': {'G0H1': 1, 'G0H0': -1},
    'input_given_G1': {'G1H1': 1, 'G1H0': -1},
    'joint': {'G1H1': 1, 'G0H0': -1},
    'interaction': {'G1H1': 1, 'G1H0': -1, 'G0H1': -1, 'G0H0': 1}}
checks = 0
def require(ok, text):
    global checks
    checks += 1
    if not ok:
        raise AssertionError(text)

def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def encode(row):
    result = {}
    for k, v in row.items():
        if isinstance(v, F):
            result[k] = str(v)
            result[k + '_decimal'] = float(v)
        else:
            result[k] = v
    return result

def save_csv(name, rows):
    rows = [encode(r) for r in rows]
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with (OUT / name).open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)

require(not OUT.is_relative_to(Z), 'output must be outside the archived input tree')
require(not OUT.is_relative_to(args.reference_tables.resolve()),
        'output must be outside the accepted reference-table tree')
require(not OUT.exists() or not any(OUT.iterdir()), 'use a new or empty output directory')
references = json.loads(args.input_references.read_text(encoding='utf-8'))
for binding in references['inputs']:
    base = Z if binding['root'] == 'archive' else args.reference_tables.resolve()
    path = base / binding['path']
    require(hashlib.sha256(path.read_bytes()).hexdigest() == binding['sha256'],
            'input hash differs from the source record: ' + binding['path'])
models, seed_allocations, means, degrees, calls, policies = load_crc(Z, CONFIG)
mapping = load_categories(Z / 'data/crc/categories.json', models, 'CRC')
require(len(models) == 85, '85 models')
counts = Counter(m['donor_id'] for m in models)
require(len(counts) == 83, '83 donors')
by_mid = {m['model_id']: m for m in models}
alpha = {m['model_id']: F(1, 83 * counts[m['donor_id']]) for m in models}
require(sum(alpha.values()) == 1, 'equal donor weights')
def cat_member(g, c):
    return c == 'ALL' or (c == 'BOTH_CORE' and mapping['core'][g] == c) or (
        c == 'ONCOGENE' and mapping['role'][g] == c)
def label_group(sample, gene):
    s, y = calls[(sample, gene)]
    return ('positive' if y == 1 else 'zero') if s == 'MEASURED' else 'unknown'

levels, contrast_rows = [], []
for method in METHODS:
    for model in models:
        mid, sample = model['model_id'], model['sample_id']
        identity = {k: model[k] for k in ('model_id', 'sample_id', 'donor_id')}
        for cell in CELLS:
            vector = means[(method, mid, cell)]
            for category in CATS:
                mass = dict.fromkeys(FIELDS, F())
                for gene, weight in vector.items():
                    if not cat_member(gene, category):
                        continue
                    state, value = calls[(sample, gene)]
                    mass['returned'] += weight
                    group = label_group(sample, gene)
                    mass[group] += weight
                    if group == 'unknown':
                        field = {'ABSENT': 'absent', 'EXPLICIT_NULL': 'explicit_null',
                                 'IDENTITY_CONFLICT': 'identity_conflict',
                                 'CONFLICTING_CALLS': 'call_conflict'}[state]
                        mass[field] += weight
                require(mass['returned'] == mass['positive'] + mass['zero'] + mass['unknown'],
                        'model category mass conservation')
                levels.append({'method': method, **identity, 'cell': cell, 'category': category,
                               **mass})
        for name, multipliers in CONTRASTS.items():
            vector = defaultdict(F)
            for cell, sign in multipliers.items():
                for gene, w in means[(method, mid, cell)].items():
                    vector[gene] += sign * w
            for category in CATS:
                vals = {f: F() for f in ('returned', 'positive', 'zero', 'unknown',
                                         'positive_lower', 'positive_upper')}
                for gene, coefficient in vector.items():
                    if not cat_member(gene, category):
                        continue
                    group = label_group(sample, gene)
                    vals['returned'] += coefficient
                    vals[group] += coefficient
                    if group == 'unknown':
                        vals['positive_lower'] += min(coefficient, 0)
                        vals['positive_upper'] += max(coefficient, 0)
                vals['positive_lower'] += vals['positive']
                vals['positive_upper'] += vals['positive']
                require(vals['returned'] == vals['positive'] + vals['zero'] + vals['unknown'],
                        'model contrast mass conservation')
                contrast_rows.append({'method': method, **identity, 'contrast': name,
                                      'category': category, **vals})

def aggregate(rows, fields, state_key, unit):
    grouped = {}
    for r in rows:
        key = (r['method'], r[state_key], r['category'])
        if unit == 'donor':
            key += (r['donor_id'],)
        if key not in grouped:
            grouped[key] = {'method': key[0], state_key: key[1], 'category': key[2],
                            **({'donor_id': key[3]} if unit == 'donor' else {}),
                            **dict.fromkeys(fields, F())}
        weight = F(1, counts[r['donor_id']]) / (83 if unit == 'population' else 1)
        for field in fields:
            grouped[key][field] += weight * r[field]
    return list(grouped.values())

donor_levels = aggregate(levels, FIELDS, 'cell', 'donor')
population_levels = aggregate(levels, FIELDS, 'cell', 'population')
contrast_fields = ('returned', 'positive', 'zero', 'unknown', 'positive_lower', 'positive_upper')
donor_contrasts = aggregate(contrast_rows, contrast_fields, 'contrast', 'donor')
population_contrasts = aggregate(contrast_rows, contrast_fields, 'contrast', 'population')
for rows in (contrast_rows, donor_contrasts, population_contrasts):
    for r in rows:
        lo, hi = r['positive_lower'], r['positive_upper']
        r['positive_direction'] = 'negative' if hi < 0 else 'positive' if lo > 0 else (
            'zero' if lo == hi == 0 else 'unresolved')

# Full old/new oncogene identity accounting uses seed-mean weight within model.
# A retained weight is min(old,new) at literal sample--gene identity, before
# averaging models/donors; gained and lost weights are separately nonnegative.
identities = []
oncogenes = sorted({g for m in models for g in m['query'] if mapping['role'][g] == 'ONCOGENE'})
for method in METHODS:
    for model in models:
        mid, sample = model['model_id'], model['sample_id']
        old, new = (means[(method, mid, c)] for c in ('G0H0', 'G1H1'))
        for gene in sorted(set(old) | set(new)):
            if mapping['role'][gene] != 'ONCOGENE':
                continue
            a, b = old.get(gene, F()), new.get(gene, F())
            flow = {'old': a, 'retained': min(a, b), 'gain': max(b-a, F()),
                    'loss': max(a-b, F()), 'new': b, 'net': b-a}
            require(flow['retained'] + flow['loss'] == a and flow['retained'] + flow['gain'] == b,
                    'identity flow conservation')
            state, value = calls[(sample, gene)]
            identities.append({'method': method, **{k: model[k] for k in
                               ('model_id', 'sample_id', 'donor_id')}, 'gene': gene,
                               'core': mapping['core'][gene], 'label_group': label_group(sample, gene),
                               'label_state': state, 'label_value': value,
                               'donor_weight': alpha[mid], **flow})

def flows(unit):
    # Explicit complete grids retain zero-flow genes/donors and all three statuses.
    units = oncogenes if unit == 'gene' else sorted(counts) if unit == 'donor' else ['ALL_83_DONORS']
    buckets = {}
    for method in METHODS:
        for name in units:
            for group in GROUPS:
                r = {'method': method, 'label_group': group, **dict.fromkeys(FLOW, F()),
                     'selected_identities': 0, 'gained_identities': 0, 'lost_identities': 0}
                if unit == 'gene':
                    r.update(gene=name, core=mapping['core'][name])
                elif unit == 'donor':
                    r.update(donor_id=name)
                buckets[(method, name, group)] = r
    for r in identities:
        name = r['gene'] if unit == 'gene' else r['donor_id'] if unit == 'donor' else 'ALL_83_DONORS'
        out = buckets[(r['method'], name, r['label_group'])]
        weight = r['donor_weight'] * (83 if unit == 'donor' else 1)
        for field in FLOW:
            out[field] += weight * r[field]
        out['selected_identities'] += 1
        out['gained_identities'] += r['gain'] > 0
        out['lost_identities'] += r['loss'] > 0
    for r in buckets.values():
        require(r['old'] == r['retained'] + r['loss'] and r['new'] == r['retained'] + r['gain']
                and r['net'] == r['gain'] - r['loss'], 'aggregated flow conservation')
    return list(buckets.values())
gene_flows, donor_flows, population_flows = flows('gene'), flows('donor'), flows('population')

lookup = {(r['method'], r['cell'], r['category']): r for r in population_levels}
shares = []
for r in population_levels:
    if r['category'] == 'ALL':
        continue
    total = lookup[(r['method'], r['cell'], 'ALL')]['returned']
    require(total > 0, 'pooled returned denominator is positive')
    shares.append({**r, 'all_returned': total,
                   'share_of_all_returned': r['returned'] / total,
                   'pooled_positive_fraction': r['positive'] / r['returned'] if r['returned'] else None})

# Exact reconciliation against accepted old/joint-new source tables and identities.
old_tables = args.reference_tables.resolve()
for r in read_csv(old_tables / 'population_levels.csv'):
    if (r['partition'], r['category']) not in [('ALL','ALL'), ('core','BOTH_CORE'), ('role','ONCOGENE')]:
        continue
    method = 'DawnRank' if r['method'] == 'DawnRank_strict' else r['method']
    got = lookup[(method, r['cell'], r['category'])]
    for old_field, field in [('returned','returned'), ('known_positive','positive'),
                             ('known_zero','zero'), ('unknown','unknown')]:
        require(F(r[old_field]) == got[field], 'accepted absolute level reconciliation')
existing_ledger = read_csv(old_tables / 'ledger.csv')
expected = {(r['method'].replace('DawnRank_strict','DawnRank'),r['model_id'],r['gene']): r
            for r in existing_ledger if r['role'] == 'ONCOGENE'}
require(set(expected) == {(r['method'],r['model_id'],r['gene']) for r in identities},
        'oncogene old/new support equals accepted ledger')
for r in identities:
    e = expected[(r['method'], r['model_id'], r['gene'])]
    require(F(e['old_weight']) == r['old'] and F(e['new_weight']) == r['new']
            and e['core'] == r['core'] and e['state'] == r['label_state']
            and e['sample_id'] == r['sample_id'] and e['donor_id'] == r['donor_id'],
            'oncogene identity reconciliation')
for r in population_flows:
    for col, cell in [('old','G0H0'), ('new','G1H1')]:
        require(r[col] == lookup[(r['method'],cell,'ONCOGENE')][r['label_group']],
                'flow equals category absolute level')

OUT.mkdir(parents=True, exist_ok=True)
outputs = {'four_arm_model_levels.csv': levels, 'four_arm_donor_levels.csv': donor_levels,
           'four_arm_population_levels.csv': population_levels,
           'four_arm_model_contrasts.csv': contrast_rows,
           'four_arm_donor_contrasts.csv': donor_contrasts,
           'four_arm_population_contrasts.csv': population_contrasts,
           'oncogene_identity_flows.csv': identities, 'oncogene_gene_flows.csv': gene_flows,
           'oncogene_donor_flows.csv': donor_flows, 'oncogene_population_flows.csv': population_flows,
           'category_returned_shares.csv': shares}
for name, rows in outputs.items():
    save_csv(name, rows)
source_files = [Z / name for name in FILES] + [old_tables/'population_levels.csv', old_tables/'ledger.csv']
report = {'status': 'PASS', 'models': 85, 'donors': 83, 'budget': 10,
          'prodigy_seed_count': 5, 'dawnrank_seed_semantics': 'deterministic aliases; strict certified allocations',
          'checks': checks, 'native_fits': 0, 'rerankings': 0,
          'contrasts': CONTRASTS, 'oncogene_query_gene_universe': len(oncogenes),
          'weighting': 'mean five PRODIGY seed allocations within model, equal model weights within donor, equal 83 donors',
          'flow_semantics': 'min and positive weight differences computed at literal sample-gene identity AFTER seed averaging',
          'share_semantics': 'ratio of equal-donor mean category returned mass to equal-donor mean all returned mass; not mean donor-specific ratios',
          'uncertainty': 'finite descriptive allocation and exact shared-label bounds only; no new sampling intervals or p-values',
          'sources': [{'path': str(p), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in source_files],
          'outputs': {name: {'rows': len(rows), 'sha256': hashlib.sha256((OUT/name).read_bytes()).hexdigest()}
                      for name, rows in outputs.items()},
          'focal_unknown_all_four_cells': [encode(r) for r in population_levels
                                          if r['category'] != 'ALL' and r['unknown'] != 0]}
(OUT/'biology_derivation_manifest.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps({'checks': checks, 'outputs': report['outputs'],
                  'focal_unknown': report['focal_unknown_all_four_cells'],
                  'contrasts': [encode(r) for r in population_contrasts if r['category'] != 'ALL'],
                  'flows': [encode(r) for r in population_flows]}, indent=2))
