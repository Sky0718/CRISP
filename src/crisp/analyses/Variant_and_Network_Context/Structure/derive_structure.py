import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
from statistics import median

PARSER = argparse.ArgumentParser(description='Fixed-protocol structural aggregation; no fitting or functional-label access.')
PARSER.add_argument('--source', type=Path, required=True, help='Directory containing the accepted Graph_Controls inputs.')
PARSER.add_argument('--output', type=Path, required=True, help='Existing output directory without accepted output files.')
PARSER.add_argument('--protocol', type=Path, default=Path(__file__).resolve().with_name('PROTOCOL.txt'))
ARGS = PARSER.parse_args()
ROOT = ARGS.output.resolve()
SOURCE = ARGS.source.resolve()
PROTOCOL = ARGS.protocol.resolve()
GRAPHS = ('G0_T700', 'G1_T700')
CORE = ('BOTH_CORE', 'ORGANOID_ONLY', 'CELL_LINE_ONLY', 'NOT_LISTED_IN_EITHER')
ROLE = ('ONCOGENE', 'TSG', 'DUAL_ONCOGENE_TSG', 'OTHER_OR_UNSPECIFIED', 'UNLISTED')
CATS = ('ALL',) + CORE + ROLE
STRATA = ('decreased', 'unchanged', 'increased')
CONTRASTS = {
    'DawnRank_joint': ('G0_T700_H0', 'G1_T700_H1'),
    'DawnRank_graph_H0': ('G0_T700_H0', 'G1_T700_H0'),
    'DawnRank_graph_H1': ('G0_T700_H1', 'G1_T700_H1'),
    'Degree': ('G0_T700_DEGREE', 'G1_T700_DEGREE'),
}


def read_csv(path, delimiter=','):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def load(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    return value


def save(name, obj):
    path = ROOT / name
    assert not path.exists(), f'Preserve existing output: {path}'
    path.write_text(json.dumps(serial(obj), indent=2) + '\n', encoding='utf-8')


def binding(path):
    return {'path': str(path), 'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def midpercentile(degrees):
    n = len(degrees)
    assert n > 1
    counts = Counter(degrees.values())
    lower = 0
    by_degree = {}
    for degree in sorted(counts):
        ties = counts[degree]
        by_degree[degree] = F(2 * lower + ties - 1, 2 * (n - 1))
        lower += ties
    return {gene: by_degree[degree] for gene, degree in degrees.items()}


def fixtures():
    assert midpercentile({'a': 0, 'b': 0, 'c': 2, 'd': 4}) == {
        'a': F(1, 6), 'b': F(1, 6), 'c': F(2, 3), 'd': F(1)}
    assert midpercentile({'a': 0, 'b': 0}) == {'a': F(1, 2), 'b': F(1, 2)}
    assert F(1, 4) + F(1, 4) + F(1, 2) == 1
    differences = [F(1), F(-1)]
    assert sum(max(v, F()) for v in differences) == 1
    assert sum(max(-v, F()) for v in differences) == 1
    assert sum(differences) == 0


def main():
    fixtures()
    assert PROTOCOL.exists() and ROOT.is_dir() and SOURCE.is_dir()
    inputs = [SOURCE / 'ALLOCATION_FREEZE.json', SOURCE / 'gene_degrees.csv',
              SOURCE / 'population_state_levels.csv', SOURCE / 'graph_metadata.csv',
              *(SOURCE / 'graphs' / f'{g}.tsv' for g in GRAPHS)]
    input_bindings = [dict(binding(p), relative_path=p.relative_to(SOURCE).as_posix()) for p in inputs]
    freeze = load(inputs[0])
    models = freeze['models']
    queries = freeze['queries']
    annotations = {g: tuple(c) for g, c in freeze['annotations'].items()}
    genes = sorted(annotations)
    donors = Counter(m['donor_id'] for m in models)
    by_sample = {m['sample_id']: m for m in models}
    assert len(genes) == 7399
    assert len(models) == len(by_sample) == 85 and len(donors) == 83
    assert set(queries) == set(by_sample)
    assert tuple(freeze['categories']) == CATS
    alpha = {m['sample_id']: F(m['alpha']) for m in models}
    assert all(alpha[m['sample_id']] == F(1, 83 * donors[m['donor_id']]) for m in models)
    assert sum(alpha.values(), F()) == 1
    assert all(len(q) > 1 and len(q) == len(set(q)) and set(q) <= set(genes) for q in queries.values())
    assert all(len(c) == 3 and c[0] == 'ALL' and c[1] in CORE and c[2] in ROLE for c in annotations.values())

    degree = {g: {} for g in GRAPHS}
    for row in read_csv(inputs[1]):
        if row['graph'] in GRAPHS:
            assert row['gene'] not in degree[row['graph']]
            assert annotations[row['gene']] == ('ALL', row['core'], row['role'])
            degree[row['graph']][row['gene']] = int(row['degree'])
    metadata = {r['graph']: r for r in read_csv(SOURCE / 'graph_metadata.csv')}
    edge_checks = {}
    for graph, expected_edges in zip(GRAPHS, (119423, 65783)):
        assert set(degree[graph]) == set(genes)
        edges = read_csv(SOURCE / 'graphs' / f'{graph}.tsv', '\t')
        rebuilt = dict.fromkeys(genes, 0)
        seen = set()
        for row in edges:
            a, b = row['gene1'], row['gene2']
            assert a != b and a in rebuilt and b in rebuilt
            pair = tuple(sorted((a, b)))
            assert pair not in seen
            seen.add(pair)
            rebuilt[a] += 1
            rebuilt[b] += 1
        assert rebuilt == degree[graph]
        assert len(edges) == expected_edges and sum(rebuilt.values()) == 2 * expected_edges
        assert int(metadata[graph]['edges']) == expected_edges
        assert int(metadata[graph]['isolates']) == sum(v == 0 for v in rebuilt.values())
        edge_checks[graph] = dict(edges=len(edges), isolates=sum(v == 0 for v in rebuilt.values()), degree_sum=sum(rebuilt.values()))

    domain_pct = {g: midpercentile(degree[g]) for g in GRAPHS}
    query_pct = {s: {g: midpercentile({gene: degree[g][gene] for gene in q}) for g in GRAPHS} for s, q in queries.items()}
    domain_detail = [dict(gene=gene, core=annotations[gene][1], role=annotations[gene][2],
                          old_degree=degree[GRAPHS[0]][gene], new_degree=degree[GRAPHS[1]][gene],
                          old_domain_percentile=domain_pct[GRAPHS[0]][gene],
                          new_domain_percentile=domain_pct[GRAPHS[1]][gene]) for gene in genes]
    domain_summary, opportunity_summary = [], []
    for cat in CATS:
        members = [gene for gene in genes if cat in annotations[gene]]
        row = dict(category=cat, genes=len(members))
        for prefix, graph in zip(('old', 'new'), GRAPHS):
            values = [degree[graph][gene] for gene in members]
            row[prefix + '_mean_degree'] = F(sum(values), len(values))
            row[prefix + '_median_degree'] = F(median(values))
            row[prefix + '_isolates'] = sum(v == 0 for v in values)
            row[prefix + '_mean_domain_percentile'] = sum((domain_pct[graph][gene] for gene in members), F()) / len(members)
        domain_summary.append(row)
        opportunities = {s: [gene for gene in q if cat in annotations[gene]] for s, q in queries.items()}
        mass = sum((alpha[s] * len(q) for s, q in opportunities.items()), F())
        row = dict(category=cat, donor_weighted_eligible_genes=mass,
                   uniform_query_ten_slot_reference=sum((alpha[s] * 10 * F(len(q), len(queries[s])) for s, q in opportunities.items()), F()))
        for prefix, graph in zip(('old', 'new'), GRAPHS):
            row[prefix + '_connected_eligible_genes'] = sum((alpha[s] * sum(degree[graph][gene] > 0 for gene in q) for s, q in opportunities.items()), F())
            total_pct = sum((alpha[s] * sum((query_pct[s][graph][gene] for gene in q), F()) for s, q in opportunities.items()), F())
            row[prefix + '_mean_query_percentile'] = total_pct / mass if mass else None
        opportunity_summary.append(row)

    states_needed = tuple(dict.fromkeys(s for pair in CONTRASTS.values() for s in pair))
    allocations = {st: {s: {g: F(a) for g, a in weights.items()} for s, weights in freeze['states'][st].items()} for st in states_needed}
    for state, ss in allocations.items():
        assert set(ss) == set(queries)
        for s, aa in ss.items():
            assert set(aa) <= set(queries[s]) and all(0 < a <= 1 for a in aa.values())
            assert sum(aa.values(), F()) == min(10, len(queries[s]))
    levels = {(r['state'], r['category']): F(r['R']) for r in read_csv(SOURCE / 'population_state_levels.csv') if r['state'] in states_needed}
    for state in states_needed:
        for cat in CATS:
            total = sum((alpha[s] * a for s, aa in allocations[state].items() for gene, a in aa.items() if cat in annotations[gene]), F())
            assert total == levels[(state, cat)], (state, cat, total)

    flows = defaultdict(lambda: [F(), F(), F()])
    donor_flows = defaultdict(lambda: [F(), F(), F()])
    query_detail = []
    for s, q in queries.items():
        model = by_sample[s]
        donor = model['donor_id']
        for gene in q:
            p0, p1 = (query_pct[s][graph][gene] for graph in GRAPHS)
            change = p1 - p0
            stratum = 'increased' if change > 0 else 'decreased' if change < 0 else 'unchanged'
            query_detail.append(dict(sample_id=s, model_id=model['model_id'], donor_id=donor, alpha=alpha[s],
                                     gene=gene, core=annotations[gene][1], role=annotations[gene][2],
                                     old_degree=degree[GRAPHS[0]][gene], new_degree=degree[GRAPHS[1]][gene],
                                     old_query_percentile=p0, new_query_percentile=p1, percentile_change=change,
                                     stratum=stratum, allocations={st: allocations[st][s].get(gene, F()) for st in states_needed}))
            for contrast, (old_state, new_state) in CONTRASTS.items():
                delta = allocations[new_state][s].get(gene, F()) - allocations[old_state][s].get(gene, F())
                amounts = (max(delta, F()), max(-delta, F()), delta)
                for cat in annotations[gene]:
                    key = contrast, cat, stratum
                    dkey = contrast, donor, cat, stratum
                    for i, value in enumerate(amounts):
                        flows[key][i] += alpha[s] * value
                        donor_flows[dkey][i] += F(1, donors[donor]) * value

    for contrast, pair in CONTRASTS.items():
        for cat in CATS:
            sums = [sum((flows[(contrast, cat, st)][i] for st in STRATA), F()) for i in range(3)]
            assert sums[0] - sums[1] == sums[2] == levels[(pair[1], cat)] - levels[(pair[0], cat)]
            for st in STRATA:
                vals = flows[(contrast, cat, st)]
                assert vals[0] - vals[1] == vals[2]
                for i in range(3):
                    assert sum((donor_flows[(contrast, d, cat, st)][i] for d in donors), F()) / 83 == vals[i]
        for partition in (CORE, ROLE):
            for st in STRATA:
                for i in range(3):
                    assert sum((flows[(contrast, cat, st)][i] for cat in partition), F()) == flows[(contrast, 'ALL', st)][i]

    flow_rows = [dict(contrast=contrast, category=cat, stratum=st, gain=flows[(contrast, cat, st)][0],
                      loss=flows[(contrast, cat, st)][1], net=flows[(contrast, cat, st)][2])
                 for contrast in CONTRASTS for cat in CATS for st in STRATA]
    donor_rows = [dict(contrast=contrast, donor_id=d, category=cat, stratum=st,
                      gain=donor_flows[(contrast, d, cat, st)][0], loss=donor_flows[(contrast, d, cat, st)][1],
                      net=donor_flows[(contrast, d, cat, st)][2])
                  for contrast in CONTRASTS for d in sorted(donors) for cat in CATS for st in STRATA]
    save('domain_gene_records.json', domain_detail)
    save('domain_summary.json', domain_summary)
    save('query_summary.json', opportunity_summary)
    save('query_gene_records.json', query_detail)
    save('allocation_flow_summary.json', flow_rows)
    save('donor_flow_records.json', donor_rows)
    save('CHECKS.json', dict(status='PASS', utc=datetime.now(timezone.utc).isoformat(),
                            protocol=binding(PROTOCOL), code=binding(Path(__file__)), inputs=input_bindings,
                            models=85, donors=83, genes=7399, eligible_model_gene_records=len(query_detail),
                            categories=list(CATS), contrasts=CONTRASTS, graph_checks=edge_checks,
                            checks=['synthetic_ties_isolates_gross_flows', 'exact_source_degree_reconstruction',
                                    'exact_donor_weights', 'complete_queries_categories', 'accepted_state_R_reconciliation',
                                    'gross_gain_loss_identity', 'donor_population_reconciliation', 'both_partition_and_stratum_closure'],
                            model_fits=0, functional_labels_read=0, p_values=0))
    print(json.dumps(dict(status='PASS', models=85, donors=83, genes=7399, eligible_records=len(query_detail),
                          flow_rows=len(flow_rows), donor_rows=len(donor_rows))))
    for cat in CATS:
        d = next(r for r in domain_summary if r['category'] == cat)
        q = next(r for r in opportunity_summary if r['category'] == cat)
        print(cat, 'degree', float(d['old_mean_degree']), float(d['new_mean_degree']),
              'query_pct', float(q['old_mean_query_percentile']), float(q['new_mean_query_percentile']),
              'eligible', float(q['donor_weighted_eligible_genes']))
    for cat in ('BOTH_CORE', 'ONCOGENE'):
        for contrast in CONTRASTS:
            print(contrast, cat, [(st, [float(v) for v in flows[(contrast, cat, st)]]) for st in STRATA])


if __name__ == '__main__':
    main()
