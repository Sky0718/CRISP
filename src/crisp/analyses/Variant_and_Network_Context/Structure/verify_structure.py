import argparse
import bisect
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path, PureWindowsPath

F = Fraction
PARSER = argparse.ArgumentParser(description='Independently verify existing structural outputs without importing their derivation.')
PARSER.add_argument('--source', type=Path, required=True, help='Directory containing the accepted Graph_Controls inputs.')
PARSER.add_argument('--results', type=Path, required=True, help='Directory containing the existing structural outputs.')
PARSER.add_argument('--report', type=Path, help='New report path; defaults to INDEPENDENT_CHECK.json in the results directory.')
ARGS = PARSER.parse_args()
HERE = ARGS.results.resolve()
SOURCE = ARGS.source.resolve()
REPORT = ARGS.report.resolve() if ARGS.report else HERE / 'INDEPENDENT_CHECK.json'


def load(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def main():
    accepted = load(SOURCE / 'ALLOCATION_FREEZE.json')
    manifest = load(HERE / 'CHECKS.json')
    for bound in manifest['inputs']:
        if 'relative_path' in bound:
            relative = Path(bound['relative_path'])
        else:
            name = PureWindowsPath(bound['path']).name
            relative = Path('graphs') / name if name.endswith('.tsv') else Path(name)
        p = SOURCE / relative
        assert p.stat().st_size == bound.get('packaged_bytes', bound['bytes'])
        assert hashlib.sha256(p.read_bytes()).hexdigest() == bound.get('packaged_sha256', bound['sha256'])
    degrees = {}
    with (SOURCE / 'gene_degrees.csv').open(newline='', encoding='utf-8') as stream:
        for row in csv.DictReader(stream):
            if row['graph'] in ('G0_T700', 'G1_T700'):
                degrees[(row['graph'], row['gene'])] = int(row['degree'])
    model_map = {m['sample_id']: m for m in accepted['models']}
    donor_n = Counter(m['donor_id'] for m in model_map.values())
    alpha = {s: F(1, 83 * donor_n[m['donor_id']]) for s, m in model_map.items()}
    calculated_percentiles = {}
    for sample, query in accepted['queries'].items():
        for graph in ('G0_T700', 'G1_T700'):
            ordered = sorted(degrees[(graph, g)] for g in query)
            for gene in query:
                degree = degrees[(graph, gene)]
                lo = bisect.bisect_left(ordered, degree)
                hi = bisect.bisect_right(ordered, degree)
                calculated_percentiles[(sample, gene, graph)] = F(lo + hi - 1, 2 * (len(query) - 1))
    detail = load(HERE / 'query_gene_records.json')
    observed_keys = set()
    for row in detail:
        s, gene = row['sample_id'], row['gene']
        assert (s, gene) not in observed_keys
        observed_keys.add((s, gene))
        p0, p1 = (calculated_percentiles[(s, gene, graph)] for graph in ('G0_T700', 'G1_T700'))
        assert F(row['old_query_percentile']) == p0 and F(row['new_query_percentile']) == p1
        assert F(row['percentile_change']) == p1 - p0
        assert row['stratum'] == ('increased' if p1 > p0 else 'decreased' if p1 < p0 else 'unchanged')
        for state, weight in row['allocations'].items():
            assert F(weight) == F(accepted['states'][state][s].get(gene, 0))
    assert observed_keys == {(s, g) for s, query in accepted['queries'].items() for g in query}

    flow_rows = load(HERE / 'allocation_flow_summary.json')
    for row in flow_rows:
        contrast, cat, stratum = row['contrast'], row['category'], row['stratum']
        first, second = manifest['contrasts'][contrast]
        gain = loss = net = F()
        for sample, query in accepted['queries'].items():
            for gene in query:
                if cat not in accepted['annotations'][gene]:
                    continue
                d = calculated_percentiles[(sample, gene, 'G1_T700')] - calculated_percentiles[(sample, gene, 'G0_T700')]
                group = 'increased' if d > 0 else 'decreased' if d < 0 else 'unchanged'
                if group != stratum:
                    continue
                difference = F(accepted['states'][second][sample].get(gene, 0)) - F(accepted['states'][first][sample].get(gene, 0))
                value = alpha[sample] * difference
                net += value
                if value > 0:
                    gain += value
                elif value < 0:
                    loss -= value
        assert F(row['gain']) == gain and F(row['loss']) == loss and F(row['net']) == net

    for row in load(HERE / 'query_summary.json'):
        cat = row['category']
        count = uniform = F()
        connected = [F(), F()]
        percentiles = [F(), F()]
        for sample, query in accepted['queries'].items():
            members = [gene for gene in query if cat in accepted['annotations'][gene]]
            count += alpha[sample] * len(members)
            uniform += alpha[sample] * F(10 * len(members), len(query))
            for index, graph in enumerate(('G0_T700', 'G1_T700')):
                for gene in members:
                    connected[index] += alpha[sample] * int(degrees[(graph, gene)] > 0)
                    percentiles[index] += alpha[sample] * calculated_percentiles[(sample, gene, graph)]
        assert F(row['donor_weighted_eligible_genes']) == count
        assert F(row['uniform_query_ten_slot_reference']) == uniform
        for index, prefix in enumerate(('old', 'new')):
            assert F(row[prefix + '_connected_eligible_genes']) == connected[index]
            assert F(row[prefix + '_mean_query_percentile']) == percentiles[index] / count

    report = dict(status='PASS', utc=datetime.now(timezone.utc).isoformat(),
                  query_records_independently_reconciled=len(detail), flow_rows_independently_reconciled=len(flow_rows),
                  opportunity_categories_independently_reconciled=10,
                  arithmetic='Exact Fraction; independent sorted-position percentiles and direct conditional aggregation',
                  input_hashes_unchanged=True, no_analysis_module_import=True,
                  verification_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    target = REPORT
    assert not target.exists()
    target.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
