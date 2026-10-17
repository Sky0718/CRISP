import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = None
ARGS = None
SEEDS = (271828, 314159, 173205)


def save(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def count_matched(edges, count, seed):
    boundary = sorted(edges.values(), reverse=True)[count - 1]
    high = {key: score for key, score in edges.items() if score > boundary}
    tied = [key for key, score in edges.items() if score == boundary]
    ordered = sorted(tied, key=lambda key: hashlib.sha256(f'{seed}\t{key[0]}\t{key[1]}'.encode()).digest())
    needed = count - len(high)
    high.update({key: boundary for key in ordered[:needed]})
    assert len(high) == count
    return high, dict(boundary_score=boundary, above_boundary=len(high) - needed, tied_at_boundary=len(tied), ties_retained=needed, tie_seed=seed)


def source_edges(links, info, common):
    with gzip.open(info, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = list(reader)
        identifier = reader.fieldnames[0]
    mapping = {row[identifier]: row['preferred_name'] for row in rows}
    multiplicity = Counter(mapping.values())
    valid = {key: value for key, value in mapping.items() if multiplicity[value] == 1 and value in common}
    edges = {}
    directions = {}
    counts = Counter()
    with gzip.open(links, 'rt', encoding='utf-8') as handle:
        fields = next(handle).split()
        assert fields[:2] == ['protein1', 'protein2'] and fields[-1] == 'combined_score'
        for line in handle:
            values = line.split()
            counts['raw_rows'] += 1
            left, right = values[:2]
            score = int(values[-1])
            assert 0 <= score <= 1000
            if left not in valid or right not in valid:
                continue
            a, b = valid[left], valid[right]
            if a == b:
                continue
            key = (a, b) if a < b else (b, a)
            bit = 1 if a < b else 2
            assert key not in edges or edges[key] == score
            edges[key] = score
            directions[key] = directions.get(key, 0) | bit
    assert all(value == 3 for value in directions.values())
    counts['common_raw_undirected'] = len(edges)
    return edges, dict(counts)


def main():
    assert (HERE / 'PROTOCOL.txt').exists()
    fixture = {('A', 'B'): 900, ('A', 'C'): 800, ('A', 'D'): 800, ('B', 'C'): 700}
    first, receipt = count_matched(fixture, 2, 271828)
    assert len(first) == 2 and first[('A', 'B')] == 900 and receipt['ties_retained'] == 1
    assert first == count_matched(fixture, 2, 271828)[0]
    assert len({key: score for key, score in fixture.items() if score > 700}) == 3
    save('GRAPH_FIXTURES.json', {'status': 'PASS', 'fixtures': ['strict_threshold', 'boundary_tie', 'deterministic_hash', 'exact_count']})
    genes = ARGS.genes.read_text().splitlines()
    assert len(genes) == len(set(genes)) == 7399
    common = set(genes)
    sources = {
        'G0': (ARGS.old_links, ARGS.old_info),
        'G1': (ARGS.new_links, ARGS.new_info)
    }
    accepted = {'G0': {}, 'G1': {}}
    with ARGS.edge_reference.open(newline='', encoding='utf-8') as handle:
        for row in csv.DictReader(handle):
            key = (row['gene1'], row['gene2'])
            for label, column in [('G0', 'score_old'), ('G1', 'score_new')]:
                if row[column]:
                    accepted[label][key] = int(row[column])
    manifest = []
    source_log = []
    (ROOT / 'graphs').mkdir(exist_ok=True)
    for label, (links, info) in sources.items():
        edges, counts = source_edges(links, info, common)
        baseline = {key: score for key, score in edges.items() if score > 700}
        assert baseline == accepted[label]
        expected = 119423 if label == 'G0' else 65783
        assert len(baseline) == expected
        raw_path = ROOT / 'graphs' / f'{label}_raw_common.tsv.gz'
        with gzip.open(raw_path, 'wt', encoding='utf-8', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['gene1', 'gene2', 'score'])
            writer.writerows((a, b, score) for (a, b), score in sorted(edges.items()))
        variants = [(f'{label}_T700', baseline, {'kind': 'reused_baseline', 'threshold': 700, 'new_fit': False})]
        for threshold in (650, 750):
            variants.append((f'{label}_T{threshold}', {key: score for key, score in edges.items() if score > threshold}, {'kind': 'threshold', 'threshold': threshold, 'new_fit': True}))
        target = 65783 if label == 'G0' else 119423
        for seed in SEEDS:
            selected, details = count_matched(edges, target, seed)
            variants.append((f'{label}_M{target}_S{seed}', selected, {'kind': 'matched_count', 'target_edges': target, 'new_fit': True, **details}))
        for name, selected, details in variants:
            path = ROOT / 'graphs' / f'{name}.tsv'
            with path.open('w', encoding='utf-8', newline='') as handle:
                writer = csv.writer(handle, delimiter='\t')
                writer.writerow(['gene1', 'gene2', 'score'])
                writer.writerows((a, b, score) for (a, b), score in sorted(selected.items()))
            degrees = Counter(gene for key in selected for gene in key)
            manifest.append({'graph': name, 'release': label, 'path': str(path), 'edges': len(selected), 'isolates': sum(gene not in degrees for gene in genes), 'maximum_degree': max(degrees.values()), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), **details})
        source_log.append({'release': label, 'links': str(links), 'links_bytes': links.stat().st_size, 'info': str(info), 'info_bytes': info.stat().st_size, 'accepted_baseline_exact': True, **counts})
    assert len(manifest) == 12 and sum(item['new_fit'] for item in manifest) == 10
    save('GRAPH_MANIFEST.json', {'status': 'PASS', 'utc': datetime.now(timezone.utc).isoformat(), 'protocol_sha256': hashlib.sha256((HERE / 'PROTOCOL.txt').read_bytes()).hexdigest(), 'source_logs': source_log, 'graphs': manifest, 'new_fit_states': 1700, 'baseline_fit_states_reused': 340})
    with (ROOT / 'graph_metadata.csv').open('w', newline='', encoding='utf-8') as handle:
        columns = ['graph', 'release', 'kind', 'new_fit', 'edges', 'isolates', 'maximum_degree', 'threshold', 'boundary_score', 'above_boundary', 'tied_at_boundary', 'ties_retained', 'tie_seed']
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(manifest)
    print(json.dumps({'status': 'PASS', 'graphs': len(manifest), 'new_states': 1700, 'counts': [{k: row.get(k) for k in ['graph', 'edges', 'boundary_score', 'tied_at_boundary', 'ties_retained']} for row in manifest]}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rebuild threshold and count-matched graphs from specified STRING archives.')
    parser.add_argument('--genes', type=Path, required=True, help='Common-domain genes.txt')
    parser.add_argument('--old-links', type=Path, required=True)
    parser.add_argument('--old-info', type=Path, required=True)
    parser.add_argument('--new-links', type=Path, required=True)
    parser.add_argument('--new-info', type=Path, required=True)
    parser.add_argument('--edge-reference', type=Path, required=True, help='Common edge ledger with gene1, gene2, score_old and score_new columns')
    parser.add_argument('--output', type=Path, required=True, help='New, nonexistent output directory')
    ARGS = parser.parse_args()
    ROOT = ARGS.output.resolve()
    if ROOT.exists():
        raise FileExistsError('Choose a new output directory.')
    ROOT.mkdir(parents=True)
    main()

