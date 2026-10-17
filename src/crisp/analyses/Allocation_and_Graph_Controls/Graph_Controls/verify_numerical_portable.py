import argparse
import csv
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy as np
import scipy
from scipy.sparse import csc_matrix

ROOT = Path(__file__).resolve().parent
OLD = ROOT / 'runtime/reference_allocations/numerical'


def read_csv(path):
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def independent_bound(operator, damping, q, x):
    u = 2.0 ** -53
    gamma = lambda n: n * u / (1 - n * u)
    row_nnz = np.diff(operator.indptr)
    col_nnz = np.diff(operator.tocsc().indptr)
    signal = damping * operator.dot(x)
    restart = (1 - damping) * q
    residual = signal + restart - x
    allowance = gamma(row_nnz + 3) / (1 - gamma(row_nnz + 3)) * abs(signal)
    allowance += gamma(2) / (1 - gamma(2)) * abs(restart)
    allowance += gamma(3) / (1 - gamma(3)) * (abs(signal) + abs(restart) + abs(x))
    rounding = np.sum(allowance) / (1 - gamma(len(x)))
    norm = float(np.max(operator.multiply(damping[:, None]).sum(axis=0))) / (1 - gamma(max(col_nnz) + 3))
    assert norm < 1 and np.all(np.isfinite(x)) and np.all(x >= 0)
    bound = (np.sum(abs(residual)) / (1 - gamma(len(x))) + rounding) / (1 - norm)
    return float(bound), norm


def main():
    global ROOT, OLD
    parser = argparse.ArgumentParser(description='Check sparse fixed-point values and ranking bounds.')
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    ROOT = args.directory.resolve()
    OLD = ROOT / 'runtime/reference_allocations/numerical'
    if args.output.exists() or args.output.resolve().is_relative_to(ROOT):
        raise ValueError('Choose a new output path outside the input directory.')
    native = json.loads((ROOT / 'NATIVE_EXIT.json').read_text())
    assert native['states'] == 1700
    manifest = json.loads((ROOT / 'GRAPH_MANIFEST.json').read_text())['graphs']
    models = (OLD / 'models.txt').read_text().splitlines()
    mutations = read_csv(OLD / 'mutation.csv')
    eligible = np.array([[int(row[sample]) for sample in models] for row in mutations], dtype=bool)
    diagnostics = {(row['state'], row['sample_id']): row for row in read_csv(ROOT / 'cutoff_certificates.csv')}
    checks = []
    for graph in manifest:
        if not graph['new_fit']:
            continue
        name = graph['graph']
        matrix = csc_matrix((np.fromfile(ROOT / 'native' / f'{name}_x.bin', dtype='<f8'), np.fromfile(ROOT / 'native' / f'{name}_i.bin', dtype='<i4'), np.fromfile(ROOT / 'native' / f'{name}_p.bin', dtype='<i4')), shape=(7399, 7399)).tocsr()
        damping = np.fromfile(ROOT / 'native' / f'{name}_d.bin', dtype='<f8')
        for h in (0, 1):
            state = f'{name}_H{h}'
            q = np.fromfile(ROOT / 'native' / f'H{h}_q.bin', dtype='<f8').reshape((7399, 85), order='F')
            expected_q = np.fromfile(OLD / f'H{h}_q.bin', dtype='<f8').reshape((7399, 85), order='F')
            assert np.array_equal(q, expected_q)
            native_scores = np.fromfile(ROOT / 'native' / f'{state}_scores.bin', dtype='<f8').reshape((7399, 85), order='F')
            current = q.copy()
            constant = (1 - damping[:, None]) * q
            active = np.ones(85, dtype=bool)
            iterations = np.zeros(85, dtype=int)
            for iteration in range(1, 10001):
                ids = np.flatnonzero(active)
                if not len(ids):
                    break
                update = damping[:, None] * matrix.dot(current[:, ids]) + constant[:, ids]
                movement = np.sqrt(np.sum((update - current[:, ids]) ** 2, axis=0))
                current[:, ids] = update
                stopped = ids[movement < 1e-12]
                iterations[stopped] = iteration
                active[stopped] = False
            assert not np.any(active)
            for index, sample in enumerate(models):
                record = diagnostics[(state, sample)]
                delta = abs(current[:, index] - native_scores[:, index])
                assert np.all(delta <= 1e-12 + 1e-10 * abs(native_scores[:, index]))
                assert iterations[index] == int(record['iterations'])
                ids = np.flatnonzero(eligible[:, index])
                rank1 = ids[np.argsort(-current[ids, index], kind='stable')]
                rank2 = ids[np.argsort(-native_scores[ids, index], kind='stable')]
                assert np.array_equal(rank1, rank2)
                bound, norm = independent_bound(matrix, damping, q[:, index], native_scores[:, index])
                stored = float(record['error_bound'])
                assert abs(bound - stored) <= max(1e-18, 1e-4 * stored)
                assert abs(norm - float(record['operator_l1_upper'])) < 1e-14
                gap = float(native_scores[rank1[9], index] - native_scores[rank1[10], index]) if len(rank1) > 10 else None
                certified = gap is None or gap > 2 * max(bound, stored)
                assert certified == (record['certified'] == 'TRUE')
                checks.append(dict(state=state, sample_id=sample, max_absolute_difference=float(delta.max()), iterations=int(iterations[index]), bound=bound, full_query_order_identical=True, cutoff_certified=bool(certified)))
            print(state, len(checks), flush=True)
    assert len(checks) == 1700
    report = dict(status='PASS', utc=datetime.now(timezone.utc).isoformat(), python=sys.version, numpy=np.__version__, scipy=scipy.__version__, vectors=len(checks), score_count=7399 * len(checks), maximum_absolute_difference=max(row['max_absolute_difference'] for row in checks), certified=sum(row['cutoff_certified'] for row in checks), complete_query_orders_identical=True, checks=checks)
    args.output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({key: value for key, value in report.items() if key != 'checks'}))


if __name__ == '__main__':
    main()
