"""Descriptive reaggregation of fixed query records; no fitting or outcome access."""
import argparse
import hashlib
import json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[2] / 'Variant_and_Network_Context/Structure/query_gene_records.json')
parser.add_argument('--output', type=Path, required=True, help='New JSON output path')
args = parser.parse_args()
source = args.source.resolve()
assert hashlib.sha256(source.read_bytes()).hexdigest() == '83bae1782a9d02a6defeece2b9c066649c68e1263947583252e9717d918f2514'
if args.output.exists():
    raise FileExistsError(args.output)
rows = json.loads(source.read_text(encoding='utf-8-sig'))
assert len(rows) == 14919
result = {}
for category, field, selected_stratum, direction in [
    ('BOTH_CORE', 'core', 'decreased', 'loss'),
    ('ONCOGENE', 'role', 'increased', 'gain'),
]:
    opportunities, flows = defaultdict(F), defaultdict(F)
    for row in rows:
        if row[field] != category:
            continue
        weight = F(row['alpha'])
        change = F(row['new_query_percentile']) - F(row['old_query_percentile'])
        stratum = 'increased' if change > 0 else 'decreased' if change < 0 else 'unchanged'
        assert stratum == row['stratum']
        opportunities[stratum] += weight
        allocation = row['allocations']
        delta = F(allocation['G1_T700_H1']) - F(allocation['G0_T700_H0'])
        flows[stratum] += weight * max(delta if direction == 'gain' else -delta, F())
    total_o, total_f = sum(opportunities.values()), sum(flows.values())
    result[category] = {
        'opportunities': {k: str(v) for k,v in opportunities.items()},
        'flow': {k: str(v) for k,v in flows.items()},
        'opportunity_share': str(opportunities[selected_stratum] / total_o),
        'flow_share': str(flows[selected_stratum] / total_f),
        'opportunity_percent': float(100 * opportunities[selected_stratum] / total_o),
        'flow_percent': float(100 * flows[selected_stratum] / total_f),
    }
assert result['BOTH_CORE']['opportunity_share'] == '199/430'
assert result['ONCOGENE']['opportunity_share'] == '692/931'
assert result['BOTH_CORE']['flow_share'] == '17/18'
assert result['ONCOGENE']['flow_share'] == '79/81'
result['source_sha256'] = hashlib.sha256(source.read_bytes()).hexdigest()
result['definition'] = 'Donor-adjusted model-query opportunities are weighted by alpha; selected flow is alpha times positive allocation difference. Descriptive only; strata are not independent observations.'
output = args.output
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))
