import argparse
import csv
import hashlib
import json
import re
import urllib.parse
from collections import Counter
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--work-dir', type=Path, default=Path(__file__).resolve().parent)
args = p.parse_args()
b = args.work_dir.resolve()

def read(name, sep=','):
    with (b / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter=sep))

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

raw = read('inputs/accepted_qualifying_variant_records.csv')
queries = read('allele_queries.csv')
alleles = read('allele_identity_coverage.csv')
records = read('record_identity_coverage.csv')
matches = read('exact_CIViC_matches.csv')
output_functional = read('accepted_Functional_evidence.csv')
source_vcfs = read('source_vcf_assembly.csv')
variant = read('sources/01-Oct-2026-VariantSummaries.tsv', '\t')
profile = read('sources/01-Oct-2026-MolecularProfileSummaries.tsv', '\t')
evidence = read('sources/01-Oct-2026-AcceptedClinicalEvidenceSummaries.tsv', '\t')
profiles = {r['molecular_profile_id']: r for r in profile}
variants = {r['variant_id']: r for r in variant}
assert len(raw) == len(records) == 184
assert len({r['gene'] for r in raw}) == 57
assert len(queries) == len(alleles) == 144
assert sum(int(r['source_record_count']) for r in alleles) == 184
assert {r['allele_key'] for r in queries} == {r['allele_key'] for r in alleles}
assert [int(r['source_row']) for r in records] == list(range(1, 185))
assert all(r['final_functional_direction'] == 'not_assigned' for r in records + alleles)
assert len(source_vcfs) == 170
assert all(r['assembly_status'] == 'verified_GRCh38' and r['matches_accepted_source_receipt'] == 'True' for r in source_vcfs)
assert Counter(r['assemblies'] for r in source_vcfs) == {'NCBI38': 93, 'GRCh38': 77}

verified_receipts = []
for folder in ['sources', 'registry_cache']:
    for receipt_path in (b / folder).glob('*.receipt.json'):
        receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
        source_path = receipt_path.with_name(receipt_path.name.removesuffix('.receipt.json'))
        assert digest(source_path) == receipt['sha256']
        assert source_path.stat().st_size == receipt['bytes']
        verified_receipts.append(str(source_path.relative_to(b)))

independent_matches = set()
registry_states = Counter()
query_times = []
for q in queries:
    cache_name = hashlib.sha256(q['public_hgvs_query'].encode()).hexdigest() + '.json'
    cache = b / 'registry_cache' / cache_name
    receipt = json.loads(cache.with_suffix('.json.receipt.json').read_text())
    url = urllib.parse.urlparse(receipt['url'])
    assert url.netloc == 'reg.clinicalgenome.org' and url.path == '/allele'
    assert urllib.parse.parse_qs(url.query) == {'hgvs': [q['public_hgvs_query']]}
    assert receipt['http_status'] == 200
    query_times.append(receipt['retrieved_utc'])
    payload = json.loads(cache.read_text(encoding='utf-8'))
    ident = payload.get('@id', '')
    registered = bool(re.fullmatch(r'http://reg.genome.network/allele/CA[0-9]+', ident))
    registry_states['registered' if registered else 'unregistered_temporary'] += 1
    r = next(r for r in alleles if r['allele_key'] == q['allele_key'])
    assert bool(r['registry_ca']) == registered
    if r['GRCh38_coordinates_verified'] != 'True':
        assert r['coverage_status'] == 'identity_unresolved'
        assert len(q['reference']) != len(q['alternative'])
        continue
    pos, ref, alt = int(q['position']), q['reference'], q['alternative']
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    ref, alt, pos = ref[prefix:], alt[prefix:], pos + prefix
    suffix = 0
    while suffix < min(len(ref), len(alt)) and ref[-suffix-1] == alt[-suffix-1]:
        suffix += 1
    if suffix:
        ref, alt = ref[:-suffix], alt[:-suffix]
    observed = [(x['start'], x['end'], x['referenceAllele'], x['allele'])
                for a in payload['genomicAlleles']
                if a.get('referenceGenome') == 'GRCh38' and a.get('chromosome') == q['chromosome'][3:]
                for x in a.get('coordinates', [])]
    assert (pos - 1, pos - 1 + len(ref), ref, alt) in observed
    if registered:
        ca = ident.split('/')[-1]
        independent_matches.update((q['allele_key'], v['variant_id']) for v in variant
                                   if v['allele_registry_id'] == ca and v['gene'] in q['genes'].split(';'))
assert independent_matches == {(r['allele_key'], r['civic_variant_id']) for r in matches}
assert all(r['match_method'] == 'G_canonical_allele' for r in matches)
assert registry_states == {'registered': 136, 'unregistered_temporary': 8}

assert len(evidence) == 4940
assert all(e['evidence_status'] == 'accepted' for e in evidence)
assert all(e['molecular_profile_id'] in profiles for e in evidence)
functional = [e for e in evidence if e['evidence_type'] == 'Functional']
assert len(functional) == 174
assert all(profiles[e['molecular_profile_id']]['variant_ids'].isdigit() for e in functional)
functional_vids = {profiles[e['molecular_profile_id']]['variant_ids'] for e in functional}
matched_vids = {v for _, v in independent_matches}
assert not (functional_vids & matched_vids)
matched_mpids = {variants[v]['single_variant_molecular_profile_id'] for v in matched_vids}
assert not [e for e in functional if e['molecular_profile_id'] in matched_mpids]
assert output_functional == []
positive_control = next(e for e in functional if e['evidence_id'] == '7114')
assert positive_control['molecular_profile'] == 'TP53 G245S'
control_profile = profiles[positive_control['molecular_profile_id']]
assert control_profile['variant_ids'] in variants
assert variants[control_profile['variant_ids']]['variant'] == 'G245S'
assert variants[control_profile['variant_ids']]['gene'] == 'TP53'

summary = {
    'status': 'PASS', 'offline_only': True, 'outcome_join_performed': False,
    'source_records': 184, 'unique_alleles': 144, 'source_genes': 57, 'source_vcf_headers_and_hashes_verified': 170,
    'cached_source_files_verified': len(verified_receipts), 'registry_request_count': len(queries),
    'registry_query_time_first_utc': min(query_times), 'registry_query_time_last_utc': max(query_times),
    'registry_identifier_status': dict(registry_states),
    'registry_verified_unique_alleles': sum(r['GRCh38_coordinates_verified'] == 'True' for r in alleles),
    'unresolved_normalized_indels': 10, 'exact_CIViC_alleles': len(matched_vids),
    'exact_CIViC_records': sum(bool(r['exact_CIViC_variant_ids']) for r in records),
    'exact_CIViC_genes': len({r['gene'] for r in records if r['exact_CIViC_variant_ids']}),
    'CIViC_global_accepted_evidence_count': len(evidence), 'CIViC_global_accepted_Functional_count': len(functional),
    'CIViC_global_Functional_with_valid_simple_profile': len(functional),
    'candidate_Functional_profile_intersection': 0,
    'CIViC_positive_control_EID': '7114', 'CIViC_positive_control_gene_variant': 'TP53 G245S',
    'CIViC_positive_control_profile_id': positive_control['molecular_profile_id'],
    'CIViC_positive_control_variant_id': control_profile['variant_ids'],
    'CIViC_accepted_evidence_latest_review_utc': max(r['last_review_date'] for r in evidence),
    'candidate_variants_without_a_profile_summary': sum(r['single_profile_id'] not in profiles for r in matches),
    'new_functional_assignments': 0,
    'code_sha256': {name: digest(b / name) for name in ['annotate_coverage.py', 'verify_coverage.py']},
    'output_sha256': {name: digest(b / name) for name in ['source_vcf_assembly.csv', 'allele_queries.csv',
        'allele_identity_coverage.csv', 'exact_CIViC_matches.csv', 'accepted_Functional_evidence.csv',
        'record_identity_coverage.csv']},
}
(b / 'independent_qc.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
print(json.dumps(summary))
