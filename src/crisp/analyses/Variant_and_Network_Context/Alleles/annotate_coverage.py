import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent
SOURCE = None
SNAPSHOT = '01-Oct-2026'
API = 'https://reg.clinicalgenome.org/allele'
ASSEMBLY_URL = ('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/'
                'GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt')

def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def sha(data):
    return hashlib.sha256(data).hexdigest()

def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

def readcsv(path, delimiter=','):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))

def writecsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        out = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        out.writeheader()
        out.writerows(rows)

def fetch(url, target):
    receipt = target.with_suffix(target.suffix + '.receipt.json')
    if target.exists() and receipt.exists():
        meta = json.loads(receipt.read_text())
        if meta['url'] != url or sha(target.read_bytes()) != meta['sha256']:
            raise ValueError('Cache receipt mismatch: ' + str(target))
        return target.read_bytes(), meta
    request = urllib.request.Request(url, headers={'User-Agent': 'Research-allele-coverage/1.0'})
    try:
        with urllib.request.urlopen(request, timeout=25) as response:
            data = response.read()
            status = response.status
            headers = dict(response.headers)
    except urllib.error.HTTPError as e:
        data, status, headers = e.read(), e.code, dict(e.headers)
    except Exception as e:
        data, status, headers = str(e).encode(), 'transport_error', {}
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
    meta = {'url': url, 'retrieved_utc': now(), 'http_status': status,
            'bytes': len(data), 'sha256': sha(data), 'headers': headers}
    dump(receipt, meta)
    return data, meta

def allele_key(row):
    return ':'.join(str(row[k]) for k in ['chromosome', 'position', 'reference', 'alternative'])

def make_hgvs(row, accessions):
    chrom, pos = row['chromosome'], int(row['position'])
    ref, alt = row['reference'], row['alternative']
    if chrom not in accessions or not re.fullmatch('[ACGT]+', ref) or not re.fullmatch('[ACGT]+', alt):
        return '', 'unsupported_contig_or_allele'
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    prefix = accessions[chrom] + ':g.'
    if len(ref) == len(alt) == 1:
        return f'{prefix}{pos}{ref}>{alt}', 'snv_reference_in_query'
    interval = str(pos) if len(ref) == 1 else f'{pos}_{pos+len(ref)-1}'
    if ref and not alt:
        return f'{prefix}{interval}del{ref}', 'deletion_reference_in_query'
    if alt and not ref:
        return f'{prefix}{pos-1}_{pos}ins{alt}', 'insertion_anchor_from_accepted_vcf'
    if ref and alt:
        return f'{prefix}{interval}del{ref}ins{alt}', 'delins_reference_in_query'
    return '', 'no_sequence_change'

def prepare():
    if SOURCE is None:
        raise ValueError('prepare requires --source-dir pointing to the Candidates folder')
    records = readcsv(SOURCE / 'inputs/accepted_qualifying_variant_records.csv')
    assert len(records) == 184
    assert len({r['gene'] for r in records}) == 57
    assert len({allele_key(r) for r in records}) == 144
    src_bytes = (SOURCE / 'inputs/accepted_qualifying_variant_records.csv').read_bytes()
    (BASE / 'inputs').mkdir(exist_ok=True)
    (BASE / 'inputs/accepted_qualifying_variant_records.csv').write_bytes(src_bytes)
    dump(BASE / 'inputs/source_receipt.json', {
        'path': str(SOURCE / 'inputs/accepted_qualifying_variant_records.csv'),
        'sha256': sha(src_bytes), 'bytes': len(src_bytes), 'record_count': len(records),
        'unique_alleles': 144, 'genes': 57,
        'protocol_sha256': sha((BASE / 'PROTOCOL.txt').read_bytes()), 'prepared_utc': now()})

    sources = json.loads((SOURCE / 'candidate_context_manifest.json').read_text())['accepted_variant_sources']
    assert len(sources) == 170
    assembly_rows = []
    contig_assemblies = {}
    for source in sources:
        original = source['original_vcf']
        path = Path(original['path'])
        header, assemblies, contigs = [], set(), {}
        with gzip.open(path, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                header.append(line)
                if line.startswith('##contig='):
                    identifier = re.search(r'ID=([^,>]+)', line)
                    assembly = re.search(r'assembly=([^,>]+)', line)
                    if identifier and assembly:
                        contigs[identifier.group(1)] = assembly.group(1)
                        assemblies.add(assembly.group(1))
        actual_hash = sha(path.read_bytes())
        provenance_ok = actual_hash == original['sha256'] and path.stat().st_size == original['bytes']
        build_ok = bool(assemblies) and assemblies <= {'NCBI38', 'GRCh38'}
        entry = {'model_id': source['model_id'], 'source_kind': source['kind'],
                 'source_vcf': str(path), 'source_bytes': path.stat().st_size,
                 'sha256': actual_hash, 'matches_accepted_source_receipt': provenance_ok,
                 'assemblies': ';'.join(sorted(assemblies)), 'contig_count': len(contigs),
                 'header_sha256': sha(''.join(header).encode()),
                 'assembly_status': 'verified_GRCh38' if build_ok and provenance_ok else 'unresolved'}
        assembly_rows.append(entry)
        contig_assemblies[(source['model_id'], source['kind'])] = contigs
    writecsv(BASE / 'source_vcf_assembly.csv', assembly_rows)
    assert all(r['assembly_status'] == 'verified_GRCh38' for r in assembly_rows)
    for r in records:
        assert contig_assemblies[(r['model_id'], r['source_kind'])].get(r['chromosome']) in {'NCBI38', 'GRCh38'}

    data, receipt = fetch(ASSEMBLY_URL, BASE / 'sources/GRCh38.p14_assembly_report.txt')
    assert receipt['http_status'] == 200, receipt
    accessions = {}
    for line in data.decode().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        fields = line.split('\t')
        if fields[1] == 'assembled-molecule' and fields[9].startswith('chr'):
            accessions[fields[9]] = fields[6]
    assert accessions['chr12'] == 'NC_000012.12'
    dump(BASE / 'sources/GRCh38_chromosome_accessions.json', accessions)

    for name in ['VariantSummaries.tsv', 'MolecularProfileSummaries.tsv',
                 'AcceptedClinicalEvidenceSummaries.tsv', 'AcceptedAssertionSummaries.tsv',
                 'civic_accepted.vcf']:
        filename = SNAPSHOT + '-' + name
        _, meta = fetch('https://civicdb.org/downloads/' + SNAPSHOT + '/' + filename,
                        BASE / 'sources' / filename)
        assert meta['http_status'] == 200, meta
    for slug, url in {
        'CIViC_licensing.html': 'https://civicdb.org/releases/licensing',
        'CIViC_functional_rules.html': 'https://docs.civicdb.org/en/latest/curating/evidence.html',
        'CIViC_significance.html': 'https://docs.civicdb.org/en/latest/model/evidence/significance.html',
        'CIViC_direction.html': 'https://docs.civicdb.org/en/latest/model/evidence/direction.html',
    }.items():
        fetch(url, BASE / 'sources' / slug)

    unique = defaultdict(list)
    for i, r in enumerate(records, 1):
        unique[allele_key(r)].append((i, r))
    queries = []
    for key, members in sorted(unique.items()):
        r = members[0][1]
        hgvs, status = make_hgvs(r, accessions)
        queries.append({'allele_key': key, 'chromosome': r['chromosome'], 'position': r['position'],
                        'reference': r['reference'], 'alternative': r['alternative'],
                        'genes': ';'.join(sorted({x['gene'] for _, x in members})),
                        'transcripts': ';'.join(sorted({x['transcript'] for _, x in members})),
                        'protein_changes': ';'.join(sorted({x['source_protein_change'] for _, x in members})),
                        'source_record_count': len(members),
                        'source_row_numbers': ';'.join(str(i) for i, _ in members),
                        'assembly': 'GRCh38', 'public_hgvs_query': hgvs, 'hgvs_conversion_status': status})
    writecsv(BASE / 'allele_queries.csv', queries)
    print(json.dumps({'prepared': len(queries), 'source_vcfs_verified': len(assembly_rows)}), flush=True)

def query(limit=None):
    queries = readcsv(BASE / 'allele_queries.csv')
    for idx, q in enumerate(queries[:limit], 1):
        hgvs = q['public_hgvs_query']
        if not hgvs:
            continue
        key = sha(hgvs.encode())
        cache = BASE / 'registry_cache' / (key + '.json')
        existing = cache.exists()
        url = API + '?' + urllib.parse.urlencode({'hgvs': hgvs})
        _, meta = fetch(url, cache)
        print(json.dumps({'query': idx, 'total': len(queries), 'http_status': meta['http_status'],
                          'cache_reused': existing}), flush=True)
        if not existing:
            time.sleep(0.55)

def summarize():
    for directory in ['sources', 'registry_cache']:
        for receipt_path in (BASE / directory).glob('*.receipt.json'):
            receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
            path = receipt_path.with_name(receipt_path.name.removesuffix('.receipt.json'))
            assert sha(path.read_bytes()) == receipt['sha256'], str(path)
    source_receipt = json.loads((BASE / 'inputs/source_receipt.json').read_text(encoding='utf-8'))
    assert sha((BASE / 'inputs/accepted_qualifying_variant_records.csv').read_bytes()) == source_receipt['sha256']
    assert sha((BASE / 'PROTOCOL.txt').read_bytes()) == source_receipt.get('packaged_protocol_sha256', source_receipt['protocol_sha256'])
    queries = readcsv(BASE / 'allele_queries.csv')
    records = readcsv(BASE / 'inputs/accepted_qualifying_variant_records.csv')
    variants = readcsv(BASE / 'sources' / (SNAPSHOT + '-VariantSummaries.tsv'), '\t')
    profiles = {m['molecular_profile_id']: m for m in readcsv(
        BASE / 'sources' / (SNAPSHOT + '-MolecularProfileSummaries.tsv'), '\t')}
    evidence = readcsv(BASE / 'sources' / (SNAPSHOT + '-AcceptedClinicalEvidenceSummaries.tsv'), '\t')
    assert all(e['evidence_status'] == 'accepted' for e in evidence)
    evidence_by_profile = defaultdict(list)
    for e in evidence:
        evidence_by_profile[e['molecular_profile_id']].append(e)
    variants_by_ca = defaultdict(list)
    for v in variants:
        if re.fullmatch(r'CA[0-9]+', v['allele_registry_id']):
            variants_by_ca[v['allele_registry_id']].append(v)

    allele_rows, evidence_rows, match_rows = [], [], []
    for q in queries:
        path = BASE / 'registry_cache' / (sha(q['public_hgvs_query'].encode()) + '.json')
        receipt = json.loads(path.with_suffix('.json.receipt.json').read_text()) if path.exists() else {}
        try:
            payload = json.loads(path.read_text()) if path.exists() else {}
        except (ValueError, UnicodeDecodeError):
            payload = {}
        identity_url = payload.get('@id', '')
        ca = identity_url.rsplit('/', 1)[-1] if re.search(r'/CA[0-9]+$', identity_url) else ''
        pos = int(q['position'])
        ref, alt = q['reference'], q['alternative']
        while ref and alt and ref[0] == alt[0]:
            ref, alt, pos = ref[1:], alt[1:], pos + 1
        while ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        expected = {'start': pos - 1, 'end': pos - 1 + len(ref), 'referenceAllele': ref, 'allele': alt}
        build38 = [a for a in payload.get('genomicAlleles', [])
                   if a.get('referenceGenome') == 'GRCh38' and
                   a.get('chromosome') == q['chromosome'].removeprefix('chr')]
        identity_verified = any(a.get('coordinates') == [expected] for a in build38)
        matches = {}
        if identity_verified:
            if ca:
                for v in variants_by_ca.get(ca, []):
                    if v['gene'] in q['genes'].split(';'):
                        matches[v['variant_id']] = (v, 'G_canonical_allele')
            for a in payload.get('genomicAlleles', []):
                if a.get('referenceGenome') not in {'GRCh37', 'GRCh38'} or len(a.get('coordinates', [])) != 1:
                    continue
                c = a['coordinates'][0]
                if len(c.get('referenceAllele', '')) != 1 or len(c.get('allele', '')) != 1:
                    continue
                for v in variants:
                    if v['variant_id'] in matches or v['gene'] not in q['genes'].split(';'):
                        continue
                    if (v['reference_build'] == a['referenceGenome'] and
                        v['chromosome'].removeprefix('chr') == a.get('chromosome') and
                        v['start'] == str(c['start'] + 1) and v['stop'] == str(c['end']) and
                        v['reference_bases'] == c['referenceAllele'] and v['variant_bases'] == c['allele']):
                        matches[v['variant_id']] = (v, 'G_exact_registry_mapped_SNV')
        all_evidence, functional = {}, {}
        for vid, (v, method) in matches.items():
            mpid = v['single_variant_molecular_profile_id']
            profile = profiles.get(mpid, {})
            simple = profile.get('variant_ids', '').strip() == vid
            match_rows.append({'allele_key': q['allele_key'], 'source_genes': q['genes'],
                'source_protein_changes': q['protein_changes'], 'match_method': method,
                'civic_variant_id': vid, 'civic_gene': v['gene'], 'civic_variant': v['variant'],
                'civic_url': v['variant_civic_url'], 'registry_ca': ca,
                'civic_ca': v['allele_registry_id'], 'civic_build': v['reference_build'],
                'civic_chromosome': v['chromosome'], 'civic_start': v['start'], 'civic_stop': v['stop'],
                'civic_ref': v['reference_bases'], 'civic_alt': v['variant_bases'],
                'civic_transcript': v['representative_transcript'], 'single_profile_id': mpid,
                'profile_variant_ids': profile.get('variant_ids', ''), 'verified_simple_profile': simple,
                'variant_flagged': v['is_flagged'], 'profile_flagged': profile.get('is_flagged', '')})
            if not simple:
                continue
            for e in evidence_by_profile.get(mpid, []):
                all_evidence[e['evidence_id']] = e
                if e['evidence_type'] != 'Functional':
                    continue
                functional[e['evidence_id']] = e
                evidence_rows.append({'allele_key': q['allele_key'], 'genes': q['genes'],
                    'source_protein_changes': q['protein_changes'], 'civic_variant_id': vid,
                    'civic_variant': v['variant'], 'match_method': method,
                    **e, 'variant_flagged': v['is_flagged'], 'profile_flagged': profile.get('is_flagged', ''),
                    'primary_experiment_verification': 'pending',
                    'functional_assignment_status': 'database_curated_coverage_only'})
        if not path.exists():
            status = 'query_pending'
        elif not identity_verified:
            status = 'identity_unresolved'
        elif not matches:
            status = 'exact_genomic_identity_no_CIViC_match'
        elif not functional:
            status = 'exact_CIViC_match_no_accepted_Functional'
        else:
            status = 'accepted_Functional_coverage_pending_primary_verification'
        allele_rows.append({**q, 'registry_http_status': receipt.get('http_status', ''),
            'registry_ca': ca, 'registry_returned_id': identity_url,
            'GRCh38_coordinates_verified': identity_verified,
            'identity_unresolved_reason': '' if identity_verified else ('registry_normalized_representation_not_locally_validated' if build38 else 'no_verified_GRCh38_representation'),
            'registry_error_type': payload.get('errorType', ''),
            'registry_error_description': payload.get('description', ''),
            'exact_CIViC_variant_ids': ';'.join(sorted(matches)),
            'exact_CIViC_variant_names': ';'.join(v['gene']+' '+v['variant'] for v, _ in matches.values()),
            'match_methods': ';'.join(sorted({method for _, method in matches.values()})),
            'accepted_all_evidence_ids': ';'.join(sorted(all_evidence)),
            'accepted_Functional_evidence_ids': ';'.join(sorted(functional)),
            'accepted_Functional_evidence_count': len(functional),
            'Functional_significance_and_direction': ';'.join(sorted({e['significance']+' / '+e['evidence_direction'] for e in functional.values()})),
            'coverage_status': status, 'final_functional_direction': 'not_assigned'})
    writecsv(BASE / 'allele_identity_coverage.csv', allele_rows)
    writecsv(BASE / 'exact_CIViC_matches.csv', match_rows, list(match_rows[0]) if match_rows else ['allele_key','civic_variant_id'])
    fields = list(evidence_rows[0]) if evidence_rows else ['allele_key','genes','evidence_id','primary_experiment_verification']
    writecsv(BASE / 'accepted_Functional_evidence.csv', evidence_rows, fields)
    bykey = {r['allele_key']: r for r in allele_rows}
    record_rows = []
    for i, r in enumerate(records, 1):
        x = bykey[allele_key(r)]
        record_rows.append({'source_row': i, 'model_id': r['model_id'], 'donor_id': r['donor_id'],
            'gene': r['gene'], 'source_kind': r['source_kind'], 'source_record': r['source_record'],
            'allele_key': allele_key(r), 'source_transcript': r['transcript'],
            'source_protein_change': r['source_protein_change'],
            'coverage_status': x['coverage_status'], 'registry_ca': x['registry_ca'],
            'exact_CIViC_variant_ids': x['exact_CIViC_variant_ids'],
            'accepted_Functional_evidence_ids': x['accepted_Functional_evidence_ids'],
            'final_functional_direction': 'not_assigned'})
    writecsv(BASE / 'record_identity_coverage.csv', record_rows)
    summary = {'created_utc': now(), 'snapshot': SNAPSHOT, 'outcome_join_performed': False,
        'source_records': len(records), 'unique_alleles': len(allele_rows), 'genes': len({r['gene'] for r in records}),
        'unique_allele_status': dict(Counter(r['coverage_status'] for r in allele_rows)),
        'source_record_status': dict(Counter(r['coverage_status'] for r in record_rows)),
        'registered_CA_identifiers': sum(bool(r['registry_ca']) for r in allele_rows),
        'exact_CIViC_matched_alleles': sum(bool(r['exact_CIViC_variant_ids']) for r in allele_rows),
        'exact_CIViC_matched_records': sum(bool(r['exact_CIViC_variant_ids']) for r in record_rows),
        'exact_CIViC_matched_genes': len({r['gene'] for r in record_rows if r['exact_CIViC_variant_ids']}),
        'Functional_covered_alleles': sum(bool(r['accepted_Functional_evidence_ids']) for r in allele_rows),
        'Functional_covered_records': sum(bool(r['accepted_Functional_evidence_ids']) for r in record_rows),
        'Functional_covered_genes': len({r['gene'] for r in record_rows if r['accepted_Functional_evidence_ids']}),
        'distinct_Functional_evidence_ids': len({r['evidence_id'] for r in evidence_rows}),
        'primary_experiments_verified': 0, 'final_functional_directions_assigned': 0,
        'CIViC_release_variant_rows': len(variants), 'CIViC_release_accepted_evidence_rows': len(evidence),
        'CIViC_global_builds': dict(Counter(r['reference_build'] for r in variants))}
    dump(BASE / 'coverage_summary.json', summary)
    assert len(allele_rows) == 144 and len(record_rows) == 184 and len(bykey) == 144
    assert sum(Counter(r['coverage_status'] for r in record_rows).values()) == 184
    print(json.dumps(summary), flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('command', choices=['prepare', 'query', 'summarize'])
    p.add_argument('--limit', type=int)
    p.add_argument('--source-dir', type=Path)
    p.add_argument('--output-dir', type=Path)
    args = p.parse_args()
    if args.output_dir is not None:
        BASE = args.output_dir.resolve()
    if args.source_dir is not None:
        SOURCE = args.source_dir.resolve()
    if args.command == 'prepare': prepare()
    elif args.command == 'query': query(args.limit)
    else: summarize()
