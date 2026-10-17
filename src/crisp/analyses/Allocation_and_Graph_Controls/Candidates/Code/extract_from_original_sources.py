"""Fixed-ledger candidate context for analysis. No fitting or reranking.

Reuses stored allocations and accepted input-only VCF decisions. A consequence or
gene-role annotation is never converted into a gain/loss-of-function assertion.
Python standard library only. The source extraction pass records every relevant
qualifying allele and retains all original source identifiers.
"""
import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path


def read_csv(p):
    with p.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def read_json(p):
    op = gzip.open if str(p).endswith('.gz') else open
    with op(p, 'rt', encoding='utf-8-sig') as f:
        return json.load(f, parse_float=str)


def write_csv(p, rows):
    assert rows, p
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with p.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fields, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def ident(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return {'path': p.as_posix(), 'bytes': p.stat().st_size, 'sha256': h.hexdigest()}


def weightdict(w):
    return {g: F(x) for g, x in (w.items() if isinstance(w, dict) else w)}


def allocations(archive):
    result = {}
    for method, path in [('PRODIGY', 'results/crc_prodigy_allocations.json'),
                         ('DawnRank', 'results/crc_precision/allocations.json')]:
        tmp = defaultdict(dict)
        for r in read_json(archive/path)['allocations']:
            if r['budget'] != 10 or r.get('method', method) != method:
                continue
            cell = r.get('cell', r.get('state')).replace('_', '')
            if cell not in ('G0H0', 'G1H1'):
                continue
            seed = str(r['seed']) if method == 'PRODIGY' else 'deterministic'
            key = method, r['model_id'], cell
            v = weightdict(r['weights'])
            assert seed not in tmp[key] or tmp[key][seed] == v
            tmp[key][seed] = v
        for key, values in tmp.items():
            assert len(values) == (5 if method == 'PRODIGY' else 1)
            genes = set().union(*(v.keys() for v in values.values()))
            result[key] = {g: sum((v.get(g, F()) for v in values.values()), F())/len(values) for g in genes}
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True, help='Directory containing the colorectal functional analysis CSVs')
    ap.add_argument('--archive', type=Path, required=True, help='Extracted base research-data archive')
    ap.add_argument('--variants', type=Path, required=True)
    ap.add_argument('--assay', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    source = args.source
    archive = args.archive
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    flows = read_csv(source/'oncogene_identity_flows.csv')
    lfcrows = read_csv(source/'current_continuous_selected_identities.csv')
    lfcmap = {(r['model_id'], r['gene']): r for r in lfcrows}
    oldpositive = read_csv(source/'current_continuous_changed_positive_context.csv')
    assert len(oldpositive) == 34
    degree = read_json(archive/'results/crc_precision/degree_allocations.json')
    models = degree['models']
    assert len(models) == 85 and len({r['donor_id'] for r in models}) == 83
    byid = {r['model_id']: r for r in models}
    query = {r['model_id']: set(r['query']) for r in models}
    vectors = allocations(archive)
    selected = {(r['model_id'], r['gene']) for r in flows}
    assert selected <= set(lfcmap)
    assert all(g in query[m] for m, g in selected)
    nodes = ('EGFR', 'KRAS', 'PTPN11')
    target = selected | {(m['model_id'], g) for m in models for g in nodes}
    degrees = {}
    common = set((archive/'data/crc/first/common_genes.txt').read_text().splitlines())
    assert len(common) == 7399
    for state in ('G0', 'G1'):
        c = Counter()
        with (archive/f'data/crc/first/edges_{state}.tsv').open(newline='') as f:
            seen = set()
            for r in csv.DictReader(f, delimiter='\t'):
                if r['gene1'] not in common or r['gene2'] not in common:
                    continue
                edge = tuple(sorted((r['gene1'], r['gene2'])))
                assert edge not in seen and edge[0] != edge[1]
                seen.add(edge)
                c.update(edge)
            assert len(seen) == {'G0':119423, 'G1':65783}[state]
        degrees[state] = c
    variantrows = []
    ledger_sources = []
    for m in models:
        mid = m['model_id']
        wanted = {g for mm,g in target if mm == mid}
        for kind in ('snp','indel'):
            path = args.variants/f'{mid}.{kind}.csv.gz'
            meta = read_json(args.variants/f'{mid}.{kind}.json')
            receipt = meta['ledger']
            assert path.stat().st_size == receipt['bytes']
            ledger_sources.append({'model_id':mid,'kind':kind, 'accepted_ledger':receipt,
                                   'original_vcf':meta['source']['source_identity_reused']})
            with gzip.open(path,'rt',newline='') as f:
                for r in csv.DictReader(f):
                    if r['source_symbol'] not in wanted or r['eligible'] != '1':
                        continue
                    assert not r['reasons']
                    fmt = dict(zip(r['format'].split(':'),r['organoid_values'].split(':')))
                    vd = r['VD'].split('|')
                    protein = [x for x in vd if x.startswith('p.')]
                    variantrows.append({'model_id':mid,'sample_id':m['sample_id'],
                        'donor_id':m['donor_id'],'source_kind':kind,'gene':r['source_symbol'],
                        **r, 'source_protein_change':'|'.join(protein),
                        'organoid_VAF':fmt.get('VAF',''),
                        'organoid_sample':meta['source']['organoid'],
                        'normal_sample':meta['source']['normal'],
                        'functional_direction':'not_assigned',
                        'function_basis':'Source sequence consequence only; no allele-level functional classification'})
        print('variant ledgers read', mid, flush=True)
    byvar = defaultdict(list)
    for r in variantrows:
        byvar[(r['model_id'],r['gene'])].append(r)
    assert all(byvar[key] for key in selected), 'Selected mutation-root gene lacks source allele'
    continuous = read_json(args.assay/'table6_MAPPED_VALUES.json.gz')
    binary = read_json(args.assay/'table5_MAPPED_VALUES.json.gz')

    def context(mid,g):
        m = byid[mid]
        vv = byvar[(mid,g)]
        return {'query_eligible':int(g in query[mid]),
            'qualifying_source_records':len(vv),
            'source_consequences':'|'.join(sorted({v['VC'] for v in vv})),
            'source_protein_changes':'|'.join(sorted({v['source_protein_change'] for v in vv if v['source_protein_change']})),
            'source_VD_annotations':json.dumps([v['VD'] for v in vv],separators=(',',':')),
            'source_alleles':json.dumps([f"{v['chromosome']}:{v['position']}:{v['reference']}>{v['alternative']}" for v in vv],separators=(',',':')),
            'source_VAFs':'|'.join(v['organoid_VAF'] for v in vv),
            'variant_level_function':'not_assigned',
            'common_domain_degree_G0':degrees['G0'][g],
            'common_domain_degree_G1':degrees['G1'][g],
            'degree_old_allocation':str(F(m['full_degree']['G0'].get(g,'0'))),
            'degree_new_allocation':str(F(m['full_degree']['G1'].get(g,'0')))}

    identityrows=[]
    for r in flows:
        mid,g = r['model_id'],r['gene']
        lc = lfcmap[(mid,g)]
        assert lc['binary_state']==r['label_state']
        assert lc['binary_call']==r['label_value']
        for state,col in [('G0H0','old'),('G1H1','new')]:
            assert vectors[(r['method'],mid,state)].get(g,F()) == F(r[col])
        identityrows.append({**r, **{k:lc[k] for k in ['continuous_state','supplied_LFC','negative_LFC','retained_library','source_rows','raw_representations']}, **context(mid,g)})
    positivecontext=[{**r,**context(r['model_id'],r['gene'])} for r in oldpositive]
    zeroids=[r for r in identityrows if r['label_group']=='zero']
    assert all(r['continuous_state']=='measured' for r in identityrows)
    genes=[]
    allgenes=sorted({r['gene'] for r in identityrows})
    for method in ('PRODIGY','DawnRank'):
        for g in allgenes:
            rr=[r for r in zeroids if r['method']==method and r['gene']==g]
            row={'method':method,'gene':g,'common_domain_degree_G0':degrees['G0'][g],
                 'common_domain_degree_G1':degrees['G1'][g], 'selected_zero_models':len(rr),
                 'gain_models':sum(F(r['gain'])>0 for r in rr),
                 'loss_models':sum(F(r['loss'])>0 for r in rr)}
            for k in ('old','new','gain','loss','net'):
                value=sum((F(r['donor_weight'])*F(r[k]) for r in rr),F())
                row[k]=str(value)
                row[k+'_decimal']=float(value)
            for k in ('gain','loss'):
                A=sum((F(r['donor_weight'])*F(r[k]) for r in rr),F())
                S=sum((F(r['donor_weight'])*F(r[k])*F(r['negative_LFC']) for r in rr),F())
                row[k+'_negative_LFC_sum']=str(S)
                row[k+'_mean_negative_LFC']=float(S/A) if A else ''
            genes.append(row)
    generef=read_csv(source/'oncogene_gene_flows.csv')
    for row in genes:
        matches=[r for r in generef if r['method']==row['method'] and r['gene']==row['gene'] and r['label_group']=='zero']
        assert len(matches)==1
        for k in ('old','new','gain','loss','net'):
            assert F(row[k])==F(matches[0][k]), (row,k)
    tops={}
    for method in ('PRODIGY','DawnRank'):
        rr=sorted([r for r in genes if r['method']==method and F(r['net'])>0],key=lambda r:(-F(r['net']),r['gene']))
        tops[method]=[r['gene'] for r in rr[:5]]
    union=sorted(set().union(*map(set,tops.values())))
    toprows=[r for r in genes if r['gene'] in union]
    nodecontext=[]
    for m in models:
        for g in nodes:
            mid,sample=m['model_id'],m['sample_id']
            b = binary.get(sample,{}).get(g,'absent')
            val = continuous.get(sample,{}).get(g)
            row={'model_id':mid,'sample_id':sample,'donor_id':m['donor_id'],'gene':g,
                 'binary_source_value':b if b is not None else 'null',
                 'supplied_LFC':val if val is not None else '', **context(mid,g)}
            for method in ('PRODIGY','DawnRank'):
                for state,label in [('G0H0','old'),('G1H1','new')]:
                    row[f'{method}_{label}']=str(vectors[(method,mid,state)].get(g,F()))
            nodecontext.append(row)
    assert len(nodecontext)==255
    summary={}
    for method in ('PRODIGY','DawnRank'):
        rr=[r for r in genes if r['method']==method]
        net=sum((F(r['net']) for r in rr),F())
        positive_net=sum((max(F(),F(r['net'])) for r in rr),F())
        kept=sum((F(r['net']) for r in rr if r['gene'] in union),F())
        summary[method]={'zero_net':str(net),'positive_gene_net':str(positive_net),
            'table_union_net':str(kept),'table_union_share_total_net':float(kept/net),
            'table_union_share_positive_gene_net':float(kept/positive_net),'top5':tops[method],
            'selected_zero_identities':sum(r['method']==method for r in zeroids)}
    nodestats=[]
    for g in nodes:
        rr=[r for r in nodecontext if r['gene']==g]
        nr={'gene':g,'models':85,'query_eligible_models':sum(r['query_eligible'] for r in rr),
            'measured_models':sum(r['supplied_LFC']!='' for r in rr),
            'source_positive_models':sum(str(r['binary_source_value'])=='1' for r in rr),
            'positive_but_query_ineligible':sum(str(r['binary_source_value'])=='1' and not r['query_eligible'] for r in rr)}
        for method in ('PRODIGY','DawnRank'):
            nr[method+'_gained_models']=sum(F(r[method+'_new'])>F(r[method+'_old']) for r in rr)
            nr[method+'_lost_models']=sum(F(r[method+'_new'])<F(r[method+'_old']) for r in rr)
        nodestats.append(nr)
    outputs={'oncogene_all_identity_context.csv':identityrows,
             'oncogene_zero_identity_context.csv':zeroids,
             'oncogene_zero_gene_contributors.csv':genes,
             'oncogene_zero_main_table.csv':toprows,
             'changed_positive_34_context.csv':positivecontext,
             'qualifying_variant_context.csv':variantrows,
             'three_node_full_cohort_context.csv':nodecontext,
             'three_node_summary.csv':nodestats}
    for name,rr in outputs.items():
        write_csv(out/name,rr)
    sources=[ident(source/p) for p in ['oncogene_identity_flows.csv','oncogene_gene_flows.csv',
        'current_continuous_selected_identities.csv','current_continuous_changed_positive_context.csv']]
    sources += [ident(archive/p) for p in ['results/crc_precision/degree_allocations.json','results/crc_prodigy_allocations.json','results/crc_precision/allocations.json','data/crc/first/edges_G0.tsv','data/crc/first/edges_G1.tsv']]
    report={'status':'COMPLETE','models':85,'donors':83,'changed_positive_identities':34,
        'selected_oncogene_identities':len(selected),'source_qualifying_records':len(variantrows),
        'identity_context_rows':len(identityrows),'zero_identity_rows':len(zeroids),
        'main_table_union':union,'method_summary':summary,'three_node_summary':nodestats,
        'source_inputs':sources,'accepted_variant_sources':ledger_sources,
        'assay_inputs':[ident(args.assay/p) for p in ['table5_MAPPED_VALUES.json.gz','table6_MAPPED_VALUES.json.gz']],
        'notes':['No new ranking, label, query or variant eligibility was calculated.',
            'All variant functional directions remain unassigned.',
            'Common-domain degree is descriptive; PRODIGY uses a wider weighted graph.',
            'Three nodes are an explanatory subset, not a complete pathway or independent validation.'],
        'outputs':[ident(out/name) for name in outputs]}
    (out/'candidate_context_manifest.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['status','selected_oncogene_identities','source_qualifying_records','main_table_union','method_summary','three_node_summary']},indent=2))


if __name__=='__main__':
    main()

