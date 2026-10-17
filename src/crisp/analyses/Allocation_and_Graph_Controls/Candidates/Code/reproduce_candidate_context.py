"""Rebuild all eight candidate-context CSVs using packaged inputs only.

Python >=3.9 standard library. Run from any directory:
  python reproduce_candidate_context.py --source-data PATH --output NEW_DIRECTORY
Source-data contains portable_manifest.json and inputs/. All divisions use exact
rational arithmetic.
"""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path


def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path,data):
    fields=list(dict.fromkeys(k for r in data for k in r))
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,lineterminator='\n')
        w.writeheader()
        w.writerows(data)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-data',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    src=a.source_data.resolve()
    out=a.output.resolve()
    assert out!=src and out!=src/'inputs','Use a separate output directory'
    manifest=json.loads((src/'portable_manifest.json').read_text())
    for r in manifest['inputs']:
        p=src/r['path']
        assert p.is_file() and p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],p
    inp=src/'inputs'
    s=json.loads((inp/'source_context.json').read_text())
    models=s['models']
    assert len(models)==85 and len({r['donor_id'] for r in models})==83
    byid={m['model_id']:m for m in models}
    query={m['model_id']:set(m['query']) for m in models}
    degrees=s['common_domain_degree']
    vectors={(r['method'],r['model_id'],r['cell']):{g:F(v) for g,v in r['weights'].items()} for r in s['allocations']}
    flows=rows(inp/'oncogene_identity_flows.csv')
    lfcmap={(r['model_id'],r['gene']):r for r in rows(inp/'current_continuous_selected_identities.csv')}
    positive=rows(inp/'current_continuous_changed_positive_context.csv')
    assert len(positive)==34
    variantrows=rows(inp/'accepted_qualifying_variant_records.csv')
    byvar=defaultdict(list)
    for r in variantrows:
        assert r['eligible']=='1' and not r['reasons']
        assert r['functional_direction']=='not_assigned'
        byvar[(r['model_id'],r['gene'])].append(r)
    assert all(r['gene'] in query[r['model_id']] and byvar[(r['model_id'],r['gene'])] for r in flows)

    def context(mid,g):
        vv=byvar[(mid,g)]
        m=byid[mid]
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

    identities=[]
    for r in flows:
        mid,g=r['model_id'],r['gene']
        lc=lfcmap[(mid,g)]
        assert lc['binary_state']==r['label_state'] and lc['binary_call']==r['label_value']
        assert lc['continuous_state']=='measured'
        for state,col in [('G0H0','old'),('G1H1','new')]:
            assert vectors[(r['method'],mid,state)].get(g,F())==F(r[col])
        identities.append({**r,**{k:lc[k] for k in ['continuous_state','supplied_LFC','negative_LFC','retained_library','source_rows','raw_representations']},**context(mid,g)})
    positivecontext=[{**r,**context(r['model_id'],r['gene'])} for r in positive]
    zeroids=[r for r in identities if r['label_group']=='zero']
    genes=[]
    for method in ('PRODIGY','DawnRank'):
        for g in sorted({r['gene'] for r in identities}):
            rr=[r for r in zeroids if r['method']==method and r['gene']==g]
            r={'method':method,'gene':g,'common_domain_degree_G0':degrees['G0'][g],
               'common_domain_degree_G1':degrees['G1'][g],'selected_zero_models':len(rr),
               'gain_models':sum(F(x['gain'])>0 for x in rr),'loss_models':sum(F(x['loss'])>0 for x in rr)}
            for k in ('old','new','gain','loss','net'):
                value=sum((F(x['donor_weight'])*F(x[k]) for x in rr),F())
                r[k]=str(value)
                r[k+'_decimal']=float(value)
            for k in ('gain','loss'):
                amount=sum((F(x['donor_weight'])*F(x[k]) for x in rr),F())
                total=sum((F(x['donor_weight'])*F(x[k])*F(x['negative_LFC']) for x in rr),F())
                r[k+'_negative_LFC_sum']=str(total)
                r[k+'_mean_negative_LFC']=float(total/amount) if amount else ''
            genes.append(r)
    references={(r['method'],r['gene']):r for r in rows(inp/'oncogene_gene_flows.csv') if r['label_group']=='zero'}
    for r in genes:
        for k in ('old','new','gain','loss','net'):
            assert F(r[k])==F(references[(r['method'],r['gene'])][k])
    tops={method:[r['gene'] for r in sorted([r for r in genes if r['method']==method and F(r['net'])>0],key=lambda r:(-F(r['net']),r['gene']))[:5]] for method in ('PRODIGY','DawnRank')}
    union=sorted(set().union(*map(set,tops.values())))
    toprows=[r for r in genes if r['gene'] in union]
    nodecontext=[]
    assays={(r['model_id'],r['gene']):r for r in rows(inp/'three_node_source_assay.csv')}
    for m in models:
        for g in s['nodes']:
            mid=m['model_id']
            source=assays[(mid,g)]
            assert source['sample_id']==m['sample_id'] and source['donor_id']==m['donor_id']
            r={**source,**context(mid,g)}
            for method in ('PRODIGY','DawnRank'):
                for state,col in [('G0H0','old'),('G1H1','new')]:
                    r[f'{method}_{col}']=str(vectors[(method,mid,state)].get(g,F()))
            nodecontext.append(r)
    assert len(nodecontext)==255
    nodestats=[]
    for g in s['nodes']:
        rr=[r for r in nodecontext if r['gene']==g]
        r={'gene':g,'models':85,'query_eligible_models':sum(x['query_eligible'] for x in rr),
           'measured_models':sum(x['supplied_LFC']!='' for x in rr),
           'source_positive_models':sum(x['binary_source_value']=='1' for x in rr),
           'positive_but_query_ineligible':sum(x['binary_source_value']=='1' and not x['query_eligible'] for x in rr)}
        for method in ('PRODIGY','DawnRank'):
            r[method+'_gained_models']=sum(F(x[method+'_new'])>F(x[method+'_old']) for x in rr)
            r[method+'_lost_models']=sum(F(x[method+'_new'])<F(x[method+'_old']) for x in rr)
        nodestats.append(r)
    products={'oncogene_all_identity_context.csv':identities,
        'oncogene_zero_identity_context.csv':zeroids,'oncogene_zero_gene_contributors.csv':genes,
        'oncogene_zero_main_table.csv':toprows,'changed_positive_34_context.csv':positivecontext,
        'qualifying_variant_context.csv':variantrows,'three_node_full_cohort_context.csv':nodecontext,
        'three_node_summary.csv':nodestats}
    out.mkdir(parents=True,exist_ok=True)
    results=[]
    for name,data in products.items():
        write(out/name,data)
        digest=sha(out/name)
        assert digest==manifest['expected_output_sha256'][name],name
        results.append({'file':name,'rows':len(data),'sha256':digest,'exact_match':True})
    report={'status':'PASS','input_files_verified':len(manifest['inputs']),
        'exact_outputs':results,'main_table_union':union,'nodes':nodestats,
        'functional_direction_assigned':False,
        'exposure':'Post-exposure, source-anchored descriptive diagnostic; not independent validation.'}
    (out/'REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
