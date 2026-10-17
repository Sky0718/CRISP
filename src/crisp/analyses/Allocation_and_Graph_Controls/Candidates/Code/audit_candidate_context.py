"""Independent algebraic and identity checks of the finite context products."""
import argparse
import csv
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source-data',type=Path,required=True)
    ap.add_argument('--report',type=Path,required=True)
    a=ap.parse_args()
    p=a.source_data
    rr=read(p/'oncogene_all_identity_context.csv')
    assert len(rr)==292
    assert len({(r['model_id'],r['gene']) for r in rr})==177
    assert len({(r['method'],r['model_id'],r['gene']) for r in rr})==len(rr)
    vr=read(p/'qualifying_variant_context.csv')
    assert len(vr)==184
    source_counts=Counter((v['model_id'],v['gene']) for v in vr)
    assert all(v['functional_direction']=='not_assigned' and v['eligible']=='1' for v in vr)
    expected={'PRODIGY':F(118,415),'DawnRank':F(55,166)}
    checks={}
    for m in expected:
        sub=[r for r in rr if r['method']==m and r['label_group']=='zero']
        net=sum((F(r['donor_weight'])*(F(r['new'])-F(r['old'])) for r in sub),F())
        gross=sum((F(r['donor_weight'])*(F(r['gain'])-F(r['loss'])) for r in sub),F())
        assert net==gross==expected[m]
        checks[m]={'net_zero_change':str(net),'zero_identity_records':len(sub)}
    for r in rr:
        assert F(r['new'])-F(r['old'])==F(r['net'])==F(r['gain'])-F(r['loss'])
        assert r['query_eligible']=='1' and r['variant_level_function']=='not_assigned'
        assert int(r['qualifying_source_records'])==source_counts[(r['model_id'],r['gene'])]
        assert F(r['negative_LFC'])==-F(r['supplied_LFC'])
    positive=read(p/'changed_positive_34_context.csv')
    assert len(positive)==len({(r['model_id'],r['gene']) for r in positive})==34
    assert all(r['binary_call']=='1' and any(F(r[m+'_old_weight'])!=F(r[m+'_new_weight']) for m in expected) for r in positive)
    nodes=read(p/'three_node_full_cohort_context.csv')
    state=json.loads((p/'inputs/source_context.json').read_text())
    queries={r['model_id']:set(r['query']) for r in state['models']}
    node_checks={}
    for g in ['EGFR','KRAS','PTPN11']:
        x=[r for r in nodes if r['gene']==g]
        assert len(x)==len({r['model_id'] for r in x})==85
        assert all(int(r['query_eligible'])==int(g in queries[r['model_id']]) for r in x)
        assert all(r['supplied_LFC']!='' and r['binary_source_value'] in ['0','1'] for r in x)
        c=Counter((int(r['query_eligible']),int(r['binary_source_value'])) for r in x)
        node_checks[g]={f'eligible{q}_positive{b}':c[(q,b)] for q in [0,1] for b in [0,1]}
        for r in x:
            if r['query_eligible']=='0':
                assert all(F(r[m+'_'+v])==0 for m in expected for v in ['old','new'])
    assert node_checks['EGFR']=={'eligible0_positive0':15,'eligible0_positive1':69,'eligible1_positive0':0,'eligible1_positive1':1}
    assert node_checks['PTPN11']=={'eligible0_positive0':24,'eligible0_positive1':60,'eligible1_positive0':0,'eligible1_positive1':1}
    assert node_checks['KRAS']=={'eligible0_positive0':38,'eligible0_positive1':13,'eligible1_positive0':6,'eligible1_positive1':28}
    report={'status':'PASS','allocation_checks':checks,'exact_query_source_call_cross_tabs':node_checks,
        'all_changed_positive_identities':34,'all_qualifying_variant_records':184,
        'scope':'Independent algebraic/counting implementation on immutable packaged inputs; not independent experimental validation.'}
    a.report.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
