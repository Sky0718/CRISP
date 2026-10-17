import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path


def load(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def main(folder, output):
    freeze = load(folder/'ALLOCATION_FREEZE.json')
    contrast_plan = load(folder/'CONTRAST_FREEZE.json')['contrasts']
    labels = {(r['sample_id'],r['gene']):(r['state'],r['value']) for r in load(folder/'ALL_SELECTED_LABELS.json')['calls']}
    grouped = defaultdict(list)
    with gzip.open(folder/'full_query_rankings.csv.gz','rt',encoding='utf-8',newline='') as stream:
        for row in csv.DictReader(stream):
            grouped[(row['state'],row['sample_id'])].append(row)
    allocations = {}
    tied_cutoffs = 0
    for key, ranks in grouped.items():
        assert [int(row['rank_order']) for row in ranks] == list(range(1,len(ranks)+1))
        assert {row['gene'] for row in ranks} == set(freeze['queries'][key[1]])
        scores = [(float.fromhex(row['binary64_score_hex']),row['gene']) for row in ranks]
        assert scores == sorted(scores,key=lambda item:(-item[0],item[1]))
        assert all(float(row['score']) == score for row,(score,_) in zip(ranks,scores))
        buckets = defaultdict(list)
        for score,gene in scores:
            buckets[score].append(gene)
        remaining = 10
        weights = {}
        for score in sorted(buckets,reverse=True):
            members = buckets[score]
            if not remaining:
                break
            taken = min(remaining,len(members))
            if taken < len(members):
                tied_cutoffs += 1
            weights.update({gene:Fraction(taken,len(members)) for gene in members})
            remaining -= taken
        assert sum(weights.values(),Fraction()) == min(10,len(scores))
        assert weights == {gene:Fraction(value) for gene,value in freeze['states'][key[0]][key[1]].items()}
        assert all(Fraction(row['allocation']) == weights.get(row['gene'],Fraction()) for row in ranks)
        allocations[key] = weights
    assert len(allocations) == 3060
    models = freeze['models']
    donor_counts = Counter(row['donor_id'] for row in models)
    assert len(models)==85 and len(donor_counts)==83
    assert all(Fraction(row['alpha']) == Fraction(1,83*donor_counts[row['donor_id']]) for row in models)
    units = {('population','all'):[(row['sample_id'],Fraction(1,83*donor_counts[row['donor_id']])) for row in models]}
    units.update({('donor',donor):[(row['sample_id'],Fraction(1,donor_counts[donor])) for row in models if row['donor_id']==donor] for donor in donor_counts})
    units.update({('model',row['model_id']):[(row['sample_id'],Fraction(1))] for row in models})
    for plan in contrast_plan.values():
        for sample,stored in plan['model_coefficients'].items():
            old = allocations[(plan['old_state'],sample)]
            new = allocations[(plan['new_state'],sample)]
            reconstructed = {gene:new.get(gene,Fraction())-old.get(gene,Fraction()) for gene in set(old)|set(new) if new.get(gene,Fraction())!=old.get(gene,Fraction())}
            assert reconstructed == {gene:Fraction(value) for gene,value in stored.items()}
    def selected(state,unit,identity,category):
        values = {}
        for sample,weight in units[(unit,identity)]:
            for gene,amount in allocations[(state,sample)].items():
                if category in freeze['annotations'][gene]:
                    values[(sample,gene)] = weight*amount
        return values
    def totals(values):
        result = dict.fromkeys(('R','H','Z','ABSENT','EXPLICIT_NULL'),Fraction())
        for key,weight in values.items():
            state,positive = labels[key]
            result['R'] += weight
            if state=='MEASURED':
                result['H' if positive else 'Z'] += weight
            else:
                result[state] += weight
        result['U'] = result['ABSENT']+result['EXPLICIT_NULL']
        result['positive_fraction'] = result['H']/result['R'] if result['R'] and not result['U'] else None
        result['fraction_lower'] = result['H']/result['R'] if result['R'] else None
        result['fraction_upper'] = (result['H']+result['U'])/result['R'] if result['R'] else None
        return result
    def interval(values):
        measured = Fraction()
        low_unknown,high_unknown = Fraction(),Fraction()
        for key,coefficient in values.items():
            state,positive = labels[key]
            if state=='MEASURED':
                measured += coefficient*positive
            elif coefficient<0:
                low_unknown += coefficient
            else:
                high_unknown += coefficient
        lo,hi=measured+low_unknown,measured+high_unknown
        direction='NEGATIVE' if hi<0 else 'POSITIVE' if lo>0 else 'ZERO' if lo==hi==0 else 'UNRESOLVED'
        return dict(known=measured,unknown_negative=low_unknown,unknown_positive=high_unknown,lower=lo,upper=hi,width=hi-lo,sign=direction,nonzero_shared_identities=sum(value!=0 for value in values.values()))
    def compare(record,expected):
        for key,value in expected.items():
            if value is None:
                assert record[key]=='',(key,record[key],value)
            elif isinstance(value,str):
                assert record[key]==value,(key,record[key],value)
            else:
                assert Fraction(record[key])==value,(key,record[key],value)
    level_count=0
    with (folder/'state_levels.csv').open(encoding='utf-8',newline='') as stream:
        for row in csv.DictReader(stream):
            compare(row,totals(selected(row['state'],row['unit'],row['id'],row['category'])))
            level_count+=1
    contrast_count=0
    with (folder/'contrast_endpoints.csv').open(encoding='utf-8',newline='') as stream:
        for row in csv.DictReader(stream):
            plan=contrast_plan[row['contrast']]
            old=selected(plan['old_state'],row['unit'],row['id'],row['category'])
            new=selected(plan['new_state'],row['unit'],row['id'],row['category'])
            a,b=totals(old),totals(new)
            coeff={key:new.get(key,Fraction())-old.get(key,Fraction()) for key in set(old)|set(new)}
            found=interval(coeff)
            found.update(old_R=a['R'],new_R=b['R'],delta_R=b['R']-a['R'],old_H=a['H'],new_H=b['H'],old_U=a['U'],new_U=b['U'],old_fraction=a['positive_fraction'],new_fraction=b['positive_fraction'])
            if a['R'] and b['R']:
                normalised={key:new.get(key,Fraction())/b['R']-old.get(key,Fraction())/a['R'] for key in set(old)|set(new)}
                ratio=interval(normalised)
                found.update(delta_fraction_lower=ratio['lower'],delta_fraction_upper=ratio['upper'],fraction_sign=ratio['sign'])
            else:
                found.update(delta_fraction_lower=None,delta_fraction_upper=None,fraction_sign='UNDEFINED_EMPTY_CATEGORY')
            compare(row,found)
            contrast_count+=1
    assert (level_count,contrast_count)==(60840,45630)
    report=dict(status='PASS',utc=datetime.now(timezone.utc).isoformat(),complete_allocation_reconstructions=len(allocations),fractional_cutoff_states=tied_cutoffs,model_contrast_coefficient_maps=27*85,level_rows=level_count,contrast_rows=contrast_count,arithmetic='Exact fractions from complete query scores; independent grouping implementation, no analysis-module imports',source_hashes={name:hashlib.sha256((folder/name).read_bytes()).hexdigest() for name in ('full_query_rankings.csv.gz','ALLOCATION_FREEZE.json','CONTRAST_FREEZE.json','ALL_SELECTED_LABELS.json','state_levels.csv','contrast_endpoints.csv')})
    if output:
        assert not output.exists(),str(output)
        output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--directory',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    main(args.directory,args.output)
