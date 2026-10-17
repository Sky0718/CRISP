"""Independent verification of analysis continuous and label-flip outputs.

This file does not import analyse.py. It rebuilds CRC allocation weights from
the accepted original JSONs, uses a different bootstrap aggregation route,
and verifies flip minimality by dynamic programming over all measured gains.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction as Q
from pathlib import Path

import numpy as np

N = 0


def check(ok, message):
    global N
    N += 1
    if not ok:
        raise AssertionError(message)


def rows(p):
    with p.open(encoding='utf-8-sig', newline='') as h:
        return list(csv.DictReader(h))


def jread(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))


def weights(obj):
    return {g:Q(w) for g,w in (obj.items() if isinstance(obj,dict) else obj)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    a=p.parse_args()
    check(not a.report.exists(),'Verification report must be new')
    inp=a.results/'inputs'
    ledger=rows(inp/'union_ledger.csv')
    sources={(r['sample_id'],r['gene']):r for r in rows(inp/'continuous_identities.csv')}
    degree=jread(a.archive/'results/crc_precision/degree_allocations.json')
    models={r['model_id']:r for r in degree['models']}
    counts=Counter(r['donor_id'] for r in models.values())
    check(len(models)==85 and len(counts)==83,'Full85/83 from original allocation JSON')
    vectors={}
    archive_paths=[]
    for policy,rel in [('PRODIGY','results/crc_prodigy_allocations.json'),('DawnRank','results/crc_precision/allocations.json')]:
        path=a.archive/rel;archive_paths.append(path)
        data=jread(path)
        temp=defaultdict(dict)
        for r in data['allocations']:
            if r['budget'] != 10:continue
            cell=r.get('cell',r.get('state')).replace('_','')
            if cell not in ('G0H0','G1H1'):continue
            t=int(cell=='G1H1');seed=str(r['seed'])
            check(r.get('numerically_certified',True),'Certified current DawnRank state')
            temp[r['model_id'],t][seed]=weights(r['weights'])
        for mid in models:
            for t in (0,1):
                sv=temp[mid,t]
                expected={'104729','130363','155921','196613','228017'} if policy=='PRODIGY' else {'deterministic'}
                check(set(sv)==expected,'Original complete seed set')
                genes=set().union(*(set(v) for v in sv.values()))
                vectors[policy,mid,t]={g:sum((v.get(g,Q()) for v in sv.values()),Q())/len(sv) for g in genes}
    archive_paths.append(a.archive/'results/crc_precision/degree_allocations.json')
    for mid,m in models.items():
        for t,cell in enumerate(('G0','G1')):
            vectors['degree',mid,t]=weights(m['full_degree'][cell])
    old_union={(r['policy'],r['model_id'],r['gene']):r for r in ledger if r['context']=='CRC'}
    accepted_models={(r['policy'],r['model_id'],int(r['state'])):r for r in rows(a.results/'continuous_model_levels.csv')}
    # Rebuild union and continuous values from complete original selected sets.
    independently_rebuilt={}
    for (policy,mid,t),vector in vectors.items():
        sample=models[mid]['sample_id']
        values=[Q(),Q(),Q()]
        for gene,w in vector.items():
            if not w:continue
            sr=sources[sample,gene]
            if sr['role']=='UNLISTED':continue
            check((policy,mid,gene) in old_union,'Every selected cancer identity present in union ledger')
            check(Q(old_union[policy,mid,gene]['new' if t else 'old'])==w,'Union ledger matches original seed-averaged allocation')
            values[0]+=w
            if sr['supplied_LFC']:
                raw=json.loads(sr['raw_representations'])
                val=Q(raw[0][5])
                check(val==Q(sr['supplied_LFC']),'Raw source decimal agrees')
                values[1]+=w;values[2]-=w*val
        expected=accepted_models[policy,mid,t]
        check(all(values[i]==Q(expected[k]) for i,k in enumerate(('R','A','S'))),'Exact independently reconstructed model totals')
        independently_rebuilt[policy,mid,t]=values
    for r in ledger:
        if r['context']=='CRC':
            for t,field in enumerate(('old','new')):
                check(Q(r[field])==vectors[r['policy'],r['model_id'],t].get(r['gene'],Q()),'No spurious union-ledger allocation')
    manifest=jread(a.results/'ANALYSIS_MANIFEST.json')
    donors=manifest['donor_ids'];policies=('PRODIGY','DawnRank','degree')
    donor_array=np.zeros((83,3,2,3))
    for di,donor in enumerate(donors):
        mids=[m for m in models if models[m]['donor_id']==donor]
        for pi,policy in enumerate(policies):
            for t in (0,1):
                donor_array[di,pi,t,:]=[float(sum((independently_rebuilt[policy,mid,t][k] for mid in mids),Q())/len(mids)) for k in range(3)]
    indices=np.load(a.results/'donor_draw_indices.npy',allow_pickle=False)
    regenerated=np.random.Generator(np.random.PCG64(470016)).integers(0,83,size=(10000,83))
    check(np.array_equal(indices,regenerated),'Frozen random generator and complete donor draws')
    multiplicity=np.asarray([np.bincount(draw,minlength=83) for draw in indices])
    drawn=(multiplicity@donor_array.reshape(83,-1)/83).reshape(10000,3,2,3)
    observed=donor_array.mean(axis=0)
    def contrasts(array):
        # Build individual family columns independently, not the production stack operation.
        out=[]
        for pi in range(3):
            v=array[...,pi,:,:]
            for metric in ('R','A','coverage','S','mean_depletion'):
                if metric=='R':out.append(v[...,1,0]-v[...,0,0])
                elif metric=='A':out.append(v[...,1,1]-v[...,0,1])
                elif metric=='coverage':out.append(v[...,1,1]/v[...,1,0]-v[...,0,1]/v[...,0,0])
                elif metric=='S':out.append(v[...,1,2]-v[...,0,2])
                else:out.append(v[...,1,2]/v[...,1,1]-v[...,0,2]/v[...,0,1])
        return np.stack(out,axis=-1)
    boot=contrasts(drawn); point=contrasts(observed)
    expected_boot=np.load(a.results/'continuous_resampled_endpoints.npy',allow_pickle=False)
    max_boot_error=float(np.max(np.abs(boot-expected_boot)))
    check(max_boot_error<1e-12,'Independent count-matrix bootstrap agrees for all150000 endpoints')
    sd=np.std(boot,axis=0,ddof=1)
    standard=np.divide(boot-point,sd,out=np.zeros_like(boot),where=sd>0)
    maxima=np.max(np.abs(standard),axis=1)
    q=float(np.quantile(maxima,.95,method='linear'))
    check(abs(q-manifest['multiplier'])<1e-12,'Independent envelope multiplier')
    for i,r in enumerate(rows(a.results/'continuous_endpoints.csv')):
        check(abs(float(r['estimate'])-point[i])<1e-12,'Envelope point')
        check(abs(float(r['lower'])-(point[i]-q*sd[i]))<1e-12,'Envelope lower')
        check(abs(float(r['upper'])-(point[i]+q*sd[i]))<1e-12,'Envelope upper')
    # B: independent exact subset maximisation, including all nonbeneficial measured flips.
    verified=[]
    witness=rows(a.results/'flip_minimal_witness.csv')
    curves=rows(a.results/'flip_curves.csv')
    for summary in rows(a.results/'flip_summary.csv'):
        context,policy=summary['context'],summary['policy']
        rr=[r for r in ledger if (r['context'],r['policy'])==(context,policy)]
        denom=[sum((Q(r['population_weight'])*Q(r[f]) for r in rr),Q()) for f in ('old','new')]
        identities={}
        for r in rr:
            key=r['sample_id'],r['gene']
            check(key not in identities,'One model-gene row per policy in original ledger')
            identities[key]=r
        def endpoint(flipped):
            sums=[Q(),Q()]
            for key,r in identities.items():
                c=Q(r['population_weight'])*(Q(r['new'])/denom[1]-Q(r['old'])/denom[0])
                y=(1-int(r['label_value']) if key in flipped else int(r['label_value'])) if r['label_state']=='MEASURED' else int(c>0)
                for t,f in enumerate(('old','new')):sums[t]+=Q(r['population_weight'])*Q(r[f])*y
            return sums[1]/denom[1]-sums[0]/denom[0]
        base=endpoint(set());check(base==Q(summary['upper_before']),'Direct two-state baseline bound')
        k=int(summary['minimal_flips'])
        gains=[]
        for key,r in identities.items():
            if r['label_state']=='MEASURED':
                c=Q(r['population_weight'])*(Q(r['new'])/denom[1]-Q(r['old'])/denom[0])
                gains.append(c*(1-2*int(r['label_value'])))
        # Exact DP stores best gain for exactly j flips; no greedy ordering is imported.
        dp=[Q()]+[None]*k
        for gain in gains:
            for j in range(k,0,-1):
                if dp[j-1] is not None:
                    proposed=dp[j-1]+gain
                    if dp[j] is None or proposed>dp[j]:dp[j]=proposed
        if k:
            check(base+max(x for x in dp[:k] if x is not None)<0,'DP proves no at-most-k-1 witness exists')
            check(base+dp[k]>=0,'DP proves k sufficient')
        else:check(base>=0,'Zero-flip case already nonnegative')
        ww=[r for r in witness if (r['context'],r['policy'])==(context,policy)]
        ids={(r['sample_id'],r['gene']) for r in ww}
        check(len(ww)==len(ids)==k,'Exactly k distinct measured identities')
        after=endpoint(ids)
        check(after==Q(summary['upper_at_k']),'Direct two-state witness attainment')
        if k:check(after==base+dp[k],'Deterministic witness is globally optimal at k')
        check(len(gains)==int(summary['measured_identities']),'Measured identity denominator')
        verified.append(dict(context=context,policy=policy,k=k,measured=len(gains),upper_before=str(base),upper_at_k=str(after),dp_verified=True))
    # C: compare every role contribution with the previously accepted exact table.
    prior={(r['context'],r['policy'],r['unit']):r for r in rows(inp/'accepted_role_decomposition.csv') if r['partition']=='role'}
    for r in rows(a.results/'four_role_decomposition.csv'):
        old=prior[r['context'],r['policy'],r['role']]
        for field in ('within','reweight','entry','exit'):
            check(Q(r[field])==Q(old[field]),'Four-role contribution equals independently accepted prior decomposition')
    report=dict(status='PASS',checks=N,independent_crc_allocations=len(vectors),
                continuous_bootstrap_endpoint_count=boot.size,maximum_absolute_bootstrap_discrepancy=max_boot_error,
                envelope_multiplier=q,flips=verified,new_native_fits=0,
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                original_json_inputs=[dict(path=str(p.relative_to(a.archive)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in archive_paths])
    a.report.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
