"""Reproduce Supplementary Data 11 using only this folder and Python >=3.10."""
import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True,help='New output folder outside the source-data directory')
    ap.add_argument('--data',type=Path,required=True,help='External Union_Sensitivity_Analysis data directory')
    a=ap.parse_args();base=a.data.resolve();code=Path(__file__).resolve().parent
    out=a.output.resolve()
    if out.exists() or out.is_relative_to(base) or out.is_relative_to(code):
        raise SystemExit('Choose a new output directory outside the input and code directories.')
    with (base/'INPUT_OUTPUT_MAP.csv').open(encoding='utf-8',newline='') as h: provenance=list(csv.DictReader(h))
    for r in provenance:
        path=base/r['package_path']
        if sha(path)!=r['sha256']: raise SystemExit('Checksum mismatch: '+r['package_path'])
    out.mkdir(parents=True)
    def run(script,*args):
        subprocess.run([sys.executable,'-B',str(code/script),*map(str,args)],check=True)
    run('analyse_fraction.py','--union',base/'inputs/union','--output',out/'results')
    run('analyse_budget_controls.py','--archive',base/'inputs/archive','--controls',base/'inputs/controls',
        '--union',base/'inputs/union','--output',out/'sensitivity_results')
    run('verify_fraction.py','--union',base/'inputs/union','--category',base/'inputs/category',
        '--analysis',out,'--controls',base/'inputs/controls','--archive',base/'inputs/archive',
        '--output',out/'INDEPENDENT_VALIDATION.json')
    matches=[]
    for folder in ('results','sensitivity_results'):
        for expected in sorted((base/'expected'/folder).glob('*.csv')):
            actual=out/folder/expected.name
            if sha(actual)!=sha(expected): raise SystemExit('Reproduction mismatch: '+folder+'/'+expected.name)
            matches.append(dict(file=folder+'/'+expected.name,sha256=sha(actual),byte_identical=True))
    expected=json.loads((base/'expected/INDEPENDENT_VALIDATION.json').read_text(encoding='utf-8'))
    actual=json.loads((out/'INDEPENDENT_VALIDATION.json').read_text(encoding='utf-8'))
    if expected!=actual: raise SystemExit('Independent validation mismatch')
    # The generated manifests describe the user's chosen location; remove absolute path prefixes.
    def clean(value):
        if isinstance(value,dict): return {k:clean(v) for k,v in value.items()}
        if isinstance(value,list): return [clean(v) for v in value]
        if isinstance(value,str):
            try:
                p=Path(value)
                if p.is_absolute():
                    if p.is_relative_to(base): return 'package/'+p.relative_to(base).as_posix()
                    if p.is_relative_to(out): return 'output/'+p.relative_to(out).as_posix()
            except (ValueError,OSError): pass
        return value
    for p in out.glob('*/manifest.json'):
        p.write_text(json.dumps(clean(json.loads(p.read_text())),indent=2)+'\n',encoding='utf-8')
    report=dict(status='PASS',python=sys.version.split()[0],input_files_verified=len(provenance),
                byte_identical_csv_outputs=matches,independent_validation_identical=True,
                dependencies='Python standard library only; no native fits or external downloads')
    (out/'REPRODUCTION_CHECK.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
