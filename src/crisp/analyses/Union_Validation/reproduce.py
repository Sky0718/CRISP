"""Reproduce all finite outputs from the small, exact supplied input bundle."""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--archive',type=Path)
    p.add_argument('--data',type=Path,required=True,help='External Union_Validation data directory')
    a=p.parse_args()
    here=Path(__file__).resolve().parent
    expected=a.data.resolve()/'results'
    manifest=json.loads((expected/'ANALYSIS_MANIFEST.json').read_text(encoding='utf-8'))
    work=a.work.resolve()
    if work.exists() or work.is_relative_to(a.data.resolve()) or work.is_relative_to(here):
        raise ValueError('Choose a new work directory outside the input and code directories')
    work.mkdir(parents=True)
    source=work/'source'
    for item in manifest['source_inputs']:
        original=expected/'inputs'/item['file']
        if hashlib.sha256(original.read_bytes()).hexdigest()!=item['sha256']:
            raise ValueError('Input hash mismatch: '+item['file'])
        target=source/item['source_relative']
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(original,target)
    out=work/'results'
    subprocess.run([sys.executable,'-B',str(here/'analyse.py'),'--source',str(source),'--output',str(out)],check=True)
    csvs=[]
    for file in sorted(expected.glob('*.csv')):
        if file.read_bytes() != (out/file.name).read_bytes():
            raise AssertionError('CSV reproduction mismatch: '+file.name)
        csvs.append(file.name)
    arrays=[]
    for file in sorted(expected.glob('*.npy')):
        lhs=np.load(file,allow_pickle=False);rhs=np.load(out/file.name,allow_pickle=False)
        if not np.array_equal(lhs,rhs):raise AssertionError('NumPy reproduction mismatch: '+file.name)
        arrays.append(file.name)
    independent='not requested; original allocation archive was not supplied'
    if a.archive:
        subprocess.run([sys.executable,'-B',str(here/'verify.py'),'--results',str(out),
                        '--archive',str(a.archive.resolve()),'--report',str(out/'INDEPENDENT_VALIDATION.json')],check=True)
        independent='PASS against supplied original allocation JSONs'
    report=dict(status='PASS',csv_files_exact=csvs,arrays_exact=arrays,independent_validation=independent,
                numpy=np.__version__,new_native_fits=0)
    (work/'REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
