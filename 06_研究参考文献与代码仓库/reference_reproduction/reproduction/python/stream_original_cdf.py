"""Bounded-memory original Python WMMSE CDF runner.

Executes original main_DNN.py's pre-DNN section one setup at a time in
memory, with 100 realizations, fixed AP locations, original algorithms,
and resumable RNG state. DNN curves need unavailable trained models and
are deliberately excluded. Default target is the paper's 2000 setups.
"""
import argparse
import os
from pathlib import Path
import pickle
import sys
import time

import numpy as np
from scipy.io import savemat

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import APLocation_Generation

LABELS = ('SE_MR_equal','SE_RZF_equal','SE_MR_Giovanni19',
          'SE_RZF_Giovanni19','SE_MR_WMMSE_ADMM','SE_RZF_WMMSE_ADMM',
          'SE_MR_WMMSE_PF_ADMM','SE_RZF_WMMSE_PF_ADMM')
CHECKPOINT = ROOT/'reproduction'/'python'/'stream_checkpoint.pkl'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--setups',type=int,default=2000)
    parser.add_argument('--batch-size',type=int,default=25)
    args=parser.parse_args()
    if not 1<=args.setups<=2000:
        parser.error('setups must be 1..2000')
    if not 1<=args.batch_size<=2000:
        parser.error('batch-size must be 1..2000')
    old_int=getattr(np,'int',None)
    old_load=np.load
    np.int=int
    np.random.seed(4873256)
    ap=APLocation_Generation.RandomAPLocations(16,1000)
    np.random.seed(4873256)
    source=(ROOT/'main_DNN.py').read_text(encoding='utf-8')
    source=source.split('dataset_B_MR = np.real(dataset_B_MR)',1)[0]
    source=source.replace('np.random.seed(4873256)#int(time.time()))',
                          '# seeded by streaming wrapper',1)
    source=source.replace('nbrOfSetups = 2000','nbrOfSetups = 1',1)
    compiled=compile(source,str(ROOT/'main_DNN.py'),'exec')
    def supplied_load(path,*a,**kw):
        if str(path).replace('\\','/').endswith('/new_storage/APpositions.npy'):
            return ap.copy()
        return old_load(path,*a,**kw)
    np.load=supplied_load
    samples={name:np.full((20,2000),np.nan) for name in LABELS}
    start=0
    if CHECKPOINT.exists():
        with CHECKPOINT.open('rb') as f: state=pickle.load(f)
        start=state['completed']; samples=state['samples']; ap=state['ap']
        np.random.set_state(state['rng'])
        print('resumed',start,flush=True)
    started=time.perf_counter()
    try:
        end=min(args.setups,start+args.batch_size)
        for setup in range(start,end):
            scope={'__name__':'__main__','__file__':str(ROOT/'main_DNN.py')}
            exec(compiled,scope)
            for name in LABELS:
                samples[name][:,setup]=scope[name][:,0]
            if (setup+1)%25==0 or setup+1==end:
                state=dict(completed=setup+1,samples=samples,ap=ap,
                           rng=np.random.get_state())
                temporary=CHECKPOINT.with_suffix('.tmp')
                with temporary.open('wb') as f: pickle.dump(state,f)
                os.replace(temporary,CHECKPOINT)
                print('completed',setup+1,'elapsed_s',round(time.perf_counter()-started,2),flush=True)
    finally:
        np.load=old_load
        if old_int is None: del np.int
        else: np.int=old_int
    if args.setups==2000 and end==2000:
        out=ROOT/'reproduction'/'python'/'original_cdf_2000.mat'
        savemat(out,{name:values for name,values in samples.items()})
        print('FORMAL_WMMSE_CDF_INPUT',out,flush=True)
    else:
        print('BATCH_COMPLETE',end,'TARGET',args.setups,flush=True)


if __name__=='__main__': main()
