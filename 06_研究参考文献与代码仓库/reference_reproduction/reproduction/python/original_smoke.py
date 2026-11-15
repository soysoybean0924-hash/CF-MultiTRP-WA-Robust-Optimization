"""Run original source in memory with reduced setup/realization counts.

This is a smoke test, never the formal 2000 x 100 paper reproduction.
Algorithm source is not edited on disk. DNN section is not executed.
"""
from pathlib import Path
import os
import sys
import numpy as np
from scipy.io import savemat

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import APLocation_Generation
import WMMSE_ADMM_Func
import WMMSE_PF_ADMM


def main():
    original_int = np.int if hasattr(np, 'int') else None
    original_load = np.load
    np.int = int
    np.random.seed(4873256)
    positions = APLocation_Generation.RandomAPLocations(16, 1000)

    def supplied_load(path, *args, **kwargs):
        if str(path).replace('\\', '/').endswith('/new_storage/APpositions.npy'):
            return positions
        return original_load(path, *args, **kwargs)

    source = (ROOT / 'main_DNN.py').read_text(encoding='utf-8')
    source = source.replace('nbrOfSetups = 2000', 'nbrOfSetups = 1', 1)
    realizations = int(os.environ.get('SMOKE_REALIZATIONS', '100'))
    assert 1 <= realizations <= 100
    source = source.replace('nbrOfRealizations = 100', f'nbrOfRealizations = {realizations}', 1)
    marker = 'dataset_B_MR = np.real(dataset_B_MR)'
    assert marker in source
    source = source.split(marker, 1)[0]
    namespace = {'__name__': '__main__', '__file__': str(ROOT / 'main_DNN.py')}
    np.load = supplied_load
    original_rand = np.random.rand
    originals = [WMMSE_ADMM_Func.WMMSE_ADMM_timing,
                 WMMSE_PF_ADMM.WMMSE_ADMM_timing]
    draws = {'sumSE_MR': [], 'sumSE_RZF': [], 'PF_MR': [], 'PF_RZF': []}
    current = [None]
    calls = [0, 0]
    def recorded_rand(*shape):
        value = original_rand(*shape)
        if current[0] is not None:
            draws[current[0]].append(value.copy())
        return value
    def wrap(kind):
        def call(*args, **kwargs):
            names = (('sumSE_MR','sumSE_RZF'),('PF_MR','PF_RZF'))[kind]
            name = names[calls[kind]]
            calls[kind] += 1
            current[0] = name
            try:
                return originals[kind](*args, **kwargs)
            finally:
                current[0] = None
        return call
    np.random.rand = recorded_rand
    WMMSE_ADMM_Func.WMMSE_ADMM_timing = wrap(0)
    WMMSE_PF_ADMM.WMMSE_ADMM_timing = wrap(1)
    try:
        exec(compile(source, str(ROOT / 'main_DNN.py'), 'exec'), namespace)
    finally:
        np.load = original_load
        np.random.rand = original_rand
        WMMSE_ADMM_Func.WMMSE_ADMM_timing = originals[0]
        WMMSE_PF_ADMM.WMMSE_ADMM_timing = originals[1]
        if original_int is None:
            del np.int
        else:
            np.int = original_int
    out = ROOT / 'reproduction' / 'inputs' / f'original_smoke_nonideal_{realizations}.mat'
    keys = ['APpositions', 'H', 'Hhat', 'Np', 'CorrR', 'betas',
            'a_MR', 'B_MR', 'a_RZF', 'B_RZF',
            'mu_MR_WMMSE_ADMM', 'mu_RZF_WMMSE_ADMM', 'SE_MR_WMMSE_ADMM',
            'SE_RZF_WMMSE_ADMM', 'mu_MR_WMMSE_PF_ADMM',
            'mu_RZF_WMMSE_PF_ADMM', 'pilotIndex', 'prelogFactor', 'Pmax']
    data = {key: namespace[key] for key in keys}
    for name, values in draws.items():
        data['draws_' + name] = np.stack(values, axis=2)
        data['outer_' + name] = len(values)
    savemat(out, data)
    print('saved', out)
    print('sumSE_MR', float(namespace['SE_MR_WMMSE_ADMM'].sum()))
    print('sumSE_RZF', float(namespace['SE_RZF_WMMSE_ADMM'].sum()))
    print('outer_counts', {key: len(value) for key,value in draws.items()})


if __name__ == '__main__':
    main()
