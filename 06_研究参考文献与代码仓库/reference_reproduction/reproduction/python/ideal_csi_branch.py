"""One-setup ideal-CSI *diagnostic* on the original exported true H.

Ideal CSI here changes the MR/RZF beam input from Hhat to H and recomputes
Monte Carlo a/B on that beam input. Geometry, H, optimization, and prelog stay
the same. This is separate from the original pilot-MMSE reproduction.
"""
from pathlib import Path
import sys
import numpy as np
from scipy.io import loadmat, savemat

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import WMMSE_ADMM_Func
import WMMSE_PF_ADMM
import SpectralEfficiencyDownlink


def statistics(H, p=100):
    M, N, K, L = H.shape
    result = {}
    for name in ('MR', 'RZF'):
        a = np.zeros((L, K), complex)
        mean = np.zeros((K, K, L), complex)
        second = np.zeros((K, K, L))
        for n in range(N):
            for l in range(L):
                V = H[:, n, :, l]
                if name == 'MR':
                    W = V
                else:
                    W = np.linalg.solve(p * V @ V.conj().T + np.eye(M), V)
                W = W / np.linalg.norm(W, axis=0)
                effective = H[:, n, :, l].conj().T @ W
                a[l, :] += np.diag(effective) / N
                mean[:, :, l] += effective / N
                second[:, :, l] += np.abs(effective)**2 / N
        B = np.zeros((L, L, K, K), complex)
        for k in range(K):
            for i in range(K):
                z = mean[k, i, :]
                B[:, :, k, i] = np.outer(z, z.conj())
        for l in range(L):
            B[l, l, :, :] = second[:, :, l]
        result['a_' + name] = np.abs(a)
        result['B_' + name] = B.real
    return result


def main():
    source = ROOT / 'reproduction/inputs/original_smoke_nonideal_100.mat'
    original = loadmat(source)
    H = original['H']
    assert H.shape == (4, 100, 20, 16)
    data = statistics(H)
    data['Pmax'] = original['Pmax']
    data['prelogFactor'] = original['prelogFactor']
    rng = np.random.get_state()
    np.random.seed(4873256)
    for mode in ('sumSE', 'PF'):
        for beam in ('MR', 'RZF'):
            label = mode + '_' + beam
            draws = []
            old_rand = np.random.rand
            def recorded_rand(*shape):
                value = old_rand(*shape)
                draws.append(value.copy())
                return value
            np.random.rand = recorded_rand
            try:
                algorithm = (WMMSE_ADMM_Func if mode == 'sumSE'
                             else WMMSE_PF_ADMM)
                mu = algorithm.WMMSE_ADMM_timing(16, 20, 1000,
                                                 data['a_' + beam],
                                                 data['B_' + beam])
            finally:
                np.random.rand = old_rand
            se = SpectralEfficiencyDownlink.Calculate_SINR_and_SE_DL(
                data['a_' + beam], data['B_' + beam], 0.95, mu, 1000)
            data['mu_' + label] = mu
            data['SE_' + label] = se
            data['draws_' + label] = np.stack(draws, axis=2)
            print(label, 'outer', len(draws), 'sum_SE', float(np.sum(se)),
                  'max_power', float(np.max(np.sum(mu**2, axis=0))), flush=True)
    np.random.set_state(rng)
    target = ROOT / 'reproduction/inputs/original_smoke_ideal_100.mat'
    savemat(target, data)
    print('IDEAL_CSI_DIAGNOSTIC_SAVED', target, flush=True)


if __name__ == '__main__':
    main()
