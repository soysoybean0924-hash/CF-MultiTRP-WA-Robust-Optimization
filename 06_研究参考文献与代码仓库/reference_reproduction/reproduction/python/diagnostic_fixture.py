"""Small synthetic a/B fixture; algorithm diagnostic, not paper reproduction."""
from pathlib import Path
import numpy as np
from scipy.io import savemat
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from WMMSE_ADMM_Func import WMMSE_ADMM_timing


def objective(mu, a, b):
    numerator = np.sum(mu * a, axis=0) ** 2
    denominator = np.ones(a.shape[1])
    for k in range(a.shape[1]):
        for i in range(a.shape[1]):
            denominator[k] += mu[:, i] @ b[:, :, k, i] @ mu[:, i]
    return float(np.log2(1 + numerator / (denominator - numerator)).sum())


def main():
    rng = np.random.default_rng(27)
    l, k, pmax = 3, 2, 1.0
    # PSD B is formed from deterministic effective-channel vectors.
    h = rng.normal(size=(l, k, k))
    for user in range(k):
        h[:, user, user] = np.abs(h[:, user, user])
    b = np.empty((l, l, k, k))
    for user in range(k):
        for stream in range(k):
            v = h[:, user, stream]
            b[:, :, user, stream] = np.outer(v, v)
    a = np.abs(h[:, np.arange(k), np.arange(k)])
    mu = 0.1 * np.sqrt(pmax / k) * np.ones((l, k))
    recorded = []
    norm_counts = []
    original_rand = np.random.rand
    original_norm = np.linalg.norm
    original_log2 = np.log2
    objective_trace = []
    def capture_rand(*shape):
        value = rng.random(shape)
        recorded.append(value)
        norm_counts.append(0)
        return value
    def capture_norm(*args, **kwargs):
        if recorded and (kwargs.get('ord') == 'fro' or
                         (len(args) >= 2 and args[1] == 'fro')):
            norm_counts[-1] += 1
        return original_norm(*args, **kwargs)
    def capture_log2(*args, **kwargs):
        value = original_log2(*args, **kwargs)
        objective_trace.append(float(np.sum(value)))
        return value
    np.random.rand = capture_rand
    np.linalg.norm = capture_norm
    np.log2 = capture_log2
    try:
        original_mu = WMMSE_ADMM_timing(l, k, pmax, a, b)
    finally:
        np.random.rand = original_rand
        np.linalg.norm = original_norm
        np.log2 = original_log2
    draws = np.zeros((l, k, max(20, len(recorded))))
    for t, draw in enumerate(recorded):
        draws[:, :, t] = draw
    numerator = np.sum(mu*a,axis=0)**2
    denominator = np.ones(k)
    for user in range(k):
        for stream in range(k):
            denominator[user] += mu[:,stream] @ b[:,:,user,stream] @ mu[:,stream]
    v = np.sqrt(numerator)/denominator
    e = 1-numerator/denominator
    w = 1/e
    first_q = mu*(1+draws[:,:,0])
    first_mu = np.empty_like(mu)
    for stream in range(k):
        A = 0.0005*np.eye(l)
        for user in range(k):
            A += w[user]*v[user]**2*b[:,:,user,stream]
        c = w[stream]*v[stream]*a[:,stream]
        first_mu[:,stream] = np.linalg.solve(A,c+0.0005*first_q[:,stream])
    first_q_projected = first_mu.copy()
    for ap in range(l):
        norm = np.linalg.norm(first_q_projected[ap,:])
        if norm>np.sqrt(pmax):
            first_q_projected[ap,:] *= np.sqrt(pmax)/norm
    out = Path(__file__).resolve().parents[1] / 'inputs' / 'diagnostic_fixture.mat'
    savemat(out, dict(a=a, B=b, Pmax=pmax, draws=draws,
                      initial_mu=mu, initial_sumSE=objective(mu, a, b),
                      python_mu=original_mu, python_final_sumSE=objective(original_mu.T,a,b),
                      python_outer_count=len(recorded), python_norm_counts=norm_counts,
                      python_objective_trace=objective_trace,
                      first_v=v, first_e=e,
                      first_w=w, first_q=first_q, first_mu=first_mu,
                      first_q_projected=first_q_projected))
    print(f'file={out} initial_sumSE={objective(mu, a, b):.15g} '
          f'final_sumSE={objective(original_mu.T,a,b):.15g} outer={len(recorded)} '
          f'max_power={np.max(np.sum(original_mu**2,axis=0)):.15g} '
          f'norm_calls={norm_counts}')


if __name__ == '__main__':
    main()
