"""Explicitly corrected feasible-output diagnostic, separate from original code.

Differences: test both ADMM residuals, return projected consensus q, and
record caps. This is NOT used to claim exact reproduction of the paper.
"""
from pathlib import Path
import numpy as np
from scipy.io import loadmat, savemat


def objective(mu, a, b):
    l, k = a.shape
    numerator = np.sum(mu * a, axis=0) ** 2
    denominator = np.ones(k)
    for user in range(k):
        for stream in range(k):
            denominator[user] += mu[:, stream] @ b[:, :, user, stream] @ mu[:, stream]
    return float(np.sum(np.log2(1 + numerator / (denominator - numerator))))


def solve(a, b, pmax, draws, max_outer=20, max_inner=100000,
          primal_tol=1e-5, dual_tol=1e-5):
    l, k = a.shape
    mu = 0.1 * np.sqrt(pmax / k) * np.ones((l, k))
    rho = 0.001
    trajectory = [objective(mu, a, b)]
    histories = []
    reason = 'outer_limit'
    previous = 0.0
    for outer in range(max_outer):
        if (trajectory[-1] - previous) ** 2 <= 0.01:
            reason = 'original_outer_stop_with_feasible_inner'
            break
        previous = trajectory[-1]
        numerator = np.sum(mu * a, axis=0) ** 2
        denominator = np.ones(k)
        for user in range(k):
            for stream in range(k):
                denominator[user] += mu[:, stream] @ b[:, :, user, stream] @ mu[:, stream]
        v = np.sqrt(numerator) / denominator
        e = 1 - numerator / denominator
        w = 1 / e
        ainverse = np.empty((k, l, l))
        c = np.empty((l, k))
        for stream in range(k):
            A = rho / 2 * np.eye(l)
            for user in range(k):
                A += w[user] * v[user] ** 2 * b[:, :, user, stream]
            ainverse[stream] = np.linalg.inv(A)
            c[:, stream] = w[stream] * v[stream] * a[:, stream]
        q = mu * (1 + draws[:, :, outer])
        g = np.zeros_like(mu)
        history = []
        for inner in range(1, max_inner + 1):
            old_q = q.copy()
            for stream in range(k):
                mu[:, stream] = ainverse[stream] @ (
                    c[:, stream] + rho / 2 * (q[:, stream] + g[:, stream]))
            q = mu - g
            norm = np.linalg.norm(q, axis=1)
            q *= np.minimum(1, np.sqrt(pmax) / np.maximum(norm, np.finfo(float).eps))[:, None]
            g += q - mu
            primal = np.linalg.norm(mu - q)
            dual = rho * np.linalg.norm(q - old_q)
            history.append((primal, dual))
            if inner >= 6 and primal <= primal_tol and dual <= dual_tol:
                break
        histories.append(np.asarray(history))
        mu = q.copy()  # consensus projection, not post-hoc normalization
        trajectory.append(objective(mu, a, b))
        if inner == max_inner:
            reason = 'inner_limit'
            break
    return dict(mu=mu.T, objective=np.asarray(trajectory),
                histories=histories, reason=reason,
                max_power=float(np.max(np.sum(mu ** 2, axis=1))))


def main():
    root = Path(__file__).resolve().parents[1]
    s = loadmat(root / 'inputs' / 'diagnostic_fixture.mat')
    result = solve(s['a'], s['B'], float(s['Pmax'][0, 0]), s['draws'])
    output = root / 'python' / 'corrected_sumse_diagnostic.mat'
    savemat(output, dict(mu=result['mu'], objective=result['objective'],
                         max_power=result['max_power'],
                         inner_counts=[len(x) for x in result['histories']],
                         final_primal=[x[-1, 0] for x in result['histories']],
                         final_dual=[x[-1, 1] for x in result['histories']]))
    print(f"corrected_final={result['objective'][-1]:.15g} "
          f"max_power={result['max_power']:.15g} "
          f"inner_counts={[len(x) for x in result['histories']]} "
          f"stop={result['reason']} output={output}")


if __name__ == '__main__':
    main()
