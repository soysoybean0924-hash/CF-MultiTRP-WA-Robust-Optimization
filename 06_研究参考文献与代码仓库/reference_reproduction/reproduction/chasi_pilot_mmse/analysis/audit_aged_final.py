"""Independent true-channel evaluation of completed aged-v1 WA checkpoints."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BRANCHES = ("ideal_sumse", "pilot_sumse", "statistical_sumse",
            "ideal_pf", "pilot_pf", "statistical_pf")


def complex_value(raw):
    return raw["real"].astype("f8") + 1j * raw["imag"].astype("f8")


def text16(raw):
    return "".join(chr(int(c)) for c in raw.ravel())


def scalar_at(handle, history, field, index):
    value = history[field][0, index]
    if isinstance(value, h5py.Reference):
        value = handle[value][()]
    return float(np.asarray(value).ravel()[0])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--checkpoints", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    digest = hashlib.sha256(args.input.read_bytes()).hexdigest()
    assert digest == args.input.with_suffix(args.input.suffix + ".sha256").read_text().strip().lower()
    args.output.mkdir(parents=True, exist_ok=True)
    groups = []
    with h5py.File(args.input) as f:
        resources = f["data/resources"]
        for g in range(17):
            h = complex_value(f[resources["H_true"][g, 0]][()]).transpose(3, 2, 1, 0)
            ids = f[resources["ueIds"][g, 0]][()].astype(int).ravel() - 1
            nu = float(f[resources["resourceWeightNu"][g, 0]][()].ravel()[0])
            assert h.shape == (len(ids), 21, 4, 64)
            groups.append((h, ids, nu))
    summary = {}
    user_rates = {}
    trajectories = []
    for branch in BRANCHES:
        with h5py.File(args.checkpoints / f"{branch}.mat") as f:
            state = f["state"]
            assert text16(state["modelId"][()]) == "target_gaussian_aged_pilot_reuse_v1"
            assert text16(state["inputSha256"][()]) == digest
            assert text16(state["branchId"][()]) == branch
            assert bool(state["complete"][()].ravel()[0])
            assert text16(state["stopReason"][()]) == "target_native_tolerance"
            outer = int(state["iter"][()].ravel()[0])
            history = state["history"]
            assert history["trueSumSE"].shape[1] == outer
            for i in range(outer):
                trajectories.append((branch, i + 1,
                    scalar_at(f, history, "estimatedObjective", i),
                    scalar_at(f, history, "trueSumSE", i),
                    scalar_at(f, history, "truePF", i),
                    scalar_at(f, history, "relativeObjectiveChange", i),
                    scalar_at(f, history, "precoderChange", i),
                    int(scalar_at(f, history, "innerSolved", i))))
            saved_sum = scalar_at(f, history, "trueSumSE", outer - 1)
            saved_pf = scalar_at(f, history, "truePF", outer - 1)
            saved_p5 = scalar_at(f, history, "trueRate5", outer - 1)
            precoder = state["precoder"]
            rates = np.zeros(210)
            power = np.zeros(21)
            scheduled = np.zeros(210, dtype=bool)
            for g, (h, ids, nu) in enumerate(groups):
                x = complex_value(f[precoder[g, 0]][()]).transpose(3, 2, 1, 0)
                assert x.shape == (21, len(ids), 64, 2)
                scheduled[ids] = True
                power += nu * np.sum(abs(x) ** 2, axis=(1, 2, 3))
                eff = np.einsum("klrt,litd->kird", h, x, optimize=True)
                for k in range(len(ids)):
                    all_cov = np.eye(4) + np.einsum("ird,isd->rs", eff[k], eff[k].conj())
                    wanted = eff[k, k] @ eff[k, k].conj().T
                    interference = all_cov - wanted
                    sa, la = np.linalg.slogdet(all_cov)
                    si, li = np.linalg.slogdet(interference)
                    assert np.isclose(sa, 1) and np.isclose(si, 1)
                    rates[ids[k]] += nu * max(0.0, (la - li) / np.log(2))
        total = float(rates.sum())
        pf = float(np.log1p(rates[scheduled] / 0.01).sum())
        p5 = float(np.percentile(rates, 5))
        assert abs(total - saved_sum) < 1e-5, (branch, total, saved_sum)
        assert abs(pf - saved_pf) < 1e-5, (branch, pf, saved_pf)
        assert abs(p5 - saved_p5) < 1e-5, (branch, p5, saved_p5)
        assert power.max() <= 20 + 1e-7, (branch, power.max())
        summary[branch] = {"outer": outer, "trueSumSE": total, "truePF": pf,
                           "p5SE": p5, "maxTrpPowerW": float(power.max()),
                           "scheduledUsers": int(scheduled.sum()),
                           "zeroSEUsers": int(np.count_nonzero(rates == 0))}
        user_rates[branch] = rates
        with (args.output / f"{branch}_trp_power.csv").open("w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(("trpId", "powerW", "budgetW"))
            writer.writerows((l + 1, power[l], 20) for l in range(21))
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    with (args.output / "convergence.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("branch", "outer", "designObjective", "trueSumSE", "truePF",
                         "relativeObjectiveChange", "precoderChange", "innerSolved"))
        writer.writerows(trajectories)
    with (args.output / "per_user_se.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("userId", *BRANCHES))
        writer.writerows((k + 1, *(user_rates[b][k] for b in BRANCHES)) for k in range(210))
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    for axis, objective in zip(axes, ("sumse", "pf")):
        for prefix in ("ideal", "pilot", "statistical"):
            branch = f"{prefix}_{objective}"
            values = np.sort(user_rates[branch])
            axis.step(values, np.arange(1, 211) / 210, where="post", label=branch)
        axis.set(xlabel="True user SE (bit/s/Hz)", ylabel="Cumulative proportion",
                 title=objective, ylim=(0, 1))
        axis.grid(True, alpha=0.3)
        axis.legend()
    fig.savefig(args.output / "true_user_se_cdf.png", dpi=200)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    for axis, objective in zip(axes, ("sumse", "pf")):
        for prefix in ("ideal", "pilot", "statistical"):
            branch = f"{prefix}_{objective}"
            rows = [r for r in trajectories if r[0] == branch]
            axis.plot([r[1] for r in rows], [r[2] for r in rows], label=branch)
        axis.set(xlabel="Completed outer iteration", ylabel="Design objective", title=objective)
        axis.grid(True, alpha=0.3)
        axis.legend()
    fig.savefig(args.output / "design_objective_trace.png", dpi=200)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
