"""Independently recompute true log-det UE rates and TRP power from saved pilot checkpoints.

Run on four completed, frozen checkpoint copies; never point it at a live writer.
The default input and checkpoint directory are the original baseline files.
"""
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

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "reproduction/chasi_pilot_mmse"
INPUT = BASE / "inputs/target_static_gaussian_v1_seed19.mat"
DEFAULT_OUTPUT = BASE / "analysis/checkpoint_truth_20261027"
EXPECTED = "a55c8c16646d16ecc1150304f5676ff5c857d63d9c4c1d8c34ccd041dfc601da"
BRANCHES = ("ideal_sumse", "pilot_sumse", "ideal_pf", "pilot_pf")


def complex_value(raw):
    return raw["real"].astype("f8") + 1j * raw["imag"].astype("f8")


def last_scalar(file, history, field):
    value = history[field][0, -1]
    if isinstance(value, h5py.Reference):
        value = file[value][()]
    return float(np.asarray(value).ravel()[0])


def text16(raw):
    return "".join(chr(int(c)) for c in raw.ravel())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint-dir", type=Path, default=BASE / "checkpoints")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--expected-sha", default=EXPECTED)
    parser.add_argument("--nmse-csv", type=Path,
                        default=BASE / "analysis/input_nmse_20261027/per_user.csv")
    args = parser.parse_args()
    output = args.output_dir
    digest = hashlib.sha256(args.input.read_bytes()).hexdigest()
    assert digest == args.expected_sha.lower()
    output.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.input) as f:
        r = f["data/resources"]
        groups = []
        for g in range(17):
            h = complex_value(f[r["H_true"][g, 0]][()]).transpose(3, 2, 1, 0)
            ids = f[r["ueIds"][g, 0]][()].astype(int).ravel() - 1
            nu = float(f[r["resourceWeightNu"][g, 0]][()].ravel()[0])
            assert h.shape == (len(ids), 21, 4, 64)
            groups.append((h, ids, nu))
    summary = {}
    rates_by_branch = {}
    power_by_branch = {}
    for branch in BRANCHES:
        with h5py.File(args.checkpoint_dir / f"{branch}.mat") as f:
            state = f["state"]
            assert text16(state["inputSha256"][()]) == digest
            assert text16(state["branchId"][()]) == branch
            iteration = int(state["iter"][()].ravel()[0])
            history = state["history"]
            assert history["trueSumSE"].shape[1] == iteration
            saved_sum = last_scalar(f, history, "trueSumSE")
            saved_pf = last_scalar(f, history, "truePF")
            saved_p5 = last_scalar(f, history, "trueRate5")
            precoder = state["precoder"]
            assert precoder.shape == (17, 1)
            rates = np.zeros(210, dtype=np.float64)
            power = np.zeros(21, dtype=np.float64)
            scheduled = np.zeros(210, dtype=bool)
            for g, (h, ids, nu) in enumerate(groups):
                x = complex_value(f[precoder[g, 0]][()]).transpose(3, 2, 1, 0)
                assert x.shape == (21, len(ids), 64, 2)
                scheduled[ids] = True
                power += nu * np.sum(abs(x) ** 2, axis=(1, 2, 3))
                eff = np.einsum("klrt,litd->kird", h, x, optimize=True)
                for k in range(len(ids)):
                    all_cov = np.eye(4, dtype=complex) + np.einsum("ird,isd->rs", eff[k], eff[k].conj())
                    wanted = eff[k, k] @ eff[k, k].conj().T
                    interf = all_cov - wanted
                    sign_all, log_all = np.linalg.slogdet(all_cov)
                    sign_int, log_int = np.linalg.slogdet(interf)
                    assert np.isclose(sign_all, 1) and np.isclose(sign_int, 1)
                    rates[ids[k]] += nu * max(0.0, (log_all - log_int) / np.log(2))
        total = float(rates.sum())
        pf = float(np.log1p(rates[scheduled] / 0.01).sum())
        p5 = float(np.percentile(rates, 5))
        assert abs(total - saved_sum) < 1e-5, (branch, total, saved_sum)
        assert abs(pf - saved_pf) < 1e-5, (branch, pf, saved_pf)
        assert abs(p5 - saved_p5) < 1e-5, (branch, p5, saved_p5)
        assert power.max() <= 20 + 1e-7, (branch, power.max())
        summary[branch] = {"outer": iteration, "users": 210, "scheduledUsers": int(scheduled.sum()),
                           "trueSumSE": total, "savedTrueSumSE": saved_sum,
                           "truePF": pf, "savedTruePF": saved_pf,
                           "p5SE": p5, "savedP5SE": saved_p5,
                           "zeroSEUsers": int(np.count_nonzero(rates == 0)),
                           "bottom5MeanSE": float(np.sort(rates)[:11].mean()),
                           "maxTrpPowerW": float(power.max()),
                           "powerViolations": int(np.count_nonzero(power > 20 + 1e-7))}
        rates_by_branch[branch] = rates
        power_by_branch[branch] = power
        with (output / f"{branch}_trp_power.csv").open("w", newline="") as file:
            w = csv.writer(file); w.writerow(["trpId", "powerW", "budgetW"])
            w.writerows((l + 1, power[l], 20) for l in range(21))
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    nmse_path = args.nmse_csv
    with nmse_path.open(newline="") as file:
        nmse = np.array([float(row["sampleNmse"]) for row in csv.DictReader(file)])
    assert nmse.shape == (210,)
    with (output / "per_user_se.csv").open("w", newline="") as file:
        w = csv.writer(file)
        w.writerow(["userId", "sampleNmse", "idealSumSE", "pilotSumSE", "idealPF", "pilotPF"])
        for k in range(210):
            w.writerow([k + 1, nmse[k], *(rates_by_branch[b][k] for b in BRANCHES)])
    fig, axs = plt.subplots(2, 2, figsize=(12, 9), layout="constrained")
    for ax, branch in zip(axs.ravel(), BRANCHES):
        ax.scatter(nmse, rates_by_branch[branch], s=14, alpha=0.65)
        ax.set(xlabel="Per UE sample NMSE (dimensionless)", ylabel="True user SE (bit/s/Hz)",
               title=f"{branch} | completed outer {summary[branch]['outer']}")
        ax.grid(True, alpha=0.3)
    outer_note = (f"common completed outer {summary[BRANCHES[0]]['outer']}"
                  if len({row['outer'] for row in summary.values()}) == 1
                  else "unequal outer counts")
    fig.suptitle(f"Static pilot LMMSE hypothesis v1 | seed 19 | {outer_note}, diagnostic only")
    fig.savefig(output / "user_nmse_vs_true_se.png", dpi=300, facecolor="white")
    fig.savefig(output / "user_nmse_vs_true_se.pdf", facecolor="white")
    plt.close(fig)
    fig, axs = plt.subplots(2, 2, figsize=(13, 8), layout="constrained")
    for ax, branch in zip(axs.ravel(), BRANCHES):
        ax.bar(np.arange(1, 22), power_by_branch[branch], color="#126eae")
        ax.axhline(20, color="#d45b30", linestyle="--", label="20 W budget")
        ax.set(xlabel="TRP index", ylabel="Weighted transmit power (W)",
               title=f"{branch} | completed outer {summary[branch]['outer']}",
               xlim=(0.2, 21.8), ylim=(0, 22))
        ax.legend(); ax.grid(axis="y", alpha=0.3)
    fig.suptitle(f"Static pilot LMMSE hypothesis v1 | seed 19 | power feasibility, {outer_note}")
    fig.savefig(output / "trp_power.png", dpi=300, facecolor="white")
    fig.savefig(output / "trp_power.pdf", facecolor="white")
    plt.close(fig)
    if len({row["outer"] for row in summary.values()}) == 1:
        fig, axs = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
        for ax, branches, title in zip(axs,
                (("ideal_sumse", "pilot_sumse"), ("ideal_pf", "pilot_pf")),
                ("sum-SE design", "PF design")):
            for branch, color in zip(branches, ("#126eae", "#d45b30")):
                values = np.sort(rates_by_branch[branch])
                ax.step(values, np.arange(1, 211) / 210, where="post", color=color,
                        label=f"{branch} | p5={summary[branch]['p5SE']:.3f}")
            ax.set(xlabel="True user SE (bit/s/Hz)", ylabel="Cumulative proportion",
                   ylim=(0, 1), title=title)
            ax.grid(True, alpha=0.3); ax.legend()
        fig.suptitle(f"Static pilot LMMSE v1 | seed 19 | common Htrue | completed outer {summary[BRANCHES[0]]['outer']} | 210 UEs")
        fig.savefig(output / "user_cdf.png", dpi=300, facecolor="white")
        fig.savefig(output / "user_cdf.pdf", facecolor="white")
        plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
