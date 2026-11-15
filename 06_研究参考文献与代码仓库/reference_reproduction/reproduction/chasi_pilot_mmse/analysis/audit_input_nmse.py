"""Read-only NMSE audit of the frozen static-pilot MATLAB v7.3 input.

Run at repository root with: python reproduction/chasi_pilot_mmse/analysis/audit_input_nmse.py
"""
import csv
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "reproduction/chasi_pilot_mmse/inputs/target_static_gaussian_v1_seed19.mat"
OUTPUT = ROOT / "reproduction/chasi_pilot_mmse/analysis/input_nmse_20261027"
EXPECTED = "a55c8c16646d16ecc1150304f5676ff5c857d63d9c4c1d8c34ccd041dfc601da"


def z(x):
    return x["real"].astype("f8") + 1j * x["imag"].astype("f8")


def cdf(x):
    x = np.sort(x[np.isfinite(x) & (x >= 0)])
    return x, np.arange(1, len(x) + 1) / len(x)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=SOURCE)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    parser.add_argument("--expected-sha", default=EXPECTED)
    parser.add_argument("--experiment-id", default="baseline_srs_m10_p8")
    args = parser.parse_args()
    digest = hashlib.sha256(args.input.read_bytes()).hexdigest()
    assert digest == args.expected_sha.lower(), digest
    args.output_dir.mkdir(parents=True, exist_ok=True)
    true_e = np.zeros((210, 21)); err_e = true_e.copy()
    post_e = true_e.copy(); prior_e = true_e.copy()
    serving = np.zeros((210, 21), bool); counts = np.zeros((210, 21), int)
    checks = []
    with h5py.File(args.input) as f:
        r = f["data/resources"]
        for g in range(17):
            get = lambda key: f[r[key][g, 0]][()]
            h = z(get("H_true")); hhat = z(get("H"))
            ce = get("C").astype("f8"); gain = get("prior_gain").astype("f8")
            users = get("ueIds").astype(int).ravel() - 1
            service = get("servingMask").astype(bool)
            weight = float(get("resourceWeightNu").ravel()[0])
            assert h.shape == (64, 4, 21, len(users))
            assert ce.shape == (4, 4, 21, len(users))
            assert gain.shape == (21, len(users)) and service.shape == (len(users), 21)
            t = np.sum(abs(h) ** 2, axis=(0, 1))
            e = np.sum(abs(h - hhat) ** 2, axis=(0, 1))
            p = 64 * np.trace(ce, axis1=0, axis2=1)
            q = 64 * 4 * gain
            actual = float(e.sum() / t.sum())
            saved = float(f[f["data/diagnostics/nmse"][g, 0]][()].ravel()[0])
            assert abs(actual - saved) < 5e-7, (g, actual, saved)
            checks.append({"group": g + 1, "saved": saved, "recomputed": actual})
            for j, k in enumerate(users):
                true_e[k] += weight * t[:, j]
                err_e[k] += weight * e[:, j]
                post_e[k] += weight * p[:, j]
                prior_e[k] += weight * q[:, j]
                serving[k] |= service[j]
                counts[k] += 1
    valid = true_e > 0
    sample = np.divide(err_e, true_e, out=np.full_like(true_e, np.nan), where=valid)
    theory = np.divide(post_e, prior_e, out=np.full_like(true_e, np.nan), where=prior_e > 0)
    masks = {"all": valid, "serving": valid & serving, "nonserving": valid & ~serving}
    summary = {"modelId": "target_gaussian_static_pilot_lmmse_hypothesis_v1",
               "experimentId": args.experiment_id,
               "inputSha256": digest, "sceneSeed": 19, "measurementSeedBase": 1701,
               "groups": 17, "users": 210, "trps": 21,
               "weight": "rb_count/273 per scheduled group occurrence",
               "sample": "weighted sum ||Htrue-Hhat||_F^2 / weighted sum ||Htrue||_F^2",
               "theory": "weighted sum 64*trace(Ce) / weighted sum 64*4*prior_gain",
               "subsets": {}, "groupChecks": checks}
    for name, mask in masks.items():
        summary["subsets"][name] = {
            "links": int(mask.sum()), "sampleNmse": float(err_e[mask].sum() / true_e[mask].sum()),
            "theoryNmse": float(post_e[mask].sum() / prior_e[mask].sum()),
            "trueEnergy": float(true_e[mask].sum()), "errorEnergy": float(err_e[mask].sum()),
            "posteriorEnergy": float(post_e[mask].sum()), "priorEnergy": float(prior_e[mask].sum())}
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    with (args.output_dir / "per_user.csv").open("w", newline="") as file:
        w = csv.writer(file)
        w.writerow(["userId", "scheduledGroups", "sampleNmse", "theoryNmse", "trueEnergy", "errorEnergy", "posteriorEnergy", "priorEnergy"])
        for k in range(210):
            w.writerow([k + 1, int(counts[k].max()), err_e[k].sum()/true_e[k].sum(), post_e[k].sum()/prior_e[k].sum(), true_e[k].sum(), err_e[k].sum(), post_e[k].sum(), prior_e[k].sum()])
    with (args.output_dir / "per_link.csv").open("w", newline="") as file:
        w = csv.writer(file)
        w.writerow(["userId", "trpId", "serving", "scheduledGroups", "sampleNmse", "theoryNmse", "trueEnergy", "errorEnergy", "posteriorEnergy", "priorEnergy"])
        for k in range(210):
            for l in range(21):
                w.writerow([k+1, l+1, int(serving[k,l]), counts[k,l], sample[k,l], theory[k,l], true_e[k,l], err_e[k,l], post_e[k,l], prior_e[k,l]])
    plt.rcParams.update({"font.size": 11, "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white"})
    fig, axs = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    for x, label, color in [(err_e.sum(1)/true_e.sum(1), "Sample, 210 UEs", "#126eae"), (post_e.sum(1)/prior_e.sum(1), "Theory, 210 UEs", "#d45b30")]:
        a, b = cdf(x); axs[0].step(a, b, where="post", label=label, color=color)
    for mask, label, color in [(masks["all"], "All 4410 links", "#126eae"), (masks["serving"], "Serving links", "#d45b30")]:
        a, b = cdf(sample[mask]); axs[1].step(a, b, where="post", label=label, color=color)
    axs[0].set(xlabel="Per UE NMSE (dimensionless)", ylabel="Cumulative proportion", ylim=(0,1), title="Per UE estimation quality")
    axs[1].set(xlabel="Per link sample NMSE (dimensionless)", ylabel="Cumulative proportion", ylim=(0,1), title="All and serving links")
    axs[1].set_xscale("log")
    for ax in axs:
        ax.grid(True, alpha=0.3); ax.legend()
    fig.suptitle(f"Static pilot LMMSE hypothesis v1 | {args.experiment_id} | seed 19 | 17 groups | RB-weighted energy")
    fig.savefig(args.output_dir / "nmse_distribution.png", dpi=300)
    fig.savefig(args.output_dir / "nmse_distribution.pdf")
    plt.close(fig)
    print(json.dumps(summary["subsets"], indent=2))


if __name__ == "__main__":
    main()
