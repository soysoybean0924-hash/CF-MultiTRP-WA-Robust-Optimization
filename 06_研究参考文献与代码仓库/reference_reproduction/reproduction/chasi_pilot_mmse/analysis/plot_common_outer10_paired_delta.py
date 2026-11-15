"""Plot paired UE rate differences from the frozen common-outer-10 audit.

The plot is a within-model progress diagnostic, not a converged CSI loss or
an isolated estimate of pilot contamination.
"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE / "analysis/common_outer10_eval/per_user_se.csv"
SUMMARY = BASE / "analysis/common_outer10_eval/summary.json"
MANIFEST = BASE / "analysis/common_outer10_manifest.json"
OUTPUT = BASE / "figures/common_outer10_paired_delta"


def main():
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert set(summary) == set(manifest["branches"])
    assert all(item["outer"] == 10 for item in summary.values())
    assert all(item["outer"] == 10 for item in manifest["branches"].values())
    with SOURCE.open(newline="") as file:
        rows = list(csv.DictReader(file))
    assert len(rows) == 210 and [int(r["userId"]) for r in rows] == list(range(1, 211))
    nmse = np.array([float(r["sampleNmse"]) for r in rows])
    ideal_sumse = np.array([float(r["idealSumSE"]) for r in rows])
    pilot_sumse = np.array([float(r["pilotSumSE"]) for r in rows])
    ideal_pf = np.array([float(r["idealPF"]) for r in rows])
    pilot_pf = np.array([float(r["pilotPF"]) for r in rows])
    assert np.all(np.isfinite(nmse)) and np.all(nmse >= 0)
    assert np.all(np.isfinite([ideal_sumse, pilot_sumse, ideal_pf, pilot_pf]))
    for name, values in (("ideal_sumse", ideal_sumse), ("pilot_sumse", pilot_sumse),
                         ("ideal_pf", ideal_pf), ("pilot_pf", pilot_pf)):
        assert abs(values.sum() - summary[name]["trueSumSE"]) < 1e-8
    delta = {"sum-SE objective": ideal_sumse - pilot_sumse,
             "PF objective": ideal_pf - pilot_pf}
    OUTPUT.mkdir(parents=True, exist_ok=True)
    with (OUTPUT / "paired_user_delta.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("userId", "sampleNmse", "idealMinusPilotSumSE", "idealMinusPilotPF"))
        writer.writerows((i + 1, nmse[i], delta["sum-SE objective"][i],
                          delta["PF objective"][i]) for i in range(210))
    metrics = {name: {
        "users": 210, "outer": 10,
        "meanIdealMinusPilotBitPerSecondPerHz": float(value.mean()),
        "medianIdealMinusPilotBitPerSecondPerHz": float(np.median(value)),
        "usersPilotHigher": int(np.count_nonzero(value < 0)),
        "usersEqual": int(np.count_nonzero(value == 0)),
        "usersIdealHigher": int(np.count_nonzero(value > 0)),
    } for name, value in delta.items()}
    (OUTPUT / "summary.json").write_text(json.dumps({
        "modelId": "target_gaussian_static_pilot_lmmse_hypothesis_v1",
        "sceneSeed": 19, "commonCompletedOuter": 10,
        "deltaDefinition": "ideal design UE SE minus pilot design UE SE on common Htrue",
        "checkpointSha256": {k: v["sha256"] for k, v in manifest["branches"].items()},
        "evidenceBoundary": "Unconverged, single scene; observation noise and reuse are combined.",
        "objectives": metrics,
    }, indent=2), encoding="utf-8")
    plt.rcParams.update({"font.size": 10, "figure.facecolor": "white",
                         "axes.facecolor": "white", "savefig.facecolor": "white"})
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    for name, color in (("sum-SE objective", "#126eae"), ("PF objective", "#d45b30")):
        value = delta[name]
        ordered = np.sort(value)
        axes[0].step(ordered, np.arange(1, 211) / 210, where="post", color=color,
                     label=f"{name} (210 UEs)")
        axes[1].scatter(nmse, value, s=17, alpha=0.65, color=color, label=name)
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    axes[0].axvline(0, color="black", linestyle="--", linewidth=1)
    axes[0].set(xlabel="Ideal minus pilot UE SE (bit/s/Hz)",
                ylabel="Cumulative proportion", ylim=(0, 1))
    axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
    axes[1].set(xlabel="Per UE sample NMSE (dimensionless)",
                ylabel="Ideal minus pilot UE SE (bit/s/Hz)")
    fig.suptitle("Static pilot LMMSE hypothesis v1 | seed 19 | common true channel | completed outer 10 (unconverged)")
    fig.savefig(OUTPUT / "paired_user_delta.png", dpi=300)
    fig.savefig(OUTPUT / "paired_user_delta.pdf")
    plt.close(fig)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
