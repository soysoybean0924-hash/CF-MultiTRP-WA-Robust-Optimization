"""Plot only saved complete outer histories from one static-pilot campaign.

Run from repo root with default night checkpoints, or pass --checkpoint-dir.
Every curve labels its own completed outer count. Incomplete runs are not
interpolated and an ADMM final-step ratio is not an inner-step trajectory.
"""
import argparse
import csv
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "reproduction/chasi_pilot_mmse"
BRANCHES = ("ideal_sumse", "pilot_sumse", "ideal_pf", "pilot_pf")
COLORS = ("#126eae", "#d45b30", "#18845e", "#8d4fb3")
EXPECTED = "a55c8c16646d16ecc1150304f5676ff5c857d63d9c4c1d8c34ccd041dfc601da"
FIELDS = ("outer", "trueSumSE", "truePF", "trueRate5", "estimatedObjective",
          "relativeObjectiveChange", "precoderChange", "innerCount", "admmPrimal",
          "epsPrimal", "admmDual", "epsDual")


def text16(raw):
    return "".join(chr(int(c)) for c in raw.ravel())


def history(file, state, field):
    raw = state["history"][field][0]
    result = []
    for value in raw:
        if isinstance(value, h5py.Reference):
            value = file[value][()]
        result.append(float(np.asarray(value).ravel()[0]))
    return np.array(result)


def export(fig, output, stem):
    fig.savefig(output / f"{stem}.png", dpi=300, facecolor="white")
    fig.savefig(output / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint-dir", type=Path, default=BASE / "checkpoints/night_20261027")
    parser.add_argument("--output-dir", type=Path, default=BASE / "figures/night_20261027")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = {}
    for branch in BRANCHES:
        with h5py.File(args.checkpoint_dir / f"{branch}.mat") as file:
            state = file["state"]
            count = int(state["iter"][()].ravel()[0])
            assert text16(state["branchId"][()]) == branch
            assert text16(state["inputSha256"][()]) == EXPECTED
            fields = {field: history(file, state, field) for field in FIELDS}
            assert all(len(x) == count for x in fields.values())
            assert np.array_equal(fields["outer"], np.arange(1, count + 1))
            data[branch] = fields
    with (args.output_dir / "completed_outer_records.csv").open("w", newline="") as file:
        writer = csv.writer(file); writer.writerow(("branch",) + FIELDS)
        for branch in BRANCHES:
            for j in range(len(data[branch]["outer"])):
                writer.writerow((branch,) + tuple(data[branch][f][j] for f in FIELDS))
    plt.rcParams.update({"font.size": 10, "figure.facecolor": "white", "axes.facecolor": "white"})
    specs = (("trueSumSE", "True sum-SE (bit/s/Hz)", None),
             ("truePF", "True PF utility (natural log)", None),
             ("relativeObjectiveChange", "Relative outer objective change", 1e-4),
             ("precoderChange", "Relative precoder change", 5e-3),
             ("innerCount", "ADMM inner iterations", None),
             ("residualRatio", "Last ADMM residual / threshold", 1))
    fig, axs = plt.subplots(2, 3, figsize=(16, 9), layout="constrained")
    for ax, (field, label, threshold) in zip(axs.ravel(), specs):
        for branch, color in zip(BRANCHES, COLORS):
            row = data[branch]
            y = (np.maximum(row["admmPrimal"] / row["epsPrimal"],
                            row["admmDual"] / row["epsDual"])
                 if field == "residualRatio" else row[field])
            ax.plot(row["outer"], y, "o-", color=color, markersize=3,
                    label=f"{branch} (n={len(y)})")
        if threshold is not None:
            ax.axhline(threshold, color="black", linestyle="--", label=f"threshold={threshold:g}")
        if field in ("relativeObjectiveChange", "precoderChange", "residualRatio"):
            ax.set_yscale("log")
        ax.set(xlabel="Completed outer iteration", ylabel=label)
        ax.grid(True, alpha=0.3)
    axs[0, 0].legend(fontsize=8)
    counts = {len(data[branch]["outer"]) for branch in BRANCHES}
    coverage = (f"common completed outer {next(iter(counts))}"
                if len(counts) == 1 else "completed outer counts differ")
    fig.suptitle(f"Static pilot LMMSE hypothesis v1 | seed 19 | progress only; {coverage}")
    export(fig, args.output_dir, "outer_progress")
    fig, axs = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    for ax, subset, label in zip(axs,
            (BRANCHES[:2], BRANCHES[2:]),
            ("Design sum-SE (bit/s/Hz)", "Design PF utility (natural log)")):
        for branch in subset:
            ax.plot(data[branch]["outer"], data[branch]["estimatedObjective"],
                    "o-", label=f"{branch} (n={len(data[branch]['outer'])})")
        ax.set(xlabel="Completed outer iteration", ylabel=label)
        ax.grid(True, alpha=0.3); ax.legend()
    fig.suptitle("Design-channel objective; distinct from true performance")
    export(fig, args.output_dir, "design_objective")
    (args.output_dir / "identity.json").write_text(json.dumps({
        "modelId": "target_gaussian_static_pilot_lmmse_hypothesis_v1",
        "inputSha256": EXPECTED,
        "checkpointDir": str(args.checkpoint_dir),
        "completedOuter": {b: len(data[b]["outer"]) for b in BRANCHES},
        "note": "ADMM ratio is the last inner step of each outer only; no inner trajectory is implied"
    }, indent=2), encoding="utf-8")
    print({b: len(data[b]["outer"]) for b in BRANCHES})


if __name__ == "__main__":
    main()
