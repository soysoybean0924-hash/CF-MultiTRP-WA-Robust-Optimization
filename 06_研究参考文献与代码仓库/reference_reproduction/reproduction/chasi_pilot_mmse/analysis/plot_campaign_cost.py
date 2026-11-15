"""Plot measured batch wall time for the independent static-pilot campaign.

Only rows whose MATLAB process exited successfully are included. The first
checkpoint rounds predate this campaign, so cumulative time starts at the
first recorded resumed batch and is not total time from initialization.
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE / "logs/night_20261027/campaign_progress.csv"
OUT = BASE / "figures/night_20261027"
BRANCHES = ("ideal_sumse", "pilot_sumse", "ideal_pf", "pilot_pf")
COLORS = ("#126eae", "#d45b30", "#18845e", "#8d4fb3")


def main():
    rows = defaultdict(list)
    with SOURCE.open(newline="", encoding="utf-8-sig") as handle:
        for item in csv.DictReader(handle):
            if item["exitCode"] != "0":
                raise ValueError(f"Incomplete batch in campaign manifest: {item}")
            if item["status"] not in ("batch_limit", "target_native_tolerance"):
                raise ValueError(f"Unexpected stop reason: {item}")
            branch = item["branch"]
            if branch not in BRANCHES:
                raise ValueError(f"Unexpected branch: {branch}")
            rows[branch].append((int(item["outer"]), float(item["elapsedSeconds"])))
    for branch in BRANCHES:
        values = rows[branch]
        if not values or len({outer for outer, _ in values}) != len(values):
            raise ValueError(f"Missing or duplicated completed batch for {branch}")
        rows[branch] = sorted(values)
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "figure.facecolor": "white", "axes.facecolor": "white"})
    fig, axs = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    for branch, color in zip(BRANCHES, COLORS):
        values = rows[branch]
        x = [outer for outer, _ in values]
        minutes = [seconds / 60 for _, seconds in values]
        cumulative = []
        total = 0.0
        for minute in minutes:
            total += minute
            cumulative.append(total / 60)
        axs[0].plot(x, minutes, "o-", color=color, label=f"{branch} ({len(x)} batches)")
        axs[1].plot(x, cumulative, "o-", color=color, label=branch)
    axs[0].set(xlabel="Completed outer iteration", ylabel="MATLAB batch wall time (min)")
    axs[1].set(xlabel="Completed outer iteration", ylabel="Cumulative recorded batch time per branch (h)")
    for ax in axs:
        ax.grid(True, alpha=0.3)
    axs[0].legend(fontsize=8)
    fig.suptitle("Static pilot LMMSE v1 | seed 19 | serial resume; startup included, original rounds excluded")
    fig.savefig(OUT / "batch_cost.png", dpi=300, facecolor="white")
    fig.savefig(OUT / "batch_cost.pdf", facecolor="white")
    plt.close(fig)
    print({branch: (len(rows[branch]), rows[branch][-1][0]) for branch in BRANCHES})


if __name__ == "__main__":
    main()
