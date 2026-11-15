"""Export p5 history from the verified four-branch outer-10 snapshot CSV."""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE / "figures/common_outer10_progress/completed_outer_records.csv"
OUT = BASE / "figures/common_outer10_p5"
BRANCHES = ("ideal_sumse", "pilot_sumse", "ideal_pf", "pilot_pf")
COLORS = ("#126eae", "#d45b30", "#18845e", "#8d4fb3")


def main():
    with SOURCE.open(newline="", encoding="utf-8-sig") as handle:
        records = list(csv.DictReader(handle))
    by_branch = {name: [] for name in BRANCHES}
    for row in records:
        by_branch[row["branch"]].append(row)
    for name in BRANCHES:
        rows = by_branch[name]
        assert len(rows) == 10
        assert [int(float(row["outer"])) for row in rows] == list(range(1, 11))
        assert all(float(row["trueRate5"]) >= 0 for row in rows)
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "p5_history.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("branch", "outer", "p5UserTrueSEBitPerSecondPerHz"))
        for name in BRANCHES:
            writer.writerows((name, int(float(row["outer"])), row["trueRate5"]) for row in by_branch[name])
    plt.rcParams.update({"font.family": "Microsoft YaHei", "font.size": 12,
                         "axes.unicode_minus": False,
                         "figure.facecolor": "white", "axes.facecolor": "white"})
    fig, ax = plt.subplots(figsize=(9, 5.5), layout="constrained")
    for name, color in zip(BRANCHES, COLORS):
        rows = by_branch[name]
        ax.plot([int(float(row["outer"])) for row in rows],
                [float(row["trueRate5"]) for row in rows],
                linestyle="--" if "pilot" in name else "-", marker="o",
                markersize=3, color=color, label=name)
    ax.set(xlabel="完整外层迭代次数", ylabel="第5百分位用户真实SE（bit/s/Hz）",
           xlim=(1, 10), title="边缘用户频谱效率随迭代变化")
    ax.grid(alpha=0.22)
    ax.legend(frameon=False, fontsize=10)
    fig.savefig(OUT / "p5_history.png", dpi=300, facecolor="white")
    fig.savefig(OUT / "p5_history.pdf", facecolor="white")
    plt.close(fig)
    print({name: float(by_branch[name][-1]["trueRate5"]) for name in BRANCHES})


if __name__ == "__main__":
    main()
