"""Eight Chinese stage figures from frozen, independently audited outer-10 data."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "figures/common_outer10_delivery_v1"
BRANCHES = ("ideal_sumse", "pilot_sumse", "ideal_pf", "pilot_pf")
COLORS = dict(zip(BRANCHES, ("#126eae", "#d45b30", "#18845e", "#8d4fb3")))
LABELS = {"ideal_sumse": "理想CSI／sum-SE", "pilot_sumse": "导频CSI／sum-SE",
          "ideal_pf": "理想CSI／PF", "pilot_pf": "导频CSI／PF"}


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as file:
        return list(csv.DictReader(file))


def write_csv(stem, header, rows):
    with (OUT / f"{stem}.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(header)
        writer.writerows(rows)


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=300, facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def series(ax, history, field, branches=BRANCHES, scale=1):
    for branch in branches:
        rows = history[branch]
        ax.plot([int(float(row["outer"])) for row in rows],
                [float(row[field]) / scale for row in rows],
                color=COLORS[branch], linestyle="--" if "pilot" in branch else "-",
                linewidth=1.8, label=LABELS[branch])
    ax.grid(alpha=0.2)
    ax.set_xlabel("完整外层迭代次数")


def main():
    manifest = json.loads((BASE / "analysis/common_outer10_manifest.json").read_text(encoding="utf-8"))
    for branch in BRANCHES:
        item = manifest["branches"][branch]
        path = Path(item["path"])
        assert item["outer"] == 10 and path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"]
    rows = read_csv(BASE / "figures/common_outer10_progress/completed_outer_records.csv")
    history = {branch: [r for r in rows if r["branch"] == branch] for branch in BRANCHES}
    assert all(len(history[b]) == 10 and
               [int(float(r["outer"])) for r in history[b]] == list(range(1, 11)) for b in BRANCHES)
    users = read_csv(BASE / "analysis/common_outer10_eval/per_user_se.csv")
    assert len(users) == 210 and [int(r["userId"]) for r in users] == list(range(1, 211))
    rates = {b: np.array([float(r[{"ideal_sumse": "idealSumSE", "pilot_sumse": "pilotSumSE",
                                     "ideal_pf": "idealPF", "pilot_pf": "pilotPF"}[b]]) for r in users])
             for b in BRANCHES}
    nmse = np.array([float(r["sampleNmse"]) for r in users])
    assert np.all(np.isfinite(nmse)) and np.all(nmse >= 0)
    audit = json.loads((BASE / "analysis/common_outer10_eval/summary.json").read_text(encoding="utf-8"))
    for branch in BRANCHES:
        assert abs(rates[branch].sum() - audit[branch]["trueSumSE"]) < 1e-8
        assert abs(float(history[branch][-1]["trueRate5"]) - np.percentile(rates[branch], 5)) < 1e-10
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Microsoft YaHei", "font.size": 10,
                         "axes.unicode_minus": False, "figure.facecolor": "white",
                         "axes.facecolor": "white", "savefig.facecolor": "white"})

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7), layout="constrained")
    for ax, branches, title in zip(axes, (BRANCHES[:2], BRANCHES[2:]), ("sum-SE任务", "PF任务")):
        for branch in branches:
            sorted_rates = np.sort(rates[branch])
            ax.step(sorted_rates, np.arange(1, 211) / 210, where="post",
                    color=COLORS[branch], linestyle="--" if "pilot" in branch else "-",
                    label=LABELS[branch])
        ax.set(xlabel="用户真实SE（bit/s/Hz）", ylabel="累计比例（CDF）", ylim=(0, 1), title=title)
        ax.grid(alpha=0.2); ax.legend(frameon=False)
    save(fig, "01_user_true_se_cdf")
    write_csv("01_user_true_se_cdf", ("userId",) + BRANCHES,
              ((int(users[i]["userId"]), *(rates[b][i] for b in BRANCHES)) for i in range(210)))

    for stem, field, title, ylabel in (
            ("02_true_sumse", "trueSumSE", "真实总SE随迭代变化", "真实总SE（bit/s/Hz）"),
            ("03_true_pf", "truePF", "真实PF效用随迭代变化", "真实PF效用（无量纲）"),
            ("04_p5_user_se", "trueRate5", "边缘用户频谱效率随迭代变化", "第5百分位用户SE（bit/s/Hz）")):
        fig, ax = plt.subplots(figsize=(8.5, 4.8), layout="constrained")
        series(ax, history, field)
        ax.set(ylabel=ylabel, title=title, xlim=(1, 10)); ax.legend(frameon=False, fontsize=9)
        save(fig, stem)
        write_csv(stem, ("branch", "outer", field),
                  ((b, int(float(r["outer"])), float(r[field])) for b in BRANCHES for r in history[b]))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7), layout="constrained")
    for ax, branches, title, ylabel in zip(axes, (BRANCHES[:2], BRANCHES[2:]),
                                          ("sum-SE设计目标", "PF设计目标"),
                                          ("设计总SE（bit/s/Hz）", "设计PF效用（无量纲）")):
        series(ax, history, "estimatedObjective", branches)
        ax.set(title=title, ylabel=ylabel, xlim=(1, 10)); ax.legend(frameon=False)
    save(fig, "05_design_objective")
    write_csv("05_design_objective", ("branch", "outer", "estimatedObjective"),
              ((b, int(float(r["outer"])), float(r["estimatedObjective"]))
               for b in BRANCHES for r in history[b]))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7), layout="constrained")
    stop_rows = []
    for ax, field, tol, title in zip(axes, ("relativeObjectiveChange", "precoderChange"),
                                     (1e-4, 0.005), ("目标变化／阈值", "预编码变化／阈值")):
        for branch in BRANCHES:
            values = [float(r[field]) / tol for r in history[branch]]
            ax.plot(range(1, 11), values, color=COLORS[branch],
                    linestyle="--" if "pilot" in branch else "-", label=LABELS[branch])
            stop_rows.extend((branch, i, field, values[i-1], tol) for i in range(1, 11))
        ax.axhline(1, color="black", linewidth=1)
        ax.set_yscale("log")
        ax.set(xlabel="完整外层迭代次数", ylabel="变化量／对应阈值", title=title, xlim=(1, 10))
        ax.grid(alpha=0.2); ax.legend(frameon=False, fontsize=8)
    save(fig, "06_outer_stop_ratios")
    write_csv("06_outer_stop_ratios", ("branch", "outer", "field", "ratio", "threshold"), stop_rows)

    fig, ax = plt.subplots(figsize=(8.5, 4.8), layout="constrained")
    ax.step(np.sort(nmse) * 100, np.arange(1, 211) / 210, where="post", color="#8d4fb3")
    ax.set(xlabel="用户NMSE（%）", ylabel="累计比例（CDF）", ylim=(0, 1),
           title="用户信道估计误差分布")
    ax.grid(alpha=0.2)
    save(fig, "07_user_nmse_cdf")
    write_csv("07_user_nmse_cdf", ("userId", "sampleNmse", "displayNmsePercent"),
              ((i+1, nmse[i], 100*nmse[i]) for i in range(210)))

    delta_sum = rates["ideal_sumse"] - rates["pilot_sumse"]
    delta_pf = rates["ideal_pf"] - rates["pilot_pf"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7), layout="constrained")
    for ax, delta, title, color in zip(axes, (delta_sum, delta_pf), ("sum-SE任务", "PF任务"),
                                       ("#126eae", "#18845e")):
        ax.scatter(100*nmse, delta, s=18, alpha=0.7, color=color)
        ax.axhline(0, color="black", linewidth=1)
        ax.set(xlabel="用户NMSE（%）", ylabel="理想设计－导频设计的真实用户SE（bit/s/Hz）", title=title)
        ax.grid(alpha=0.2)
    save(fig, "08_nmse_paired_se_loss")
    write_csv("08_nmse_paired_se_loss",
              ("userId", "sampleNmse", "displayNmsePercent", "sumseIdealMinusPilot", "pfIdealMinusPilot"),
              ((i+1, nmse[i], 100*nmse[i], delta_sum[i], delta_pf[i]) for i in range(210)))

    (OUT / "identity.json").write_text(json.dumps({
        "modelId": "target_gaussian_static_pilot_lmmse_hypothesis_v1", "sceneSeed": 19,
        "commonCompletedOuter": 10, "userCount": 210,
        "snapshotHashes": {b: manifest["branches"][b]["sha256"] for b in BRANCHES},
        "stopping": "batch_limit; none met original outer stopping conditions",
        "source": "saved complete outer histories and independent common-true-channel audit",
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print("PILOT_OUTER10_DELIVERY_OK", len(list(OUT.glob("*.png"))))


if __name__ == "__main__":
    main()
