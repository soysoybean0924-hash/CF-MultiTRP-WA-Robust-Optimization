"""Audit one saved paired-CSI stage without starting another MATLAB process.

Run with a Python environment containing NumPy, SciPy, h5py, and Matplotlib.
The power command must run while the checkpoint still has the requested outer
index; it does not change the checkpoint or the experiment.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


BASE = Path(__file__).resolve().parent
RESULTS = BASE / "results"
CHECKPOINTS = BASE / "checkpoints"
BRANCHES = ("ideal_sumse", "nominal_sumse", "ideal_pf", "nominal_pf")


def one(pattern: str, folder: Path = RESULTS) -> Path:
    paths = list(folder.glob(pattern))
    if len(paths) != 1:
        raise RuntimeError(f"Expected one {pattern} in {folder}; found {len(paths)}")
    return paths[0]


def read_csv(path: Path) -> list[list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle))


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def audit_power(outer: int, branch: str) -> Path:
    import h5py
    import numpy as np

    if branch not in BRANCHES[:2]:
        raise ValueError("Power audit is for the two sum-SE branches")
    rate_path = one(f"*paired_formal_{branch}_outer{outer}_true_user_rates.csv")
    rates = read_csv(rate_path)
    checkpoint = one(f"*paired_formal_{branch}_checkpoint.mat", CHECKPOINTS)
    summary = json.loads(
        one(f"*paired_formal_{branch}_outer{outer}_summary.json").read_text(
            encoding="utf-8"
        )
    )
    rows: list[list[object]] = []
    with h5py.File(checkpoint, "r") as handle:
        state = handle["state"]
        saved_outer = int(state["outer"][()].item())
        if saved_outer != outer:
            raise RuntimeError(
                f"{branch}: checkpoint is outer {saved_outer}, expected {outer}"
            )
        cells = state["X"]
        service = state["serviceTRP"][()]
        users = int(state["K"][()].item())
        groups = int(state["G"][()].item())
        streams = int(state["streams"][()].item())
        if cells.shape != (groups, users) or service.shape != (2, users):
            raise RuntimeError("MATLAB cell or service-link dimensions changed")
        if len(rates) != users:
            raise RuntimeError("Endpoint user count differs from checkpoint")
        trp_power = np.zeros(int(state["L"][()].item()))
        for k in range(users):
            user_id, rate_text = rates[k]
            if int(user_id) != k + 1:
                raise RuntimeError("User IDs are not in checkpoint order")
            rate = float(rate_text)
            powers = np.zeros(2)
            for group in range(groups):
                block = handle[cells[group, k]][()]
                if block.shape[0] != streams or block.shape[1] % 2:
                    raise RuntimeError("Unexpected saved precoder dimensions")
                half = block.shape[1] // 2
                magnitude2 = block["real"] ** 2 + block["imag"] ** 2
                powers[0] += magnitude2[:, :half].sum()
                powers[1] += magnitude2[:, half:].sum()
            trp_power[int(service[0, k]) - 1] += powers[0]
            trp_power[int(service[1, k]) - 1] += powers[1]
            rows.append(
                [
                    branch,
                    outer,
                    k + 1,
                    rate,
                    int(rate == 0.0),
                    int(service[0, k]),
                    int(service[1, k]),
                    float(powers[0]),
                    float(powers[1]),
                    float(powers.sum()),
                ]
            )
        if not np.isclose(
            trp_power.max(), summary["evaluation"]["maxPowerW"], atol=1e-8
        ) or trp_power.max() > float(state["Pmax"][()].item()) + 1e-8:
            raise RuntimeError("Per-user link powers disagree with the TRP power audit")
    output = BASE / "analysis" / f"outer{outer}" / f"{branch}_user_power.csv"
    write_csv(
        output,
        [
            "branch",
            "outer",
            "userId",
            "trueRate_bit_s_Hz",
            "isExactZero",
            "serviceTRP1",
            "serviceTRP2",
            "powerTRP1_W",
            "powerTRP2_W",
            "totalUserPower_W",
        ],
        rows,
    )
    print(f"{branch} outer {outer} exact-zero users:")
    for row in rows:
        if row[4]:
            print(f"  UE {row[2]}: {row[7]:.6g} + {row[8]:.6g} = {row[9]:.6g} W")
    return output


def compare(outer: int) -> tuple[Path, Path]:
    import numpy as np
    from scipy.io import loadmat

    output_dir = BASE / "analysis" / f"outer{outer}"
    comparison: list[list[object]] = []
    trends: list[list[object]] = []
    for branch in BRANCHES:
        stem = f"*paired_formal_{branch}_outer{outer}"
        summary = json.loads(one(stem + "_summary.json").read_text(encoding="utf-8"))
        if summary["evaluation"]["outerCount"] != outer:
            raise RuntimeError(f"{branch} has not reached outer {outer}")
        objectives = np.asarray(read_csv(one(stem + "_objectives.csv")), dtype=float)
        with one(stem + "_snapshots.csv").open(newline="", encoding="utf-8") as handle:
            snapshots = list(csv.DictReader(handle))
        diagnostic = loadmat(
            one(stem + ".mat"),
            squeeze_me=True,
            struct_as_record=False,
            variable_names=["diagnostic"],
        )["diagnostic"]
        scales = np.atleast_1d(diagnostic.objectiveScale)
        if len(objectives) != outer + 1 or len(snapshots) != outer + 1:
            raise RuntimeError(f"{branch}: missing an outer snapshot")
        if len(scales) != outer + 1:
            raise RuntimeError(f"{branch}: missing an objective scale")
        aligned_key = "trueSumSE_bit_s_Hz" if branch.endswith("sumse") else "truePF"
        first = len(comparison)
        for n, row in enumerate(snapshots):
            if int(row["outer"]) != n or int(objectives[n, 0]) != n:
                raise RuntimeError(f"{branch}: outer indices disagree")
            design = objectives[n, 1]
            actual = float(row[aligned_key])
            design_change = design - objectives[n - 1, 1] if n else math.nan
            actual_change = actual - float(snapshots[n - 1][aligned_key]) if n else math.nan
            normalized = abs(design_change) / scales[n - 1] if n else math.nan
            comparison.append(
                [
                    branch,
                    n,
                    design,
                    actual,
                    design - actual,
                    design_change,
                    actual_change,
                    normalized,
                    float(row["precoderChange"]) if n else math.nan,
                    float(row["covarianceChange"]) if n else math.nan,
                    float(row["trueSumSE_bit_s_Hz"]),
                    float(row["truePF"]),
                ]
            )
        aligned = comparison[first:]
        if branch.startswith("ideal") and max(abs(r[4]) for r in aligned) > 1e-8:
            raise RuntimeError(f"{branch}: ideal design and Htrue evaluation disagree")
        latest = aligned[-10:]
        # Ten transitions, outer N-9 through N, use the original objectiveScale.
        trends.append(
            [
                branch,
                outer - 9,
                outer,
                summary["evaluation"]["stopReason"],
                summary["evaluation"]["optimizerSeconds"],
                latest[-1][2] - aligned[-11][2],
                latest[-1][3] - aligned[-11][3],
                latest[0][7],
                latest[-1][7],
                max(r[7] for r in latest),
                latest[0][8],
                latest[-1][8],
                latest[0][9],
                latest[-1][9],
                max(r[9] for r in latest),
                sum(r[7] <= 1e-5 and r[9] <= 1e-4 for r in latest),
            ]
        )
    comparison_path = output_dir / "design_vs_true.csv"
    write_csv(
        comparison_path,
        [
            "branch",
            "outer",
            "designObjective",
            "trueAlignedMetric",
            "designMinusTrue",
            "designStepChange",
            "trueStepChange",
            "normalizedDesignChange",
            "precoderChange",
            "covarianceChange",
            "trueSumSE_bit_s_Hz",
            "truePF",
        ],
        comparison,
    )
    trend_path = output_dir / "last10_trends.csv"
    write_csv(
        trend_path,
        [
            "branch",
            "fromOuter",
            "toOuter",
            "stopReason",
            "batchOptimizerSeconds",
            "designGain10",
            "trueGain10",
            "normalizedDesignChangeFirst",
            "normalizedDesignChangeLast",
            "normalizedDesignChangeMax",
            "precoderChangeFirst",
            "precoderChangeLast",
            "covarianceChangeFirst",
            "covarianceChangeLast",
            "covarianceChangeMax",
            "roundsMeetingBothOuterThresholds",
        ],
        trends,
    )
    audit_zero_history(outer, output_dir)
    plot_comparison(comparison, output_dir)
    return comparison_path, trend_path


def audit_zero_history(outer: int, output_dir: Path) -> None:
    history: list[list[object]] = []
    counts: list[list[object]] = []
    for branch in BRANCHES[:2]:
        for summary_path in RESULTS.glob(
            f"*paired_formal_{branch}_outer*_summary.json"
        ):
            saved = json.loads(summary_path.read_text(encoding="utf-8"))
            endpoint = int(saved["evaluation"]["outerCount"])
            if endpoint > outer:
                continue
            user_path = one(
                f"*paired_formal_{branch}_outer{endpoint}_true_user_rates.csv"
            )
            rates = read_csv(user_path)
            if len(rates) != 210:
                raise RuntimeError(f"{branch} outer {endpoint}: expected 210 users")
            zeros = [(int(r[0]), float(r[1])) for r in rates if float(r[1]) == 0.0]
            counts.append([branch, endpoint, len(zeros), 210])
            history.extend([branch, endpoint, user_id, rate] for user_id, rate in zeros)
    counts.sort(key=lambda row: (BRANCHES.index(row[0]), row[1]))
    history.sort(key=lambda row: (BRANCHES.index(row[0]), row[1], row[2]))
    write_csv(
        output_dir / "exact_zero_user_history.csv",
        ["branch", "outer", "userId", "trueRate_bit_s_Hz"],
        history,
    )
    write_csv(
        output_dir / "exact_zero_counts.csv",
        ["branch", "outer", "exactZeroCount", "totalUsers"],
        counts,
    )


def plot_comparison(rows: list[list[object]], output_dir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "Microsoft YaHei",
            "axes.unicode_minus": False,
            "font.size": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), constrained_layout=True)
    for axis, kind, label in zip(
        axes,
        ("sumse", "pf"),
        ("总频谱效率 (bit/s/Hz)", "PF 效用 (无量纲)"),
    ):
        ideal = [r for r in rows if r[0] == "ideal_" + kind]
        nominal = [r for r in rows if r[0] == "nominal_" + kind]
        axis.plot(
            [r[1] for r in ideal],
            [r[3] for r in ideal],
            color="#607080",
            linewidth=1.8,
            label="理想 CSI：真实评价",
        )
        axis.plot(
            [r[1] for r in nominal],
            [r[2] for r in nominal],
            color="#c3542c",
            linestyle="--",
            linewidth=1.8,
            label="非理想 CSI：设计目标",
        )
        axis.plot(
            [r[1] for r in nominal],
            [r[3] for r in nominal],
            color="#1769a3",
            linewidth=1.8,
            label="非理想 CSI：真实评价",
        )
        axis.set_title("sum-SE 目标" if kind == "sumse" else "PF 目标")
        axis.set_xlabel("WMMSE 外层迭代次数")
        axis.set_ylabel(label)
        axis.grid(color="#d9dfe3", alpha=0.65, linewidth=0.6)
        axis.spines[["top", "right"]].set_visible(False)
        axis.legend(frameon=False, fontsize=8, loc="best")
    for extension in ("png", "pdf"):
        fig.savefig(output_dir / f"design_vs_true.{extension}", dpi=250)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer", type=int, required=True)
    parser.add_argument("--power-branch", choices=BRANCHES[:2])
    parser.add_argument("--compare", action="store_true")
    args = parser.parse_args()
    if args.outer < 10 or not (args.power_branch or args.compare):
        parser.error("Choose --compare or --power-branch and an outer >= 10")
    if args.power_branch:
        print(audit_power(args.outer, args.power_branch))
    if args.compare:
        print(*compare(args.outer), sep="\n")


if __name__ == "__main__":
    main()
