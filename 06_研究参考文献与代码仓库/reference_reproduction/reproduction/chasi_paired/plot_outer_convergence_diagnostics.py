"""Plot saved outer convergence diagnostics without touching the live solver."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUTPUT = ROOT / "figures" / "solver_diagnostics"
SCENE = "umi19p4_s27_u10_t2_fc3.5_g17"
RUN_TAG = "paired_formal"
BRANCHES = (
    ("ideal_sumse", "Ideal CSI / sum-SE", "#14519e", "-"),
    ("nominal_sumse", "Nonideal CSI / sum-SE", "#d45e1f", "--"),
    ("ideal_pf", "Ideal CSI / PF", "#1f8050", "-."),
    ("nominal_pf", "Nonideal CSI / PF", "#9257a8", ":"),
)
COLUMNS = (
    "branchIndex",
    "outer",
    "normalizedDesignObjectiveChange",
    "transmitCovarianceChange",
    "primalResidualRatio",
    "dualResidualRatio",
    "innerSolved",
)


def latest_rows(branch_index: int, branch: str) -> np.ndarray:
    stem = f"{SCENE}_{RUN_TAG}_{branch}"
    summaries = list(RESULTS.glob(f"{stem}_outer*_summary.json"))
    if not summaries:
        raise FileNotFoundError(f"No completed batch for {branch}")
    items = []
    for path in summaries:
        summary = json.loads(path.read_text(encoding="utf-8"))
        items.append((int(summary["evaluation"]["outerCount"]), path))
    outer, summary_path = max(items, key=lambda item: item[0])
    batch_stem = summary_path.name.removesuffix("_summary.json")
    saved = loadmat(RESULTS / f"{batch_stem}.mat", simplify_cells=True)[
        "diagnostic"
    ]
    objective = np.asarray(saved["objective"], dtype=float).reshape(-1)
    scale = np.asarray(saved["objectiveScale"], dtype=float).reshape(-1)
    diagnostic = np.loadtxt(RESULTS / f"{batch_stem}_diagnostics.csv", delimiter=",")
    objective_csv = np.loadtxt(RESULTS / f"{batch_stem}_objectives.csv", delimiter=",")
    assert len(objective) == len(scale) == outer + 1
    assert diagnostic.shape == (outer, 9)
    np.testing.assert_array_equal(diagnostic[:, 0], np.arange(1, outer + 1))
    np.testing.assert_array_equal(objective_csv[:, 0], np.arange(outer + 1))
    np.testing.assert_allclose(objective_csv[:, 1], objective, rtol=0, atol=1e-8)
    change = np.abs(np.diff(objective)) / scale[:-1]
    assert np.isfinite(change).all() and (change >= 0).all()
    rows = np.column_stack(
        (
            np.full(outer, branch_index),
            diagnostic[:, 0],
            change,
            diagnostic[:, 5],
            diagnostic[:, 1] / diagnostic[:, 2],
            diagnostic[:, 3] / diagnostic[:, 4],
            diagnostic[:, 7],
        )
    )
    return rows


def write_csv(path: Path, rows: np.ndarray) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(COLUMNS)
        for row in rows:
            writer.writerow(
                [int(row[0]), int(row[1]), *(f"{value:.17g}" for value in row[2:6]), int(row[6])]
            )
    os.replace(temporary, path)


def plot(rows_by_branch: list[np.ndarray], zoom: bool, stem: str) -> None:
    plt.rcParams.update(
        {"font.family": "Arial", "font.size": 11, "figure.facecolor": "white"}
    )
    fig, axes = plt.subplots(2, 1, figsize=(10, 7.4), layout="constrained")
    for panel, (column, threshold, ylabel) in enumerate(
        (
            (2, 1e-5, "Normalized design objective change"),
            (3, 1e-4, "Relative transmit covariance change"),
        )
    ):
        ax = axes[panel]
        for rows, (_, label, color, style) in zip(rows_by_branch, BRANCHES):
            ax.plot(
                rows[:, 1],
                rows[:, column],
                color=color,
                linestyle=style,
                linewidth=1.7,
                label=label,
            )
        ax.axhline(threshold, color="0.4", linestyle=":", linewidth=1)
        ax.set_yscale("log", nonpositive="mask")
        ax.set_xlabel("WMMSE outer iteration")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.14)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, loc="best", fontsize=9)
    fig.suptitle("Late outer iterations" if zoom else "Outer convergence diagnostics")
    for suffix in ("png", "pdf"):
        target = OUTPUT / f"{stem}.{suffix}"
        temporary = target.with_suffix(target.suffix + ".tmp")
        fig.savefig(temporary, format=suffix, dpi=300, facecolor="white")
        os.replace(temporary, target)
    plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    all_rows = [latest_rows(index, branch) for index, (branch, *_rest) in enumerate(BRANCHES, 1)]
    full_stem = f"{RUN_TAG}_outer_convergence_log"
    write_csv(OUTPUT / f"{full_stem}.csv", np.vstack(all_rows))
    plot(all_rows, False, full_stem)
    late_rows = [rows[-40:] for rows in all_rows]
    late_stem = f"{RUN_TAG}_late_outer_zoom"
    write_csv(OUTPUT / f"{late_stem}.csv", np.vstack(late_rows))
    plot(late_rows, True, late_stem)
    trend_path = OUTPUT / f"{RUN_TAG}_late_outer_trends.csv"
    temporary = trend_path.with_suffix(trend_path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "branchIndex", "startOuter", "endOuter",
                "objectiveChangeStart", "objectiveChangeEnd", "objectiveDecreasePercent",
                "covarianceChangeStart", "covarianceChangeEnd", "covarianceDecreasePercent",
                "objectiveDescendingSteps", "covarianceDescendingSteps", "stepCount",
                "objectiveEndOverThreshold", "covarianceEndOverThreshold",
            )
        )
        for rows in late_rows:
            objective_decrease = 100 * (1 - rows[-1, 2] / rows[0, 2])
            covariance_decrease = 100 * (1 - rows[-1, 3] / rows[0, 3])
            writer.writerow(
                (
                    int(rows[0, 0]), int(rows[0, 1]), int(rows[-1, 1]),
                    f"{rows[0, 2]:.17g}", f"{rows[-1, 2]:.17g}",
                    f"{objective_decrease:.8g}",
                    f"{rows[0, 3]:.17g}", f"{rows[-1, 3]:.17g}",
                    f"{covariance_decrease:.8g}",
                    int(np.count_nonzero(np.diff(rows[:, 2]) < 0)),
                    int(np.count_nonzero(np.diff(rows[:, 3]) < 0)),
                    len(rows) - 1,
                    f"{rows[-1, 2] / 1e-5:.8g}",
                    f"{rows[-1, 3] / 1e-4:.8g}",
                )
            )
    os.replace(temporary, trend_path)
    print("diagnostics", [int(rows[-1, 1]) for rows in all_rows], OUTPUT)


if __name__ == "__main__":
    main()
