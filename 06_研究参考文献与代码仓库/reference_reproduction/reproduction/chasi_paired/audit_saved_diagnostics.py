"""Audit the latest completed paired batches without reading the live checkpoint."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUTPUT = ROOT / "audit"
SCENE = "umi19p4_s27_u10_t2_fc3.5_g17"
RUN = "paired_formal"
BRANCHES = ("ideal_sumse", "nominal_sumse", "ideal_pf", "nominal_pf")


def write_csv(path: Path, header: tuple[str, ...], rows: list[tuple]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    batch_endpoints: dict[str, set[int]] = {branch: set() for branch in BRANCHES}
    with (ROOT / f"{RUN}_convergence_batches.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            batch_endpoints[row["branch"]].add(int(row["actualOuter"]))

    residual_rows: list[tuple] = []
    design_rows: list[tuple] = []
    covariance: dict[str, np.ndarray] = {}
    for branch in BRANCHES:
        stem = f"{SCENE}_{RUN}_{branch}"
        choices = []
        for path in RESULTS.glob(f"{stem}_outer*_summary.json"):
            summary = json.loads(path.read_text(encoding="utf-8"))
            choices.append((int(summary["evaluation"]["outerCount"]), path))
        outer, summary_path = max(choices, key=lambda item: item[0])
        batch_stem = summary_path.name.removesuffix("_summary.json")
        saved = loadmat(RESULTS / f"{batch_stem}.mat", simplify_cells=True)[
            "diagnostic"
        ]
        d = np.loadtxt(RESULTS / f"{batch_stem}_diagnostics.csv", delimiter=",")
        objective = np.loadtxt(RESULTS / f"{batch_stem}_objectives.csv", delimiter=",")
        with (RESULTS / f"{batch_stem}_snapshots.csv").open(newline="") as stream:
            snapshots = list(csv.DictReader(stream))
        assert d.shape == (outer, 9) and objective.shape == (outer + 1, 3)
        assert len(snapshots) == outer + 1
        assert np.array_equal(d[:, 0], np.arange(1, outer + 1))
        assert np.array_equal(objective[:, 0], np.arange(outer + 1))
        rho = np.asarray(saved["rhoTrace"]).reshape(-1)
        rho_changes = np.asarray(saved["rhoChanges"]).reshape(-1)
        inner_solved = np.asarray(saved["innerSolved"]).reshape(-1)
        covariance[branch] = d[:, 5]
        for idx in range(outer):
            n = idx + 1
            primal_ratio = d[idx, 1] / d[idx, 2]
            dual_ratio = d[idx, 3] / d[idx, 4]
            prev_primal = d[idx - 1, 1] / d[idx - 1, 2] if idx else float("nan")
            prev_dual = d[idx - 1, 3] / d[idx - 1, 4] if idx else float("nan")
            residual_rows.append(
                (
                    branch, n, int(n - 1 in batch_endpoints[branch]),
                    int(d[idx, 6]), int(bool(inner_solved[idx])),
                    "residual_thresholds" if inner_solved[idx] else "inner_limit",
                    rho[idx], int(rho_changes[idx]),
                    d[idx, 1], d[idx, 2], primal_ratio,
                    d[idx, 3], d[idx, 4], dual_ratio,
                    prev_primal, primal_ratio - prev_primal,
                    prev_dual, dual_ratio - prev_dual,
                )
            )
        for idx, row in enumerate(snapshots):
            assert int(row["outer"]) == idx
            design_rows.append(
                (
                    branch, idx, objective[idx, 1], objective[idx, 2],
                    float(row["trueSumSE_bit_s_Hz"]), float(row["truePF"]),
                    float(row["p5_bit_s_Hz"]),
                )
            )

    write_csv(
        OUTPUT / "paired_formal_residual_audit_all.csv",
        (
            "branch", "outer", "batchResumeStep", "innerCount", "innerSolved",
            "innerStopReason", "finalRho", "rhoChanges", "primal", "epsPrimal",
            "primalRatio", "dual", "epsDual", "dualRatio", "priorPrimalRatio",
            "primalRatioJump", "priorDualRatio", "dualRatioJump",
        ),
        residual_rows,
    )
    write_csv(
        OUTPUT / "paired_formal_design_vs_true_all.csv",
        (
            "branch", "outer", "designObjective", "designSumSE", "trueSumSE",
            "truePF", "trueP5SE",
        ),
        design_rows,
    )
    overlap_rows = []
    for first in range(len(BRANCHES)):
        for second in range(first + 1, len(BRANCHES)):
            a = covariance[BRANCHES[first]]
            b = covariance[BRANCHES[second]]
            n = min(len(a), len(b))
            overlap_rows.append(
                (
                    BRANCHES[first], BRANCHES[second], n,
                    float(np.corrcoef(a[:n], b[:n])[0, 1]),
                    float(np.max(np.abs(a[:n] - b[:n]))),
                    float(np.mean(np.abs(a[:n] - b[:n]))),
                    int(np.count_nonzero(a[:n] == b[:n])),
                )
            )
    write_csv(
        OUTPUT / "paired_formal_covariance_overlap.csv",
        ("firstBranch", "secondBranch", "commonOuters", "pearsonCorrelation",
         "maximumAbsoluteDifference", "meanAbsoluteDifference", "exactlyEqualCount"),
        overlap_rows,
    )
    print("latest_outer", {branch: len(covariance[branch]) for branch in BRANCHES})
    print("residual_ratio_max", max(row[10] for row in residual_rows), max(row[13] for row in residual_rows))
    for branch in BRANCHES:
        subset = [row for row in residual_rows if row[0] == branch and row[1] > 1]
        top_primal = sorted(subset, key=lambda row: abs(row[15]), reverse=True)[:5]
        top_dual = sorted(subset, key=lambda row: abs(row[17]), reverse=True)[:5]
        print(branch, "primal_jumps", [(r[1], round(r[15], 4), r[2], r[3], r[6]) for r in top_primal])
        print(branch, "dual_jumps", [(r[1], round(r[17], 4), r[2], r[3], r[6]) for r in top_dual])


if __name__ == "__main__":
    main()
