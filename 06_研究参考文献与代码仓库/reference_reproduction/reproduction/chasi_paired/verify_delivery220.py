"""Read-only numerical audit of the eight outer-220 delivery figures."""

from pathlib import Path
import json

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent
FIG = BASE / "figures" / "delivery220_v1"
RESULTS = BASE / "results"
PREFIX = "umi19p4_s27_u10_t2_fc3.5_g17_paired_formal_"
BRANCHES = ["ideal_sumse", "nominal_sumse", "ideal_pf", "nominal_pf"]
NAMES = [
    "01_user_true_se_cdf",
    "02_mean_true_se",
    "03_p5_true_se",
    "04_mean_true_pf",
    "05_mean_design_objective",
    "06_outer_stop_ratios",
    "11_user_nmse_cdf",
    "12_user_nmse_vs_paired_loss",
]


def same(actual, expected, label, atol=1e-9):
    if not np.allclose(actual, expected, rtol=0, atol=atol):
        raise AssertionError(f"{label}: maximum difference {np.max(np.abs(actual-expected))}")


def main():
    for name in NAMES:
        for ext in ("png", "pdf", "fig", "csv"):
            path = FIG / f"{name}.{ext}"
            assert path.is_file() and path.stat().st_size > 100, path
    snap, objective, rates = [], [], []
    for branch in BRANCHES:
        stem = RESULTS / f"{PREFIX}{branch}_outer220"
        summary = json.loads(Path(f"{stem}_summary.json").read_text())
        assert summary["evaluation"]["outerCount"] == 220
        assert summary["evaluation"]["stopReason"] == "outer_limit"
        s = pd.read_csv(f"{stem}_snapshots.csv")
        o = pd.read_csv(f"{stem}_objectives.csv", header=None)
        r = pd.read_csv(f"{stem}_true_user_rates.csv", header=None)
        assert np.array_equal(s.outer, np.arange(221))
        assert np.array_equal(o.iloc[:, 0], np.arange(221))
        assert np.array_equal(r.iloc[:, 0], np.arange(1, 211))
        same([r.iloc[:, 1].sum()], [summary["evaluation"]["trueSumSE"]], branch)
        snap.append(s)
        objective.append(o.iloc[:, 1].to_numpy())
        rates.append(r.iloc[:, 1].to_numpy())
    rates = np.column_stack(rates)
    wide_columns = ["idealSumSE", "nominalSumSE", "idealPF", "nominalPF"]
    expected = {
        "02_mean_true_se": np.column_stack([s.trueSumSE_bit_s_Hz / 210 for s in snap]),
        "03_p5_true_se": np.column_stack([s.p5_bit_s_Hz for s in snap]),
        "04_mean_true_pf": np.column_stack([s.truePF / 210 for s in snap]),
        "05_mean_design_objective": np.column_stack(objective) / 210,
    }
    for name, values in expected.items():
        plotted = pd.read_csv(FIG / f"{name}.csv")
        same(plotted[wide_columns].to_numpy(), values, name)
    cdf = pd.read_csv(FIG / "01_user_true_se_cdf.csv")
    for j in range(4):
        part = cdf[cdf.branchIndex == j + 1]
        assert len(part) == 211
        sorted_rates = np.sort(rates[:, j])
        same(part.trueSE_bit_s_Hz.to_numpy(), np.r_[0, sorted_rates], "CDF SE")
        same(part.cdf.to_numpy(), np.arange(211) / 210, "CDF fraction")
        same(part.p5_bit_s_Hz.to_numpy(), np.quantile(sorted_rates, .05), "CDF p5")
    existing = pd.read_csv(BASE / "figures" / "solver_diagnostics" /
                           "paired_formal_outer_convergence_log.csv")
    stop = pd.read_csv(FIG / "06_outer_stop_ratios.csv")
    for j, suffix in enumerate(["IdealSumSE", "NominalSumSE", "IdealPF", "NominalPF"], 1):
        part = existing[existing.branchIndex == j]
        same(stop[f"objectiveRatio{suffix}"], part.normalizedDesignObjectiveChange / 1e-5,
             f"objective stop {suffix}")
        same(stop[f"covarianceRatio{suffix}"], part.transmitCovarianceChange / 1e-4,
             f"covariance stop {suffix}")
    energy = pd.read_csv(BASE / "analysis" / "energy_audit_20261026" /
                         "per_user_energy.csv")
    nmse = energy.actualErrorEnergy / energy.actualTrueEnergy
    same(nmse, energy.sampleNmse, "sample NMSE")
    nmse_cdf = pd.read_csv(FIG / "11_user_nmse_cdf.csv")
    same(nmse_cdf.userNmse, np.r_[min(nmse), np.sort(nmse)], "NMSE CDF")
    same(nmse_cdf.cdf, np.arange(211) / 210, "NMSE CDF fraction")
    loss = pd.read_csv(FIG / "12_user_nmse_vs_paired_loss.csv")
    assert np.array_equal(loss.userId, np.arange(1, 211))
    same(loss.sampleNmseAllLinks, nmse, "paired NMSE")
    same(loss.sumSETaskLoss_bit_s_Hz, rates[:, 0] - rates[:, 1], "sum-SE loss")
    same(loss.pfTaskLoss_bit_s_Hz, rates[:, 2] - rates[:, 3], "PF loss")
    assert np.count_nonzero(loss.sumSETaskLoss_bit_s_Hz < 0) == 4
    power = np.loadtxt(BASE / "figures" / "audit220" / "trp_power_W.csv", delimiter=",")
    assert power.shape == (21, 5) and np.all(power[:, 1:] <= 20 + 1e-7)
    for name in ("trp_power", "outer_cost"):
        for ext in ("png", "pdf", "fig"):
            assert (BASE / "figures" / "audit220" / f"{name}.{ext}").is_file()
    print("DELIVERY220_VERIFY_OK: 8 new figure groups, 2 indexed originals; data match source")


if __name__ == "__main__":
    main()
