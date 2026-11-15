"""Plot two explicitly labelled Htrue user-SE CDF comparisons from complete batches.

This observer never reads or writes the live checkpoint or true-metric trace.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
SCENE = "umi19p4_s27_u10_t2_fc3.5_g17_paired_formal"
BRANCHES = ("ideal_sumse", "nominal_sumse", "ideal_pf", "nominal_pf")
LABELS = ("Ideal CSI · sum-SE", "Nonideal CSI · sum-SE", "Ideal CSI · PF", "Nonideal CSI · PF")
COLORS = ("#164b83", "#d36b31", "#28794a", "#8654a1")
STYLES = ("-", "--", "-.", ":")


def latest_complete(branch: str) -> int:
    prefix = f"{SCENE}_{branch}_outer"
    found = []
    for path in (ROOT / "results").glob(f"{prefix}*_summary.json"):
        summary = json.loads(path.read_text(encoding="utf-8"))
        outer = int(summary["evaluation"]["outerCount"])
        if summary["evaluation"]["stopReason"] not in ("outer_limit", "paired_objective_covariance_admm_stop"):
            continue
        found.append(outer)
    return max(found)


def read_rates(branch: str, outer: int) -> np.ndarray:
    stem = f"{SCENE}_{branch}_outer{outer}"
    result = ROOT / "results"
    summary = json.loads((result / f"{stem}_summary.json").read_text(encoding="utf-8"))
    assert int(summary["evaluation"]["outerCount"]) == outer
    values = np.loadtxt(result / f"{stem}_true_user_rates.csv", delimiter=",")
    assert values.shape == (210, 2)
    assert np.array_equal(values[:, 0], np.arange(1, 211))
    assert np.all(np.isfinite(values[:, 1])) and np.all(values[:, 1] >= -1e-9)
    return np.maximum(values[:, 1], 0)


def draw(name: str, endpoints: dict[str, int]) -> None:
    output = ROOT / "figures" / "audit_20261025" / name
    output.mkdir(parents=True, exist_ok=True)
    sample_rows = []
    cdf_rows = []
    fig, (ax, zoom) = plt.subplots(1, 2, figsize=(10.2, 4.5), facecolor="white", gridspec_kw={"width_ratios": (3, 2)})
    for index, branch in enumerate(BRANCHES):
        rates = read_rates(branch, endpoints[branch])
        sample_rows.extend((branch, endpoints[branch], k + 1, f"{rate:.17g}") for k, rate in enumerate(rates))
        distinct, counts = np.unique(rates, return_counts=True)
        cumulative = np.cumsum(counts) / len(rates)
        assert cumulative[-1] == 1 and np.all((cumulative >= 0) & (cumulative <= 1))
        x = np.r_[distinct[0], distinct]
        y = np.r_[0, cumulative]
        cdf_rows.extend((branch, endpoints[branch], f"{xx:.17g}", f"{yy:.17g}") for xx, yy in zip(x, y))
        for axes in (ax, zoom):
            axes.step(x, y, where="post", color=COLORS[index], linestyle=STYLES[index], linewidth=1.7, label=LABELS[index])
    for axes in (ax, zoom):
        axes.set_facecolor("white")
        axes.set_ylim(0, 1)
        axes.set_xlabel("SE per UE [bit/s/Hz]")
        axes.set_ylabel("CDF")
        axes.grid(alpha=0.16)
    zoom.set_ylim(0, 0.12)
    zoom.set_xlim(left=0)
    ax.set_title("User spectral efficiency")
    zoom.set_title("Lower tail")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "user_se_cdf.png", dpi=300, facecolor="white")
    fig.savefig(output / "user_se_cdf.pdf", facecolor="white")
    plt.close(fig)
    for filename, header, rows in (
        ("user_se_samples.csv", ("branch", "outer", "ueId", "trueSE_bit_s_Hz"), sample_rows),
        ("user_se_cdf.csv", ("branch", "outer", "rate_bit_s_Hz", "cumulativeProbability"), cdf_rows),
    ):
        with (output / filename).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(header)
            writer.writerows(rows)
    print(name, endpoints)


if __name__ == "__main__":
    latest = {branch: latest_complete(branch) for branch in BRANCHES}
    common = min(latest.values())
    draw(f"common_outer{common}", {branch: common for branch in BRANCHES})
    draw("latest_completed_endpoints", latest)
