"""Audit the presentation revision without running either optimizer."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import whosmat

BASE = Path(__file__).resolve().parent
V1 = BASE / "figures" / "delivery220_v1"
V2 = BASE / "figures" / "delivery220_v2"
ENERGY = BASE / "analysis" / "energy_audit_20261026"
NAMES = (
    "01_user_true_se_cdf", "02_mean_true_se", "03_p5_true_se",
    "04_mean_true_pf", "05_mean_design_objective", "06_outer_stop_ratios",
    "11_user_nmse_cdf", "12_user_nmse_vs_paired_loss",
)
BRANCHES = ("ideal_sumse", "nominal_sumse", "ideal_pf", "nominal_pf")
STEM = "umi19p4_s27_u10_t2_fc3.5_g17_paired_formal_"


def same(a: object, b: object, label: str, atol: float = 1e-10) -> None:
    x, y = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    assert x.shape == y.shape and np.allclose(x, y, rtol=0, atol=atol,
                                             equal_nan=True), label


def main() -> None:
    assert V2.is_dir() and V1.is_dir()
    for name in NAMES:
        for ext in ("png", "pdf", "fig", "csv"):
            path = V2 / f"{name}.{ext}"
            assert path.is_file() and path.stat().st_size > 1024, path
        assert (V2 / f"{name}.png").read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        assert (V2 / f"{name}.pdf").read_bytes()[:4] == b"%PDF"
        assert any(name == "hgS_080000" for name, *_ in whosmat(V2 / f"{name}.fig"))
        old = pd.read_csv(V1 / f"{name}.csv")
        new = pd.read_csv(V2 / f"{name}.csv")
        pd.testing.assert_frame_equal(new[old.columns], old, check_exact=True)
        if name in ("11_user_nmse_cdf", "12_user_nmse_vs_paired_loss"):
            source = (old.userNmse if name.startswith("11")
                      else old.sampleNmseAllLinks)
            same(new.userNmse_percent, source * 100, name)

    endpoint = pd.read_csv(V2 / "endpoint_metrics.csv")
    stop = pd.read_csv(V2 / "endpoint_stop.csv")
    p5 = pd.read_csv(V2 / "endpoint_p5.csv")
    assert endpoint.branch.tolist() == list(BRANCHES)
    assert stop.branch.tolist() == list(BRANCHES)
    assert len(p5) == 4
    rates = []
    for j, branch in enumerate(BRANCHES):
        prefix = BASE / "results" / f"{STEM}{branch}_outer220"
        summary = json.loads(Path(f"{prefix}_summary.json").read_text())
        ev = summary["evaluation"]
        assert ev["outerCount"] == 220 and ev["stopReason"] == "outer_limit"
        user = pd.read_csv(f"{prefix}_true_user_rates.csv", header=None)
        assert user.iloc[:, 0].tolist() == list(range(1, 211))
        rates.append(user.iloc[:, 1].to_numpy())
        same(endpoint.trueMeanSE.iloc[j], ev["trueSumSE"] / 210,
             f"true SE {branch}")
        same(endpoint.trueMeanPF.iloc[j], ev["truePF"] / 210,
             f"true PF {branch}")
        same(endpoint.p5TrueSE.iloc[j], ev["fifthPercentileSE"],
             f"p5 {branch}")
        same(p5.p5_bit_s_Hz.iloc[j], np.quantile(rates[-1], .05,
             method="linear"), f"CDF p5 {branch}")
        same(stop.objectiveChange.iloc[j], ev["lastNormalizedObjectiveChange"],
             f"objective change {branch}")
        same(stop.covarianceChange.iloc[j], ev["lastCovarianceChange"],
             f"covariance change {branch}")
        assert stop.consecutivePassed.iloc[j] == ev["consecutivePassed"]
        assert bool(stop.innerAllSolved.iloc[j]) == ev["innerAllSolved"]
        same(endpoint.maxPowerW.iloc[j], ev["maxPowerW"], f"power {branch}")
    rates = np.stack(rates, axis=1)
    same(pd.read_csv(V2 / "03_p5_true_se.csv").iloc[-1, 1:],
         p5.p5_bit_s_Hz, "#1/#3 endpoint p5")
    loss = pd.read_csv(V2 / "12_user_nmse_vs_paired_loss.csv")
    same(loss.sumSETaskLoss_bit_s_Hz, rates[:, 0] - rates[:, 1], "sum-SE loss")
    same(loss.pfTaskLoss_bit_s_Hz, rates[:, 2] - rates[:, 3], "PF loss")
    assert len(loss) == 210

    links = pd.read_csv(ENERGY / "per_link_energy.csv")
    users = pd.read_csv(ENERGY / "per_user_energy.csv")
    assert len(links) == 210 * 21 and len(users) == 210
    assert int(links.serving.sum()) == 420
    for col in ("actualTrueEnergy", "actualErrorEnergy"):
        same(links.groupby("ue")[col].sum().to_numpy(), users[col],
             f"per-user {col}", atol=1e-2)
    global_nmse = links.actualErrorEnergy.sum() / links.actualTrueEnergy.sum()
    serving = links[links.serving == 1]
    serving_nmse = (serving.actualErrorEnergy.sum() /
                    serving.actualTrueEnergy.sum())
    serving_user = serving.groupby("ue")[["actualErrorEnergy",
                                               "actualTrueEnergy"]].sum()
    serving_user_nmse = (serving_user.actualErrorEnergy /
                         serving_user.actualTrueEnergy).to_numpy()
    all_user_nmse = users.actualErrorEnergy / users.actualTrueEnergy
    same(loss.sampleNmseAllLinks, all_user_nmse, "loss NMSE source")
    energy_summary = json.loads((ENERGY / "summary.json").read_text())
    same(global_nmse, energy_summary["sampleNmse"], "global NMSE")
    same(serving_nmse, energy_summary["servingSampleNmse"], "serving NMSE")
    audit = {
        "model": "chasi_paired synthetic additive CSI error",
        "outer": 220, "stopReason": "outer_limit",
        "userCount": 210, "allLinkCount": len(links),
        "servingLinkCount": len(serving),
        "globalEnergyRatioNmse": global_nmse,
        "servingEnergyRatioNmse": serving_nmse,
        "allLinkUserNmseMean": float(all_user_nmse.mean()),
        "allLinkUserNmseMedian": float(all_user_nmse.median()),
        "allLinkUserNmseP05": float(all_user_nmse.quantile(.05)),
        "allLinkUserNmseP95": float(all_user_nmse.quantile(.95)),
        "servingUserNmseMean": float(serving_user_nmse.mean()),
        "servingUserNmseMedian": float(np.median(serving_user_nmse)),
        "servingUserNmseP05": float(np.quantile(serving_user_nmse, .05)),
        "servingUserNmseP95": float(np.quantile(serving_user_nmse, .95)),
        "sumSETaskNegativePairedLossCount": int((loss.sumSETaskLoss_bit_s_Hz < 0).sum()),
        "pfTaskNegativePairedLossCount": int((loss.pfTaskLoss_bit_s_Hz < 0).sum()),
        "trueSEEvaluation": "run_chasi_paired_branch.m calls joint_matrix_rate(Htrue,out.X,serviceTRP,rbPerGroup)",
        "signalInterferenceStatistic": "not saved in existing figure/energy/summary CSV; no sensitivity conclusion",
    }
    (V2 / "data_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("DELIVERY220_V2_VERIFY_OK", json.dumps(audit, ensure_ascii=False))


if __name__ == "__main__":
    main()
