"""Export unaltered late Htrue metrics and signed iteration increments."""

from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUTPUT = ROOT / "figures" / "solver_diagnostics"
SCENE = "umi19p4_s27_u10_t2_fc3.5_g17"
RUN = "paired_formal"
PAIRS = (
    ("true_sumse", ("ideal_sumse", "nominal_sumse"), "trueSumSE_bit_s_Hz"),
    ("true_pf", ("ideal_pf", "nominal_pf"), "truePF"),
)
WINDOW = 40
HEADER = ("branchIndex", "outer", "endpointOuter", "metricValue", "priorMetricValue", "signedDelta")


def latest_snapshots(branch: str, metric: str) -> tuple[int, list[float]]:
    stem = f"{SCENE}_{RUN}_{branch}"
    choices = []
    for path in RESULTS.glob(f"{stem}_outer*_summary.json"):
        summary = json.loads(path.read_text(encoding="utf-8"))
        choices.append((int(summary["evaluation"]["outerCount"]), path))
    endpoint, summary_path = max(choices, key=lambda item: item[0])
    batch_stem = summary_path.name.removesuffix("_summary.json")
    with (RESULTS / f"{batch_stem}_snapshots.csv").open(newline="") as stream:
        records = list(csv.DictReader(stream))
    assert len(records) == endpoint + 1
    assert [int(record["outer"]) for record in records] == list(range(endpoint + 1))
    values = [float(record[metric]) for record in records]
    return endpoint, values


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, branches, metric in PAIRS:
        rows = []
        for branch_index, branch in enumerate(branches, 1):
            endpoint, values = latest_snapshots(branch, metric)
            for outer in range(max(1, endpoint - WINDOW + 1), endpoint + 1):
                rows.append(
                    (
                        branch_index, outer, endpoint,
                        f"{values[outer]:.17g}",
                        f"{values[outer - 1]:.17g}",
                        f"{values[outer] - values[outer - 1]:.17g}",
                    )
                )
        for suffix in ("late_zoom", "signed_delta"):
            output = OUTPUT / f"{RUN}_{name}_{suffix}.csv"
            with output.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(HEADER)
                writer.writerows(rows)
        print(name, [(branch, latest_snapshots(branch, metric)[0]) for branch in branches], len(rows))


if __name__ == "__main__":
    main()
