"""Fingerprint the four independent, current paired precoder checkpoints."""

from __future__ import annotations

import csv
import hashlib
import time
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parent
SCENE = "umi19p4_s27_u10_t2_fc3.5_g17_paired_formal"
BRANCHES = ("ideal_sumse", "nominal_sumse", "ideal_pf", "nominal_pf")
OUTPUT = ROOT / "audit" / "paired_formal_checkpoint_precoder_fingerprints.csv"


def fingerprint(branch: str) -> tuple:
    path = ROOT / "checkpoints" / f"{SCENE}_{branch}_checkpoint.mat"
    for attempt in range(5):
        try:
            with h5py.File(path, "r") as file:
                state = file["state"]
                outer = int(np.asarray(state["outer"]).item())
                refs = np.asarray(state["X"])
                digest = hashlib.sha256()
                power = 0.0
                for ref in refs.flat:
                    record = np.asarray(file[ref])
                    array = np.ascontiguousarray((record["real"] + 1j * record["imag"]).T)
                    digest.update(array.view(np.float64).tobytes())
                    power += float(np.vdot(array, array).real)
                return branch, outer, refs.shape[1], refs.shape[0], digest.hexdigest(), power
        except (OSError, KeyError, ValueError):
            if attempt == 4:
                raise
            time.sleep(2)
    raise AssertionError("unreachable")


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    rows = [fingerprint(branch) for branch in BRANCHES]
    assert len({row[4] for row in rows}) == 4
    with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("branch", "outer", "users", "resourceGroups", "precoderSha256", "sumFrobeniusSquared"))
        writer.writerows(rows)
    print(rows)


if __name__ == "__main__":
    main()
