"""Independently recompute adjacent live transmit-covariance changes.

Read only. The running MATLAB checkpoint is never changed or copied.
"""

from __future__ import annotations

import csv
import hashlib
import time
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parent
BRANCH = "ideal_sumse"
STEM = f"umi19p4_s27_u10_t2_fc3.5_g17_paired_formal_{BRANCH}"
CHECKPOINT = ROOT / "checkpoints" / f"{STEM}_checkpoint.mat"
OUTPUT = ROOT / "audit" / f"{STEM}_adjacent_covariance_check.csv"
PAIRS = 3
POLL_SECONDS = 10
DEADLINE_SECONDS = 2400


def read_state() -> tuple[int, list[np.ndarray], float, float, str]:
    with h5py.File(CHECKPOINT, "r") as file:
        state = file["state"]
        outer = int(np.asarray(state["outer"]).item())
        refs = np.asarray(state["X"])
        cells = []
        digest = hashlib.sha256()
        for ref in refs.flat:
            record = np.asarray(file[ref])
            array = np.ascontiguousarray((record["real"] + 1j * record["imag"]).T)
            digest.update(array.view(np.float64).tobytes())
            cells.append(array)
        covariance = float(np.asarray(state["covarianceChange"]).reshape(-1)[-1])
        raw = float(np.asarray(state["variableChange"]).reshape(-1)[-1])
    return outer, cells, covariance, raw, digest.hexdigest()


def recompute(old: list[np.ndarray], new: list[np.ndarray]) -> tuple[float, float]:
    assert len(old) == len(new)
    numerator_squared = 0.0
    denominator_squared = 0.0
    raw_squared = 0.0
    base_power = 0.0
    for a, b in zip(old, new):
        assert a.shape == b.shape
        gram_a = a.conj().T @ a
        gram_b = b.conj().T @ b
        cross = a.conj().T @ b
        numerator_squared += max(
            0.0,
            np.vdot(gram_a, gram_a).real
            + np.vdot(gram_b, gram_b).real
            - 2 * np.vdot(cross, cross).real,
        )
        denominator_squared += np.vdot(gram_a, gram_a).real
        raw_squared += np.vdot(b - a, b - a).real
        base_power += np.vdot(a, a).real
    covariance = np.sqrt(numerator_squared) / max(1.0, np.sqrt(denominator_squared))
    raw = np.sqrt(raw_squared) / max(1.0, np.sqrt(base_power))
    return float(covariance), float(raw)


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    begin = time.monotonic()
    old = read_state()
    print("initial", old[0], old[4], flush=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ("branch", "fromOuter", "toOuter", "oldXsha256", "newXsha256",
             "independentCovarianceChange", "savedCovarianceChange", "absoluteError",
             "independentRawChange", "savedRawChange", "rawAbsoluteError")
        )
        completed = 0
        while completed < PAIRS and time.monotonic() - begin < DEADLINE_SECONDS:
            time.sleep(POLL_SECONDS)
            try:
                new = read_state()
            except (OSError, KeyError, ValueError):
                # A checkpoint may briefly be between writes; retry read-only.
                continue
            if new[0] <= old[0]:
                continue
            if new[0] != old[0] + 1:
                print("skipped", old[0], new[0], flush=True)
                old = new
                continue
            covariance, raw = recompute(old[1], new[1])
            row = (
                BRANCH, old[0], new[0], old[4], new[4],
                f"{covariance:.17g}", f"{new[2]:.17g}",
                f"{abs(covariance - new[2]):.17g}",
                f"{raw:.17g}", f"{new[3]:.17g}",
                f"{abs(raw - new[3]):.17g}",
            )
            writer.writerow(row)
            stream.flush()
            print("pair", row[1], row[2], row[7], row[10], flush=True)
            old = new
            completed += 1
    if completed < PAIRS:
        raise TimeoutError(f"Captured {completed} of {PAIRS} adjacent pairs")


if __name__ == "__main__":
    main()
