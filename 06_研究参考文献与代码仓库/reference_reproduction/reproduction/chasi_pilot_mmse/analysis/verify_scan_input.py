"""Verify a frozen scan input keeps the baseline physical scene and start point.

The scan may alter Hhat/C through measurement quality or pilot reuse. It must
not alter Htrue, scheduling, resource weights, or initial precoding.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

BASE = Path(__file__).resolve().parents[1]
SCAN = BASE / "scans/static_v1"
REFERENCE = BASE / "inputs/target_static_gaussian_v1_seed19.mat"
FIELDS = ("H_true", "ueIds", "servingMask", "resourceWeightNu",
          "prior_gain", "home_trp", "tx_rb_count")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entry(file, field, group):
    raw = file[f"data/resources/{field}"][group, 0]
    return file[raw][()] if isinstance(raw, h5py.Reference) else np.asarray(raw)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("scan_id", choices=("srs_m20_p8", "srs_0_p8", "srs_m10_p26"))
    args = parser.parse_args()
    manifest = json.loads((SCAN / "config_manifest.json").read_text(encoding="utf-8"))
    point = next(row for row in manifest["points"] if row["scanId"] == args.scan_id)
    config = BASE / point["config"]
    target = SCAN / "inputs" / f"{args.scan_id}.mat"
    assert sha(config) == point["sha256"]
    assert sha(REFERENCE) == manifest["baselineInputSha256"]
    assert target.is_file() and (target.with_suffix(".mat.sha256")).is_file()
    target_hash = sha(target)
    assert target_hash == target.with_suffix(".mat.sha256").read_text().strip().lower()
    config_data = json.loads(config.read_text(encoding="utf-8"))
    assert config_data["scanId"] == args.scan_id
    assert config_data[point["control"]] == point["value"]
    changed_estimates = 0
    with h5py.File(REFERENCE) as old, h5py.File(target) as new:
        assert old["data/resources/H_true"].shape == new["data/resources/H_true"].shape == (17, 1)
        for group in range(17):
            for field in FIELDS:
                if not np.array_equal(entry(old, field, group), entry(new, field, group)):
                    raise AssertionError(f"Physical scene differs: group {group + 1}, {field}")
            if not np.array_equal(entry(old, "H", group), entry(new, "H", group)):
                changed_estimates += 1
            old_precoder = old[old["data/initialPrecoder"][group, 0]][()]
            new_precoder = new[new["data/initialPrecoder"][group, 0]][()]
            assert np.array_equal(old_precoder, new_precoder), f"Initial precoder differs, group {group + 1}"
    assert changed_estimates > 0, "Measurement control changed no estimates"
    result = {"scanId": args.scan_id, "configSha256": point["sha256"],
              "inputSha256": target_hash, "baselineInputSha256": sha(REFERENCE),
              "physicalGroupsIdentical": 17, "initialPrecoderGroupsIdentical": 17,
              "groupsWithDifferentHhat": changed_estimates,
              "solverStatus": "not_started"}
    output = SCAN / "verification" / f"{args.scan_id}.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
