"""Preserve one equal-iteration, read-only view of the pilot campaign.

The formal runner remains the sole writer of its night checkpoints. This
watcher copies each branch only while its saved state is exactly outer 2.
"""
import hashlib
import json
import os
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "reproduction/chasi_pilot_mmse"
NIGHT = BASE / "checkpoints/night_20261027"
DEST = BASE / "checkpoints/common_outer2_20261027"
MANIFEST = BASE / "analysis/common_outer2_manifest.json"
EXPECTED_INPUT = "a55c8c16646d16ecc1150304f5676ff5c857d63d9c4c1d8c34ccd041dfc601da"
BRANCHES = ("ideal_sumse", "pilot_sumse", "ideal_pf", "pilot_pf")


def chars(data):
    return "".join(map(chr, data.ravel()))


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def identity(path):
    with h5py.File(path) as f:
        s = f["state"]
        return int(s["iter"][()].ravel()[0]), chars(s["branchId"][()]), chars(s["inputSha256"][()])


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    deadline = time.monotonic() + 3 * 3600
    captured = {}
    while time.monotonic() < deadline:
        for branch in BRANCHES:
            target = DEST / f"{branch}.mat"
            if target.is_file():
                outer, name, input_sha = identity(target)
                if (outer, name, input_sha) != (2, branch, EXPECTED_INPUT):
                    raise RuntimeError(f"Invalid existing snapshot {target}")
                captured[branch] = {"outer": 2, "sha256": sha(target), "path": str(target)}
                continue
            source = (BASE / "checkpoints/ideal_sumse.mat") if branch == "ideal_sumse" else (NIGHT / f"{branch}.mat")
            if not source.is_file() or Path(str(source) + ".partial.mat").exists():
                continue
            try:
                outer, name, input_sha = identity(source)
            except OSError:
                continue  # An atomic rename may have occurred during this read.
            if name != branch or input_sha != EXPECTED_INPUT:
                raise RuntimeError(f"Source identity mismatch: {source}")
            if outer > 2:
                raise RuntimeError(f"Missed outer 2 for {branch}: source already {outer}")
            if outer < 2:
                continue
            before = sha(source)
            temporary = Path(str(target) + ".partial-copy")
            shutil.copy2(source, temporary)
            if sha(temporary) != before or sha(source) != before:
                temporary.unlink(missing_ok=True)
                continue
            os.replace(temporary, target)
            captured[branch] = {"outer": 2, "sha256": before, "path": str(target)}
            print(f"SNAPSHOT {branch} outer=2 sha256={before}", flush=True)
        MANIFEST.write_text(json.dumps({"updatedUtc": datetime.now(timezone.utc).isoformat(), "branches": captured}, indent=2), encoding="utf-8")
        if len(captured) == 4:
            print("COMMON_OUTER2_COMPLETE", flush=True)
            return
        time.sleep(20)
    raise TimeoutError(f"Did not capture four outer-2 states within 3 hours: {list(captured)}")


if __name__ == "__main__":
    main()
