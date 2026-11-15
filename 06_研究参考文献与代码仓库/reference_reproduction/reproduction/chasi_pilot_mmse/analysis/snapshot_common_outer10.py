"""Capture immutable outer-10 views without touching serial writer checkpoints."""
import hashlib
import json
import os
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE / "checkpoints/night_20261027"
DEST = BASE / "checkpoints/common_outer10_20261027"
MANIFEST = BASE / "analysis/common_outer10_manifest.json"
INPUT_SHA = "a55c8c16646d16ecc1150304f5676ff5c857d63d9c4c1d8c34ccd041dfc601da"
BRANCHES = ("ideal_sumse", "pilot_sumse", "ideal_pf", "pilot_pf")
OUTER = 10


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path):
    with h5py.File(path) as handle:
        state = handle["state"]
        number = int(state["iter"][()].ravel()[0])
        branch = "".join(map(chr, state["branchId"][()].ravel()))
        input_hash = "".join(map(chr, state["inputSha256"][()].ravel()))
    return number, branch, input_hash


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    found = {}
    deadline = time.monotonic() + 3 * 3600
    while time.monotonic() < deadline:
        for branch in BRANCHES:
            target = DEST / f"{branch}.mat"
            if target.is_file():
                if identity(target) != (OUTER, branch, INPUT_SHA):
                    raise RuntimeError(f"Existing snapshot identity mismatch: {target}")
                found[branch] = {"outer": OUTER, "sha256": sha(target), "path": str(target)}
                continue
            source = SOURCE / f"{branch}.mat"
            if not source.is_file() or Path(str(source) + ".partial.mat").exists():
                continue
            try:
                number, name, input_hash = identity(source)
            except OSError:
                continue
            if (name, input_hash) != (branch, INPUT_SHA):
                raise RuntimeError(f"Source identity mismatch: {source}")
            if number > OUTER:
                raise RuntimeError(f"Missed exact outer {OUTER} snapshot: {branch} at {number}")
            if number < OUTER:
                continue
            before = sha(source)
            temporary = Path(str(target) + ".partial-copy")
            shutil.copy2(source, temporary)
            if sha(temporary) != before or sha(source) != before:
                temporary.unlink(missing_ok=True)
                continue
            os.replace(temporary, target)
            found[branch] = {"outer": OUTER, "sha256": before, "path": str(target)}
            print(f"SNAPSHOT {branch} outer={OUTER} sha256={before}", flush=True)
        MANIFEST.write_text(json.dumps({"updatedUtc": datetime.now(timezone.utc).isoformat(),
                                        "branches": found}, indent=2) + "\n", encoding="utf-8")
        if len(found) == len(BRANCHES):
            print("COMMON_OUTER10_COMPLETE", flush=True)
            return
        time.sleep(20)
    raise TimeoutError(f"Only captured {list(found)} at outer {OUTER}")


if __name__ == "__main__":
    main()
