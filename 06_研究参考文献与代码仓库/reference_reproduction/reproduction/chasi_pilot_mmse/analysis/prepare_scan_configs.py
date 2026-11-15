"""Freeze one-variable static-pilot scan configurations without running MATLAB.

The baseline configuration is copied exactly except for the named measurement
control, scanId, and label. Input generation and solver runs are separate steps.
"""
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE / "configs/target_static_gaussian_v1.json"
DEST = BASE / "scans/static_v1/configs"
MANIFEST = BASE / "scans/static_v1/config_manifest.json"
POINTS = (
    ("srs_m20_p8", "pilotSnrDb", -20, "SRS observation SNR -20 dB, 8 reused pilots"),
    ("srs_0_p8", "pilotSnrDb", 0, "SRS observation SNR 0 dB, 8 reused pilots"),
    ("srs_m10_p26", "pilotSequences", 26, "SRS observation SNR -10 dB, 26 orthogonal user pilots"),
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    config = json.loads(SOURCE.read_text(encoding="utf-8"))
    assert sha(SOURCE) == "03125afba774a0bffffa584d371003edf8b841a838cac98db545ff757447d833"
    assert config["sceneSeed"] == 19 and config["measurementSeed"] == 1701
    assert config["initialSeed"] == 901 and config["pilotSnrDb"] == -10
    assert config["pilotSequences"] == 8
    DEST.mkdir(parents=True, exist_ok=True)
    manifest = {"baselineConfigSha256": sha(SOURCE), "baselineInputSha256":
                "a55c8c16646d16ecc1150304f5676ff5c857d63d9c4c1d8c34ccd041dfc601da",
                "sceneSeed": 19, "measurementSeed": 1701, "initialSeed": 901,
                "points": []}
    for scan_id, parameter, value, label in POINTS:
        candidate = dict(config)
        candidate[parameter] = value
        candidate["scanId"] = scan_id
        candidate["label"] = label
        changed = {key for key in config if candidate[key] != config[key]}
        assert changed == {parameter, "label"}
        path = DEST / f"{scan_id}.json"
        content = json.dumps(candidate, indent=2, ensure_ascii=False) + "\n"
        if path.exists() and path.read_text(encoding="utf-8") != content:
            raise FileExistsError(f"Refusing to overwrite altered scan config: {path}")
        path.write_text(content, encoding="utf-8")
        manifest["points"].append({"scanId": scan_id, "control": parameter,
                                   "value": value, "config": str(path.relative_to(BASE)),
                                   "sha256": sha(path), "inputStatus": "not_generated",
                                   "solverStatus": "not_started"})
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
