"""Read-only host gate for the serial pilot campaign (no WMI dependency)."""
import json
import sys

import psutil


def main():
    memory = psutil.virtual_memory()
    swap = psutil.swap_memory()
    battery = psutil.sensors_battery()
    matlab = []
    for process in psutil.process_iter(["pid", "name", "create_time"]):
        if "matlab" in (process.info["name"] or "").lower():
            try:
                command = process.cmdline()
            except (psutil.AccessDenied, psutil.NoSuchProcess):
                command = ["command line unavailable; conservatively block"]
            matlab.append({"pid": process.pid, "startedUnix": process.info["create_time"], "command": command})
    state = {
        "availablePhysicalGiB": round(memory.available / 2**30, 3),
        "availableVirtualApproxGiB": round((memory.available + swap.free) / 2**30, 3),
        "batteryPercent": None if battery is None else battery.percent,
        "powerPlugged": None if battery is None else battery.power_plugged,
        "matlab": matlab,
    }
    print(json.dumps(state, ensure_ascii=False))
    if matlab:
        return 4
    if battery is not None and not battery.power_plugged and battery.percent < 50:
        return 3
    if state["availablePhysicalGiB"] < 2 or state["availableVirtualApproxGiB"] < 4:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
