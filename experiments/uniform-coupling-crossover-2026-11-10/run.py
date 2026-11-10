"""Locate changes of the largest pair difference as all bonds are scaled."""

import csv
import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
PREVIOUS = OUTPUT.parent / "uniform-coupling-sweep-2026-11-10"
MODEL = run_path(str(PREVIOUS / "run.py"))
PAIRS, CHANNELS = MODEL["PAIRS"], MODEL["CHANNELS"]
T, DT = MODEL["T"], MODEL["FINE_DT"]
J_STEP = 0.01
ROOT_TOLERANCE = 5e-8
# Reflection gives 01 = 23 and 02 = 13; rank representative pairs individually.
REPRESENTATIVES = ("01", "02", "03", "12")


def sample(coupling, dt):
    time = np.arange(round(T / dt) + 1) * dt
    opened, ring, _, audit = MODEL["measure"](time, coupling)
    delta = ring - opened
    power = np.trapezoid(delta**2, x=time, axis=0) / T
    scores = dict(zip(PAIRS, np.sqrt(power.sum(axis=1)).tolist()))
    channel_rms = {pair: dict(zip(CHANNELS, np.sqrt(row).tolist()))
                   for pair, row in zip(PAIRS, power)}
    i, k = PAIRS.index("01"), CHANNELS.index("YY")
    yy = delta[:, i, k]

    def peak(index):
        return {"t": float(time[index]), "delta": float(yy[index]),
                "open": float(opened[index, i, k]), "ring": float(ring[index, i, k])}

    return {
        "J": float(coupling),
        "pair_rms": scores,
        "pair_channel_rms": channel_rms,
        "leader": max(REPRESENTATIVES, key=scores.__getitem__),
        "outer_YY": {
            "rms": channel_rms["01"]["YY"],
            "positive_peak": peak(np.argmax(yy)),
            "negative_peak": peak(np.argmin(yy)),
            "abs_peak": peak(np.argmax(np.abs(yy))),
        },
        "audit": audit,
    }


def crossing(before, after, left, right, dt):
    def gap(row):
        return row["pair_rms"][before] - row["pair_rms"][after]

    assert gap(sample(left, dt)) > 0 > gap(sample(right, dt))
    while right - left > ROOT_TOLERANCE:
        middle = (left + right) / 2
        if gap(sample(middle, dt)) > 0:
            left = middle
        else:
            right = middle
    return sample((left + right) / 2, dt)


def checked_crossing(before, after, left, right):
    grids = {str(dt): crossing(before, after, left, right, dt) for dt in (2 * DT, DT, DT / 2)}
    main = grids[str(DT)]
    return {"from": before, "to": after, "bracket": [left, right], "J": main["J"],
            "by_dt": {step: row["J"] for step, row in grids.items()},
            "pair_rms": main["pair_rms"], "pair_channel_rms": main["pair_channel_rms"],
            "outer_YY": main["outer_YY"]}


def write_csv(rows):
    with (OUTPUT / "sweep.csv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["J", *(f"pair_{pair}" for pair in PAIRS), "outer_YY_rms",
                         "outer_YY_peak_signed", "outer_YY_peak_t", "outer_YY_positive_peak",
                         "outer_YY_negative_peak", "leader",
                         *(f"{pair}_{channel}_rms" for pair in PAIRS for channel in CHANNELS)])
        for row in rows:
            yy = row["outer_YY"]
            writer.writerow([f"{row['J']:.2f}", *row["pair_rms"].values(), yy["rms"],
                             yy["abs_peak"]["delta"], yy["abs_peak"]["t"],
                             yy["positive_peak"]["delta"], yy["negative_peak"]["delta"],
                             row["leader"],
                             *(row["pair_channel_rms"][p][a] for p in PAIRS for a in CHANNELS)])


def main():
    grid = np.linspace(1, 2, round(1 / J_STEP) + 1)
    rows = [sample(j, DT) for j in grid]
    adjacent = list(zip(rows, rows[1:]))
    leaders = [checked_crossing(a["leader"], b["leader"], a["J"], b["J"])
               for a, b in adjacent if a["leader"] != b["leader"]]
    seam_outer = [checked_crossing("03", "01", a["J"], b["J"]) for a, b in adjacent
                  if a["pair_rms"]["03"] > a["pair_rms"]["01"]
                  and b["pair_rms"]["03"] < b["pair_rms"]["01"]]
    yy_min = min(rows, key=lambda row: row["outer_YY"]["rms"])
    yy_max = max(rows, key=lambda row: row["outer_YY"]["rms"])

    previous = json.loads((PREVIOUS / "results.json").read_text(encoding="utf-8"))
    archive_error = max(abs(row["pair_channel_rms"][p][a]
                            - previous["cases"][str(row["J"])]["fine"]["pair_channels"][p][a]["rms"])
                        for row in (rows[0], rows[-1]) for p in PAIRS for a in CHANNELS)
    assert archive_error < 1e-12
    crossing_dt_shift = max(abs(value - root["J"]) for root in leaders + seam_outer
                            for value in root["by_dt"].values())
    assert crossing_dt_shift < 1e-6
    audit = {"archive_endpoint_rms_error": archive_error, "max_crossing_dt_shift": crossing_dt_shift,
             "max_state_norm_error": max(row["audit"]["max_norm_error"] for row in rows),
             "max_symmetry_error": max(row["audit"]["max_symmetry_error"] for row in rows)}
    results = {"N": 4, "Omega": previous["Omega"], "initial": previous["initial"], "T": T,
               "dt": DT, "J_step": J_STEP, "delta": "ring - open", "leader_crossings": leaders,
               "seam_vs_outer_crossings": seam_outer, "outer_YY_grid_min": yy_min,
               "outer_YY_grid_max": yy_max, "audit": audit}

    write_csv(rows)
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    for root in leaders + seam_outer:
        print(f"{root['from']} -> {root['to']} at J={root['J']:.7f}, dt audit {root['by_dt']}")
    for name, row in (("min", yy_min), ("max", yy_max)):
        print(f"Outer YY grid {name}: J={row['J']:.2f}, RMS={row['outer_YY']['rms']:.9f}")


if __name__ == "__main__":
    main()
