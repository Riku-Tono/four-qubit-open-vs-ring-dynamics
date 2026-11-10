"""Decompose the outer–central crossing by channel, time, and observation window."""

import csv
import json
from functools import cache
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
SWEEP = OUTPUT.parent / "uniform-coupling-sweep-2026-11-10"
PREVIOUS = OUTPUT.parent / "uniform-coupling-crossover-2026-11-10"
MODEL = run_path(str(SWEEP / "run.py"))
PAIRS, CHANNELS = MODEL["PAIRS"], MODEL["CHANNELS"]
OUTER, CENTER = PAIRS.index("01"), PAIRS.index("12")
REPRESENTATIVES = [PAIRS.index(pair) for pair in ("01", "02", "03", "12")]
T, DT = MODEL["T"], MODEL["FINE_DT"]
COUPLINGS = np.linspace(1.0, 2.0, 101)
WINDOWS = np.linspace(2.688, 4.688, 101).round(3)
WINDOW_CHECKPOINTS = (2.688, 3.188, 3.588, T, 3.788, 4.188, 4.688)
TIME_BANDS = (0.0, 1.0, 2.0, 3.0, T)
J_BEFORE, J_AFTER = 1.55, 1.60


def time_grid(end, dt=DT):
    return np.arange(round(end / dt) + 1) * dt


def trace(coupling, time):
    opened, ring, _, audit = MODEL["measure"](time, coupling)
    delta = ring - opened
    squares = delta**2
    increments = (squares[:-1] + squares[1:]) * (time[1] / 2)
    integral = np.concatenate([np.zeros_like(squares[:1]), np.cumsum(increments, axis=0)])
    return opened, ring, delta, integral, audit


def bisect(gap, left, right):
    left_sign = np.signbit(gap(left))
    assert left_sign != np.signbit(gap(right))
    while right - left > 5e-8:
        middle = (left + right) / 2
        if np.signbit(gap(middle)) == left_sign:
            left = middle
        else:
            right = middle
    return float((left + right) / 2)


def crossings(values, gap):
    indices = np.flatnonzero(np.signbit(values[:-1]) != np.signbit(values[1:]))
    return [bisect(gap, COUPLINGS[i], COUPLINGS[i + 1]) for i in indices]


def window_result(window, power):
    time = time_grid(window)
    pairs = power.sum(axis=-1)
    leaders = np.asarray(REPRESENTATIVES)[np.argmax(pairs[:, REPRESENTATIVES], axis=1)]
    changes = np.flatnonzero(leaders[:-1] != leaders[1:])

    @cache
    def pair_power_at(j):
        return trace(j, time)[3][-1].sum(axis=-1)

    def gap(j, a, b):
        end = pair_power_at(j)
        return float(end[a] - end[b])

    switches = []
    for i in changes:
        a, b = leaders[i:i + 2]
        j = bisect(lambda j: gap(j, a, b), COUPLINGS[i], COUPLINGS[i + 1])
        switches.append({"from": PAIRS[a], "to": PAIRS[b], "J": j})
    equalities = crossings(pairs[:, OUTER] - pairs[:, CENTER], lambda j: gap(j, OUTER, CENTER))
    return {"T": window, "first_leader": PAIRS[leaders[0]], "leader_crossings": switches,
            "outer_center_equalities": equalities}


def yy_crossings(power, dt):
    time = time_grid(T, dt)

    def gap(j):
        end = trace(j, time)[3][-1]
        return float(end[OUTER, 1] - end[CENTER, 1])

    return crossings(power[:, OUTER, 1] - power[:, CENTER, 1], gap)


def write_window_grid(powers):
    with (OUTPUT / "windows.csv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["T", "J", "leader", *(f"pair_{PAIRS[p]}" for p in REPRESENTATIVES),
                         "outer_minus_center"])
        for t, window in enumerate(WINDOWS):
            for j, coupling in enumerate(COUPLINGS):
                pair_power = powers[j, t].sum(axis=-1) / window
                leader = max(REPRESENTATIVES, key=pair_power.__getitem__)
                writer.writerow([f"{window:.3f}", f"{coupling:.2f}", PAIRS[leader],
                                 *np.sqrt(pair_power[REPRESENTATIVES]),
                                 pair_power[OUTER] - pair_power[CENTER]])


def timing(crossing_j, dt, writer=None):
    time = time_grid(T, dt)
    band_indices = np.rint(np.asarray(TIME_BANDS) / dt).astype(int)
    output = {}
    for case, j in (("before", J_BEFORE), ("crossing", crossing_j), ("after", J_AFTER)):
        opened, ring, delta, integral, audit = trace(j, time)
        instant = delta[:, OUTER]**2 - delta[:, CENTER]**2
        cumulative = (integral[:, OUTER] - integral[:, CENTER]) / T
        points = {}
        for t in (2.25, 2.437, T):
            i = round(t / dt)
            points[f"{t:.3f}"] = {
                "open_outer_ZZ": float(opened[i, OUTER, 2]), "ring_outer_ZZ": float(ring[i, OUTER, 2]),
                "open_center_ZZ": float(opened[i, CENTER, 2]), "ring_center_ZZ": float(ring[i, CENTER, 2]),
            }
        output[case] = {"J": j, "gap_by_band": np.diff(cumulative[band_indices], axis=0).tolist(),
                        "full_gap_by_channel": cumulative[-1].tolist(), "ZZ_snapshots": points,
                        "audit": audit}
        if writer is not None:
            pair_indices = [OUTER, CENTER]
            rows = zip(time, delta[:, pair_indices, 2], instant, cumulative,
                       opened[:, pair_indices, 2], ring[:, pair_indices, 2])
            writer.writerows([j, case, f"{t:.3f}", *zz, *g, *acc, *opened_zz, *ring_zz]
                             for t, zz, g, acc, opened_zz, ring_zz in rows)
    return output


def main():
    time = time_grid(WINDOWS[-1])
    time_indices = np.rint(WINDOWS / DT).astype(int)
    powers = np.stack([trace(j, time)[3][time_indices] for j in COUPLINGS])
    fixed = {str(t): window_result(t, powers[:, np.where(WINDOWS == t)[0][0]])
             for t in WINDOW_CHECKPOINTS}
    reference = fixed[str(T)]
    crossing_j, = reference["outer_center_equalities"]
    prior = json.loads((PREVIOUS / "results.json").read_text(encoding="utf-8"))
    archive_error = abs(crossing_j - prior["leader_crossings"][-1]["J"])
    assert archive_error < 1e-7

    baseline_power = powers[:, np.where(WINDOWS == T)[0][0]]
    yy_j, = yy_crossings(baseline_power, DT)
    yy_refined, = yy_crossings(baseline_power, DT / 2)
    yy_rms = float(np.sqrt(trace(yy_j, time_grid(T))[3][-1, OUTER, 1] / T))

    with (OUTPUT / "waveforms.csv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["J", "case", "t", "delta_outer_ZZ", "delta_center_ZZ",
                         *(f"g_{axis}" for axis in CHANNELS), *(f"cumulative_{axis}" for axis in CHANNELS),
                         "open_outer_ZZ", "open_center_ZZ", "ring_outer_ZZ", "ring_center_ZZ"])
        cases = timing(crossing_j, DT, writer)
    fine_cases = timing(crossing_j, DT / 2)
    error = max(float(np.max(np.abs(np.asarray(cases[key]["gap_by_band"])
                                   - fine_cases[key]["gap_by_band"]))) for key in cases)
    assert error < 1e-7 and abs(yy_j - yy_refined) < 1e-6
    write_window_grid(powers)
    results = {
        "N": 4, "Omega": MODEL["REFERENCE"]["OMEGA"], "initial": "0000", "reference_T": T, "dt": DT,
        "J_step": 0.01, "T_step": 0.02, "channels": CHANNELS, "time_bands": TIME_BANDS,
        "cases": cases, "window_checkpoints": fixed,
        "YY_crossing": {"J": yy_j, "common_rms": yy_rms, "half_step_J": yy_refined},
        "audit": {"archive_crossing_error": archive_error, "max_time_band_change_dt_halving": error},
    }
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"window_checkpoints": fixed, "YY_crossing": results["YY_crossing"],
                      "audit": results["audit"]}, indent=2))


if __name__ == "__main__":
    main()
