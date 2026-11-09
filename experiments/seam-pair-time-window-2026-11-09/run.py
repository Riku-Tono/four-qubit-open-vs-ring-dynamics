"""Decompose the N=4 to N=6 change in S at zero endpoint field."""

import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
PREVIOUS = OUTPUT.parent / "endpoint-field-six-qubit-2026-11-09"
MODEL = run_path(str(PREVIOUS / "run.py"))
BASELINE, DT = MODEL["END"], MODEL["DT"]
LINK, REFERENCE = MODEL["LINK"], MODEL["REFERENCE"]
OFFSETS = (-1, -0.5, 0, 0.5, 1)
WINDOWS = [round(BASELINE + offset, 3) for offset in OFFSETS]
MATCH_N6 = ("01", "04", "05", "14", "15", "45")


def integrate(power, dt):
    increments = (power[:-1] + power[1:]) * (dt / 2)
    return np.concatenate([np.zeros_like(power[:1]), np.cumsum(increments, axis=0)])


def measure(size, time):
    model = MODEL["Chain"](size)
    states = [model.evolve(time, link, 0) for link in (LINK, REFERENCE)]
    norm_error = float(np.max(np.abs(np.linalg.norm(states, axis=-1) - 1)))
    linked, reference = map(model.connected, states)
    # Arrays are indexed by time, state, pair, then XX/YY/ZZ.
    power = np.swapaxes((linked - reference)**2, -1, -2)
    reflected = [model.pairs.index((size - 1 - j, size - 1 - i)) for i, j in model.pairs]
    reflection_error = float(np.max(np.abs(power[:, 0] - power[:, 1, reflected])))
    assert norm_error < 1e-12 and reflection_error < 1e-12
    coarse = integrate(power[::2, 0], DT)
    fine = integrate(power[:, 0], DT / 2)[::2]
    assert np.all(np.diff(coarse, axis=0) >= 0)
    labels = [f"{i}{j}" for i, j in model.pairs]
    return labels, coarse, fine, {"max_norm_error": norm_error, "max_reflected_pair_power_error": reflection_error}


def brackets(time, values):
    changes = np.flatnonzero(np.sign(values[1:]) != np.sign(values[:-1]))
    return [[round(float(time[i]), 3), round(float(time[i + 1]), 3)] for i in changes]


def ring_distance_two(labels, size, integrals):
    distances = [int(pair[1]) - int(pair[0]) for pair in labels]
    return integrals[:, [min(d, size - d) == 2 for d in distances]].sum(axis=1)


def main():
    fine_time = np.arange(round(WINDOWS[-1] / (DT / 2)) + 1) * (DT / 2)
    time = fine_time[::2]
    labels4, c4, f4, audit4 = measure(4, fine_time)
    labels6, c6, f6, audit6 = measure(6, fine_time)
    archived = json.loads((PREVIOUS / "results.json").read_text(encoding="utf-8"))
    anchor = round(BASELINE / DT)
    archive_error = 0.0
    for size, coarse, fine, labels, seam in ((4, c4, f4, labels4, "03"), (6, c6, f6, labels6, "05")):
        for grid, result in (("primary", coarse), ("refined", fine)):
            old = archived[f"N{size}"][grid]["direct_derivative_at_zero"]
            archive_error = max(archive_error, abs(float(result[anchor].sum()) - old["denominator_A_B"][0]),
                                abs(float(result[anchor, labels.index(seam)].sum()) - old["numerator_A_B"][0]))
    assert archive_error < 1e-12

    i4, i6 = c4.sum(axis=-1), c6.sum(axis=-1)
    matched_indices = [labels6.index(label) for label in MATCH_N6]
    extra_indices = [i for i in range(len(labels6)) if i not in matched_indices]
    drops = i4 - i6[:, matched_indices]
    extras = i6[:, extra_indices].sum(axis=1)
    s4, s6 = i4.sum(axis=1), i6.sum(axis=1)
    net = s4 - s6
    distance_two4 = ring_distance_two(labels4, 4, i4)
    distance_two6 = ring_distance_two(labels6, 6, i6)
    closure = float(np.max(np.abs(net - (drops.sum(axis=1) - extras))))
    assert closure < 1e-12

    def row(index):
        return {"T": round(float(time[index]), 3), "S4": float(s4[index]), "S6": float(s6[index]),
                "S4_minus_S6": float(net[index]), "matched_drops": drops[index].tolist(),
                "central_9_total": float(extras[index]),
                "cross_seam_drop_fraction": float((drops[index, 1] + drops[index, 4]) / net[index])}

    checkpoints = []
    for end in WINDOWS:
        index = round(end / DT)
        checkpoints.append({**row(index), "N4_pair_channels": dict(zip(labels4, c4[index].tolist())),
                            "N6_pair_channels": dict(zip(labels6, c6[index].tolist())),
                            "ring_distance_two_total_N4_N6": [float(distance_two4[index]), float(distance_two6[index])]})

    first, last = (round(end / DT) for end in (WINDOWS[0], WINDOWS[-1]))
    interval = slice(first, last + 1)
    in_window = time[interval]
    fractions = (drops[interval, 1] + drops[interval, 4]) / net[interval]
    distance_two_fractions = (distance_two4[interval] - distance_two6[interval]) / net[interval]
    top_two = np.sort(np.argsort(drops[interval], axis=1)[:, -2:], axis=1)
    summary = {
        "S6_over_S4_min_max": [float(np.min(s6[interval] / s4[interval])), float(np.max(s6[interval] / s4[interval]))],
        "cross_seam_fraction_min_max": [float(np.min(fractions)), float(np.max(fractions))],
        "ring_distance_two_fraction_min_max": [float(np.min(distance_two_fractions)),
                                               float(np.max(distance_two_fractions))],
        "both_cross_seam_pairs_are_largest_drops": bool(np.all(top_two == [1, 4])),
        "largest_drop_switch_brackets": brackets(in_window, drops[interval, 1] - drops[interval, 4]),
        "seam_drop_sign_change_brackets": brackets(in_window, drops[interval, 2]),
    }
    numerical = {"N4": audit4, "N6": audit6, "archive_baseline_max_error": archive_error,
                 "max_decomposition_closure_error": closure,
                 "max_pair_integral_change_dt_halving": float(max(np.max(np.abs((a - b).sum(axis=-1)))
                                                                   for a, b in ((c4, f4), (c6, f6)))),
                 "max_S_change_dt_halving": float(max(np.max(np.abs((a - b).sum(axis=(1, 2))))
                                                       for a, b in ((c4, f4), (c6, f6))))}
    results = {"Omega": MODEL["OMEGA"], "J": MODEL["J"], "lambda": LINK, "reference_lambda": REFERENCE,
               "epsilon": 0, "initial": {"N4": "1000", "N6": "100000"}, "baseline_T": BASELINE,
               "dt": DT, "audit_dt": DT / 2, "channels": ["XX", "YY", "ZZ"],
               "matched_pairs_N4_to_N6": list(zip(labels4, MATCH_N6)),
               "central_N6_pairs": [labels6[i] for i in extra_indices],
               "checkpoints": checkpoints, "scan": [row(i) for i in range(first, last + 1, 10)],
               "window_summary": summary, "audit": numerical}
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"baseline": row(anchor), "window_summary": summary, "audit": numerical}, indent=2))


if __name__ == "__main__":
    main()
