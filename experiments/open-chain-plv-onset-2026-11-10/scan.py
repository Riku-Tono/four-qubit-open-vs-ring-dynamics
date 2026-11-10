"""Resolve the first sampled PLV failure on the previous fixed error directions."""

import csv
import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
PREVIOUS = OUTPUT.parent / "open-chain-independent-angle-errors-2026-11-10"
OLD = run_path(str(PREVIOUS / "sweep.py"))
GATES, MAIN_GATE, MIN_PLV, MAX_LAG, LEFT, RIGHT, PAIRS = (
    OLD[key] for key in ("GATES", "MAIN_GATE", "MIN_PLV", "MAX_LAG", "LEFT", "RIGHT", "PAIRS")
)
WIDTHS = np.arange(75, 101) / 1000
DT, FINE_DT, GUARD = 0.02, 0.001, 0.001
BATCH = 48


def write_json(name, data):
    (OUTPUT / name).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(name, rows):
    with (OUTPUT / name).open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def spectral_checks(model, coefficients):
    line = int(np.argmin(np.abs(model.frequencies - OLD["TRACKED_LINE"])))
    assert abs(model.frequencies[line] - OLD["TRACKED_LINE"]) < 1e-9
    absolute = np.abs(coefficients)
    envelope_lower = 2 * absolute[..., line] - absolute.sum(axis=-1) - 2e-10
    assert envelope_lower.min() > max(GATES)

    # Factor out the shared carrier before bounding derivatives of the unit phasors.
    detuning = np.abs(model.frequencies - model.frequencies[line])
    phase_prime = (absolute @ detuning) / envelope_lower
    phase_second = (absolute @ detuning**2) / envelope_lower + phase_prime**2
    relative_second = (phase_second[..., LEFT] + phase_second[..., RIGHT]
                       + (phase_prime[..., LEFT] + phase_prime[..., RIGHT])**2)
    return line, envelope_lower, float(np.max(relative_second))


def evaluate(coefficients, kernel):
    count, sites, modes = coefficients.shape
    signal = (coefficients.reshape(-1, modes) @ kernel).reshape(count, sites, -1)
    env = np.abs(signal)
    minimum = env.min(axis=-1)
    amplitude = np.ptp(signal.real, axis=-1) / 2
    unit = signal / env
    unit[..., [0, -1]] *= np.sqrt(0.5)
    locking = (unit @ unit.swapaxes(-1, -2).conj())[:, LEFT, RIGHT] / (unit.shape[-1] - 1)
    return locking, minimum, amplitude


def evaluate_chunks(coefficients, kernel, batch):
    chunks = [evaluate(coefficients[i:i + batch], kernel) for i in range(0, len(coefficients), batch)]
    return tuple(np.concatenate(parts) for parts in zip(*chunks))


def scan(model, coefficients, windows, derivative_bound):
    coarse_kernel = [model.prepare(OLD["time_grid"](bounds, DT)).kernel for bounds in windows.values()]
    fine_kernel = [model.prepare(OLD["time_grid"](bounds, FINE_DT)).kernel for bounds in windows.values()]
    shape = (len(WIDTHS), coefficients.shape[1], len(windows))
    locking = np.empty((*shape, len(PAIRS)), dtype=complex)
    envelope, amplitude = (np.empty((*shape, 4)) for _ in range(2))
    refined = np.zeros(shape, dtype=bool)
    audit = {"refined_case_windows": 0, "refined_status_changes": 0, "max_locking_change": 0.0}
    error_bound = (DT**2 + FINE_DT**2) * derivative_bound / 12 + 1e-12
    assert error_bound < GUARD
    for w, width in enumerate(WIDTHS):
        for v, (coarse, fine) in enumerate(zip(coarse_kernel, fine_kernel)):
            z, env, amp = evaluate_chunks(coefficients[w], coarse, BATCH)
            near = np.any(np.abs(np.abs(z) - MIN_PLV) <= GUARD, axis=-1)
            if np.any(near):
                precise, env[near], amp[near] = evaluate_chunks(coefficients[w, near], fine, BATCH // 4)
                audit["refined_case_windows"] += int(near.sum())
                audit["refined_status_changes"] += int(np.count_nonzero((np.abs(z[near]) < MIN_PLV) !=
                                                                          (np.abs(precise) < MIN_PLV)))
                audit["max_locking_change"] = max(audit["max_locking_change"], float(np.max(np.abs(precise - z[near]))))
                z[near] = precise
            locking[w, :, v], envelope[w, :, v], amplitude[w, :, v], refined[w, :, v] = z, env, amp, near
        print(f"Evaluated width {width:.3f}", flush=True)
    assert audit["max_locking_change"] < error_bound
    assert min(envelope.min(), amplitude.min()) > max(GATES)
    lag_error = np.rad2deg(np.arcsin(error_bound / np.abs(locking).min()))
    assert np.max(np.abs(np.angle(locking, deg=True))) + lag_error < MAX_LAG
    return locking, envelope, amplitude, refined, audit


def audit_endpoints(pass_windows, labels):
    match = {(float(row["width"]), row["direction"]): int(row["pass_windows"])
             for row in csv.DictReader((PREVIOUS / "case_summary.csv").open(encoding="utf-8", newline=""))
             if float(row["width"]) in (WIDTHS[0], WIDTHS[-1]) and float(row["gate"]) == MAIN_GATE}
    for index in (0, len(WIDTHS) - 1):
        assert all(pass_windows[index, j] == match[WIDTHS[index], direction] for j, (direction, _) in enumerate(labels))
    return len(match)


def main():
    nominal, windows, _, _, checksum = OLD["REFERENCE"]["sources"]()
    model = OLD["Model"]("open")
    labels, directions = OLD["samples"]()
    angles = nominal + WIDTHS[:, None, None, None] * directions
    states = OLD["product_state"](angles).reshape(-1, 16)
    coefficients = OLD["spectrum"](model, states).reshape(len(WIDTHS), len(labels), 4, -1)
    line, envelope_lower, derivative_bound = spectral_checks(model, coefficients)
    print(f"Spectral envelope lower bound {envelope_lower.min():.9f}", flush=True)
    locking, envelope, amplitude, refined, audit = scan(model, coefficients, windows, derivative_bound)
    plv = np.minimum(1, np.abs(locking))
    lag = np.angle(locking, deg=True)
    failure = plv < MIN_PLV
    pass_windows = (~failure.any(axis=-1)).sum(axis=-1)
    audit["endpoint_cases_matched"] = audit_endpoints(pass_windows, labels)
    names = list(windows)
    first = int(np.flatnonzero(failure.any(axis=(1, 2, 3)))[0])
    assert first > 0

    def sample(w, d, v, p, gate):
        return dict(width=float(WIDTHS[w]), gate=gate, direction=labels[d][0], family=labels[d][1],
                    window=names[v], pair=PAIRS[p], plv=float(plv[w, d, v, p]), lag_deg=float(lag[w, d, v, p]),
                    coverage=1.0, min_envelope=float(envelope[w, d, v].min()), passes=not bool(failure[w, d, v, p]),
                    target_dt=FINE_DT if refined[w, d, v] else DT)

    witness_rows, trace, cases = [], [], []
    for d, v, p in np.argwhere(failure[first]):
        for w in range(len(WIDTHS)):
            for gate in GATES:
                trace.append(sample(w, d, v, p, gate))
        for w in (first - 1, first):
            row = sample(w, d, v, p, MAIN_GATE)
            target = model.prepare(OLD["time_grid"](windows[names[v]], FINE_DT / 4)).kernel
            accurate, env, amp = evaluate(coefficients[w, d:d + 1], target)
            row["plv_fine"] = float(abs(accurate[0, p]))
            row["lag_deg_fine"] = float(np.angle(accurate[0, p], deg=True))
            row["pass_windows"] = int(pass_windows[w, d])
            assert (row["plv_fine"] >= MIN_PLV) == row["passes"]
            assert env.min() > max(GATES) and amp.min() > max(GATES)
            witness_rows.append(row)
            cases.append(dict(case=f"w{WIDTHS[w]:.3f}-{labels[d][0]}", width=float(WIDTHS[w]), direction=labels[d][0],
                              window=names[v], pair=PAIRS[p], gate=MAIN_GATE,
                              offsets_rad=(WIDTHS[w] * directions[d]).tolist()))

    width_rows, direction_rows = [], []
    for w, width in enumerate(WIDTHS):
        d, v, p = np.unravel_index(np.argmin(plv[w]), plv[w].shape)
        extrema = dict(min_plv=float(plv[w, d, v, p]), min_plv_direction=labels[d][0], min_plv_pair=PAIRS[p],
                       min_plv_window=names[v], max_abs_lag=float(np.abs(lag[w]).max()),
                       min_envelope=float(envelope[w].min()), min_half_amplitude=float(amplitude[w].min()), min_coverage=1.0)
        for gate in GATES:
            width_rows.append(dict(width=float(width), gate=gate, conditions=len(labels),
                                   all_pass_cases=int(np.count_nonzero(pass_windows[w] == len(names))),
                                   failing_directions=int(np.count_nonzero(pass_windows[w] < len(names))),
                                   failing_direction_windows=int(np.count_nonzero(failure[w].any(axis=-1))),
                                   failing_pair_windows=int(np.count_nonzero(failure[w])), **extrema))
        for d, (direction, family) in enumerate(labels):
            v, p = np.unravel_index(np.argmin(plv[w, d]), plv[w, d].shape)
            direction_rows.append(dict(width=float(width), direction=direction, family=family,
                                       max_abs_error=float(width * np.abs(directions[d]).max()),
                                       **dict(zip(OLD["AXES"], (width * directions[d]).ravel().tolist())),
                                       pass_windows_all_gates=int(pass_windows[w, d]),
                                       min_plv=float(plv[w, d, v, p]), min_plv_window=names[v], min_plv_pair=PAIRS[p],
                                       max_abs_lag=float(np.abs(lag[w, d]).max()), min_envelope=float(envelope[w, d].min()),
                                       min_half_amplitude=float(amplitude[w, d].min()), min_coverage=1.0))
    violations = [sample(w, d, v, p, MAIN_GATE) for w, d, v, p in np.argwhere(failure)]
    for filename, rows in (("width_summary.csv", width_rows), ("direction_summary.csv", direction_rows),
                           ("witness_trace.csv", trace), ("first_witnesses.csv", witness_rows), ("violations.csv", violations)):
        write_csv(filename, rows)
    np.savez_compressed(OUTPUT / "all_windows.npz", widths=WIDTHS, directions=np.asarray([d for d, _ in labels]),
                        offsets_at_unit_width=directions, windows=np.asarray(names), window_bounds=np.asarray(list(windows.values())),
                        gates=GATES, pairs=PAIRS, locking=locking, site_min_envelope=envelope, half_amplitude=amplitude,
                        fine_grid=refined, envelope_lower=envelope_lower)

    provenance = {"selection_sha256": checksum, "windows": windows, "gates": GATES,
                  "carrier_frequency": float(model.frequencies[line]), "nominal_angles_rad": nominal.tolist()}
    write_json("diagnostic_cases.json", {**provenance, "cases": cases})
    write_json("scan.json", {**provenance, "directions": len(labels), "widths": WIDTHS.tolist(),
                            "first_failure_width": float(WIDTHS[first]), "previous_width": float(WIDTHS[first - 1]),
                            "minimum_envelope_lower": float(envelope_lower.min()), "primary_dt": DT, "fine_dt": FINE_DT,
                            "primary_locking_error_bound": DT**2 * derivative_bound / 12,
                            "fine_locking_error_bound": FINE_DT**2 * derivative_bound / 12,
                            "first_witnesses": witness_rows, "audit": audit})


if __name__ == "__main__":
    main()
