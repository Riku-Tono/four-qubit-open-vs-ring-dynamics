"""Ablate individual spectral lines at the first sampled PLV crossing."""

import argparse
import csv
import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
EXPERIMENTS = OUTPUT.parent
PREVIOUS = EXPERIMENTS / "open-chain-independent-angle-errors-2026-11-10" / "diagnostics.py"
REFERENCE = EXPERIMENTS / "open-chain-phase-spectrum-robustness-2026-11-10" / "robustness.py"
MODEL = run_path(str(PREVIOUS))
GATES, PAIRS, LEFT, RIGHT = (MODEL[name] for name in ("GATES", "PAIRS", "LEFT", "RIGHT"))
CARRIER_OMEGA, TARGET_OMEGA = 0.50239249619252, 0.93063614364892
ALL_WINDOW_VARIANTS = ("full", "omit_0930636", "carrier_plus_0930636", "restore_carrier_magnitude")


def save_csv(name, rows):
    with (OUTPUT / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def variants(coefficients, nominal, carrier, target, frequencies):
    line = np.arange(len(frequencies))
    dominant = coefficients[:, carrier]
    unit = dominant / np.abs(dominant)
    restored = coefficients.copy()
    restored[:, carrier] = unit * np.abs(nominal[:, carrier])
    restored_target = coefficients.copy()
    restored_target[:, target] = nominal[:, target] / nominal[:, carrier] * dominant
    yield "full", None, coefficients
    yield "carrier_only", frequencies[carrier], coefficients * (line == carrier)
    yield "omit_0930636", frequencies[target], coefficients * (line != target)
    yield "carrier_plus_0930636", frequencies[target], coefficients * ((line == carrier) | (line == target))
    yield "align_carrier_phase", frequencies[carrier], coefficients * unit.conj()[:, None]
    yield "restore_carrier_magnitude", frequencies[carrier], restored
    yield "restore_normalized_0930636", frequencies[target], restored_target
    for mode, omega in enumerate(frequencies):
        if mode != carrier:
            yield "omit_line", omega, coefficients * (line != mode)


def phase_metrics(relative, mask, weights, pair_index):
    if not mask.any():
        return None, None, None, None
    lock = np.average(relative[mask], axis=0, weights=weights[mask])
    plv = np.abs(lock).clip(max=1)
    lags = np.angle(lock, deg=True)
    lag = float(lags[pair_index]) if plv[pair_index] > 1e-12 else None
    all_lag = float(np.max(np.abs(lags))) if plv.min() > 1e-12 else None
    return float(plv[pair_index]), lag, float(plv.min()), all_lag


def evaluate(case, time, frequencies, coefficients, nominal, carrier, target, selected_variants=None):
    pair_index = PAIRS.index(case["pair"])
    weights = np.ones(len(time))
    weights[[0, -1]] = 0.5
    kernel = np.exp(1j * time[:, None] * frequencies)
    full = kernel @ coefficients.T
    original_minimum = np.abs(full).min(axis=1)
    rows = []
    for name, omega, altered in variants(coefficients, nominal, carrier, target, frequencies):
        if selected_variants is not None and name not in selected_variants:
            continue
        signal = full if name == "full" else kernel @ altered.T
        minimum, relative = MODEL["phase_samples"](signal)
        minimum_half_amplitude = float(np.ptp(signal.real, axis=0).min() / 2)
        for gate in GATES:
            valid = minimum >= gate
            original = original_minimum >= gate
            coverage = float(np.average(valid, weights=weights))
            plv, lag, minimum_plv, all_lag = phase_metrics(relative, valid, weights, pair_index)
            original_plv, original_lag, original_minimum_plv, _ = phase_metrics(relative, original, weights, pair_index)
            passed = (coverage >= MODEL["MIN_COVERAGE"] and minimum_half_amplitude >= gate
                      and minimum_plv is not None and minimum_plv >= MODEL["MIN_PLV"]
                      and all_lag is not None and all_lag <= MODEL["MAX_LAG"])
            rows.append(dict(case=case["case"], width=case["width"], direction=case["direction"],
                             window=case["window"], pair=case["pair"], gate=gate, variant=name,
                             frequency=None if omega is None else float(omega),
                             min_envelope=float(minimum.min()), min_half_amplitude=minimum_half_amplitude,
                             coverage=coverage, plv=plv, lag_deg=lag, min_plv=minimum_plv,
                             original_mask_plv=original_plv, original_mask_lag_deg=original_lag,
                             original_mask_min_plv=original_minimum_plv, max_abs_lag_deg=all_lag,
                             passes=bool(passed)))
    return rows, full


def check_all_windows(case, rows, time_grids, frequencies, coefficients, nominal, carrier, target):
    output = [row for row in rows if row["variant"] in ALL_WINDOW_VARIANTS]
    for window, time in time_grids.items():
        if window == case["window"]:
            continue
        evaluated, _ = evaluate({**case, "window": window}, time, frequencies, coefficients, nominal, carrier, target,
                                ALL_WINDOW_VARIANTS)
        output.extend(evaluated)
    summary = {}
    for gate in GATES:
        summary[str(gate)] = {}
        for variant in ALL_WINDOW_VARIANTS:
            subset = [row for row in output if row["gate"] == gate and row["variant"] == variant]
            worst = min(subset, key=lambda row: row["min_plv"])
            summary[str(gate)][variant] = {"pass_windows": sum(row["passes"] for row in subset),
                                           "min_plv": worst["min_plv"], "min_plv_window": worst["window"],
                                           "min_coverage": min(row["coverage"] for row in subset)}
    return output, summary


def mode_rows(case, coefficients, frequencies, carrier, target):
    i, j = (int(site) for site in case["pair"])
    normalized = coefficients / coefficients[:, carrier, None]
    for k, omega in enumerate(frequencies):
        yield dict(case=case["case"], width=case["width"], direction=case["direction"],
                   window=case["window"], pair=case["pair"], frequency=float(omega),
                   is_carrier=k == carrier, is_target=k == target,
                   amplitude_left=float(abs(coefficients[i, k])), amplitude_right=float(abs(coefficients[j, k])),
                   phase_left_deg=float(np.angle(coefficients[i, k], deg=True)),
                   phase_right_deg=float(np.angle(coefficients[j, k], deg=True)),
                   normalized_mismatch=float(abs(normalized[i, k] - normalized[j, k])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", nargs="?", type=Path, default=OUTPUT / "diagnostic_cases.json")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    nominal_angles, windows, dt, _, checksum = run_path(str(REFERENCE))["sources"]()
    assert manifest.get("selection_sha256", checksum) == checksum
    model = MODEL["Model"]("open")
    frequencies = model.frequencies
    carrier = int(np.argmin(abs(frequencies - CARRIER_OMEGA)))
    target = int(np.argmin(abs(frequencies - TARGET_OMEGA)))
    assert abs(frequencies[carrier] - CARRIER_OMEGA) < 1e-9 and abs(frequencies[target] - TARGET_OMEGA) < 1e-9
    nominal_state = MODEL["product_state"](nominal_angles)
    _, nominal = MODEL["grouped_coefficients"](model, nominal_state)
    rows, modes, cases, all_window_rows = [], [], [], []
    time_grids = {window: MODEL["time_grid"](bounds, dt / 2) for window, bounds in windows.items()}
    crossing_width = max(case["width"] for case in manifest["cases"])
    audit = {"dt": dt / 2, "coarse_dt": dt, "reconstruction_error": 0.0,
             "max_grid_plv_change": 0.0, "max_phase_alignment_plv_change": 0.0, "pass_flips": []}
    for case in manifest["cases"]:
        bounds = windows[case["window"]]
        state = MODEL["product_state"](nominal_angles + np.asarray(case["offsets_rad"]).reshape(4, 2))
        dc, coefficients = MODEL["grouped_coefficients"](model, state)
        time = time_grids[case["window"]]
        fine, full = evaluate(case, time, frequencies, coefficients, nominal, carrier, target)
        coarse, _ = evaluate(case, MODEL["time_grid"](bounds, dt), frequencies, coefficients, nominal, carrier, target)
        direct = model.direct(state, time[::317])
        audit["reconstruction_error"] = max(audit["reconstruction_error"], float(np.max(abs(dc + full[::317].real - direct))))
        for a, b in zip(fine, coarse):
            if a["plv"] is not None and b["plv"] is not None:
                audit["max_grid_plv_change"] = max(audit["max_grid_plv_change"], abs(a["plv"] - b["plv"]))
            if a["passes"] != b["passes"]:
                audit["pass_flips"].append([case["case"], a["gate"], a["variant"], a["frequency"]])
        main = {row["variant"]: row for row in fine if row["variant"] != "omit_line" and row["gate"] == MODEL["MAIN_GATE"]}
        alignment_delta = abs(main["align_carrier_phase"]["original_mask_plv"] - main["full"]["original_mask_plv"])
        audit["max_phase_alignment_plv_change"] = max(audit["max_phase_alignment_plv_change"], alignment_delta)
        omissions = [row for row in fine if row["variant"] == "omit_line" and row["gate"] == MODEL["MAIN_GATE"]]
        omissions.sort(key=lambda row: row["original_mask_plv"], reverse=True)
        summary = {**case, "window_bounds": bounds, "carrier": MODEL["spectrum_summary"](frequencies, coefficients, carrier),
                   "main_gate_variants": main, "ranked_mode_omissions": omissions}
        if case["width"] == crossing_width:
            extended, summary["all_windows"] = check_all_windows(case, fine, time_grids, frequencies, coefficients,
                                                                 nominal, carrier, target)
            all_window_rows.extend(extended)
        cases.append(summary)
        rows.extend(fine)
        modes.extend(mode_rows(case, coefficients, frequencies, carrier, target))
    assert audit["reconstruction_error"] < 1e-10 and audit["max_phase_alignment_plv_change"] < 1e-12
    save_csv("mechanism.csv", rows)
    save_csv("mechanism_modes.csv", modes)
    save_csv("mechanism_windows.csv", all_window_rows)
    result = {"selection_sha256": checksum, "carrier_omega": float(frequencies[carrier]),
              "target_omega": float(frequencies[target]), "cases": cases, "audit": audit,
              "original_mask_definition": "All four envelopes of the unmodified signal at or above the gate."}
    (OUTPUT / "mechanism.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
