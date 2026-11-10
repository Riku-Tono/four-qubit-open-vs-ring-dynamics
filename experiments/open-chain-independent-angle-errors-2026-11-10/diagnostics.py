"""Dissect frozen boundary cases into carrier, residual and amplitude contributions."""

import argparse
import csv
import json
from pathlib import Path
from runpy import run_path
import sys

import numpy as np


OUTPUT = Path(__file__).resolve().parent
EXPERIMENTS = OUTPUT.parent
sys.path.insert(0, str(EXPERIMENTS / "four-site-in-phase-2026-11-10"))
from metrics import GATES, LEFT, MAIN_GATE, MAX_LAG, MIN_COVERAGE, MIN_PLV, PAIRS, RIGHT, Model, describe, phase_samples, product_state, time_grid


REFERENCE = EXPERIMENTS / "open-chain-phase-spectrum-robustness-2026-11-10" / "robustness.py"
CARRIER = 0.50239249619252
VARIANTS = {
    "full": "Original physical state and all its frequencies.",
    "carrier_only": "Retain the 0.502392496 line at its original site amplitudes and phases.",
    "align_carrier_phase": "Rotate each site's entire analytic signal so the carrier phases coincide. This keeps envelopes and normalized residuals exactly fixed.",
    "restore_carrier_magnitude": "Restore nominal carrier magnitudes, keeping the perturbed carrier phases and unnormalized residual coefficients.",
    "restore_site_scale": "Rescale each entire analytic signal to the nominal carrier magnitude, keeping all instantaneous phases and normalized residuals.",
    "restore_normalized_residual": "Use the nominal residual coefficients relative to the carrier, keeping perturbed carrier phases and magnitudes.",
    "omit_largest_residual_mismatch": "Omit the non-carrier line with the largest normalized coefficient difference across six pairs.",
    "carrier_and_largest_residual_mismatch": "Keep the carrier and the non-carrier line with the largest normalized coefficient difference.",
}


def write_csv(name, rows):
    with (OUTPUT / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def grouped_coefficients(model, state):
    """Return the coefficients in R_i = dc_i + Re sum c_i(w) exp(iwt)."""
    energy_state = state @ model.vectors
    products = energy_state[model.a].conj() * energy_state[model.b]
    coefficients = np.add.reduceat(products * model.observed, model.group_starts, axis=-1)
    dc = np.einsum("a,abi,b->i", energy_state.conj(), model.dc, energy_state).real
    return dc, coefficients


def normalized_mismatch(coefficients, carrier):
    normalized = coefficients / coefficients[:, carrier, None]
    mismatch = np.abs(normalized[LEFT] - normalized[RIGHT])
    mismatch[:, carrier] = 0
    return mismatch


def spectrum_summary(frequency, coefficients, carrier):
    magnitude = np.abs(coefficients)
    dominant = coefficients[:, carrier]
    carrier_magnitude = magnitude[:, carrier]
    power = np.sum(magnitude**2, axis=1)
    residual = magnitude.sum(axis=1) - carrier_magnitude
    relative = np.divide(residual, carrier_magnitude, out=np.full(4, np.nan), where=carrier_magnitude > 0)
    carrier_phase = np.angle(dominant, deg=True)
    lag = np.angle(dominant[LEFT] * dominant[RIGHT].conj(), deg=True)
    pair_angles = {pair: float(value) if carrier_magnitude[i] > 0 and carrier_magnitude[j] > 0 else None
                   for pair, value, i, j in zip(PAIRS, lag, LEFT, RIGHT)}
    top_mismatch = []
    if np.all(carrier_magnitude >= 1e-14):
        mismatch = normalized_mismatch(coefficients, carrier)
        mode_order = np.argsort(mismatch.max(axis=0))[::-1][:5]
        for index in mode_order:
            pair_index = np.argmax(mismatch[:, index])
            top_mismatch.append({"omega": float(frequency[index]), "pair": PAIRS[pair_index],
                                 "abs_normalized_difference": float(mismatch[pair_index, index])})
    return {
        "carrier_amplitude": carrier_magnitude.tolist(),
        "carrier_phase_deg": [float(value) if radius > 0 else None for value, radius in zip(carrier_phase, carrier_magnitude)],
        "carrier_power_ratio": (carrier_magnitude**2 / np.where(power > 0, power, 1)).tolist(),
        "carrier_pair_lag_deg": pair_angles,
        "max_carrier_pair_lag_deg": max((abs(value) for value in pair_angles.values() if value is not None), default=None),
        "residual_l1": residual.tolist(),
        "residual_l1_over_carrier": [float(value) if np.isfinite(value) else None for value in relative],
        "envelope_bound_triangle": np.maximum(0, carrier_magnitude - residual).tolist(),
        "largest_line_is_carrier": (np.argmax(magnitude, axis=1) == carrier).tolist(),
        "largest_normalized_residual_mismatches": top_mismatch,
    }


def variants(coefficients, nominal, carrier):
    dominant = coefficients[:, carrier]
    reference = nominal[:, carrier]
    magnitude = np.abs(dominant)
    if np.any(magnitude < 1e-14):
        # A missing carrier has undefined phase and relative residual.
        return {"full": coefficients, "carrier_only": coefficients * (np.arange(coefficients.shape[1]) == carrier)}
    unit = dominant / magnitude
    scale = np.abs(reference) / magnitude
    modes = np.arange(coefficients.shape[1])
    residual_modes = modes[modes != carrier]
    strongest = residual_modes[np.argmax(normalized_mismatch(coefficients, carrier)[:, residual_modes].max(axis=0))]
    restored_carrier = coefficients.copy()
    restored_carrier[:, carrier] = unit * np.abs(reference)
    return {
        "full": coefficients,
        "carrier_only": coefficients * (np.arange(coefficients.shape[1]) == carrier),
        "align_carrier_phase": coefficients * unit.conj()[:, None],
        "restore_carrier_magnitude": restored_carrier,
        "restore_site_scale": coefficients * scale[:, None],
        "restore_normalized_residual": nominal * (dominant / reference)[:, None],
        "omit_largest_residual_mismatch": coefficients * (modes != strongest),
        "carrier_and_largest_residual_mismatch": coefficients * ((modes == carrier) | (modes == strongest)),
    }


def state_spectrum(case, frequencies, coefficients, carrier):
    rows = []
    magnitudes = np.abs(coefficients)
    powers = np.sum(magnitudes**2, axis=-1)
    for index, omega in enumerate(frequencies):
        for site, value in enumerate(coefficients[:, index]):
            magnitude = magnitudes[site, index]
            rows.append({"case": case, "mode": f"m{index + 1:02}", "omega": float(omega), "is_carrier": index == carrier,
                         "site": site, "real": float(value.real), "imag": float(value.imag),
                         "amplitude": float(magnitude), "phase_deg": float(np.angle(value, deg=True)) if magnitude else None,
                         "power_ratio": float(magnitude**2 / powers[site]) if powers[site] else None})
    return rows


def phase_on_mask(relative, mask, weights):
    if not mask.any():
        return None, None
    average = np.average(relative[mask], axis=0, weights=weights[mask])
    plv = float(np.min(np.abs(average)))
    lag = float(np.max(np.abs(np.angle(average, deg=True)))) if plv > 1e-12 else None
    return lag, min(1, plv)


def describe_variant(case, variant, window, bounds, gate, occupation, signal, ungated, envelope, matched):
    metrics = describe(occupation, signal, gate)
    minimum_amplitude = min(metrics["amplitude"])
    lag = metrics["max_abs_lag"]
    plv = metrics["min_plv"]
    flags = {
        "half_range_fail": minimum_amplitude < gate,
        "coverage_fail": metrics["coverage"] < MIN_COVERAGE,
        "lag_fail": lag is not None and lag > MAX_LAG,
        "plv_fail": plv is not None and plv < MIN_PLV,
        "phase_undefined": lag is None or plv is None,
    }
    assert metrics["passes"] == (not any(flags.values()))
    lag_pair = max((p for p in PAIRS if metrics["pairs"][p]["lag_deg"] is not None),
                   key=lambda p: abs(metrics["pairs"][p]["lag_deg"]), default=None)
    plv_pair = min((p for p in PAIRS if metrics["pairs"][p]["plv"] is not None),
                   key=lambda p: metrics["pairs"][p]["plv"], default=None)
    return {
        "case": case, "variant": variant, "gate": gate, "window": window,
        "start": bounds[0], "end": bounds[1], "passes": metrics["passes"], **flags,
        "coverage": metrics["coverage"], "min_half_range": minimum_amplitude,
        "min_envelope": float(np.min(envelope)), "max_abs_lag_deg": lag, "min_plv": plv,
        "max_abs_lag_pair": lag_pair, "min_plv_pair": plv_pair,
        "ungated_max_abs_lag_deg": ungated["max_abs_lag"], "ungated_min_plv": ungated["min_plv"],
        "original_mask_max_abs_lag_deg": matched[0], "original_mask_min_plv": matched[1],
        **{f"amplitude{i}": value for i, value in enumerate(metrics["amplitude"])},
        **{f"min_envelope{i}": float(value) for i, value in enumerate(np.min(envelope, axis=0))},
        **{f"{m}_{p}": metrics["pairs"][p][m] for m in ("lag_deg", "plv") for p in PAIRS},
    }


def margin(row):
    lag = row["max_abs_lag_deg"]
    plv = row["min_plv"]
    return min(row["min_half_range"] / row["gate"] - 1,
               (row["coverage"] - MIN_COVERAGE) / (1 - MIN_COVERAGE),
               (MAX_LAG - lag) / MAX_LAG if lag is not None else -2,
               (plv - MIN_PLV) / (1 - MIN_PLV) if plv is not None else -2)


def summarize(rows):
    result = {"pass_windows": sum(row["passes"] for row in rows), "worst_window": min(rows, key=margin)["window"]}
    for key in ("half_range_fail", "coverage_fail", "lag_fail", "plv_fail", "phase_undefined"):
        result[key + "_windows"] = [row["window"] for row in rows if row[key]]
    fields = {
        "coverage": min, "min_half_range": min, "min_envelope": min, "max_abs_lag_deg": max,
        "min_plv": min, "ungated_max_abs_lag_deg": max, "ungated_min_plv": min,
        "original_mask_max_abs_lag_deg": max, "original_mask_min_plv": min,
    }
    for field, select in fields.items():
        defined = [row for row in rows if row[field] is not None]
        worst = select(defined, key=lambda row: row[field]) if defined else None
        result[field] = worst[field] if worst is not None else None
        result[field + "_window"] = worst["window"] if worst is not None else None
    return result


def waveform_rows(case, window, bounds, frequencies, dc, coefficients, carrier):
    # Output has plotting resolution; all window metrics use the original fine grid.
    t = time_grid(bounds, 0.015)
    signal = np.exp(1j * t[:, None] * frequencies) @ coefficients.T
    dominant = np.exp(1j * t[:, None] * frequencies[carrier]) * coefficients[:, carrier]
    env = np.abs(signal)
    minimum, relative = phase_samples(signal)
    lags = np.angle(relative, deg=True)
    columns = ("case", "window", "t", *(f"R{i}" for i in range(4)), *(f"envelope{i}" for i in range(4)),
               "min_envelope", *(f"phase_{p}" for p in PAIRS), *(f"carrier_R{i}" for i in range(4)))
    return [dict(zip(columns, (case, window, float(t[i]), *r, *env[i], minimum[i], *lags[i], *dominant_r)))
            for i, (r, dominant_r) in enumerate(zip(dc + signal.real, dc + dominant.real))]


def check_finer(model, summaries, windows, nominal, carrier, dt):
    """Halve the grid spacing in each representative plot window."""
    audit = {"dt": dt / 2, "compared": 0, "pass_flips": [], "max_changes": {}}
    for alias, case in summaries.items():
        if alias == "nominal":
            continue
        time = time_grid(windows[case["plot_window"]], dt / 2)
        kernel = np.exp(1j * time[:, None] * model.frequencies)
        dc, coefficients = grouped_coefficients(model, product_state(np.asarray(case["angles_rad"])))
        for variant, altered in variants(coefficients, nominal, carrier).items():
            signal = kernel @ altered.T
            fine = describe(dc + signal.real, signal, case["gate"])
            old = case["plot_window_metrics"][variant]
            audit["compared"] += 1
            if fine["passes"] != old["passes"]:
                audit["pass_flips"].append([alias, case["plot_window"], case["gate"], variant])
            fields = {"max_abs_lag_deg": fine["max_abs_lag"], "min_plv": fine["min_plv"],
                      "coverage": fine["coverage"], "min_half_range": min(fine["amplitude"]),
                      "min_envelope": float(np.abs(signal).min())}
            for key, value in fields.items():
                if old[key] is not None and value is not None:
                    audit["max_changes"][key] = max(audit["max_changes"].get(key, 0), abs(old[key] - value))
    return audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", nargs="?", type=Path, default=OUTPUT / "diagnostic_cases.json")
    parser.add_argument("--dt", type=float, help="Override the prior grid spacing for an independent check.")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    nominal_angles, windows, previous_dt, previous, checksum = run_path(str(REFERENCE))["sources"]()
    dt = previous_dt if args.dt is None else args.dt
    assert dt > 0
    assert manifest.get("selection_sha256", checksum) == checksum
    cases = {"nominal": {"offsets_rad": np.zeros((4, 2)).tolist(), "gate": MAIN_GATE}, **manifest["cases"]}
    assert not np.any(cases["nominal"]["offsets_rad"])
    model = Model("open")
    frequency = model.frequencies
    carrier = int(np.argmin(np.abs(frequency - CARRIER)))
    assert abs(frequency[carrier] - CARRIER) < 1e-9
    nominal_dc, nominal_coeff = grouped_coefficients(model, product_state(nominal_angles))
    prepared = {}
    for window, bounds in windows.items():
        time = time_grid(bounds, dt)
        weights = np.ones(len(time))
        weights[[0, -1]] = 0.5
        prepared[window] = time, np.exp(1j * time[:, None] * frequency), weights
    spectra, waveform, window_rows, summaries = [], [], [], {}
    reconstruction_error = 0.0
    phase_alignment_error = 0.0
    phase_alignment_plv_error = 0.0
    nominal_lag_error = 0.0
    nominal_plv_error = 0.0
    requested_windows = set()
    for case, spec in cases.items():
        offsets = np.asarray(spec["offsets_rad"])
        assert offsets.shape == (4, 2)
        angles = nominal_angles + offsets
        state = product_state(angles)
        dc, coefficients = grouped_coefficients(model, state)
        spectra.extend(state_spectrum(case, frequency, coefficients, carrier))
        signals = variants(coefficients, nominal_coeff, carrier)
        results = {variant: {gate: [] for gate in GATES} for variant in signals}
        for window, bounds in windows.items():
            time, kernel, weights = prepared[window]
            full = kernel @ coefficients.T
            original_envelope = np.abs(full)
            original_minimum = original_envelope.min(axis=1)
            original_masks = {gate: original_minimum >= gate for gate in GATES}
            reference = model.direct(state, time[::317])
            reconstruction_error = max(reconstruction_error, float(np.max(np.abs(dc + full[::317].real - reference))))
            for variant, variant_coefficients in signals.items():
                signal = full if variant == "full" else kernel @ variant_coefficients.T
                env = np.abs(signal)
                if variant == "align_carrier_phase":
                    phase_alignment_error = max(phase_alignment_error, float(np.max(np.abs(env - original_envelope))))
                occupation = dc + signal.real
                ungated = describe(occupation, signal, 1e-14)
                defined, relative = phase_samples(signal)
                for gate in GATES:
                    mask = original_masks[gate]
                    matched = phase_on_mask(relative, mask, weights) if np.all(defined[mask] > 1e-14) else (None, None)
                    row = describe_variant(case, variant, window, bounds, gate, occupation, signal, ungated, env, matched)
                    window_rows.append(row)
                    results[variant][gate].append(row)
                    if variant == "align_carrier_phase" and matched[1] is not None:
                        old_plv = results["full"][gate][-1]["original_mask_min_plv"]
                        phase_alignment_plv_error = max(phase_alignment_plv_error, abs(matched[1] - old_plv))
                    if case == "nominal" and variant == "full" and dt == previous_dt:
                        old = previous[window][str(gate)]
                        assert row["passes"] == old["passes"]
                        nominal_lag_error = max(nominal_lag_error, abs(row["max_abs_lag_deg"] - old["max_abs_lag"]))
                        nominal_plv_error = max(nominal_plv_error, abs(row["min_plv"] - old["min_plv"]))
        gate = float(spec.get("gate", MAIN_GATE))
        assert gate in GATES
        worst = min(results["full"][gate], key=margin)
        plot_window = spec.get("plot_window") or worst["window"]
        assert plot_window in windows
        requested_windows.add(plot_window)
        if case != "nominal":
            waveform.extend(waveform_rows(case, plot_window, windows[plot_window], frequency, dc, coefficients, carrier))
        summaries[case] = {
            **spec, "angles_rad": angles.tolist(), "gate": gate, "plot_window": plot_window,
            "spectrum": spectrum_summary(frequency, coefficients, carrier),
            "worst_physical_window": worst,
            "plot_window_metrics": {variant: next(row for row in by_gate[gate] if row["window"] == plot_window)
                                    for variant, by_gate in results.items()},
            "variants": {variant: {str(g): summarize(rows) for g, rows in by_gate.items()}
                         for variant, by_gate in results.items()},
        }
    for window in sorted(requested_windows):
        waveform.extend(waveform_rows("nominal", window, windows[window], frequency, nominal_dc, nominal_coeff, carrier))
    assert reconstruction_error < 1e-10 and phase_alignment_error < 1e-12 and phase_alignment_plv_error < 1e-12
    assert nominal_lag_error < 1e-9 and nominal_plv_error < 1e-9
    fine_audit = check_finer(model, summaries, windows, nominal_coeff, carrier, dt)
    write_csv("diagnostic_spectrum.csv", spectra)
    write_csv("diagnostic_waveforms.csv", waveform)
    write_csv("diagnostic_windows.csv", window_rows)
    summary = {
        "prior_selection_sha256": checksum, "windows": windows, "dt": dt, "carrier_omega": float(frequency[carrier]),
        "variant_definitions": VARIANTS, "ungated_phase_epsilon": 1e-14,
        "original_mask_phase_definition": "Each diagnostic variant uses the full physical signal's common envelope mask for the original_mask phase columns. Other phase columns use each variant's natural mask.",
        "cases": summaries,
        "audit": {"direct_reconstruction_error": reconstruction_error,
                  "phase_alignment_envelope_invariance_error": phase_alignment_error,
                  "phase_alignment_original_mask_plv_invariance_error": phase_alignment_plv_error,
                  "nominal_lag_error_deg": nominal_lag_error, "nominal_plv_error": nominal_plv_error,
                  "fine_grid": fine_audit},
    }
    (OUTPUT / "diagnostics.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                                              encoding="utf-8")


if __name__ == "__main__":
    main()
