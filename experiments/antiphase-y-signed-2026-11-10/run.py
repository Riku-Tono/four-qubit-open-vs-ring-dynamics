"""Check the even and odd response to a local Y rotation."""

import csv
import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
PREVIOUS = run_path(str(OUTPUT.parent / "antiphase-y-rotation-2026-11-10" / "run.py"))
COUPLINGS, TOPOLOGIES = PREVIOUS["COUPLINGS"], PREVIOUS["TOPOLOGIES"]
WINDOWS, DT, GATE = PREVIOUS["WINDOWS"], PREVIOUS["DT"], PREVIOUS["GATE"]
OCCUPATION = PREVIOUS["OCCUPATION"]
OMEGA = PREVIOUS["REFERENCE"]["OMEGA"]
MAGNITUDES = PREVIOUS["ANGLES"][1:]
ANGLES = (0.0,) + tuple(angle for magnitude in MAGNITUDES for angle in (magnitude, -magnitude))


def majorana_coefficients(time, coupling):
    """Evaluate the 00 and 02 entries of exp(M t) for the open chain."""
    edges = np.array([OMEGA, coupling, OMEGA, coupling, OMEGA, coupling, OMEGA])
    matrix = np.diag(edges, 1) - np.diag(edges, -1)
    energy, vectors = np.linalg.eigh(1j * matrix)
    products = vectors[0] * vectors[[0, 2]].conj()
    propagator = np.exp(-1j * time[:, None] * energy) @ products.T
    return -propagator[:, 0].real, 0.5 * propagator[:, 1].real


def ring_coefficients(prepared, time, other_site):
    """Propagate |1000> and |0000> to evaluate the ring diagonal and cross terms."""
    energy, vectors, _, _ = prepared
    phases = np.exp(-1j * time[:, None] * energy)
    a = (phases * vectors[0b1000]) @ vectors.T
    b = (phases * vectors[0b0000]) @ vectors.T
    observable = OCCUPATION[:, 0] + OCCUPATION[:, other_site] - 1
    even = np.abs(b)**2 @ observable
    odd = -((a.conj() * b) @ observable).real
    return even, odd


def max_abs(values):
    return float(np.max(np.abs(values)))


def peak_data(values, time, weights):
    peak = int(np.argmax(np.abs(values)))
    return {
        "max_abs": float(abs(values[peak])),
        "signed_peak": float(values[peak]),
        "time_at_peak": float(time[peak]),
        "window_mean": float(np.average(values, weights=weights)),
    }


def target_phase(signal, weights, pair, mask=None):
    phase = PREVIOUS["phase_stats"](signal, weights, mask=mask)
    return {"coverage": phase["coverage"], "order": phase["order"], **phase["pairs"].get(pair, {})}


def summarize(time, delta, occupation, oscillation, coefficients, other_site):
    weights = np.ones(len(time))
    weights[[0, -1]] = 0.5
    target_pair = f"0{other_site}"
    dc = occupation[:, 0] - oscillation[:, 0].real
    residual_dc = dc[:, 0] + dc[:, other_site] - 1
    angles = {}
    for angle, values, signal, constant in zip(ANGLES, delta, oscillation, residual_dc):
        angles[str(angle)] = {
            **peak_data(values, time, weights),
            "initial_value": float(values[0]),
            "spectral_dc": float(constant),
            "phase": target_phase(signal, weights, target_pair),
        }

    pairs = {}
    for index, angle in enumerate(MAGNITUDES, start=1):
        plus, minus = 2 * index - 1, 2 * index
        even = (delta[plus] + delta[minus]) / 2
        odd = (delta[plus] - delta[minus]) / 2
        even_factor, odd_factor = np.sin(angle / 2)**2, np.sin(angle)
        even_peak, odd_peak = max_abs(even), max_abs(odd)
        plus_valid = np.all(np.abs(oscillation[plus]) >= GATE, axis=1)
        minus_valid = np.all(np.abs(oscillation[minus]) >= GATE, axis=1)
        pairs[str(angle)] = {
            "even_max_abs": even_peak,
            "odd_max_abs": odd_peak,
            "even_max_abs_over_sin_half_squared": even_peak / even_factor,
            "odd_max_abs_over_sin": odd_peak / odd_factor,
            "even_pointwise_error": max_abs(even - even_factor * coefficients[0]),
            "odd_pointwise_error": max_abs(odd - odd_factor * coefficients[1]),
            "phase_on_common_support": {
                "plus": target_phase(oscillation[plus], weights, target_pair, mask=minus_valid),
                "minus": target_phase(oscillation[minus], weights, target_pair, mask=plus_valid),
            },
        }

    log2_orders = {}
    for component in ("even", "odd"):
        peaks = np.array([item[f"{component}_max_abs"] for item in pairs.values()])
        log2_orders[component] = (
            np.log2(peaks[1:] / peaks[:-1]).tolist() if np.all(peaks > 1e-11) else None
        )
    return {
        "dt": float(time[1] - time[0]),
        "sample_count": len(time),
        "angles": angles,
        "angle_pairs": pairs,
        "log2_orders_between_adjacent_magnitudes": log2_orders,
        "coefficient_peaks": {"even": max_abs(coefficients[0]), "odd": max_abs(coefficients[1])},
    }


def write_waveforms(writer, topology, coupling, time, delta, window):
    if window == "extended":
        return
    indices = np.r_[np.arange(0, len(time) - 1, 5), len(time) - 1]
    if window == "long":
        indices = indices[time[indices] > WINDOWS["short"][1]]
    sampled = delta[:, indices].T
    for t, values in zip(time[indices], sampled):
        row = [topology, coupling, f"{t:.10g}", f"{values[0]:.12g}"]
        for plus, minus in values[1:].reshape(-1, 2):
            row.extend(f"{v:.12g}" for v in (plus, minus, (plus + minus) / 2, (plus - minus) / 2))
        writer.writerow(row)


def main():
    results = {
        "Omega": OMEGA,
        "angles_rad": ANGLES,
        "couplings": COUPLINGS,
        "windows": WINDOWS,
        "phase_gate_all_sites": GATE,
        "zero_gap_tolerance": PREVIOUS["ZERO_GAP"],
        "cases": {},
    }
    max_errors = {}
    with (OUTPUT / "waveforms.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["topology", "J", "t", "delta_0"] + [
            f"{name}_{angle}"
            for angle in MAGNITUDES
            for name in ("delta_plus", "delta_minus", "even", "odd")
        ])
        for topology, link, other_site in TOPOLOGIES:
            by_j = {}
            for coupling in COUPLINGS:
                prepared = PREVIOUS["prepare"](coupling, link)
                by_window = {}
                for window, bounds in WINDOWS.items():
                    time = PREVIOUS["time_grid"](bounds, DT)
                    occupation, oscillation, errors = PREVIOUS["sample"](prepared, time, angles=ANGLES)
                    delta = occupation[:, :, 0] + occupation[:, :, other_site] - 1
                    if topology == "open":
                        coefficients = majorana_coefficients(time, coupling)
                    else:
                        coefficients = ring_coefficients(prepared, time, other_site)
                    stats = summarize(time, delta, occupation, oscillation, coefficients, other_site)
                    checks = {
                        **errors,
                        "zero_rotation_residual": max_abs(delta[0]),
                        "initial_signed_residual": max_abs(delta[:, 0] + np.sin(np.asarray(ANGLES) / 2)**2),
                        f"{topology}_even_formula": max(
                            pair["even_pointwise_error"] for pair in stats["angle_pairs"].values()
                        ),
                        f"{topology}_odd_formula": max(
                            pair["odd_pointwise_error"] for pair in stats["angle_pairs"].values()
                        ),
                    }
                    if coupling == 0:
                        checks["uncoupled_odd"] = max(
                            pair["odd_max_abs"] for pair in stats["angle_pairs"].values()
                        )
                    assert max(checks.values()) < 1e-11, (topology, coupling, window, checks)
                    for name, value in checks.items():
                        max_errors[name] = max(max_errors.get(name, 0.0), value)
                    by_window[window] = stats
                    write_waveforms(writer, topology, coupling, time, delta, window)
                by_j[str(coupling)] = by_window
            results["cases"][topology] = by_j
    results["max_audit_errors"] = max_errors
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n",
                                         encoding="utf-8")
    print(json.dumps(max_errors, indent=2))


if __name__ == "__main__":
    main()
