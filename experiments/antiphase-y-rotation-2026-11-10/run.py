"""Measure how a local Y rotation changes occupations and their phases."""

import csv
import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
PREVIOUS = run_path(str(OUTPUT.parent / "phase-alignment-2026-11-10" / "run.py"))
MODEL, OCCUPATION = PREVIOUS["MODEL"], PREVIOUS["OCCUPATION"]
REFERENCE = PREVIOUS["PREVIOUS"]["REFERENCE"]
PAIRS, LEFT, RIGHT = PREVIOUS["PAIRS"], PREVIOUS["LEFT"], PREVIOUS["RIGHT"]
time_grid = PREVIOUS["time_grid"]
DT, GATE, GATES = PREVIOUS["DT"], PREVIOUS["GATE"], PREVIOUS["GATES"]
TOPOLOGIES = (("open", 0, 3), ("ring", 1, 2))
COUPLINGS = (0.0, 0.7)
ANGLES = (0.0, 0.025, 0.05, 0.1, 0.2)
WINDOWS = {"short": PREVIOUS["WINDOWS"]["short"],
           "long": PREVIOUS["WINDOWS"]["full"], "extended": (0.0, 256.0)}
ZERO_GAP = PREVIOUS["ZERO_GAP"]


def prepare(coupling, link):
    base = MODEL.hamiltonian(link, 0)
    zz = np.diag(np.diag(base))
    hamiltonian = base - zz + coupling / REFERENCE["J"] * zz
    energy, vectors = np.linalg.eigh(hamiltonian)
    gaps = energy[:, None] - energy[None, :]
    matrix = np.einsum("sa,si,sb->abi", vectors, OCCUPATION, vectors)
    return energy, vectors, gaps, matrix


def sample(prepared, time, angles=ANGLES):
    energy, vectors, gaps, matrix = prepared
    half_angles = np.asarray(angles)[:, None] / 2
    # exp(-i theta Y_0 / 2) |1000> = cos(theta/2)|1000> - sin(theta/2)|0000>.
    projection = np.cos(half_angles) * vectors[0b1000] - np.sin(half_angles) * vectors[0b0000]
    weighted = matrix[None, ...] * projection[:, :, None, None] * projection[:, None, :, None]
    dc = weighted[:, np.abs(gaps) <= ZERO_GAP].sum(axis=1)
    positive = gaps > ZERO_GAP
    frequencies = gaps[positive]
    coefficients = 2 * weighted[:, positive, :].transpose(1, 0, 2).reshape(len(frequencies), -1)
    occupation = np.empty((len(angles), len(time), 4))
    oscillation = np.empty(occupation.shape, dtype=complex)
    error = {"spectral_reconstruction": 0.0, "state_norm_squared": 0.0}
    for start in range(0, len(time), 4096):
        block = slice(start, start + 4096)
        signal = np.exp(1j * time[block, None] * frequencies) @ coefficients
        oscillation[:, block] = signal.reshape(-1, len(angles), 4).transpose(1, 0, 2)
        occupation[:, block] = dc[:, None, :] + oscillation[:, block].real
        evolution = np.exp(-1j * time[block, None] * energy)
        states = (evolution[None, ...] * projection[:, None, :]) @ vectors.T
        probabilities = np.abs(states)**2
        error["spectral_reconstruction"] = max(
            error["spectral_reconstruction"], float(np.max(np.abs(occupation[:, block] - probabilities @ OCCUPATION)))
        )
        error["state_norm_squared"] = max(
            error["state_norm_squared"], float(np.max(np.abs(probabilities.sum(axis=2) - 1)))
        )
    return occupation, oscillation, error


def phase_stats(oscillation, weights, threshold=GATE, mask=None):
    valid = np.all(np.abs(oscillation) >= threshold, axis=1)
    if mask is not None:
        valid &= mask
    coverage = float(weights[valid].sum() / weights.sum())
    if not np.any(valid):
        return {"coverage": coverage, "order": None, "pairs": {}}
    active = weights[valid]
    unit = np.exp(1j * np.angle(oscillation[valid]))
    relative = unit[:, LEFT] * unit[:, RIGHT].conj()
    locking = np.average(relative, axis=0, weights=active)
    absolute_lag = np.average(np.abs(np.angle(relative, deg=True)), axis=0, weights=active)
    return {
        "coverage": coverage,
        "order": float(np.average(np.abs(unit.mean(axis=1)), weights=active)),
        "pairs": {pair: {"plv": float(abs(value)), "lag_deg": float(np.angle(value, deg=True)),
                         "mean_abs_lag_deg": float(absolute)}
                  for pair, value, absolute in zip(PAIRS, locking, absolute_lag)},
    }


def summarize(occupation, oscillation, time, other_site, *, matched):
    weights = np.ones(len(time))
    weights[[0, -1]] = 0.5
    baseline = oscillation[0]
    baseline_valid = np.all(np.abs(baseline) >= GATE, axis=1)
    metrics = []
    for values, signal in zip(occupation, oscillation):
        delta = values[:, 0] + values[:, other_site] - 1
        peak = int(np.argmax(np.abs(delta)))
        amplitude = np.ptp(values, axis=0) / 2
        result = {
            "initial_signed_residual": float(delta[0]),
            "max_abs_residual": float(abs(delta[peak])),
            "first_time_at_max_abs_residual": float(time[peak]),
            "rms_residual": float(np.sqrt(np.average(delta**2, weights=weights))),
            "occupation_spread": float(np.sqrt(np.average(np.var(values, axis=1), weights=weights))),
            "site_amplitude": amplitude.tolist(),
            "mean_site_amplitude": float(amplitude.mean()),
            "phase_by_threshold": {str(gate): phase_stats(signal, weights, gate) for gate in GATES},
        }
        if matched:
            valid = np.all(np.abs(signal) >= GATE, axis=1)
            result["phase_matched_to_zero"] = {
                "rotated": phase_stats(signal, weights, mask=baseline_valid),
                "zero": phase_stats(baseline, weights, mask=valid),
            }
        metrics.append(result)
    return metrics


def audit_dt(coarse, fine, audit):
    for key in ("max_abs_residual", "rms_residual", "occupation_spread", "site_amplitude"):
        delta = np.asarray(coarse[key]) - np.asarray(fine[key])
        audit[key] = max(audit.get(key, 0.0), float(np.max(np.abs(delta))))
    first, second = (item["phase_by_threshold"][str(GATE)] for item in (coarse, fine))
    for key in ("coverage", "order"):
        if first[key] is not None and second[key] is not None:
            audit[key] = max(audit.get(key, 0.0), abs(first[key] - second[key]))
    for pair in first["pairs"].keys() & second["pairs"].keys():
        first_pair, second_pair = first["pairs"][pair], second["pairs"][pair]
        for key in ("plv", "lag_deg", "mean_abs_lag_deg"):
            delta = first_pair[key] - second_pair[key]
            if key == "lag_deg":
                delta = (delta + 180) % 360 - 180
            audit[key] = max(audit.get(key, 0.0), abs(delta))


def save_phases(cases):
    with (OUTPUT / "phases.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("topology", "J", "window", "theta_rad", "pair", "coverage", "phase_order",
                         "lag_deg", "plv", "mean_abs_lag_deg"))
        for topology, couplings in cases.items():
            for coupling, case in couplings.items():
                for window, data in case["windows"].items():
                    for angle, metrics in data["angles"].items():
                        phase = metrics["phase_by_threshold"][str(GATE)]
                        for pair in PAIRS:
                            entry = phase["pairs"].get(pair, {})
                            values = (phase["coverage"], phase["order"], entry.get("lag_deg"),
                                      entry.get("plv"), entry.get("mean_abs_lag_deg"))
                            writer.writerow((topology, coupling, window, angle, pair,
                                             *("" if value is None else f"{value:.10g}" for value in values)))


def main():
    cases, checks, dt_checks = {}, {}, {}
    for topology, link, other_site in TOPOLOGIES:
        cases[topology] = {}
        for coupling in COUPLINGS:
            prepared = prepare(coupling, link)
            windows = {}
            for window, bounds in WINDOWS.items():
                time = time_grid(bounds, DT)
                fine_time = np.linspace(*bounds, 2 * len(time) - 1)
                occupation, oscillation, error = sample(prepared, fine_time)
                exact_t0 = np.cos(np.asarray(ANGLES) / 2)**2
                error["initial_occupation"] = float(np.max(np.abs(occupation[:, 0, 0] - exact_t0)))
                error["zero_rotation_complement"] = float(
                    np.max(np.abs(occupation[0, :, 0] + occupation[0, :, other_site] - 1))
                )
                if topology == "open":
                    error["open_right_endpoint_invariance"] = float(
                        np.max(np.abs(occupation[:, :, 3] - occupation[0, :, 3]))
                    )
                else:
                    error["ring_reflection"] = float(np.max(np.abs(occupation[:, :, 1] - occupation[:, :, 3])))
                if coupling == 0.0:
                    factors = np.column_stack([np.cos(ANGLES), -np.ones((len(ANGLES), 3))])
                    analytic = (1 + factors[:, None, :] * np.cos(REFERENCE["OMEGA"] * fine_time[None, :, None])) / 2
                    error["uncoupled_analytic_occupation"] = float(np.max(np.abs(occupation - analytic)))
                assert max(error.values()) < 1e-11
                for key, value in error.items():
                    checks[key] = max(checks.get(key, 0.0), value)
                coarse = summarize(occupation[:, ::2], oscillation[:, ::2], time, other_site, matched=True)
                fine = summarize(occupation, oscillation, fine_time, other_site, matched=False)
                for a, b in zip(coarse, fine):
                    audit_dt(a, b, dt_checks)
                windows[window] = {"sample_count": len(time), "dt": float(time[1] - time[0]),
                                   "angles": dict(zip(map(str, ANGLES), coarse))}
            cases[topology][str(coupling)] = {"complementary_sites": [0, other_site], "windows": windows}
    results = {
        "N": 4, "Omega": REFERENCE["OMEGA"], "initial": "1000", "rotation": "exp(-i theta Y_0 / 2)",
        "angles_rad": ANGLES, "couplings": COUPLINGS, "windows": WINDOWS, "dt_target": DT,
        "phase_gate": GATE, "phase_gate_sensitivity": GATES, "cases": cases,
        "audit": {"max_errors": checks, "max_change_under_dt_halving": dt_checks},
    }
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n",
                                         encoding="utf-8")
    save_phases(cases)
    print(json.dumps(results["audit"], indent=2))


if __name__ == "__main__":
    main()
