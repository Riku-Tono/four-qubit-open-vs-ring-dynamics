"""Measure occupation spread, amplitude, and phase in the |1000> experiment."""

import csv
import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
PREVIOUS = run_path(str(OUTPUT.parent / "uniform-coupling-sweep-2026-11-10" / "run.py"))
MODEL, OCCUPATION = PREVIOUS["MODEL"], PREVIOUS["OCCUPATION"]
PAIRS = PREVIOUS["PAIRS"]
LEFT, RIGHT = np.asarray(MODEL.pairs).T
COUPLINGS = (0.0, 0.2, 0.5, 0.7, 1.0, 1.5, 2.0)
WINDOWS = {"short": (0.0, PREVIOUS["T"]), "full": (0.0, 8 * np.pi),
           "late": (4 * np.pi, 8 * np.pi)}
DT, GATE = 0.002, 0.05
GATES = (0.025, GATE, 0.1)
ZERO_GAP = 1e-10
INITIAL_INDEX = 8  # Site 0 is the leftmost bit in |1000>.


def occupation_spectrum(coupling, link):
    base = MODEL.hamiltonian(link, 0)
    zz = np.diag(np.diag(base))
    hamiltonian = base - zz + coupling / PREVIOUS["REFERENCE"]["J"] * zz
    energy, vectors = np.linalg.eigh(hamiltonian)
    gaps = energy[:, None] - energy[None, :]
    initial = vectors[INITIAL_INDEX]
    matrix = np.einsum("sa,si,sb->abi", vectors, OCCUPATION, vectors)
    weighted = matrix * initial[:, None, None] * initial[None, :, None]
    positive = gaps > ZERO_GAP
    dc = weighted[np.abs(gaps) <= ZERO_GAP].sum(axis=0)
    return energy, vectors, dc, gaps[positive], 2 * weighted[positive]


def sample(spectrum, time, topology, coupling):
    energy, vectors, dc, frequencies, coefficients = spectrum
    oscillation = np.exp(1j * time[:, None] * frequencies) @ coefficients
    occupation = dc + oscillation.real
    states = (np.exp(-1j * time[:, None] * energy) * vectors[INITIAL_INDEX]) @ vectors.T
    direct = np.abs(states)**2 @ OCCUPATION
    anti_site = 3 if topology == "open" else 2
    errors = {
        "spectral_reconstruction": float(np.max(np.abs(occupation - direct))),
        "state_norm": float(np.max(np.abs(np.linalg.norm(states, axis=1) - 1))),
        "antiphase_signal": float(np.max(np.abs(oscillation[:, 0] + oscillation[:, anti_site]))),
        "complementary_occupations": float(np.max(np.abs(occupation[:, 0] + occupation[:, anti_site] - 1))),
    }
    if topology == "ring":
        errors["reflection_signal"] = float(np.max(np.abs(oscillation[:, 1] - oscillation[:, 3])))
    if coupling == 0:
        expected = (1 + np.cos(time[:, None]) * np.array([1, -1, -1, -1])) / 2
        errors["uncoupled_occupation"] = float(np.max(np.abs(occupation - expected)))
    assert max(errors.values()) < 1e-11
    return occupation, oscillation, errors


def time_grid(bounds, dt):
    start, end = bounds
    return np.linspace(start, end, round((end - start) / dt) + 1)


def describe(occupation, oscillation):
    amplitude = np.abs(oscillation)
    unit_phase = np.exp(1j * np.angle(oscillation))
    relative_phase = unit_phase[:, LEFT] * unit_phase[:, RIGHT].conj()
    weights = np.ones(len(occupation))
    weights[[0, -1]] = 0.5
    spatial_variance = np.var(occupation, axis=1)
    peak_amplitudes = np.ptp(occupation, axis=0) / 2
    result = {
        "occupation_spread": float(np.sqrt(np.average(spatial_variance, weights=weights))),
        "mean_occupation_amplitude": float(peak_amplitudes.mean()),
        "site_occupation_amplitude": peak_amplitudes.tolist(),
        "mean_envelope": float(np.average(amplitude.mean(axis=1), weights=weights)),
        "phase_by_threshold": {},
    }
    for threshold in GATES:
        valid_weights = weights * np.all(amplitude >= threshold, axis=1)
        phases = {"coverage": float(valid_weights.sum() / weights.sum()), "order": None,
                  "min_pair_plv": None, "pairs": {}}
        if valid_weights.sum() > 0:
            order = np.abs(unit_phase.mean(axis=1))
            locking = np.average(relative_phase, axis=0, weights=valid_weights)
            phases.update(
                order=float(np.average(order, weights=valid_weights)),
                min_pair_plv=float(np.min(np.abs(locking))),
                pairs={pair: {"plv": float(abs(value)), "lag_deg": float(np.angle(value, deg=True))}
                       for pair, value in zip(PAIRS, locking)},
            )
        result["phase_by_threshold"][str(threshold)] = phases
    return result


def save_waveforms(writer, topology, coupling, time, occupation, oscillation):
    indices = np.unique(np.r_[np.arange(0, len(time), 10), len(time) - 1])
    amplitude = np.abs(oscillation[indices])
    valid = np.all(amplitude >= GATE, axis=1)
    phase = np.angle(oscillation[indices])
    unit = np.exp(1j * phase)
    pair_angle = np.abs(np.angle(unit[:, LEFT] * unit[:, RIGHT].conj(), deg=True))
    for t, values, radius, angles, gaps, keep in zip(time[indices], occupation[indices],
                                                    amplitude, phase, pair_angle, valid):
        writer.writerow([topology, coupling, f"{t:.8f}", *(f"{v:.10g}" for v in values),
                         *(f"{a:.10g}" for a in radius), *(f"{p:.10g}" for p in angles),
                         int(keep), *(f"{gap:.9g}" if keep else "" for gap in gaps)])


def main():
    results = {"N": 4, "Omega": PREVIOUS["REFERENCE"]["OMEGA"], "initial": "1000",
               "J": COUPLINGS, "windows": WINDOWS, "dt_target": DT,
               "phase_amplitude_thresholds": GATES, "main_threshold": GATE, "cases": {}}
    errors = {}
    dt_errors = {key: 0.0 for key in ("occupation_spread", "mean_occupation_amplitude",
                                      "phase_coverage", "phase_order", "min_pair_plv")}
    with (OUTPUT / "summary.csv").open("w", encoding="utf-8", newline="") as summary_file, (
        OUTPUT / "waveforms.csv"
    ).open("w", encoding="utf-8", newline="") as wave_file:
        summary, waves = csv.writer(summary_file), csv.writer(wave_file)
        summary.writerow(["topology", "J", "window", "occupation_spread", "mean_occupation_amplitude",
                          "mean_envelope", "phase_coverage", "phase_order", "min_pair_plv",
                          *(f"occupation_amplitude{i}" for i in range(4))])
        waves.writerow(["topology", "J", "t", *(f"R{i}" for i in range(4)),
                        *(f"amplitude{i}" for i in range(4)), *(f"phase{i}" for i in range(4)),
                        "phase_valid", *(f"angle_{pair}" for pair in PAIRS)])
        for topology, link in (("open", 0), ("ring", 1)):
            results["cases"][topology] = {}
            for coupling in COUPLINGS:
                spectrum = occupation_spectrum(coupling, link)
                _, _, dc, _, _ = spectrum
                case = {"dc_occupation": dc.tolist(), "windows": {}}
                results["cases"][topology][str(coupling)] = case
                for name, bounds in WINDOWS.items():
                    time = time_grid(bounds, DT)
                    occupation, oscillation, sample_errors = sample(spectrum, time, topology, coupling)
                    result = describe(occupation, oscillation)
                    case["windows"][name] = result
                    for key, value in sample_errors.items():
                        errors[key] = max(errors.get(key, 0), value)
                    fine_occupation, fine_oscillation, _ = sample(
                        spectrum, time_grid(bounds, DT / 2), topology, coupling
                    )
                    fine = describe(fine_occupation, fine_oscillation)
                    primary_phase, fine_phase = (entry["phase_by_threshold"][str(GATE)]
                                                 for entry in (result, fine))
                    checks = {
                        "occupation_spread": (result["occupation_spread"], fine["occupation_spread"]),
                        "mean_occupation_amplitude": (result["mean_occupation_amplitude"],
                                                      fine["mean_occupation_amplitude"]),
                        "phase_coverage": (primary_phase["coverage"], fine_phase["coverage"]),
                        "phase_order": (primary_phase["order"], fine_phase["order"]),
                        "min_pair_plv": (primary_phase["min_pair_plv"], fine_phase["min_pair_plv"]),
                    }
                    for key, (a, b) in checks.items():
                        if a is not None and b is not None:
                            dt_errors[key] = max(dt_errors[key], abs(a - b))
                    summary.writerow([topology, coupling, name, result["occupation_spread"],
                                      result["mean_occupation_amplitude"], result["mean_envelope"],
                                      primary_phase["coverage"], primary_phase["order"],
                                      primary_phase["min_pair_plv"], *result["site_occupation_amplitude"]])
                    if name == "full":
                        save_waveforms(waves, topology, coupling, time, occupation, oscillation)
    results["audit"] = {"max_residual": errors, "max_summary_change_dt_halving": dt_errors}
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n",
                                         encoding="utf-8")
    print(json.dumps(results["audit"], indent=2))


if __name__ == "__main__":
    main()
