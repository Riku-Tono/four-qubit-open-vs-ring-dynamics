"""Reevaluate the existing four-site preparations at the matched coupling."""

import csv
import json
from pathlib import Path

import numpy as np

from metrics import GATES, PAIRS, WINDOWS, Model, describe, time_grid


OUTPUT = Path(__file__).resolve().parent
EXPERIMENTS = OUTPUT.parent
TOPOLOGIES = ("open", "ring")
BASIS_STATES = ("0000", "1000", "0100")
ANGLES = (0.025, 0.05, 0.10, 0.20)
BASELINE_WINDOWS = ("short", "full", "holdout", "late")
SITES = ("R0", "R1", "R2", "R3")


def preparation_specs(topology):
    """Return the basis preparation and site-0 XY-axis rotation parameters."""
    source = EXPERIMENTS / "antiphase-axis-optimization-2026-11-10" / "optima.csv"
    with source.open(encoding="utf-8", newline="") as handle:
        axes = {row["window"]: float(row["beta_opt_deg"]) for row in csv.DictReader(handle)
                if row["topology"] == topology and float(row["J"]) == 0.7}
    assert set(axes) == {"short", "long"}

    specs = {bits: (bits, 0.0, None) for bits in BASIS_STATES}
    rotations = (("X", 0.0, ANGLES), ("Y", 90.0, ANGLES),
                 ("short_axis", axes["short"], ANGLES),
                 ("long_axis", axes["long"], ANGLES),
                 ("phi_pi_4_audit", -135.0, (0.20,)))
    for name, beta, magnitudes in rotations:
        for magnitude in magnitudes:
            digits = 3 if magnitude == 0.025 else 2
            for theta in (magnitude, -magnitude):
                specs[f"{name}:{theta:+.{digits}f}"] = ("1000", theta, beta)
    return specs


def preparations(topology):
    states = {}
    for name, (bits, theta, beta) in preparation_specs(topology).items():
        index = int(bits, 2)
        state = np.zeros(16, dtype=complex)
        state[index] = np.cos(theta / 2)
        if beta is not None:
            state[0b0000] = -1j * np.exp(-1j * np.deg2rad(beta)) * np.sin(theta / 2)
        states[name] = state
    return states


def audit_saved_phases(rows):
    """Check the inherited |1000> phase convention and signed X/Y rotations."""
    with (EXPERIMENTS / "phase-alignment-2026-11-10" / "results.json").open(encoding="utf-8") as handle:
        phase_reference = json.load(handle)["cases"]
    index = {(row["topology"], row["case"], row["window"], row["gate"]): row for row in rows}
    errors = dict(coverage=0.0, site_amplitude=0.0, plv=0.0, lag_deg=0.0)

    def compare_pair(row, pair, plv, lag):
        lag_error = (row[f"lag_deg_{pair}"] - lag + 180) % 360 - 180
        errors["lag_deg"] = max(errors["lag_deg"], abs(lag_error))
        errors["plv"] = max(errors["plv"], abs(row[f"plv_{pair}"] - plv))

    for topology in TOPOLOGIES:
        windows = phase_reference[topology]["0.7"]["windows"]
        for window, reference in windows.items():
            for gate in GATES:
                row = index[topology, "1000", window, gate]
                phase = reference["phase_by_threshold"][str(gate)]
                errors["coverage"] = max(errors["coverage"], abs(row["coverage"] - phase["coverage"]))
                for site, amplitude in enumerate(reference["site_occupation_amplitude"]):
                    errors["site_amplitude"] = max(errors["site_amplitude"], abs(row[f"amplitude{site}"] - amplitude))
                for pair, phase_pair in phase["pairs"].items():
                    compare_pair(row, pair, phase_pair["plv"], phase_pair["lag_deg"])

    source = EXPERIMENTS / "antiphase-xy-rotation-2026-11-10" / "phases.csv"
    with source.open(encoding="utf-8", newline="") as handle:
        for reference in csv.DictReader(handle):
            theta = float(reference["theta_rad"])
            if float(reference["J"]) != 0.7 or reference["support"] != "native" or not theta:
                continue
            digits = 3 if abs(theta) == 0.025 else 2
            name = f"{reference['axis']}:{theta:+.{digits}f}"
            window = {"short": "short", "long": "full"}[reference["window"]]
            row = index[reference["topology"], name, window, 0.05]
            errors["coverage"] = max(errors["coverage"], abs(row["coverage"] - float(reference["coverage"])))
            compare_pair(row, reference["pair"], float(reference["plv"]), float(reference["lag_deg"]))
    assert max(errors.values()) < 1e-7, errors
    return errors


def main():
    rows = []
    symmetry_error = 0.0
    with (OUTPUT / "baseline_waveforms.csv").open("w", encoding="utf-8", newline="") as handle:
        waves = csv.writer(handle)
        waves.writerow(["topology", "case", "t", *SITES])
        for topology in TOPOLOGIES:
            specs = preparation_specs(topology)
            states = preparations(topology)
            model = Model(topology)
            for window in BASELINE_WINDOWS:
                time = time_grid(WINDOWS[window])
                grid = model.prepare(time)
                names = list(states)
                for start in range(0, len(names), 16):
                    batch = names[start:start + 16]
                    occupations, signals = grid.sample(np.stack([states[name] for name in batch]))
                    for name, occupation, signal in zip(batch, occupations, signals):
                        bits, theta, beta = specs[name]
                        for gate in GATES:
                            metrics = describe(occupation, signal, gate)
                            amplitudes = metrics["amplitude"]
                            rows.append({
                                "topology": topology, "case": name, "basis": bits,
                                "theta_rad": theta, "beta_deg": beta, "window": window, "gate": gate,
                                **{key: value for key, value in metrics.items()
                                   if key not in ("amplitude", "pairs")},
                                "amp_min": min(amplitudes),
                                **{f"amplitude{site}": value for site, value in enumerate(amplitudes)},
                                **{f"{metric}_{pair}": metrics["pairs"].get(pair, {}).get(metric)
                                   for pair in PAIRS for metric in ("lag_deg", "plv", "pearson", "in_phase")},
                            })
                        if name == "0000":
                            permutation = (3, 2, 1, 0) if topology == "open" else (1, 2, 3, 0)
                            symmetry_error = max(symmetry_error, float(np.max(
                                np.abs(occupation - occupation[:, permutation])
                            )))
                        if window == "full":
                            indices = np.r_[np.arange(0, len(time) - 1, 20), len(time) - 1]
                            waves.writerows([topology, name, f"{time[i]:.8f}",
                                             *(f"{value:.10g}" for value in occupation[i])]
                                            for i in indices)

    assert symmetry_error < 1e-11
    audit = {"saved_phase_metrics_max_difference": audit_saved_phases(rows),
             "uniform_basis_occupation_symmetry_max_error": symmetry_error}
    with (OUTPUT / "baseline.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (OUTPUT / "baseline_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
