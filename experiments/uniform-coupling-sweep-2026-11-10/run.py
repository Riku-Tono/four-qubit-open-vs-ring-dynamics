"""Compare four-site open and periodic chains from |0000> at matched J."""

import csv
import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
REFERENCE = run_path(str(OUTPUT.parent / "endpoint-field-six-qubit-2026-11-09" / "run.py"))
MODEL = REFERENCE["Chain"](4)
CHANNELS = ("XX", "YY", "ZZ")
PAIRS = tuple(f"{i}{j}" for i, j in MODEL.pairs)
COUPLINGS = (0.7, 1.0, 2.0)
T, DT, FINE_DT = 3.688, 0.002, 0.001
OCCUPATION = (np.arange(16)[:, None] >> np.arange(3, -1, -1)) & 1


def evolve(time, coupling, link):
    baseline = MODEL.hamiltonian(link, 0)
    # At zero endpoint perturbation, ZZ is diagonal and the X drive is off-diagonal.
    zz = np.diag(np.diag(baseline))
    hamiltonian = (baseline - zz) + (coupling / REFERENCE["J"]) * zz
    energy, vectors = np.linalg.eigh(hamiltonian)
    return (np.exp(-1j * time[:, None] * energy) * vectors[0]) @ vectors.T


def summarize(time, opened, ring, occupations):
    delta = ring - opened
    power = np.trapezoid(delta**2, x=time, axis=0) / time[-1]
    combined = np.sqrt(np.sum(delta**2, axis=(1, 2)))
    peak = int(np.argmax(combined))
    total = power.sum()
    pair_channels = {}
    for i, pair in enumerate(PAIRS):
        pair_channels[pair] = {}
        for k, channel in enumerate(CHANNELS):
            n = int(np.argmax(np.abs(delta[:, i, k])))
            pair_channels[pair][channel] = {
                "rms": float(np.sqrt(power[i, k])),
                "peak_t": float(time[n]),
                "signed_peak": float(delta[n, i, k]),
                "open_at_peak": float(opened[n, i, k]),
                "ring_at_peak": float(ring[n, i, k]),
                "delta_at_end": float(delta[-1, i, k]),
            }
    return {
        "combined_rms": float(np.sqrt(total)),
        "combined_peak": {"t": float(time[peak]), "distance": float(combined[peak]),
                          "open_occupations": occupations[0, peak].tolist(),
                          "ring_occupations": occupations[1, peak].tolist()},
        "channel_rms": dict(zip(CHANNELS, np.sqrt(power.sum(axis=0)).tolist())),
        "channel_share": dict(zip(CHANNELS, (power.sum(axis=0) / total).tolist())),
        "pair_share": dict(zip(PAIRS, (power.sum(axis=1) / total).tolist())),
        "pair_channels": pair_channels,
    }


def measure(time, coupling):
    states = [evolve(time, coupling, link) for link in (0, 1)]
    opened, ring = [MODEL.connected(branch).swapaxes(1, 2) for branch in states]
    reflection = [MODEL.pairs.index((3 - j, 3 - i)) for i, j in MODEL.pairs]
    rotation = [MODEL.pairs.index(tuple(sorted(((i + 1) % 4, (j + 1) % 4))))
                for i, j in MODEL.pairs]
    norm_error = max(float(np.max(np.abs(np.linalg.norm(branch, axis=1) - 1)))
                     for branch in states)
    symmetry_error = float(max(np.max(np.abs(opened - opened[:, reflection])),
                               np.max(np.abs(ring - ring[:, rotation]))))
    initial_error = float(max(np.max(np.abs(opened[0])), np.max(np.abs(ring[0]))))
    assert max(norm_error, symmetry_error, initial_error) < 1e-12
    occupations = np.abs(np.asarray(states))**2 @ OCCUPATION
    audit = {"max_norm_error": norm_error, "max_symmetry_error": symmetry_error,
             "max_initial_correlation": initial_error}
    return opened, ring, occupations, audit


def main():
    time = np.arange(round(T / FINE_DT) + 1) * FINE_DT
    columns = [f"{branch}_{pair}_{channel}" for branch in ("open", "ring")
               for pair in PAIRS for channel in CHANNELS]
    results = {"N": 4, "Omega": REFERENCE["OMEGA"], "initial": "0000", "T": T,
               "dt": DT, "fine_dt": FINE_DT, "delta": "ring - open", "cases": {}}
    series = []
    for coupling in COUPLINGS:
        opened, ring, occupations, audit = measure(time, coupling)
        primary = summarize(time[::2], opened[::2], ring[::2], occupations[:, ::2])
        fine = summarize(time, opened, ring, occupations)
        rms_change = max(abs(primary["pair_channels"][p][a]["rms"]
                             - fine["pair_channels"][p][a]["rms"])
                         for p in PAIRS for a in CHANNELS)
        assert rms_change < 1e-7
        audit["max_pair_channel_rms_change_dt_halving"] = rms_change
        if coupling == 0.7:
            # The archived baseline used symmetric second-order Trotter evolution.
            index = round(3.462 / FINE_DT)
            archive_error = abs(np.linalg.norm(ring[index] - opened[index]) - 0.974637690628)
            assert archive_error < 1e-6
            audit["archive_combined_distance_error_at_t_3.462"] = float(archive_error)
        results["cases"][str(coupling)] = {"fine": fine, "primary": primary, "audit": audit}
        series.append((coupling, np.column_stack([opened[::2].reshape(len(time[::2]), -1),
                                                   ring[::2].reshape(len(time[::2]), -1)])))

    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    with (OUTPUT / "waveforms.csv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["J", "t", *columns])
        for coupling, waves in series:
            writer.writerows([coupling, f"{t:.3f}", *(f"{v:.12g}" for v in row)]
                             for t, row in zip(time[::2], waves))
    for coupling, case in results["cases"].items():
        fine = case["fine"]
        print(f"J={coupling}  RMS={fine['combined_rms']:.9f}  peak={fine['combined_peak']}")


if __name__ == "__main__":
    main()
