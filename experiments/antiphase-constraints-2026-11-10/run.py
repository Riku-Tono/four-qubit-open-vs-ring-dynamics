"""Audit the complementary occupations in the |1000> experiment (NumPy 2)."""

import json
from functools import reduce
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
SWEEP = run_path(str(OUTPUT.parent / "uniform-coupling-sweep-2026-11-10" / "run.py"))
MODEL, OCCUPATION = SWEEP["MODEL"], SWEEP["OCCUPATION"]
OMEGA, REFERENCE_J = SWEEP["REFERENCE"]["OMEGA"], SWEEP["REFERENCE"]["J"]
INITIAL = 0b1000
T, STEPS = 256.0, (0.08, 0.04, 0.02, 0.01, 0.005)
COUPLINGS = (-2.0, -0.7, 0.0, 0.2, 0.5, 0.7, 1.0, 1.5, 2.0, 4.0)
TOPOLOGIES = (("open", 0, 3), ("ring", 1, 2))


def exact_states(hamiltonian, time):
    energy, vectors = np.linalg.eigh(hamiltonian)
    return (np.exp(-1j * time[:, None] * energy) * vectors[INITIAL]) @ vectors.T


def strang_states(zz_diagonal, dt, count):
    """Apply exp(-i dt H_ZZ/2) exp(-i dt H_X) exp(-i dt H_ZZ/2)."""
    half_phase = np.exp(-0.5j * dt * zz_diagonal)
    angle = dt * OMEGA / 2
    local = np.array([[np.cos(angle), -1j * np.sin(angle)],
                      [-1j * np.sin(angle), np.cos(angle)]])
    drive = reduce(np.kron, [local] * MODEL.size)
    propagator = half_phase[:, None] * drive * half_phase[None, :]
    states = np.zeros((count, len(zz_diagonal)), dtype=complex)
    states[0, INITIAL] = 1.0
    for n in range(1, count):
        states[n] = propagator @ states[n - 1]
    return states


def occupation(states):
    probabilities = np.abs(states)**2
    raw = probabilities @ OCCUPATION
    norm_squared = probabilities.sum(axis=1)
    # Keep the raw values so drift from floating-point normalization can be audited.
    return raw / norm_squared[:, None], raw, norm_squared


def residual(occupations, other_site):
    return np.abs(occupations[:, 0] + occupations[:, other_site] - 1)


def audit(occupations, raw, norm_squared, time, other_site):
    difference = residual(occupations, other_site)
    peak = int(np.argmax(difference))
    return {
        "max_residual": float(difference[peak]),
        "first_time_at_max_residual": float(time[peak]),
        "max_raw_residual": float(residual(raw, other_site).max()),
        "max_norm_squared_error": float(np.max(np.abs(norm_squared - 1))),
    }


def main():
    results = {"N": MODEL.size, "Omega": OMEGA, "initial": "1000", "T": T,
               "J": COUPLINGS, "dt": STEPS, "cases": {}, "max_by_topology_and_dt": {}}
    fine_dt = min(STEPS)
    fine_time = np.arange(round(T / fine_dt) + 1) * fine_dt
    uncoupled_error = 0.0
    for topology, link, other_site in TOPOLOGIES:
        baseline = MODEL.hamiltonian(link, 0)
        baseline_zz = np.diag(baseline)
        drive = baseline - np.diag(baseline_zz)
        cases = {}
        for coupling in COUPLINGS:
            zz = coupling / REFERENCE_J * baseline_zz
            exact = occupation(exact_states(drive + np.diag(zz), fine_time))
            if coupling == 0:
                expected = (1 + np.cos(OMEGA * fine_time[:, None]) * np.array([1, -1, -1, -1])) / 2
                uncoupled_error = max(uncoupled_error, float(np.max(np.abs(exact[0] - expected))))
            by_dt = {}
            for dt in STEPS:
                stride = round(dt / fine_dt)
                time = fine_time[::stride]
                sampled = tuple(values[::stride] for values in exact)
                exact_audit = audit(*sampled, time, other_site)
                strang = occupation(strang_states(zz, dt, len(time)))
                strang_audit = audit(*strang, time, other_site)
                strang_audit["max_error_vs_exact_all_sites"] = float(np.max(np.abs(strang[0] - sampled[0])))
                assert exact_audit["max_residual"] < 2e-12
                assert strang_audit["max_residual"] < 2e-10
                assert strang_audit["max_norm_squared_error"] < 2e-10
                by_dt[str(dt)] = {"sample_count": len(time), "exact": exact_audit, "strang": strang_audit}
            cases[str(coupling)] = by_dt
        results["cases"][topology] = {"sites": [0, other_site], "by_J": cases}
        maxima = {}
        for dt in STEPS:
            dt_cases = [case[str(dt)] for case in cases.values()]
            maxima[str(dt)] = {
                method: {key: max(case[method][key] for case in dt_cases)
                         for key in dt_cases[0][method] if key != "first_time_at_max_residual"}
                for method in ("exact", "strang")
            }
        results["max_by_topology_and_dt"][topology] = maxima
    assert uncoupled_error < 2e-12
    results["max_uncoupled_exact_error_vs_cosine"] = uncoupled_error
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n",
                                         encoding="utf-8")
    print(json.dumps(results["max_by_topology_and_dt"], indent=2))


if __name__ == "__main__":
    main()
