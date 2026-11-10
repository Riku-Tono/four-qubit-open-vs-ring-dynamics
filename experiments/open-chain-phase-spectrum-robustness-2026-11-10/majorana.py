"""Reconstruct the frozen state's end and interior spectra with Majorana modes."""

from collections import defaultdict
from itertools import combinations, permutations, product
import json
from pathlib import Path
import sys

import numpy as np


OUTPUT = Path(__file__).resolve().parent
EXPERIMENTS = OUTPUT.parent
PRIOR = EXPERIMENTS / "open-chain-phase-persistence-2026-11-10"
sys.path.insert(0, str(EXPERIMENTS / "four-site-in-phase-2026-11-10"))
from metrics import Model, product_state


def majoranas():
    x = np.array([[0, 1], [1, 0]])
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    a, b = [], []
    prefix = np.ones((1, 1))
    for site in range(4):
        suffix = np.eye(2**(3 - site))
        a.append(np.kron(np.kron(prefix, z), suffix))
        b.append(np.kron(np.kron(prefix, y), suffix))
        prefix = np.kron(prefix, x)
    return np.array(a), np.array(b)


def cubic_coefficients(left, right, frequencies, normal_operators, state):
    """Compute R1=(1-i A0 B0 A1)/2 by the exterior cube of mode evolution."""
    triples = np.asarray(list(combinations(range(8), 3)))
    moments = np.array([np.vdot(state, 1j * (normal_operators[i] @ normal_operators[j]
                                            @ normal_operators[k] @ state)) for i, j, k in triples])
    assert np.max(np.abs(moments.imag)) < 1e-14

    # Each signed mode denotes (a_k - i sign b_k)/2, in the (a,b) basis.
    mode = np.tile(np.arange(4), 2)
    sign = np.repeat([1, -1], 4)
    unit = np.eye(4)[mode]
    basis = np.concatenate((unit, -1j * sign[:, None] * unit), axis=1) / 2
    row_a0 = left[0, mode, None] * basis
    row_b0 = (1j * sign * right[0, mode])[:, None] * basis
    row_a1 = left[1, mode, None] * basis

    # Repeated signed modes have zero wedge; group each remaining term by its
    # exact signed integer combination of the four elementary frequencies.
    signed_units = (sign[:, None] * unit).astype(int)
    coefficients = defaultdict(complex)
    for i, j, k in permutations(range(8), 3):
        rows = np.array([row_a0[i], row_b0[j], row_a1[k]])
        label = tuple(signed_units[i] + signed_units[j] + signed_units[k])
        coefficients[label] += np.linalg.det(rows[:, triples].transpose(1, 0, 2)) @ moments.real
    positive = [(label @ frequencies, label, -coefficient) for label, coefficient in coefficients.items()
                if label @ frequencies > 0]
    positive.sort(key=lambda item: item[0])
    assert len(positive) == 20
    assert sum(sum(abs(n) for n in label) == 1 for _, label, _ in positive) == 4
    assert sum(sum(abs(n) for n in label) == 3 for _, label, _ in positive) == 16
    return positive


def power_fraction(values):
    return float(abs(values[0])**2 / np.sum(np.abs(values)**2))


def spectral_row(frequency, label, coefficient):
    return {"frequency": float(frequency), "signed_mode_coefficients": [int(n) for n in label],
            "real": float(coefficient.real), "imag": float(coefficient.imag),
            "amplitude": float(abs(coefficient)), "phase_deg": float(np.angle(coefficient, deg=True))}


def reflection_transitions(model, initial):
    """Identify spectral lines inaccessible from a reflection-even initial state."""
    reflected_bits = [int(f"{i:04b}"[::-1], 2) for i in range(16)]
    reflection = np.einsum("sa,sa->a", model.vectors, model.vectors[reflected_bits])
    assert np.max(np.abs(np.abs(reflection) - 1)) < 1e-12
    even = reflection > 0
    assert np.max(np.abs(initial[~even])) < 1e-12

    allowed = even[model.a] & even[model.b] & (np.abs(model.observed[1]) > 1e-13)
    counts = np.add.reduceat(allowed.astype(int), model.group_starts)
    return {"max_odd_state_overlap": float(np.max(np.abs(initial[~even]))),
            "parity_by_energy": np.sign(reflection).astype(int).tolist(),
            "even_transition_counts": counts.tolist(),
            "forbidden_interior_frequencies": model.frequencies[counts == 0].tolist()}


def main():
    selection = json.loads((PRIOR / "selection.json").read_text(encoding="utf-8"))["selection"]["0.05"]
    angles = np.asarray(selection["families"][selection["selected_family"]]["angles_rad"])
    state = product_state(angles)
    model = Model("open")
    theta, phi = angles.T
    x, y, z = np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)
    prefix = np.r_[1, np.cumprod(x)[:-1]]
    mean_a, mean_b = z * prefix, y * prefix

    hopping = np.eye(4) - 0.7 * np.eye(4, k=-1)
    left, frequencies, right_t = np.linalg.svd(hopping)
    order = np.argsort(frequencies)
    left, frequencies, right = left[:, order], frequencies[order], right_t[order].T
    alpha, beta = left.T @ mean_a, right.T @ mean_b
    end = -left[0] * (alpha - 1j * beta) / 2

    a, b = majoranas()
    normal = np.concatenate((np.einsum("jk,jab->kab", left, a), np.einsum("jk,jab->kab", right, b)))
    inner = cubic_coefficients(left, right, frequencies, normal, state)
    reconstructed_h = 0.5j * np.einsum("jk,jab,kbc->ac", hopping, a, b)
    energies = np.sort(0.5 * np.asarray(list(product((-1, 1), repeat=4))) @ frequencies)

    initial = state @ model.vectors
    products = initial[model.a].conj() * initial[model.b]
    direct = np.add.reduceat(products * model.observed, model.group_starts, axis=-1)
    nearest = np.argmin(np.abs(frequencies[:, None] - model.frequencies), axis=1)
    dc = np.einsum("a,abi,b->i", initial.conj(), model.dc, initial)
    audit = {
        "hamiltonian_matrix_error": float(np.max(np.abs(reconstructed_h - model.hamiltonian))),
        "energy_error": float(np.max(np.abs(energies - model.energy))),
        "frequency_error": float(np.max(np.abs(np.array([item[0] for item in inner]) - model.frequencies))),
        "end_coefficient_error": float(np.max(np.abs(end - direct[0, nearest]))),
        "inner_coefficient_error": float(np.max(np.abs(np.array([item[2] for item in inner]) - direct[1]))),
        "dc_error": float(np.max(np.abs(dc - 0.5))),
    }
    assert max(audit.values()) < 1e-12, audit
    # P_n(s)=(s-1.49) P_(n-1)(s)-0.49 P_(n-2)(s), P_0=1 and P_1=s-1.
    exact_polynomial = [1000000, -5470000, 9660300, -6067849, 1000000]
    assert np.max(np.abs(np.poly(frequencies**2) - np.asarray(exact_polynomial) / 1e6)) < 1e-12
    output = {"angles_rad": angles.tolist(), "bloch_vectors": np.stack((x, y, z), axis=1).tolist(),
              "hopping": hopping.tolist(), "normal_mode_frequencies": frequencies.tolist(),
              "squared_frequency_polynomial_numerators": exact_polynomial,
              "squared_frequency_polynomial_denominator": 1000000,
              "minimum_energy_spacing": float(np.min(np.diff(energies))),
              "majorana_mean_a": mean_a.tolist(), "majorana_mean_b": mean_b.tolist(),
              "initial_mode_a": alpha.tolist(), "initial_mode_b": beta.tolist(),
              "end_site_mode_weights": left[0].tolist(),
              "lowest_mode_fraction_a": power_fraction(alpha),
              "lowest_mode_fraction_b": power_fraction(beta),
              "lowest_mode_fraction_combined": power_fraction(alpha - 1j * beta),
              "end_spectrum": [spectral_row(f, label, c) for f, label, c in zip(frequencies, np.eye(4), end)],
              "interior_spectrum": [spectral_row(f, label, c) for f, label, c in inner],
              "reflection": reflection_transitions(model, initial), "audit": audit}
    (OUTPUT / "majorana.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
