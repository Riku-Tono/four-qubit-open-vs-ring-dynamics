"""Resolve the frozen open-chain state into exact gap groups and compare phase contributions."""

import csv
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import sys

import numpy as np


OUTPUT = Path(__file__).resolve().parent
EXPERIMENTS = OUTPUT.parent
PRIOR = EXPERIMENTS / "open-chain-phase-persistence-2026-11-10"
sys.path[:0] = [str(EXPERIMENTS / "four-site-in-phase-2026-11-10"),
                str(EXPERIMENTS / "antiphase-axis-long-time-2026-11-10")]
from metrics import MAIN_GATE, OCCUPATION, Model, product_state, time_grid
from spectral_groups import certify_gap_groups

sys.path.insert(0, str(PRIOR))
from certificate import spectral_bound


def save_csv(name, rows):
    with (OUTPUT / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def frozen_state():
    source = PRIOR / "selection.json"
    frozen = json.loads(source.read_text(encoding="utf-8"))["selection"][str(MAIN_GATE)]
    angles = np.asarray(frozen["families"][frozen["selected_family"]]["angles_rad"])
    assert np.array_equal(angles, angles[::-1])
    return angles, product_state(angles), hashlib.sha256(source.read_bytes()).hexdigest()


def frequencies_from_majoranas():
    # The endpoint is linear in Majoranas; the adjacent interior Z is a product of three.
    edges = np.array([1, 0.7, 1, 0.7, 1, 0.7, 1])
    generator = np.diag(edges, 1) - np.diag(edges, -1)
    eigenvalues, vectors = np.linalg.eigh(1j * generator)
    elementary = eigenvalues[4:]
    candidates = [(frequency, "single", f"l{i + 1}") for i, frequency in enumerate(elementary)]
    for indices in combinations(range(4), 3):
        for signs in product((-1, 1), repeat=2):
            signed = (1, *signs)
            frequency = np.dot(signed, elementary[list(indices)])
            expression = "".join(f"{'+' if sign > 0 else '-'}l{i + 1}"
                                 for i, sign in zip(indices, np.sign(frequency) * np.array(signed)))
            candidates.append((abs(frequency), "triple", expression.lstrip("+")))
    return sorted(candidates), eigenvalues, vectors


def exact_gap_expansion(model, state):
    scaled = 20 * model.hamiltonian
    integer = np.rint(scaled).astype(int)
    assert np.max(np.abs(scaled - integer)) < 1e-13
    diagonal, flipped = integer[:8, :8], integer[:8, :7:-1]
    plus_block, minus_block = diagonal + flipped, diagonal - flipped
    plus_energy, plus_vectors = np.linalg.eigh(plus_block)
    minus_energy, minus_vectors = np.linalg.eigh(minus_block)
    gap_ids, audit = certify_gap_groups(plus_block, minus_block, plus_energy, minus_energy)
    plus = np.vstack((plus_vectors, plus_vectors[::-1])) / np.sqrt(2)
    minus = np.vstack((minus_vectors, -minus_vectors[::-1])) / np.sqrt(2)

    reflection = np.array([int(f"{j:04b}"[::-1], 2) for j in range(16)])
    parities = [np.sum(v * v[reflection], axis=0) for v in (plus, minus)]
    assert max(np.max(np.abs(np.abs(p) - 1)) for p in parities) < 1e-12
    even_plus, even_minus = (p > 0 for p in parities)
    initial_plus, initial_minus = state @ plus, state @ minus
    reflected_error = max(np.max(np.abs(initial_plus[~even_plus])), np.max(np.abs(initial_minus[~even_minus])))
    assert reflected_error < 1e-12
    matrix = np.einsum("sa,si,sb->abi", plus, OCCUPATION, minus)
    signed_coeffs = 2 * initial_plus[:, None, None].conj() * initial_minus[None, :, None] * matrix
    active = even_plus[:, None] & even_minus[None, :]
    signed_coeffs[~active] = 0
    positive_coeffs = np.where(gap_ids[..., None] > 0, signed_coeffs, signed_coeffs.conj())
    index = np.abs(gap_ids).ravel() - 1
    counts = np.bincount(index)
    active_counts = np.bincount(index, weights=active.ravel()).astype(int)
    gaps = np.abs(plus_energy[:, None] - minus_energy[None, :]).ravel() / 20
    frequency = np.bincount(index, weights=gaps) / counts
    summed = np.zeros((len(counts), 4), dtype=complex)
    incoherent_l1 = np.zeros(summed.shape)
    np.add.at(summed, index, positive_coeffs.reshape(-1, 4))
    np.add.at(incoherent_l1, index, np.abs(positive_coeffs).reshape(-1, 4))

    candidates, majorana_eigenvalues, majorana_vectors = frequencies_from_majoranas()
    candidate_error = np.max(np.abs(frequency - [candidate[0] for candidate in candidates]))
    assert candidate_error < 1e-12
    triple = np.array([candidate[1] == "triple" for candidate in candidates])
    forbidden = np.ix_(triple, [0, 3])
    endpoint_forbidden_max = np.max(np.abs(summed[forbidden]))
    assert endpoint_forbidden_max < 1e-12
    summed[forbidden] = 0

    rows = []
    for a, b in np.ndindex(gap_ids.shape):
        k = abs(gap_ids[a, b]) - 1
        rows.append(dict(mode=f"m{k + 1:02}", omega=float(frequency[k]), parity_plus_index=a,
                         parity_minus_index=b, energy_plus=float(plus_energy[a] / 20),
                         energy_minus=float(minus_energy[b] / 20), positive_orientation=bool(gap_ids[a, b] > 0),
                         reflection_plus=int(np.sign(parities[0][a])), reflection_minus=int(np.sign(parities[1][b])),
                         allowed_by_initial_reflection=bool(active[a, b]),
                         **{f"{part}{i}": float(getattr(c, part)) for i, c in enumerate(positive_coeffs[a, b])
                            for part in ("real", "imag")}))

    audit.update(reflection_plus_dimension=int(even_plus.sum()), reflection_minus_dimension=int(even_minus.sum()),
                 max_initial_reflection_odd_overlap=float(reflected_error),
                 max_majorana_gap_label_error=float(candidate_error),
                 max_endpoint_triple_before_structural_zero=float(endpoint_forbidden_max),
                 min_group_spacing=float(np.min(np.diff(frequency))))
    return frequency, summed.T, incoherent_l1.T, counts, active_counts, candidates, audit, rows, majorana_eigenvalues, majorana_vectors


def endpoint_projection(angles, eigenvalues, vectors, frequencies, coefficients):
    theta, phi = angles.T
    x, y, z = np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)
    prefix = np.r_[1.0, np.cumprod(x[:-1])]
    expectation = (prefix[:, None] * np.column_stack((z, y))).ravel()
    negative = np.flatnonzero(eigenvalues < 0)[::-1]
    elementary = -eigenvalues[negative]
    overlap = vectors[:, negative].conj().T @ expectation
    projected = -vectors[0, negative] * overlap
    index = np.argmin(np.abs(frequencies[None, :] - elementary[:, None]), axis=1)
    error = float(np.max(np.abs(coefficients[index] - projected)))
    assert error < 1e-12
    return dict(expectations=expectation.tolist(), overlap_power_fraction=(np.abs(overlap)**2 / np.sum(np.abs(overlap)**2)).tolist(),
                observed_weights=np.abs(vectors[0, negative]).tolist(), elementary_frequencies=elementary.tolist(),
                max_endpoint_coefficient_error=error)


def first_order_data(frequency, relative, carrier):
    offsets = frequency - frequency[carrier]
    folded = []
    for index in np.flatnonzero(np.arange(len(frequency)) != carrier):
        value = relative[index] if offsets[index] > 0 else -relative[index].conjugate()
        group = next((g for g in folded if abs(g["offset"] - abs(offsets[index])) < 1e-10), None)
        if group is None:
            folded.append(dict(offset=float(abs(offsets[index])), value=value, mode_indices=[index]))
        else:
            group["value"] += value
            group["mode_indices"].append(index)
    power = sum(abs(group["value"])**2 / 2 for group in folded)
    return power, [dict(offset=g["offset"], modes=[f"m{i + 1:02}" for i in g["mode_indices"]],
                        power_fraction=float(abs(g["value"])**2 / (2 * power)))
                   for g in sorted(folded, key=lambda g: abs(g["value"]), reverse=True)]


def ablations(windows, frequency, coefficients, carrier, triple, relative):
    grids = np.array([time_grid(bounds) for bounds in windows.values()])
    offsets = frequency - frequency[carrier]
    kernel = np.exp(1j * grids.ravel()[:, None] * offsets)
    weights = np.ones(grids.shape[1])
    weights[[0, -1]] = 0.5
    ranked = [i for i in np.argsort(-np.abs(relative)) if i != carrier]
    full = np.arange(len(frequency))
    variants = [("full", None, full), ("dominant_only", None, [carrier]),
                ("all_singles", None, np.flatnonzero(~triple)),
                ("dominant_plus_triples", None, np.r_[carrier, np.flatnonzero(triple)])]
    for position, index in enumerate(ranked):
        variants.extend((("add_one", index, [carrier, index]), ("omit_one", index, full[full != index]),
                         ("cumulative", index, [carrier, *ranked[:position + 1]])))

    carrier_lag = np.angle(coefficients[0, carrier] * coefficients[1, carrier].conj())
    rows = []
    full_lag = None
    for case, index, keep in variants:
        mask = np.isin(full, keep)
        signals = kernel @ (coefficients[:2] * mask).T
        envelope = np.abs(signals)
        # All ablations keep the carrier and all sites remain above all three gates.
        assert envelope.min() > 0.10
        lag = np.angle(signals[:, 0] * signals[:, 1].conj()).reshape(grids.shape)
        first = carrier_lag + (kernel @ (relative * mask)).imag.reshape(grids.shape)
        locking = np.average(np.exp(1j * lag), axis=-1, weights=weights)
        mean, first_mean = (np.average(values, axis=-1, weights=weights).mean() for values in (lag, first))
        rms = lambda values: float(np.rad2deg(np.sqrt(np.average(values**2, axis=-1, weights=weights).mean())))
        if case == "full":
            full_lag = lag
        rows.append(dict(case=case, mode=f"m{index + 1:02}" if index is not None else "",
                         retained_modes=" ".join(f"m{k + 1:02}" for k in keep),
                         lag_rms_deg=rms(lag - mean), peak_abs_lag_deg=float(np.rad2deg(np.max(np.abs(lag)))),
                         min_envelope=float(envelope.min()), first_order_rms_deg=rms(first - first_mean),
                         first_order_error_rms_deg=rms(lag - first), error_vs_full_rms_deg=rms(lag - full_lag),
                         worst_window_mean_lag_deg=float(np.max(np.abs(np.angle(locking, deg=True)))),
                         worst_window_min_plv=float(np.abs(locking).min()),
                         minimum_end_envelope=float(envelope[:, 0].min()),
                         minimum_interior_envelope=float(envelope[:, 1].min())))
    return rows


def main():
    angles, state, sha = frozen_state()
    model = Model("open")
    (frequency, coefficients, incoherent_l1, counts, active_counts, candidates, audit, transitions,
     eigenvalues, vectors) = exact_gap_expansion(model, state)
    triple = np.array([candidate[1] == "triple" for candidate in candidates])
    carrier = int(np.argmax(np.abs(coefficients).min(axis=0)))
    dominant = coefficients[:, carrier]
    normalized = coefficients / dominant[:, None]
    normalized[:, carrier] = 0
    relative = normalized[0] - normalized[1]
    modes = [f"m{k + 1:02}" for k in range(len(frequency))]
    fundamental = np.arange(len(frequency)) == carrier

    rows = []
    for k, (omega, (_, kind, expression)) in enumerate(zip(frequency, candidates)):
        row = dict(mode=modes[k], omega=float(omega), is_dominant=bool(fundamental[k]), frequency_class=kind,
                   gap_expression=expression, transition_count=int(counts[k]),
                   reflection_allowed_transitions=int(active_counts[k]))
        for i, c in enumerate(coefficients[:, k]):
            row.update({f"re{i}": float(c.real), f"im{i}": float(c.imag), f"abs{i}": float(abs(c)),
                        f"phase{i}_deg": float(np.angle(c, deg=True)) if abs(c) > 1e-14 else None,
                        f"power_fraction{i}": float(abs(c)**2 / np.sum(np.abs(coefficients[i])**2)),
                        f"ungrouped_l1_{i}": float(incoherent_l1[i, k])})
        row.update(normalized_diff_re=float(relative[k].real), normalized_diff_im=float(relative[k].imag),
                   normalized_diff_abs=float(abs(relative[k])),
                   first_order_rms_deg=float(np.rad2deg(abs(relative[k]) / np.sqrt(2))))
        rows.append(row)

    previous = json.loads((PRIOR / "results.json").read_text(encoding="utf-8"))
    windows = {name: data["bounds"] for name, data in previous["windows"].items() if data["role"] == "holdout"}
    variants = ablations(windows, frequency, coefficients, carrier, triple, relative)
    previous_main = previous["cases"]["selected_0.05"]
    plv_error = abs(variants[0]["worst_window_min_plv"] - min(previous_main[w]["0.05"]["min_plv"] for w in windows))
    lag_error = abs(variants[0]["worst_window_mean_lag_deg"] - max(previous_main[w]["0.05"]["max_abs_lag"] for w in windows))
    t = np.r_[0, 3.688, np.linspace(64 * np.pi, 256 * np.pi, 1001)]
    reconstructed = 0.5 + (np.exp(1j * t[:, None] * frequency) @ coefficients.T).real
    occupation_error = float(np.max(np.abs(reconstructed - model.direct(state, t))))
    assert occupation_error < 1e-10 and plv_error < 1e-10 and lag_error < 1e-7
    audit.update(occupation_reconstruction_max_error=occupation_error, previous_worst_window_plv_error=plv_error,
                 previous_worst_window_lag_error_deg=lag_error,
                 reflection_coefficient_error=float(np.max(np.abs(coefficients - coefficients[::-1]))))

    projection = endpoint_projection(angles, eigenvalues, vectors, frequency, coefficients[0])
    power, folded = first_order_data(frequency, relative, carrier)
    magnitude = np.abs(coefficients)
    residual = magnitude[:, ~fundamental].sum(axis=-1)
    individual_angle_bounds = np.arcsin(residual / np.abs(dominant))
    phase_offset = float(np.angle(dominant[0] * dominant[1].conj(), deg=True))
    summary = {
        "selection_sha256": sha, "angles_rad": angles.tolist(), "dc": [0.5] * 4,
        "exact_gap_audit": audit, "endpoint_projection": projection,
        "dominant": {"mode": modes[carrier], "omega": float(frequency[carrier]), "phases_deg": np.angle(dominant, deg=True).tolist(),
                     "magnitudes": np.abs(dominant).tolist(), "end_minus_interior_deg": phase_offset,
                     "interior_over_end_amplitude": float(abs(dominant[1] / dominant[0]))},
        "site_summaries": [{"site": i, "dominant_power_fraction": float(rows[carrier][f"power_fraction{i}"]),
                            "residual_l1": float(residual[i]), "residual_to_dominant_l1": float(residual[i] / abs(dominant[i])),
                            "triple_l1": float(magnitude[i, triple].sum()),
                            "triple_power_fraction": float(np.sum(magnitude[i, triple]**2) / np.sum(magnitude[i]**2)),
                            "independent_envelope_lower": float(abs(dominant[i]) - residual[i]),
                            "nonzero_lines": int(np.count_nonzero(magnitude[i] > 1e-14))} for i in range(4)],
        "interference": {"window_count": len(windows), "window_bounds": windows,
                         "independent_residual_lag_bound_deg": abs(phase_offset) + float(np.rad2deg(individual_angle_bounds[:2].sum())),
                         "first_order_infinite_rms_deg": float(np.rad2deg(np.sqrt(power))),
                         "first_order_beat_groups": folded,
                         "first_order_triple_l1": float(np.abs(relative[triple]).sum()),
                         "first_order_single_l1": float(np.abs(relative[~triple]).sum()),
                         "first_order_triple_variance_fraction": float(sum(g["power_fraction"] for g in folded
                             if all(triple[int(m[1:]) - 1] for m in g["modes"]))),
                         "notable_ablations": [row for row in variants if not row["mode"]]},
        "relative_residual_bound": spectral_bound(model, state),
    }
    save_csv("spectrum.csv", rows)
    save_csv("transitions.csv", transitions)
    save_csv("mode_ablation.csv", variants)
    (OUTPUT / "spectral.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"dominant": summary["dominant"], "site_summaries": summary["site_summaries"][:2],
                      "ablations": summary["interference"]["notable_ablations"], "audit_error": occupation_error}, indent=2))


if __name__ == "__main__":
    main()
