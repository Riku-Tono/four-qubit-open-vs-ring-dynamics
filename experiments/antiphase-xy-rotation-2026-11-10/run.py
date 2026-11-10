"""Compare local X and Y rotations with the previous spectral phase convention."""

import csv
import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
PREVIOUS = run_path(str(OUTPUT.parent / "antiphase-y-rotation-2026-11-10" / "run.py"))
TOPOLOGIES, COUPLINGS = PREVIOUS["TOPOLOGIES"], PREVIOUS["COUPLINGS"]
OCCUPATION = PREVIOUS["OCCUPATION"]
OMEGA = PREVIOUS["REFERENCE"]["OMEGA"]
WINDOWS = {name: PREVIOUS["WINDOWS"][name] for name in ("short", "long")}
DT, GATE, ZERO_GAP = (PREVIOUS[name] for name in ("DT", "GATE", "ZERO_GAP"))
AXES = ("X", "Y")
FACTORS = np.array([-1j, -1])  # exp(i phi) for X and Y respectively
MAGNITUDES = PREVIOUS["ANGLES"][1:]
ANGLES = np.array([0.0, *(v for a in MAGNITUDES for v in (a, -a))])
PLUS = np.arange(1, len(ANGLES), 2)
MINUS = PLUS + 1
A, B = 0b1000, 0b0000
TEST_PHASE = np.pi / 4


def max_abs(values):
    return float(np.max(np.abs(values)))


def sample(prepared, time, other):
    energy, vectors, gaps, matrix = prepared
    half = ANGLES[None, :, None] / 2
    projections = (np.cos(half) * vectors[A] +
                   FACTORS[:, None, None] * np.sin(half) * vectors[B])
    rotation_count = len(AXES) * len(ANGLES)
    test_angle = ANGLES[-2:]
    test_projections = (np.cos(test_angle[:, None] / 2) * vectors[A] +
                        np.exp(1j * TEST_PHASE) * np.sin(test_angle[:, None] / 2) * vectors[B])
    projections = np.vstack((projections.reshape(rotation_count, -1), test_projections))

    weighted = (matrix[None, ...] * projections.conj()[:, :, None, None] *
                projections[:, None, :, None])
    dc = weighted[:, np.abs(gaps) <= ZERO_GAP].sum(axis=1).real
    positive = gaps > ZERO_GAP
    frequencies = gaps[positive]
    coefficients = 2 * weighted[:, positive].transpose(1, 0, 2).reshape(len(frequencies), -1)
    signal = np.exp(1j * time[:, None] * frequencies) @ coefficients
    signal = signal.reshape(len(time), len(projections), 4).transpose(1, 0, 2)
    spectral = dc[:, None, :] + signal.real

    evolution = np.exp(-1j * time[:, None] * energy)
    states = (evolution[None, ...] * projections[:, None, :]) @ vectors.T
    probabilities = np.abs(states)**2
    direct = probabilities @ OCCUPATION
    a, b = (evolution[None, ...] * vectors[[A, B], None, :]) @ vectors.T
    s_diagonal = -2 * (OCCUPATION[:, 0] + OCCUPATION[:, other] - 1)
    s_a = np.abs(a)**2 @ s_diagonal
    s_b = np.abs(b)**2 @ s_diagonal
    cross = (a.conj() * b) @ s_diagonal
    test_residual = direct[-2:, :, 0] + direct[-2:, :, other] - 1
    test_formula = -0.5 * (
        np.cos(test_angle[:, None] / 2)**2 * s_a
        + np.sin(test_angle[:, None] / 2)**2 * s_b
        + np.sin(test_angle[:, None]) * (np.exp(1j * TEST_PHASE) * cross).real
    )
    errors = {
        "occupation_spectral_vs_state": max_abs(spectral - direct),
        "site_dc_half": max_abs(dc - 0.5),
        "norm_squared": max_abs(probabilities.sum(axis=2) - 1),
        "s_a_exact_zero": max_abs(s_a),
        "general_phi_pi_over_4": max_abs(test_residual - test_formula),
    }
    shape = (len(AXES), len(ANGLES), len(time), 4)
    return (direct[:rotation_count].reshape(shape), signal[:rotation_count].reshape(shape),
            s_b, cross, errors)


def parity(delta):
    return (delta[:, PLUS] + delta[:, MINUS]) / 2, (delta[:, PLUS] - delta[:, MINUS]) / 2


def audit_formulas(errors, topology, coupling, other, time, occupation, delta, s_b, cross):
    even, odd = parity(delta)
    half_squared = np.sin(ANGLES[PLUS, None] / 2)**2
    sin_angle = np.sin(ANGLES[PLUS, None])
    errors.update(
        initial_residual=max_abs(delta[:, :, 0] + np.sin(ANGLES / 2)**2),
        even_from_diagonal=max_abs(even + half_squared * s_b / 2),
        x_odd_from_imag=max_abs(odd[0] + sin_angle * cross.imag / 2),
        y_odd_from_real=max_abs(odd[1] - sin_angle * cross.real / 2),
    )
    if topology == "open":
        errors["open_site_3_invariance"] = max_abs(occupation[..., 3] - occupation[0, 0, :, 3])
        edges = np.array([OMEGA, coupling, OMEGA, coupling, OMEGA, coupling, OMEGA])
        generator = np.diag(edges, 1) - np.diag(edges, -1)
        energy, vectors = np.linalg.eigh(1j * generator)
        products = vectors[0] * vectors[[0, 1, 2]].conj()
        f00, f01, f02 = (np.exp(-1j * time[:, None] * energy) @ products.T).real.T
        errors["open_majorana_s_b"] = max_abs(s_b - 2 * f00)
        errors["open_majorana_cross"] = max_abs(cross - (f02 + 1j * f01))
    else:
        errors["ring_sites_1_3_reflection"] = max_abs(occupation[..., 1] - occupation[..., 3])
    if coupling == 0:
        wt, theta = OMEGA * time, ANGLES[:, None]
        exact_site0 = np.stack([(1 + np.cos(wt + theta)) / 2,
                                (1 + np.cos(theta) * np.cos(wt)) / 2])
        errors["uncoupled_site_0"] = max_abs(occupation[..., 0] - exact_site0)
        errors["uncoupled_spectator_sites"] = max_abs(
            occupation[..., 1:] - (1 - np.cos(wt[:, None])) / 2
        )
        errors["uncoupled_y_odd"] = max_abs(odd[1])
        errors["uncoupled_cross"] = max_abs(cross - 1j * np.sin(wt))
    assert max(errors.values()) < 1e-11, (topology, coupling, other, errors)


def phase_metrics(signal, weights, pair, other, mask):
    stats = PREVIOUS["phase_stats"](signal, weights, mask=mask)
    active = np.all(np.abs(signal) >= GATE, axis=1) & mask
    envelope = np.average(np.abs(signal[active]), axis=0, weights=weights[active])
    return {"coverage": stats["coverage"], **stats["pairs"][pair],
            "mean_envelope_0": envelope[0], "mean_envelope_k": envelope[other]}


def summarize(time, occupation, signal, other, s_b, cross):
    weights = np.ones(len(time))
    weights[[0, -1]] = 0.5
    delta = occupation[..., 0] + occupation[..., other] - 1
    even, odd = parity(delta)
    pair = f"0{other}"
    valid = np.all(np.abs(signal) >= GATE, axis=-1)
    angle_index = {angle: i for i, angle in enumerate(ANGLES)}
    residuals, components, phases = [], [], []
    for axis, name in enumerate(AXES):
        for i, angle in enumerate(ANGLES):
            values = delta[axis, i]
            left, right = np.argmin(values), np.argmax(values)
            amplitude = np.ptp(occupation[axis, i], axis=0) / 2
            residuals.append(dict(
                axis=name, theta_rad=angle, initial=values[0], min=values[left],
                t_min=time[left], max=values[right], t_max=time[right],
                max_abs=max_abs(values), mean=np.average(values, weights=weights),
                rms=np.sqrt(np.average(values**2, weights=weights)),
                half_range_0=amplitude[0], half_range_k=amplitude[other],
            ))
            masks = {"native": np.ones(len(time), dtype=bool),
                     "axis_common": valid[1 - axis, i],
                     "sign_common": valid[axis, angle_index[-angle]]}
            for support, mask in masks.items():
                phases.append(dict(axis=name, theta_rad=angle, support=support,
                                   **phase_metrics(signal[axis, i], weights, pair, other, mask)))
        odd_slope = -cross.imag / 2 if name == "X" else cross.real / 2
        for index, angle in enumerate(MAGNITUDES):
            components.append(dict(
                axis=name, abs_theta_rad=angle, even_max_abs=max_abs(even[axis, index]),
                odd_max_abs=max_abs(odd[axis, index]),
                odd_over_sin_max_abs=max_abs(odd[axis, index] / np.sin(angle)),
                even_formula_error=max_abs(even[axis, index] + np.sin(angle / 2)**2 * s_b / 2),
                odd_formula_error=max_abs(odd[axis, index] - np.sin(angle) * odd_slope),
            ))
    return delta, (residuals, components, phases)


def audit_sampling(changes, coarse, fine):
    for coarse_rows, fine_rows, keys in zip(coarse, fine, (
        ("max_abs", "mean", "rms", "half_range_0", "half_range_k"),
        ("even_max_abs", "odd_max_abs"),
        ("coverage", "lag_deg", "mean_abs_lag_deg", "plv", "mean_envelope_0", "mean_envelope_k"),
    )):
        for left, right in zip(coarse_rows, fine_rows):
            for key in keys:
                diff = left[key] - right[key]
                if key == "lag_deg":
                    diff = (diff + 180) % 360 - 180
                changes[key] = max(changes.get(key, 0.0), abs(float(diff)))


def write_waveforms(writer, tag, time, s_b, cross, delta):
    even, odd = parity(delta)
    selected = np.r_[np.arange(0, len(time) - 1, 5), len(time) - 1]
    for i in selected:
        values = [time[i], s_b[i], cross[i].real, cross[i].imag, delta[0, 0, i]]
        for axis in range(len(AXES)):
            for index, (plus, minus) in enumerate(zip(PLUS, MINUS)):
                values.extend((delta[axis, plus, i], delta[axis, minus, i],
                               even[axis, index, i], odd[axis, index, i]))
        writer.writerow([*tag.values(), *(f"{v:.12g}" for v in values)])


def write_summary(name, rows):
    with (OUTPUT / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    rows = ([], [], [])
    audit = dict(omega=OMEGA, dt_target=DT, phase_gate_all_sites=GATE, zero_gap_tolerance=ZERO_GAP,
                 windows=WINDOWS, max_formula_errors={}, max_change_under_dt_halving={}, cases=[])
    with (OUTPUT / "waveforms.csv").open("w", newline="", encoding="utf-8") as handle:
        waveforms = csv.writer(handle)
        waveforms.writerow(["topology", "J", "window", "t", "S_B", "C_real", "C_imag", "delta_0"] + [
            f"{axis}_{name}_{angle}" for axis in AXES for angle in MAGNITUDES
            for name in ("plus", "minus", "even", "odd")
        ])
        for topology, link, other in TOPOLOGIES:
            for coupling in COUPLINGS:
                prepared = PREVIOUS["prepare"](coupling, link)
                for window, bounds in WINDOWS.items():
                    tag = dict(topology=topology, J=coupling, window=window)
                    time = PREVIOUS["time_grid"](bounds, DT)
                    fine_time = np.linspace(*bounds, 2 * len(time) - 1)
                    occupation, signal, s_b, cross, errors = sample(prepared, fine_time, other)
                    delta_fine, fine = summarize(fine_time, occupation, signal, other, s_b, cross)
                    audit_formulas(errors, topology, coupling, other, fine_time,
                                   occupation, delta_fine, s_b, cross)
                    delta, coarse = summarize(time, occupation[:, :, ::2], signal[:, :, ::2],
                                              other, s_b[::2], cross[::2])
                    audit_sampling(audit["max_change_under_dt_halving"], coarse, fine)
                    for name, value in errors.items():
                        maxima = audit["max_formula_errors"]
                        maxima[name] = max(maxima.get(name, 0.0), value)
                    for table, entries in zip(rows, coarse):
                        table.extend({**tag, "pair": f"0{other}", **entry} for entry in entries)
                    audit["cases"].append(dict(**tag, samples=len(time), dt=time[1] - time[0],
                                               max_abs_C_real=max_abs(cross[::2].real),
                                               max_abs_C_imag=max_abs(cross[::2].imag)))
                    write_waveforms(waveforms, tag, time, s_b[::2], cross[::2], delta)
    for name, table in zip(("residuals.csv", "parity.csv", "phases.csv"), rows):
        write_summary(name, table)
    (OUTPUT / "audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: audit[key] for key in ("max_formula_errors", "max_change_under_dt_halving")}, indent=2))


if __name__ == "__main__":
    main()
