"""Compute the exact-frequency limit and finite-window slope RMS for XY axes."""

import csv
import json
from fractions import Fraction
from pathlib import Path
from runpy import run_path

import numpy as np

from spectral_groups import certify_active_support, certify_gap_groups


OUTPUT = Path(__file__).resolve().parent
HOLDOUT = run_path(str(OUTPUT.parent / "antiphase-axis-holdout-2026-11-10" / "run.py"))
OPTIMIZATION, XY = HOLDOUT["OPTIMIZATION"], HOLDOUT["XY"]
PREVIOUS = XY["PREVIOUS"]
OCCUPATION, A, B = XY["OCCUPATION"], XY["A"], XY["B"]
COUPLING, ENERGY_SCALE = 0.7, 20
CUMULATIVE_TIMES = (8 * np.pi, 16 * np.pi, 32 * np.pi, 64 * np.pi)
LATE_WINDOW = (32 * np.pi, 64 * np.pi)
CURVE_STEP_PI = 0.125
average_kernel = HOLDOUT["average_kernel"]
unit_axis, slope_rms, optimum = (OPTIMIZATION[name] for name in ("unit_axis", "slope_rms", "optimum"))


def grouped_spectrum(topology, link, other):
    """Return physical frequencies and sine-X/cosine-Y amplitudes."""
    scaled = ENERGY_SCALE * PREVIOUS["MODEL"].hamiltonian(link, 0)
    integer_h = np.rint(scaled).astype(np.int64)
    assert np.max(np.abs(scaled - integer_h)) < 1e-13
    half = len(integer_h) // 2
    direct = integer_h[:half, :half]
    flip = integer_h[:half, :half - 1:-1]
    plus_block, minus_block = direct + flip, direct - flip
    e_plus, v_plus = np.linalg.eigh(plus_block)
    e_minus, v_minus = np.linalg.eigh(minus_block)
    groups, certificate = certify_gap_groups(plus_block, minus_block, e_plus, e_minus)

    plus = np.vstack((v_plus, v_plus[::-1])) / np.sqrt(2)
    minus = np.vstack((v_minus, -v_minus[::-1])) / np.sqrt(2)
    z_pair = (2 * (1 - OCCUPATION[:, 0] - OCCUPATION[:, other])).astype(int)
    support = certify_active_support(topology, plus_block, minus_block, groups, certificate, z_pair)
    certificate["active_support"] = support
    matrix = plus.T @ (z_pair[:, None] * minus)
    w_forward = plus[A, :, None] * matrix * minus[B, None, :]
    w_reverse = plus[B, :, None] * matrix * minus[A, None, :]
    # Coefficients at E_plus-E_minus. The reverse frequencies are conjugates.
    at_signed_gap = np.stack((-(w_forward - w_reverse) / (4j), (w_forward + w_reverse) / 4), axis=-1)
    at_positive_gap = np.where(groups[..., None] > 0, at_signed_gap, at_signed_gap.conj())
    indices = np.abs(groups).ravel() - 1
    abs_gap = np.abs(e_plus[:, None] - e_minus[None, :]).ravel() / ENERGY_SCALE
    counts = np.bincount(indices)
    frequency = np.bincount(indices, weights=abs_gap) / counts
    grouped = np.zeros((len(counts), 2), dtype=complex)
    np.add.at(grouped, indices, at_positive_gap.reshape(-1, 2))
    amplitudes = np.column_stack((-2 * grouped[:, 0].imag, 2 * grouped[:, 1].real))
    inactive = np.array(support["inactive_positive_gap_ids"], dtype=int) - 1
    assert np.max(np.abs(amplitudes[inactive])) < 1e-13
    amplitudes[inactive] = 0
    assert np.max(np.abs(grouped[:, 0].real)) == 0 and np.max(np.abs(grouped[:, 1].imag)) == 0
    assert np.all(np.diff(frequency) > 0)
    incoherent = 2 * np.sum(np.abs(at_signed_gap)**2, axis=(0, 1))
    return frequency, amplitudes, counts, certificate, incoherent


def finite_gram(frequency, amplitudes, bounds):
    """Average products of sine and cosine modes analytically on the window."""
    difference = average_kernel(frequency[:, None] - frequency[None, :], bounds)
    total = average_kernel(frequency[:, None] + frequency[None, :], bounds)
    x, y = amplitudes.T
    xx = x @ (difference.real - total.real) @ x / 2
    yy = y @ (difference.real + total.real) @ y / 2
    xy = x @ (difference.imag + total.imag) @ y / 2
    return np.array([[xx, xy], [xy, yy]])


def metrics(gram):
    entry = optimum(gram)
    low, high = np.linalg.eigvalsh(gram)
    return dict(**entry, lambda_min=low, lambda_max=high, eigengap=high - low,
                gap_over_trace=(high - low) / (high + low))


def fixed_evaluation(topology, beta, gram):
    """Evaluate only the axes chosen before opening this interval."""
    axes = ("fixed_infinite", "X", "Y")
    rms = slope_rms(gram, np.vstack((unit_axis(beta), [1, 0], [0, 1])))
    best = min(rms[1:])
    return [
        dict(topology=topology, t_start=LATE_WINDOW[0], t_end=LATE_WINDOW[1],
             axis=axis, beta_deg=angle, rms=value,
             best_pure_axis="X" if rms[1] <= rms[2] else "Y", improvement_fraction=(best - value) / best)
        for axis, angle, value in zip(axes, (beta, 0.0, 90.0), rms)
    ]


def check_against_spin_model(link, other, frequency, amplitudes, windows):
    prepared = PREVIOUS["prepare"](COUPLING, link)
    maximum = 0.0
    for bounds in windows:
        direct, _, _ = HOLDOUT["residual_terms"](prepared, other, np.eye(2), average_kernel(
            prepared[2].ravel()[:, None] + prepared[2].ravel()[None, :], bounds
        ))
        error = np.max(np.abs(direct - finite_gram(frequency, amplitudes, bounds)))
        maximum = max(maximum, error)
    assert maximum < 1e-12
    return maximum


def closest_visible_modes(frequency, amplitudes):
    """Describe the nearest resolved beats without filtering the integral."""
    visible = np.flatnonzero(np.any(amplitudes != 0, axis=1))
    separations = np.diff(frequency[visible])
    smallest = np.min(separations)
    x, y = amplitudes.T
    beats = []
    # Include ties at the precision of this diagnostic.
    for index in np.flatnonzero(separations - smallest < 1e-12):
        left, right = visible[index:index + 2]
        difference = separations[index]
        mixed_sine = (x[right] * y[left] - x[left] * y[right]) / 2
        beats.append(dict(
            low_frequency=frequency[left], high_frequency=frequency[right],
            frequency_separation=difference, beat_period=2 * np.pi / difference,
            g_xy_difference_term_at_64pi=mixed_sine * average_kernel(
                difference, (0, CUMULATIVE_TIMES[-1])
            ).imag,
        ))
    return dict(
        visible_positive_count=len(visible), nearest_pairs=beats,
    )


def write_csv(name, rows):
    with (OUTPUT / (name + ".csv")).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    limits, convergence, curve, late, spectra, audits = [], [], [], [], [], {}
    for topology, link, other in XY["TOPOLOGIES"]:
        omega, amplitudes, counts, certificate, incoherent = grouped_spectrum(topology, link, other)
        support = certificate["active_support"]
        exact_xx, exact_yy = (Fraction(*support[key]) for key in ("exact_gram_xx", "exact_gram_yy"))
        exact_diagonal = np.array([float(exact_xx), float(exact_yy)])
        diagonal_error = np.max(np.abs(np.sum(amplitudes**2, axis=0) / 2 - exact_diagonal))
        assert diagonal_error < 1e-13 and exact_xx != exact_yy
        limit_gram = np.diag(exact_diagonal)
        limit = metrics(limit_gram)
        limits.append(dict(topology=topology, **limit))
        fixed_beta = limit["beta_opt_deg"]
        assert fixed_beta == (90 if exact_yy < exact_xx else 0)

        for end in CUMULATIVE_TIMES:
            gram = finite_gram(omega, amplitudes, (0, end))
            convergence.append(dict(topology=topology, T=end, **metrics(gram),
                                    gram_distance_to_limit=np.linalg.norm(gram - limit_gram, ord=2)))
        for factor in np.arange(8, 64 + CURVE_STEP_PI, CURVE_STEP_PI):
            end = factor * np.pi
            curve.append(dict(topology=topology, T=end, **metrics(finite_gram(omega, amplitudes, (0, end)))))
        late.extend(fixed_evaluation(topology, fixed_beta, finite_gram(omega, amplitudes, LATE_WINDOW)))

        active_ids = set(support["active_positive_gap_ids"])
        for i, (w, (x, y), count) in enumerate(zip(omega, amplitudes, counts), start=1):
            spectra.append(dict(topology=topology, positive_group=i, omega=w, sine_x=x,
                                cosine_y=y, parity_gap_pairs=int(count), active=i in active_ids))
        windows = [(0, end) for end in CUMULATIVE_TIMES] + [LATE_WINDOW]
        audits[topology] = dict(
            exact_gap_certificate=certificate,
            limit_diagonal_error_vs_rational=diagonal_error,
            finite_gram_max_error_vs_full_spin=check_against_spin_model(link, other, omega, amplitudes, windows),
            incoherent_diagonal_for_comparison=incoherent.tolist(),
            closest_visible_modes=closest_visible_modes(omega, amplitudes),
        )
        print(f"{topology}: limit beta={fixed_beta:.6f} deg, RMS={limit['rms_min']:.9f}")
    for name, rows in zip(("limits", "convergence", "curve", "late_holdout", "spectral_modes"),
                          (limits, convergence, curve, late, spectra)):
        write_csv(name, rows)
    metadata = dict(
        N=4, omega=XY["OMEGA"], J=COUPLING, initial="1000", energy_scale=ENERGY_SCALE,
        cumulative_end_times=CUMULATIVE_TIMES, late_fixed_evaluation=LATE_WINDOW,
        curve_step_in_pi_units=CURVE_STEP_PI, audits=audits,
    )
    (OUTPUT / "audit.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
