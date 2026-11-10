"""Minimize the window RMS of the signed residual's angular slope."""

import csv
import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
XY = run_path(str(OUTPUT.parent / "antiphase-xy-rotation-2026-11-10" / "run.py"))
PREVIOUS = XY["PREVIOUS"]
A, B = XY["A"], XY["B"]
OMEGA, DT = XY["OMEGA"], XY["DT"]
BETA_DEG = np.linspace(0, 180, 721)


def unit_axis(beta_deg):
    beta = np.deg2rad(beta_deg)
    return np.stack((np.cos(beta), np.sin(beta)), axis=-1)


def spectral_gram(prepared, other, end):
    """Integrate products of the two real slope waveforms on [0, end]."""
    _, vectors, gaps, occupation = prepared
    # The identity term has no A,B cross element at any time.
    pair_z = -2 * (occupation[:, :, 0] + occupation[:, :, other])
    cross = vectors[A, :, None] * pair_z * vectors[B, None, :]
    x = -(cross - cross.T) / (4j)
    y = (cross + cross.T) / 4
    coefficients = np.column_stack((x.ravel(), y.ravel()))
    omega = gaps.ravel()
    midpoint_phase = (omega[:, None] + omega[None, :]) * end / 2
    kernel = np.exp(1j * midpoint_phase) * np.sinc(midpoint_phase / np.pi)
    # q(t) is real. Both factors here use +i omega t, so there is no conjugate.
    gram = coefficients.T @ kernel @ coefficients
    assert np.max(np.abs(gram.imag)) < 1e-13
    return (gram.real + gram.real.T) / 2


def slope_rms(gram, directions):
    squared = np.einsum("...i,ij,...j->...", directions, gram, directions)
    return np.sqrt(np.maximum(squared, 0))


def optimum(gram):
    eigenvalues, vectors = np.linalg.eigh(gram)
    assert eigenvalues[0] >= 0 and eigenvalues[1] > eigenvalues[0]
    beta = float(np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0])) % 180)
    min_rms, max_rms = np.sqrt(eigenvalues)
    rms_x, rms_y = np.sqrt(np.diag(gram))
    best = min(rms_x, rms_y)
    improvement = best - min_rms
    return dict(
        beta_opt_deg=beta, phi_opt_deg=(-90 - beta) % 180, rms_min=min_rms,
        rms_x=rms_x, rms_y=rms_y, best_axis="X" if rms_x <= rms_y else "Y",
        rms_best_axis=best, improvement_abs=improvement,
        improvement_fraction=improvement / best if best else None,
        rms_max=max_rms, gram_xx=gram[0, 0], gram_xy=gram[0, 1], gram_yy=gram[1, 1],
    )


def trap_weights(count):
    weights = np.ones(count)
    weights[[0, -1]] = 0.5
    return weights / weights.sum()


def verify_window(prepared, other, bounds, gram, result):
    """Check the spectral integral and the finite-angle relation in spin space."""
    energy, vectors, *_ = prepared
    coarse_time = PREVIOUS["time_grid"](bounds, DT)
    time = np.linspace(*bounds, 2 * len(coarse_time) - 1)
    evolution = np.exp(-1j * time[:, None] * energy)
    a, b = (evolution[None, ...] * vectors[[A, B], None, :]) @ vectors.T
    occupied = XY["OCCUPATION"]
    diagonal = occupied[:, 0] + occupied[:, other] - 1
    cross = -2 * ((a.conj() * b) @ diagonal)
    slopes = np.column_stack((-cross.imag, cross.real)) / 2
    weights = trap_weights(len(time))
    sampled = slopes[::2]
    errors = dict(
        gram_vs_trapezoid=np.max(np.abs(
            gram - sampled.T @ (sampled * trap_weights(len(sampled))[:, None])
        )),
        gram_vs_half_step=np.max(np.abs(gram - slopes.T @ (slopes * weights[:, None]))),
    )
    assert errors["gram_vs_trapezoid"] < 1e-7 and errors["gram_vs_half_step"] < 2.5e-8

    even_coefficient = np.abs(b)**2 @ diagonal
    errors["finite_angle_odd_waveform"] = 0.0
    errors["finite_angle_paired_mean_square"] = 0.0
    for beta in (0.0, 90.0, result["beta_opt_deg"]):
        slope = slopes @ unit_axis(beta)
        phase = -1j * np.exp(-1j * np.deg2rad(beta))
        for theta in XY["MAGNITUDES"]:
            c, s = np.cos(theta / 2), np.sin(theta / 2)
            plus = np.abs(c * a + phase * s * b)**2 @ diagonal
            minus = np.abs(c * a - phase * s * b)**2 @ diagonal
            predicted_odd = np.sin(theta) * slope
            predicted_even = s**2 * even_coefficient
            odd_error = np.max(np.abs((plus - minus) / 2 - predicted_odd))
            mean_square_error = abs(weights @ (
                (plus**2 + minus**2) / 2 - predicted_even**2 - predicted_odd**2
            ))
            errors["finite_angle_odd_waveform"] = max(errors["finite_angle_odd_waveform"], odd_error)
            errors["finite_angle_paired_mean_square"] = max(
                errors["finite_angle_paired_mean_square"], mean_square_error
            )
    assert max(errors[key] for key in ("finite_angle_odd_waveform", "finite_angle_paired_mean_square")) < 1e-12
    return errors


def write_csv(name, rows):
    with (OUTPUT / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    summary, sweep, audits = [], [], []
    axes = unit_axis(BETA_DEG)
    axes[[0, len(axes) // 2, -1]] = ((1, 0), (0, 1), (-1, 0))
    for topology, link, other in XY["TOPOLOGIES"]:
        for coupling in XY["COUPLINGS"]:
            prepared = PREVIOUS["prepare"](coupling, link)
            for window, bounds in XY["WINDOWS"].items():
                tag = dict(topology=topology, J=coupling, window=window)
                end = bounds[1]
                gram = spectral_gram(prepared, other, end)
                checks = {}
                if coupling == 0:
                    wt = OMEGA * end
                    exact_gram = np.diag([(1 - np.sin(2 * wt) / (2 * wt)) / 8, 0])
                    checks["uncoupled_analytic_gram"] = np.max(np.abs(gram - exact_gram))
                    assert checks["uncoupled_analytic_gram"] < 1e-13
                    gram = exact_gram
                result = optimum(gram)
                curve = slope_rms(gram, axes)
                offset = (BETA_DEG[np.argmin(curve)] - result["beta_opt_deg"] + 90) % 180 - 90
                checks["sweep_grid_min_angle_deg"] = abs(offset)
                assert abs(offset) <= 0.125 + 1e-12
                checks.update(verify_window(prepared, other, bounds, gram, result))
                audits.append({**tag, **checks})
                summary.append({**tag, "T": end, **result})
                sweep.extend({**tag, "beta_deg": beta, "rms": rms} for beta, rms in zip(BETA_DEG, curve))
    write_csv("optima.csv", summary)
    write_csv("sweep.csv", sweep)
    metadata = dict(
        omega=OMEGA, initial="1000", axis_convention="beta from X toward Y in the XY plane",
        phase_convention="phi = -pi/2 - beta modulo 2*pi", windows=XY["WINDOWS"],
        beta_step_deg=BETA_DEG[1] - BETA_DEG[0], validation_dt_target=DT, checks=audits,
    )
    (OUTPUT / "audit.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    for row in summary:
        print(f"{row['topology']} J={row['J']} {row['window']}: "
              f"beta={row['beta_opt_deg']:.5f} deg, min RMS={row['rms_min']:.9f}")


if __name__ == "__main__":
    main()
