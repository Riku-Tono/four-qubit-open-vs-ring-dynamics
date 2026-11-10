"""Evaluate the saved short-window rotation axes on later time intervals."""

import csv
import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
SOURCE = OUTPUT.parent / "antiphase-axis-optimization-2026-11-10"
OPTIMIZATION = run_path(str(SOURCE / "run.py"))
XY, PREVIOUS = OPTIMIZATION["XY"], OPTIMIZATION["PREVIOUS"]
A, B, OCCUPATION = XY["A"], XY["B"], XY["OCCUPATION"]
COUPLING = 0.7
WINDOWS = {
    "train": XY["WINDOWS"]["short"],
    "holdout_a": (XY["WINDOWS"]["short"][1], 8 * np.pi),
    "holdout_b": (8 * np.pi, 16 * np.pi),
}
AXES = ("fixed", "X", "Y")
THETA = np.array([0.0, 0.10, -0.10, 0.20, -0.20])
ANGLE_PAIRS = ((0, 0), (1, 2), (3, 4))
BETAS = OPTIMIZATION["BETA_DEG"]
unit_axis = OPTIMIZATION["unit_axis"]
slope_rms = OPTIMIZATION["slope_rms"]
trap_weights = OPTIMIZATION["trap_weights"]


def load_training_axes():
    with (SOURCE / "optima.csv").open(encoding="utf-8", newline="") as handle:
        rows = {row["topology"]: row for row in csv.DictReader(handle)
                if row["window"] == "short" and float(row["J"]) == COUPLING}
    assert set(rows) == {name for name, _, _ in XY["TOPOLOGIES"]}
    assert all(float(row["T"]) == WINDOWS["train"][1] for row in rows.values())
    return rows


def average_kernel(frequency, bounds):
    start, end = bounds
    return np.exp(1j * frequency * (start + end) / 2) * np.sinc(
        frequency * (end - start) / (2 * np.pi)
    )


def average_product(left, right, kernel):
    # Both Fourier expansions use +i omega t and represent real waveforms.
    return np.einsum("...i,...i->...", left @ kernel, right).real


def residual_terms(prepared, other, directions, kernel):
    _, vectors, _, occupation = prepared
    diagonal = occupation[:, :, 0] + occupation[:, :, other] - np.eye(len(vectors))
    cross = vectors[A, :, None] * diagonal * vectors[B, None, :]
    slopes = np.stack(((cross - cross.T) / (2j), -(cross + cross.T) / 2))
    slopes = slopes.reshape(2, -1)
    gram = (slopes @ kernel @ slopes.T).real
    gram = (gram + gram.T) / 2
    even = (vectors[B, :, None] * diagonal * vectors[B, None, :]).ravel()
    even_square = average_product(even, even, kernel)
    even_slope = average_product(directions @ slopes, even, kernel)
    return gram, even_square, even_slope


def compare_with_pure_axes(rows, metrics):
    for key, suffix in metrics:
        x, y = rows[1][key], rows[2][key]
        best = min(x, y)
        for row in rows:
            row["best_pure_axis" + suffix] = ("X" if x <= y else "Y") if best else ""
            row["improvement_fraction" + suffix] = (best - row[key]) / best if best else None


def residual_scores(tag, betas, directions, gram, even_square, even_slope):
    linear_rms = slope_rms(gram, directions)
    linear = [{**tag, "axis": axis, "beta_deg": beta, "rms": rms}
              for axis, beta, rms in zip(AXES, betas, linear_rms)]
    compare_with_pure_axes(linear, (("rms", ""),))
    finite, signed_square = [], np.zeros((len(AXES), len(THETA)))
    for plus, minus in ANGLE_PAIRS:
        angle = THETA[plus]
        s2, sin_theta = np.sin(angle / 2)**2, np.sin(angle)
        even_sq = s2**2 * even_square
        odd_sq = (sin_theta * linear_rms)**2
        pooled_sq = even_sq + odd_sq
        cross = 2 * s2 * sin_theta * even_slope
        signed_square[:, plus] = pooled_sq + cross
        signed_square[:, minus] = pooled_sq - cross
        rows = []
        for i, (axis, beta) in enumerate(zip(AXES, betas)):
            rows.append(dict(
                **tag, axis=axis, beta_deg=beta, abs_theta_rad=angle,
                total_rms_plus=np.sqrt(signed_square[i, plus]),
                total_rms_minus=np.sqrt(signed_square[i, minus]),
                total_rms_pooled=np.sqrt(pooled_sq[i]), even_rms=np.sqrt(even_sq),
                odd_rms=np.sqrt(odd_sq[i]), twice_even_odd_inner_plus=cross[i],
            ))
        compare_with_pure_axes(rows, tuple(
            ("total_rms_" + sign, "_" + sign) for sign in ("plus", "minus", "pooled")
        ))
        finite.extend(rows)
    return linear, finite, signed_square


def occupation_scores(tag, prepared, other, bounds, betas, directions, first, kernel, gram, signed_square):
    energy, vectors, _, matrix = prepared
    half_angles = THETA[None, :, None] / 2
    phase = -directions[:, 1] - 1j * directions[:, 0]
    projection = (np.cos(half_angles) * vectors[A] +
                  phase[:, None, None] * np.sin(half_angles) * vectors[B])
    weighted = projection.conj()[..., :, None] * projection[..., None, :]
    coeff = (weighted[..., None] * matrix).reshape(len(AXES), len(THETA), -1, 4).swapaxes(-2, -1)
    mean = (coeff @ first).real
    variance = average_product(coeff, coeff, kernel) - mean**2
    std = np.sqrt(variance)
    covariance = average_product(coeff[..., 0, :], coeff[..., other, :], kernel)
    covariance -= mean[..., 0] * mean[..., other]
    std_product = std[..., 0] * std[..., other]
    correlation = np.clip(covariance / std_product, -1, 1)
    bias = mean[..., 0] + mean[..., other] - 1
    mismatch = (std[..., 0] - std[..., other])**2
    defect = 2 * std_product * (1 + correlation)

    coarse = PREVIOUS["time_grid"](bounds, XY["DT"])
    time = np.linspace(*bounds, 2 * len(coarse) - 1)
    evolution = np.exp(-1j * time[:, None] * energy)
    states = (projection[:, :, None, :] * evolution) @ vectors.T
    b_state = (evolution * vectors[B]) @ vectors.T
    diagonal = OCCUPATION[:, 0] + OCCUPATION[:, other] - 1
    cross = (states[0, 0].conj() * b_state) @ diagonal
    slopes = np.column_stack((cross.imag, -cross.real))
    direct = np.abs(states)**2 @ OCCUPATION
    del states
    minimum, maximum = np.min(direct, axis=2), np.max(direct, axis=2)
    half_range = (maximum - minimum) / 2
    coarse_range = np.ptp(direct[:, :, ::2], axis=2) / 2
    weights = trap_weights(len(time))
    direct_mean = np.einsum("ijst,s->ijt", direct, weights)
    direct_variance = np.einsum("ijst,s->ijt", (direct - direct_mean[:, :, None])**2, weights)
    delta = direct[..., 0] + direct[..., other] - 1
    checks = dict(
        time_grid_points=len(time), time_grid_dt=time[1] - time[0],
        gram_vs_trapezoid=np.max(np.abs(gram - slopes.T @ (weights[:, None] * slopes))),
        mean_vs_trapezoid=np.max(np.abs(mean - direct_mean)),
        std_vs_trapezoid=np.max(np.abs(std - np.sqrt(direct_variance))),
        residual_rms_vs_trapezoid=np.max(np.abs(np.sqrt(signed_square) - np.sqrt(delta**2 @ weights))),
        zero_angle_max_abs=np.max(np.abs(delta[:, 0])),
        pair_mse_reconstruction=np.max(np.abs(bias**2 + mismatch + defect - signed_square)),
        half_range_change_under_dt_halving=np.max(np.abs(half_range - coarse_range)),
    )
    for key in ("gram_vs_trapezoid", "mean_vs_trapezoid", "std_vs_trapezoid", "residual_rms_vs_trapezoid",
                "half_range_change_under_dt_halving"):
        assert checks[key] < 1e-6, (tag, key, checks[key])
    assert checks["zero_angle_max_abs"] < 1e-11 and checks["pair_mse_reconstruction"] < 1e-12

    sites, pairs = [], []
    for i, (axis, beta) in enumerate(zip(AXES, betas)):
        for j, theta in enumerate(THETA):
            condition = dict(**tag, axis=axis, beta_deg=beta, theta_rad=theta)
            sites.extend(dict(
                **condition, site=site, mean=mean[i, j, site], time_std=std[i, j, site],
                min_sampled=minimum[i, j, site], max_sampled=maximum[i, j, site],
                half_range=half_range[i, j, site],
            ) for site in range(4))
            pairs.append(dict(
                **condition, partner_site=other, mean_residual=bias[i, j], bias_sq=bias[i, j]**2,
                std_0=std[i, j, 0], std_k=std[i, j, other], covariance=covariance[i, j],
                correlation=correlation[i, j], std_mismatch_sq=mismatch[i, j],
                anticorrelation_defect=defect[i, j], residual_mse=signed_square[i, j],
            ))
    return sites, pairs, checks


def posthoc_diagnostics(tag, gram, fixed_beta, fixed_rms, even_square):
    """Score the diagnostic lower bound without changing the evaluation axes."""
    best = OPTIMIZATION["optimum"](gram)
    posthoc_beta, posthoc_rms = best["beta_opt_deg"], best["rms_min"]
    regret = fixed_rms - posthoc_rms
    assert regret >= -1e-13
    result = dict(
        **tag, beta_fixed_deg=fixed_beta, beta_posthoc_deg=posthoc_beta,
        rms_fixed=fixed_rms, rms_posthoc=posthoc_rms, regret=max(regret, 0.0),
        angle_distance_deg=abs((fixed_beta - posthoc_beta + 90) % 180 - 90),
        gram_xx=gram[0, 0], gram_xy=gram[0, 1], gram_yy=gram[1, 1],
    )
    for theta, label in ((0.10, "0_10"), (0.20, "0_20")):
        even_sq = np.sin(theta / 2)**4 * even_square
        best_pooled = np.hypot(np.sqrt(even_sq), np.sin(theta) * posthoc_rms)
        fixed_pooled = np.hypot(np.sqrt(even_sq), np.sin(theta) * fixed_rms)
        result["pooled_regret_theta_" + label] = max(fixed_pooled - best_pooled, 0.0)
    return result


def main():
    training = load_training_axes()
    names = ("linear", "finite", "sites", "pair_moments", "regret", "sweep")
    tables, audits = {name: [] for name in names}, []
    scan_directions = unit_axis(BETAS)
    scan_directions[[0, len(BETAS) // 2, -1]] = ((1, 0), (0, 1), (-1, 0))
    for topology, link, other in XY["TOPOLOGIES"]:
        fixed_beta = float(training[topology]["beta_opt_deg"])
        betas = (fixed_beta, 0.0, 90.0)
        directions = np.vstack((unit_axis(fixed_beta), [1, 0], [0, 1]))
        prepared = PREVIOUS["prepare"](COUPLING, link)
        frequencies = prepared[2].ravel()
        for window, bounds in WINDOWS.items():
            tag = dict(topology=topology, window=window, t_start=bounds[0], t_end=bounds[1])
            first = average_kernel(frequencies, bounds)
            kernel = average_kernel(frequencies[:, None] + frequencies[None, :], bounds)
            gram, even_square, even_slope = residual_terms(prepared, other, directions, kernel)
            linear, finite, signed_square = residual_scores(tag, betas, directions, gram, even_square, even_slope)
            sites, pairs, checks = occupation_scores(
                tag, prepared, other, bounds, betas, directions, first, kernel, gram, signed_square
            )
            for name, records in zip(names[:4], (linear, finite, sites, pairs)):
                tables[name].extend(records)
            audits.append({**tag, **checks})
            fixed_rms = linear[0]["rms"]
            if window == "train":
                assert abs(fixed_rms - float(training[topology]["rms_min"])) < 1e-13
            tables["regret"].append(posthoc_diagnostics(tag, gram, fixed_beta, fixed_rms, even_square))
            if window != "train":
                tables["sweep"].extend(dict(topology=topology, window=window, beta_deg=beta, rms=rms)
                                       for beta, rms in zip(BETAS, slope_rms(gram, scan_directions)))
    for name, rows in tables.items():
        with (OUTPUT / (name + ".csv")).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    metadata = dict(
        omega=XY["OMEGA"], J=COUPLING, initial="1000", theta_rad=THETA.tolist(), windows=WINDOWS,
        fixed_axis_source="../antiphase-axis-optimization-2026-11-10/optima.csv (J=0.7, short only)",
        fixed_beta_deg={name: float(row["beta_opt_deg"]) for name, row in training.items()},
        sweep_step_deg=BETAS[1] - BETAS[0], checks=audits,
    )
    (OUTPUT / "audit.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    for row in tables["linear"]:
        if row["axis"] == "fixed":
            print(f"{row['topology']} {row['window']}: fixed RMS {row['rms']:.9f}, "
                  f"improvement {row['improvement_fraction']:+.4%}")


if __name__ == "__main__":
    main()
