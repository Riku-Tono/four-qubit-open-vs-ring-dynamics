"""Bound occupation phase and amplitude directly from the spectrum."""

from itertools import combinations

import numpy as np


COEFFICIENT_MARGIN = 1e-10


def spectral_bound(model, state):
    """Use a common dominant line to bound all-time envelopes and phase gaps.

    The formulas are exact for exact Fourier coefficients. The small coefficient
    margin allows for floating-point roundoff in the eigendecomposition and sums.
    The coefficients have not been evaluated with interval arithmetic.
    """
    initial = state @ model.vectors
    products = initial[model.a].conj() * initial[model.b]
    coefficients = np.add.reduceat(products * model.observed, model.group_starts, axis=-1)
    magnitude = np.abs(coefficients)
    residual = magnitude.sum(axis=-1, keepdims=True) - magnitude
    lower = magnitude - residual - 2 * COEFFICIENT_MARGIN
    line = int(np.argmax(lower.min(axis=0)))
    envelope_lower = lower[:, line]
    if envelope_lower.min() <= 0:
        return {"available": False, "reason": "No common line dominates all four sites."}

    frequency = float(model.frequencies[line])
    dominant = coefficients[:, line]
    dominant_magnitude = magnitude[:, line]
    dominant_lower = dominant_magnitude - COEFFICIENT_MARGIN
    radius = residual[:, line]
    relative_radius = (radius + COEFFICIENT_MARGIN) / dominant_lower
    angular_radius = np.arcsin(relative_radius)
    dominant_uncertainty = np.arcsin(COEFFICIENT_MARGIN / dominant_magnitude)

    normalized = coefficients / dominant[:, None]
    normalized[:, line] = 0
    normalized_uncertainty = COEFFICIENT_MARGIN * (1 + radius / dominant_magnitude) / dominant_lower
    pairs = {}
    for i, j in combinations(range(4), 2):
        main_lag = abs(np.angle(dominant[i] * dominant[j].conj()))
        remainder_bound = angular_radius[i] + angular_radius[j]
        # Cancellation in the relative residual can sharpen the individual bounds.
        difference = np.abs(normalized[i] - normalized[j]).sum()
        difference += normalized_uncertainty[i] + normalized_uncertainty[j]
        ratio_radius = difference / (1 - min(relative_radius[i], relative_radius[j]))
        if ratio_radius < 1:
            remainder_bound = min(remainder_bound, np.arcsin(ratio_radius))
        lag_bound = min(np.pi, main_lag + dominant_uncertainty[i] + dominant_uncertainty[j] + remainder_bound)
        pairs[f"{i}{j}"] = {"max_abs_lag_upper_deg": float(np.rad2deg(lag_bound)),
                            "plv_lower": float(max(0, np.cos(lag_bound)))}

    period = 2 * np.pi / frequency
    # The dominant real sinusoid reaches both extrema within any full period.
    return {
        "available": True,
        "frequency": frequency,
        "period": period,
        "coefficients": [{"real": float(c.real), "imag": float(c.imag)} for c in dominant],
        "residual_l1": radius.tolist(),
        "coefficient_margin": COEFFICIENT_MARGIN,
        "envelope_lower": envelope_lower.tolist(),
        "min_envelope_lower": float(envelope_lower.min()),
        "half_range_lower_one_period": envelope_lower.tolist(),
        "max_abs_lag_upper_deg": max(item["max_abs_lag_upper_deg"] for item in pairs.values()),
        "min_plv_lower": min(item["plv_lower"] for item in pairs.values()),
        "pairs": pairs,
    }
