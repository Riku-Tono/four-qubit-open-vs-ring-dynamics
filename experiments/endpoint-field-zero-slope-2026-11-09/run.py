"""Extend the endpoint-field sweep and extrapolate the linear response."""

import json
from pathlib import Path
from runpy import run_path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


OUTPUT = Path(__file__).resolve().parent
PREVIOUS = OUTPUT.parent / "endpoint-field-scaling-2026-11-09"
MODEL = run_path(str(PREVIOUS / "run.py"))
ARCHIVE = json.loads((PREVIOUS / "results.json").read_text(encoding="utf-8"))
DT = ARCHIVE["dt"]
MAGNITUDES = [MODEL["EPS"] / 2**index for index in range(5)]
FIT_WINDOWS = [(2, 1), (3, 1), (3, 2), (5, 2)]


def measure(dt):
    time = np.arange(round(MODEL["END"] / dt) + 1) * dt
    rows = []
    for magnitude in MAGNITUDES:
        positive = MODEL["seam_shares"](time, magnitude)
        negative = MODEL["seam_shares"](time, -magnitude)
        for eps, shares, flipped in ((magnitude, positive, negative), (-magnitude, negative, positive)):
            qa, qb = map(float, shares)
            delta = qa - qb
            cross_error = qa - float(flipped[1])
            assert abs(cross_error) < 1e-12
            rows.append({"epsilon": eps, "Q_A": qa, "Q_B": qb, "same_sign_delta": delta,
                         "delta_over_epsilon": delta / eps, "cross_sign_error": cross_error})
    return rows


def fit_response(rows):
    x = np.square(MAGNITUDES)
    # The two signs give the same response by symmetry. Average their roundoff.
    y = np.array([row["delta_over_epsilon"] for row in rows]).reshape(-1, 2).mean(axis=1)
    fits = []
    for points, degree in FIT_WINDOWS:
        coefficients = np.polynomial.polynomial.polyfit(x[-points:], y[-points:], degree)
        fits.append({"points_nearest_zero": points, "degree_in_epsilon_squared": degree,
                     "intercept": float(coefficients[0]), "coefficients": coefficients.tolist()})
    return y, fits


def plot(y, local_fit):
    x = np.square(MAGNITUDES)
    fig, (all_points, zoom) = plt.subplots(1, 2, figsize=(9.3, 4.2), layout="constrained")
    fig.suptitle(r"Seam-share response at small fields   ($N=4$, $\lambda=0.5$ vs $0$)")

    all_points.scatter(x[:3], y[:3], color="#225D97", label="Previous 3 magnitudes", zorder=3)
    all_points.scatter(x[3:], y[3:], color="#B44915", marker="s", label="New 2 magnitudes", zorder=4)
    all_points.axvspan(0, 1.8e-4, color="#225D97", alpha=0.07, linewidth=0)
    all_points.set(title="Full sweep", xlim=(-1.0e-4, 2.62e-3))
    all_points.legend(frameon=False, loc="lower left")

    fit_x = np.linspace(0, 1.8e-4, 100)
    coeff = local_fit["coefficients"]
    zoom.plot(fit_x, np.polynomial.polynomial.polyval(fit_x, coeff), color="#394C60",
              linestyle="--", label="Quadratic fit, smallest 3")
    zoom.scatter(x[2], y[2], color="#225D97", zorder=3)
    zoom.scatter(x[3:], y[3:], color="#B44915", marker="s", zorder=3)
    zoom.scatter(0, coeff[0], facecolor="white", edgecolor="#394C60", marker="D", zorder=4,
                 label=rf"Extrapolated $\epsilon=0$: {coeff[0]:.5f}")
    zoom.set(title="Near the origin", xlim=(-0.9e-5, 1.8e-4))
    zoom.set_xticks([0, 5e-5, 1e-4, 1.5e-4])
    zoom.legend(frameon=False, loc="lower left")

    for axis in (all_points, zoom):
        axis.set(xlabel=r"$\epsilon^2$", ylabel=r"$\Delta Q/\epsilon$")
        axis.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useMathText=True)
        axis.grid(alpha=0.22)
    fig.savefig(OUTPUT / "q03-zero-slope.png", dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    primary, refined = measure(DT), measure(DT / 2)
    anchor_error = max(abs(a[key] - b[key]) for a, b in zip(primary, ARCHIVE["fine"])
                       for key in ("Q_A", "Q_B"))
    assert anchor_error < 1e-12
    y, fits = fit_response(primary)
    _, refined_fits = fit_response(refined)
    results = {
        "lambda": MODEL["LINK"], "reference_lambda": 0.0, "end_time": MODEL["END"],
        "dt": DT, "audit_dt": DT / 2, "primary": primary, "refined": refined,
        "fits": fits, "refined_fits": refined_fits, "previous_points_max_q03_error": anchor_error,
        "max_cross_sign_error_primary": max(abs(row["cross_sign_error"]) for row in primary),
        "max_cross_sign_error_refined": max(abs(row["cross_sign_error"]) for row in refined),
        "max_grid_difference_delta_over_epsilon": max(
            abs(a["delta_over_epsilon"] - b["delta_over_epsilon"]) for a, b in zip(primary, refined)),
    }
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    plot(y, fits[2])
    print("Wrote the zero-field extrapolation and audit to", OUTPUT)
