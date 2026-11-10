"""Render the fixed-axis RMS comparison and diagnostic holdout sweeps."""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator, StrMethodFormatter


DIRECTORY = Path(__file__).resolve().parent
TOPOLOGIES = (("open", "Open chain"), ("ring", "Ring"))
WINDOWS = (("train", "Train"), ("holdout_a", "Holdout A"), ("holdout_b", "Holdout B"))
INK, MUTED = "#172637", "#435469"
FIXED, X_REF, Y_REF = "#007c6a", "#497da8", "#ae641f"
SWEEP, POSTHOC = "#294a69", "#973d84"
AXES = (("fixed", "Fixed from Train", FIXED), ("X", "X axis", X_REF), ("Y", "Y axis", Y_REF))


def read_rows(name):
    with (DIRECTORY / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def time_label(value):
    multiple = round(value / math.pi)
    if multiple and math.isclose(value, multiple * math.pi, rel_tol=1e-9):
        return f"{multiple}\\pi" if multiple != 1 else r"\pi"
    return f"{value:.6g}"


def interval_label(bounds, *, closed_left):
    start, end = map(time_label, bounds)
    left = "[" if closed_left else "("
    return "$" + rf"{left}{start},\ {end}]$"


def style_axis(axis):
    axis.set_axisbelow(True)
    axis.grid(axis="y", color="#e1e9f0", linewidth=0.7)
    axis.tick_params(length=0, pad=5, labelsize=9, colors=MUTED)
    for spine in axis.spines.values():
        spine.set_color("#c9d4e0")


def heading(figure, title, subtitle):
    figure.text(0.08, 0.941, title, fontsize=20, color=INK, weight="bold")
    figure.text(0.08, 0.895, subtitle, fontsize=11, color=MUTED)


def draw_comparison(linear_rows, windows):
    linear = {(row["topology"], row["window"], row["axis"], None): float(row["rms"]) for row in linear_rows}
    finite = {(row["topology"], row["window"], row["axis"], float(row["abs_theta_rad"])):
              float(row["total_rms_pooled"]) for row in read_rows("finite.csv")}
    linear_limit = 1.12 * max(linear.values())
    finite_limit = 1.12 * max(finite.values())
    columns = (
        (r"Linear RMS · $\mathrm{rad}^{-1}$", linear, None, 2, linear_limit),
        (r"Total RMS at $|\theta|=0.10$ rad", finite, 0.10, 3, finite_limit),
        (r"Total RMS at $|\theta|=0.20$ rad", finite, 0.20, 3, finite_limit),
    )
    figure, axes = plt.subplots(2, 3, figsize=(12.6, 6.9))
    figure.subplots_adjust(left=0.080, right=0.976, top=0.759, bottom=0.135, wspace=0.24, hspace=0.40)
    heading(figure, "Fixed-axis RMS on training and holdout windows",
            r"$J=0.7$, $\Omega=1$ · The total RMS pools equal contributions from $+\theta$ and $-\theta$.")
    handles = [Patch(facecolor=color, label=label) for _, label, color in AXES]
    figure.legend(handles=handles, ncols=3, loc="upper left", bbox_to_anchor=(0.073, 0.885),
                  frameon=False, handlelength=1.2, columnspacing=1.7, fontsize=10, labelcolor=MUTED)

    width = 0.23
    positions = range(len(WINDOWS))
    for col, (title, values, angle, decimals, maximum) in enumerate(columns):
        figure.text(axes[0, col].get_position().x0, 0.81, title, fontsize=11.5, color=INK, weight="bold")
        for row, (topology, name) in enumerate(TOPOLOGIES):
            axis = axes[row, col]
            style_axis(axis)
            for index, (rotation_axis, _, color) in enumerate(AXES):
                heights = [values[topology, window, rotation_axis, angle] for window, _ in WINDOWS]
                axis.bar([pos + (index - 1) * width for pos in positions], heights,
                         width=0.92 * width, color=color, zorder=3)
            axis.set_xticks(list(positions), [label for _, label in WINDOWS])
            axis.set_ylim(0, maximum)
            axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
            axis.yaxis.set_major_formatter(StrMethodFormatter(f"{{x:.{decimals}f}}"))
            axis.set_title(name, loc="left", fontsize=10, color=MUTED, pad=6)

    bounds = "   ·   ".join(
        f"{label} {interval_label(windows[window], closed_left=(window == 'train'))}"
        for window, label in WINDOWS
    )
    figure.text(0.080, 0.064, rf"Time windows in $t$   ·   {bounds}", color=MUTED, fontsize=10)
    figure.savefig(DIRECTORY / "rms-comparison.png", dpi=180)


def draw_sweep(axis, points, optimum, title, top):
    beta, rms = zip(*sorted(points))
    fixed, posthoc, fixed_rms, posthoc_rms = (float(optimum[key]) for key in
                                            ("beta_fixed_deg", "beta_posthoc_deg", "rms_fixed", "rms_posthoc"))
    style_axis(axis)
    axis.grid(axis="x", color="#edf1f5", linewidth=0.7)
    axis.plot(beta, rms, color=SWEEP, linewidth=1.8)
    axis.plot(posthoc, posthoc_rms, "o", markersize=10.5, markeredgewidth=1.7,
              markeredgecolor=POSTHOC, markerfacecolor="white", zorder=4)
    axis.plot(fixed, fixed_rms, "o", markersize=6.5, markeredgewidth=0.8,
              markeredgecolor="white", markerfacecolor=FIXED, zorder=5)
    axis.set_xlim(-3, 183)
    axis.set_ylim(0, top)
    axis.set_xticks([0, 45, 90, 135, 180], ["0° (X)", "45°", "90° (Y)", "135°", "180° (X)"])
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
    axis.yaxis.set_major_formatter(StrMethodFormatter("{x:.2f}"))
    axis.set_title(title, loc="left", fontsize=11, color=INK, weight="bold", pad=8)
    axis.set_title(f"Fixed {fixed:.2f}°  ·  Post hoc {posthoc:.2f}°", loc="right",
                   fontsize=9.6, color=MUTED, pad=8)


def draw_holdouts(windows):
    sweeps = {}
    for row in read_rows("sweep.csv"):
        if row["window"] != "train":
            key = row["topology"], row["window"]
            sweeps.setdefault(key, []).append((float(row["beta_deg"]), float(row["rms"])))
    optima = {(row["topology"], row["window"]): row for row in read_rows("regret.csv")}
    top = 1.12 * max(rms for points in sweeps.values() for _, rms in points)
    figure, axes = plt.subplots(2, 2, figsize=(11.6, 6.9), sharey=True)
    figure.subplots_adjust(left=0.080, right=0.975, top=0.737, bottom=0.123, wspace=0.15, hspace=0.37)
    heading(figure, "Fixed axis and holdout-specific minima",
            r"RMS of $\partial_\theta\delta(0,t)$ · $J=0.7$, $\Omega=1$")
    handles = [Line2D([], [], color=SWEEP, linewidth=1.8, label="XY-plane sweep"),
               Line2D([], [], marker="o", linestyle="none", color=FIXED, label="Fixed from Train"),
               Line2D([], [], marker="o", linestyle="none", markeredgecolor=POSTHOC,
                      markerfacecolor="white", markeredgewidth=1.7, markersize=8, label="Post hoc minimum")]
    figure.legend(handles=handles, ncols=3, loc="upper left", bbox_to_anchor=(0.073, 0.883),
                  frameon=False, handlelength=1.65, columnspacing=1.6, fontsize=10, labelcolor=MUTED)

    for col, (window, label) in enumerate(WINDOWS[1:]):
        left = axes[0, col].get_position().x0
        figure.text(left, 0.798, f"{label}  ·  {interval_label(windows[window], closed_left=False)}",
                    fontsize=11, color=MUTED)
        for row, (topology, name) in enumerate(TOPOLOGIES):
            key = topology, window
            draw_sweep(axes[row, col], sweeps[key], optima[key], name, top)
        axes[-1, col].set_xlabel(r"Rotation-axis angle $\beta$ from X toward Y", fontsize=10, color=MUTED)
    figure.supylabel(r"RMS of linear residual coefficient ($\mathrm{rad}^{-1}$)",
                     x=0.021, y=0.43, fontsize=10.5, color=MUTED)
    figure.savefig(DIRECTORY / "holdout-axis-sweep.png", dpi=180)


def main():
    linear_rows = read_rows("linear.csv")
    windows = {row["window"]: (float(row["t_start"]), float(row["t_end"])) for row in linear_rows}
    draw_comparison(linear_rows, windows)
    draw_holdouts(windows)


if __name__ == "__main__":
    main()
