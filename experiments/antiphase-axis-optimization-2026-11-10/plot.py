"""Plot the angular sweep and exact RMS minima from the experiment CSV files."""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


DIRECTORY = Path(__file__).resolve().parent
CASES = (
    ("open", 0.7, r"Open chain · $J=0.7$"),
    ("ring", 0.7, r"Ring · $J=0.7$"),
    ("open", 0.0, r"Uncoupled · $J=0$ · both topologies"),
)
WINDOWS = (("short", r"Short window · $0\leq t\leq3.688$"),
           ("long", r"Long window · $0\leq t\leq8\pi$"))
INK, MUTED = "#172637", "#435469"
SWEEP, MINIMUM, X_REF, Y_REF = "#1264a0", "#00816f", "#4f5664", "#ab5310"


def case_key(row):
    return row["topology"], float(row["J"]), row["window"]


def read_results():
    sweeps = {(topology, coupling, window): []
              for topology, coupling, _ in CASES for window, _ in WINDOWS}
    with (DIRECTORY / "sweep.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            key = case_key(row)
            if key in sweeps:
                sweeps[key].append((float(row["beta_deg"]), float(row["rms"])))
    with (DIRECTORY / "optima.csv").open(encoding="utf-8", newline="") as handle:
        optima = {case_key(row): row for row in csv.DictReader(handle)}
    return sweeps, optima


def draw_sweep(axis, points, optimum, top, title):
    beta, rms = zip(*sorted(points))
    angle, minimum, x_ref, y_ref = (float(optimum[key]) for key in
                                   ("beta_opt_deg", "rms_min", "rms_x", "rms_y"))
    axis.set_axisbelow(True)
    axis.grid(color="#e3eaf1", linewidth=0.7)
    axis.axhline(0, color="#aab6c1", linewidth=0.8)
    axis.plot(beta, rms, color=SWEEP, linewidth=1.75)
    axis.plot(angle, minimum, "o", markersize=10, markeredgewidth=1.7,
              markeredgecolor=MINIMUM, markerfacecolor="white", zorder=4)
    axis.plot([0, 180], [x_ref, x_ref], "s", color=X_REF, markersize=4.7, zorder=5)
    axis.plot(90, y_ref, "D", color=Y_REF, markersize=4.7, zorder=5)

    axis.set_xlim(-4, 184)
    axis.set_ylim(-0.032, top + 0.008)
    axis.set_xticks([0, 45, 90, 135, 180], ["0° (X)", "45°", "90° (Y)", "135°", "180° (X)"])
    axis.set_yticks([step / 10 for step in range(round(10 * top) + 1)])
    axis.tick_params(length=0, pad=5, labelsize=9, colors=MUTED)
    for spine in axis.spines.values():
        spine.set_color("#c9d4e0")

    axis.set_title(title, loc="left", fontsize=10.5, color=INK, weight="bold", pad=9)
    axis.set_title(rf"$\beta_{{\mathrm{{min}}}}={angle:.2f}^\circ$ · RMS {minimum:.5f}",
                   loc="right", fontsize=9.4, color=MINIMUM, pad=9)


def main():
    sweeps, optima = read_results()
    top = math.ceil(10 * max(rms for points in sweeps.values() for _, rms in points)) / 10
    figure, axes = plt.subplots(3, 2, figsize=(11.7, 8.1), sharey=True)
    figure.subplots_adjust(left=0.093, right=0.976, top=0.778, bottom=0.091, wspace=0.14, hspace=0.47)
    figure.text(0.093, 0.939, "Rotation-axis sweep for the linear residual",
                fontsize=20, color=INK, weight="bold")
    figure.text(0.093, 0.899,
                r"RMS of $\partial_\theta\delta(0,t)$ · Site 0 of $|1000\rangle$ · $\Omega=1$",
                fontsize=11, color=MUTED)
    handles = [Line2D([], [], color=SWEEP, linewidth=1.75, label="XY-plane sweep"),
               Line2D([], [], marker="o", linestyle="none", markeredgecolor=MINIMUM,
                      markerfacecolor="white", markeredgewidth=1.7, markersize=8, label="Minimum"),
               Line2D([], [], marker="s", linestyle="none", color=X_REF, label="X axis"),
               Line2D([], [], marker="D", linestyle="none", color=Y_REF, label="Y axis")]
    figure.legend(handles=handles, ncols=4, loc="upper left", bbox_to_anchor=(0.086, 0.89),
                  frameon=False, handlelength=1.65, columnspacing=1.6, fontsize=10, labelcolor=MUTED)

    for column, (window, window_title) in enumerate(WINDOWS):
        left = axes[0, column].get_position().x0
        figure.text(left, 0.827, window_title, fontsize=11, color=MUTED)
        for row, (topology, coupling, title) in enumerate(CASES):
            key = topology, coupling, window
            draw_sweep(axes[row, column], sweeps[key], optima[key], top, title)
        axes[-1, column].set_xlabel(r"Rotation-axis angle $\beta$ from X toward Y",
                                   fontsize=10, color=MUTED)
    figure.supylabel(r"RMS of linear residual coefficient ($\mathrm{rad}^{-1}$)",
                     x=0.032, y=0.435, fontsize=11, color=MUTED)
    figure.savefig(DIRECTORY / "axis-rms-sweep.png", dpi=180)


if __name__ == "__main__":
    main()
