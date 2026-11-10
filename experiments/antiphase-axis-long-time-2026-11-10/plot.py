"""Plot cumulative optimal-axis convergence and its infinite-time references."""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator


DIRECTORY = Path(__file__).resolve().parent
TOPOLOGIES = (("open", "Open chain"), ("ring", "Ring"))
TIME_MARKS = (8, 16, 32, 64)
INK, MUTED = "#172637", "#435469"
CURVE, LIMIT = "#1264a0", "#a65c1e"
LIMIT_STYLE = (0, (5, 3))


def read_rows(name):
    with (DIRECTORY / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def by_topology(rows):
    cases = {topology: [] for topology, _ in TOPOLOGIES}
    for row in rows:
        cases[row["topology"]].append(row)
    for case in cases.values():
        case.sort(key=lambda row: float(row["T"]))
    return cases


def axis_angles(rows, limit):
    return [limit + (float(row["beta_opt_deg"]) - limit + 90) % 180 - 90 for row in rows]


def improvements(rows):
    return [100 * float(row["improvement_fraction"]) for row in rows]


def limits_for(values):
    lower, upper = min(values), max(values)
    span = upper - lower
    return lower - 0.068 * span, upper + 0.085 * span


def draw_panel(axis, curve, samples, limit, y_limits, topology, limit_label):
    axis.set_axisbelow(True)
    axis.grid(color="#e3eaf1", linewidth=0.7)
    axis.axhline(limit, color=LIMIT, linewidth=1.3, linestyle=LIMIT_STYLE)
    axis.plot(*curve, color=CURVE, linewidth=1.7)
    axis.plot(*samples, "o", markeredgecolor=CURVE, markerfacecolor="white",
              markeredgewidth=1.6, markersize=5.7, zorder=4)
    axis.set_xlim(6.4, 65.6)
    axis.set_ylim(*y_limits)
    axis.set_xticks(TIME_MARKS)
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
    axis.tick_params(length=0, pad=5, labelsize=9, colors=MUTED)
    for spine in axis.spines.values():
        spine.set_color("#c9d4e0")
    axis.set_title(topology, loc="left", fontsize=11, weight="bold", color=INK, pad=8)
    axis.set_title(limit_label, loc="right", fontsize=9.7, color=LIMIT, pad=8)


def main():
    curve_rows = read_rows("curve.csv")
    sample_rows = read_rows("convergence.csv")
    curves = by_topology(curve_rows)
    samples = by_topology(sample_rows)
    limits = {row["topology"]: row for row in read_rows("limits.csv")}
    gain_limits = limits_for(improvements([*curve_rows, *sample_rows, *limits.values()]))
    figure, axes = plt.subplots(2, 2, figsize=(11.7, 6.85))
    figure.subplots_adjust(left=0.078, right=0.971, top=0.754, bottom=0.158, wspace=0.16, hspace=0.40)
    figure.text(0.078, 0.942, "Optimal rotation axis as the time window grows",
                fontsize=20, color=INK, weight="bold")
    figure.text(0.078, 0.895, r"Cumulative intervals $[0,T]$ · $\Omega=1$, $J=0.7$ · RMS of the linear residual coefficient",
                fontsize=11, color=MUTED)
    handles = [Line2D([], [], color=CURVE, linewidth=1.7, label="Cumulative curve"),
               Line2D([], [], marker="o", linestyle="none", markeredgecolor=CURVE,
                      markerfacecolor="white", markeredgewidth=1.6, label="Sampled endpoints"),
               Line2D([], [], color=LIMIT, linewidth=1.3, linestyle=LIMIT_STYLE, label="Infinite-time limit")]
    figure.legend(handles=handles, ncols=3, loc="upper left", bbox_to_anchor=(0.071, 0.885),
                  frameon=False, handlelength=2, columnspacing=1.6, fontsize=10, labelcolor=MUTED)
    headings = ("Optimal axis angle β (degrees)", "RMS improvement over the better of X and Y (%)")
    for col, title in enumerate(headings):
        figure.text(axes[0, col].get_position().x0, 0.808, title, fontsize=11.2, color=INK, weight="bold")

    for row, (topology, title) in enumerate(TOPOLOGIES):
        curve, points, reference = curves[topology], samples[topology], limits[topology]
        angle_limit = float(reference["beta_opt_deg"])
        improvement_limit = improvements([reference])[0]
        curve_t = [float(item["T"]) / math.pi for item in curve]
        sample_t = [float(item["T"]) / math.pi for item in points]
        curve_angles = axis_angles(curve, angle_limit)
        sample_angles = axis_angles(points, angle_limit)
        angle_limits = limits_for([*curve_angles, *sample_angles, angle_limit])
        draw_panel(axes[row, 0], (curve_t, curve_angles), (sample_t, sample_angles),
                   angle_limit, angle_limits, title, f"Limit {angle_limit:g}° ({reference['best_axis']})")
        draw_panel(axes[row, 1], (curve_t, improvements(curve)), (sample_t, improvements(points)),
                   improvement_limit, gain_limits, title, f"Limit {improvement_limit:g}%")
    for axis in axes[-1]:
        axis.set_xlabel(r"Window endpoint $T/\pi$", fontsize=10, color=MUTED)

    figure.text(0.078, 0.055,
                "Angles are shown modulo 180° using the representative nearest the limit. The angle scales differ by topology.",
                fontsize=9.3, color=MUTED)
    figure.savefig(DIRECTORY / "long-time-convergence.png", dpi=180)


if __name__ == "__main__":
    main()
