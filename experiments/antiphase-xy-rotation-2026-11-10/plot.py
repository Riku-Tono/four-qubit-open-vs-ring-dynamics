"""Plot signed occupation residuals at |theta| = 0.20 rad from waveforms.csv."""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator


DIRECTORY = Path(__file__).resolve().parent
ANGLE = "0.2"
CASES = (
    ("open", 0.7, r"Open chain · $J=0.7$", r"$R_0+R_3-1$"),
    ("ring", 0.7, r"Ring · $J=0.7$", r"$R_0+R_2-1$"),
    ("open", 0.0, r"Uncoupled · $J=0$ · both topologies", r"$R_0+R_k-1$"),
)
WINDOWS = (("short", 3.688, "Short window"), ("long", 8 * math.pi, "Long window"))
AXES = (("X", "#1264a0"), ("Y", "#aa5310"))
SIGNS = (("plus", "-", "+"), ("minus", (0, (4, 2.4)), "−"))
CURVES = [(f"{axis}_{sign}_{ANGLE}", color, style, f"{axis}, {label}0.20")
          for axis, color in AXES for sign, style, label in SIGNS]
INK, MUTED = "#172637", "#435469"


def read_waveforms():
    cases = {(topology, coupling, window): []
             for topology, coupling, _, _ in CASES for window, _, _ in WINDOWS}
    with (DIRECTORY / "waveforms.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            key = row["topology"], float(row["J"]), row["window"]
            if key in cases:
                cases[key].append((float(row["t"]), *(float(row[name]) for name, *_ in CURVES)))
    return cases


def draw_window(axis, samples, end):
    time, *curves = zip(*samples)
    axis.set_axisbelow(True)
    axis.grid(color="#e2e9f0", linewidth=0.7)
    axis.axhline(0, color="#9aa8b5", linewidth=0.8)
    for values, (_, color, style, _) in zip(curves, CURVES):
        axis.plot(time, values, color=color, linestyle=style, linewidth=1.45)

    axis.set_xlim(0, end)
    if end == 3.688:
        axis.set_xticks([0, 1, 2, 3, end], ["0", "1", "2", "3", "3.688"])
    else:
        axis.set_xticks([i * 2 * math.pi for i in range(5)],
                        ["0", r"$2\pi$", r"$4\pi$", r"$6\pi$", r"$8\pi$"])
    axis.yaxis.set_major_locator(MultipleLocator(0.05))
    axis.tick_params(length=0, pad=5, labelsize=9, colors=MUTED)
    for spine in axis.spines.values():
        spine.set_color("#c9d4e0")


def main():
    samples = read_waveforms()
    limit = 1.1 * max(abs(value) for rows in samples.values() for _, *values in rows for value in values)
    figure, axes = plt.subplots(3, 2, figsize=(11.4, 7.8), sharey=True)
    figure.subplots_adjust(left=0.10, right=0.975, top=0.776, bottom=0.085, wspace=0.16, hspace=0.45)
    figure.text(0.10, 0.935, "Signed occupation residuals after X and Y rotations",
                fontsize=20, color=INK, weight="bold")
    figure.text(0.10, 0.893,
                r"Site 0 of $|1000\rangle$ · $\Omega=1$ · rotation angles in radians",
                fontsize=11, color=MUTED)
    legend = [Line2D([], [], color=color, linestyle=style, linewidth=1.9, label=label)
              for _, color, style, label in CURVES]
    figure.legend(handles=legend, ncols=4, loc="upper left", bbox_to_anchor=(0.092, 0.885),
                  frameon=False, handlelength=2.5, columnspacing=1.7, fontsize=10, labelcolor=MUTED)

    for row, (topology, coupling, label, ylabel) in enumerate(CASES):
        for col, (window, end, title) in enumerate(WINDOWS):
            axis = axes[row, col]
            draw_window(axis, samples[topology, coupling, window], end)
            axis.set_ylim(-limit, limit)
            axis.set_title(f"{label} · {title}" if row == 0 else label,
                           loc="left", color=INK, fontsize=10.5, weight="bold", pad=8)
            if row == len(CASES) - 1:
                axis.set_xlabel(r"Time $t$", fontsize=10, color=MUTED)
        axes[row, 0].set_ylabel(ylabel, fontsize=11, color=MUTED)

    figure.savefig(DIRECTORY / "xy-residuals.png", dpi=180)


if __name__ == "__main__":
    main()
