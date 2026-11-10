"""Plot occupation variation and selected phase differences for the Ising runs."""

import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


DIRECTORY = Path(__file__).resolve().parent
INK = "#172637"
MUTED = "#435469"
OPEN_COLOR = "#22619c"
RING_COLOR = "#ad4c21"
OPPOSITE_COLOR = "#7a438d"
OTHER_COLOR = "#117366"


def read_csv(name):
    with (DIRECTORY / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def style_axis(axis, title):
    axis.set_title(title, loc="left", color=INK, fontsize=12, weight="bold", pad=10)
    axis.set_axisbelow(True)
    axis.grid(color="#e6edf4", linewidth=0.7)
    axis.tick_params(length=0, pad=6, labelsize=9, colors=MUTED)
    for spine in axis.spines.values():
        spine.set_color("#c9d4e0")


def draw_occupation(axis, rows, metric, title, y_label):
    style_axis(axis, title)
    for topology, label, color in (("open", "Open chain", OPEN_COLOR), ("ring", "Ring", RING_COLOR)):
        samples = sorted((row for row in rows if row["topology"] == topology), key=lambda row: float(row["J"]))
        axis.plot([float(row["J"]) for row in samples], [float(row[metric]) for row in samples],
                  label=label, color=color, linewidth=2, marker="o", markersize=3.6)
    strengths = sorted({float(row["J"]) for row in rows})
    axis.set_xticks(strengths, [f"{strength:g}" for strength in strengths])
    axis.margins(y=0.08)
    axis.set_ylim(bottom=0)
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
    axis.set_xlabel(r"Uniform bond strength $J$", fontsize=10, color=MUTED)
    axis.set_ylabel(y_label, fontsize=10, color=MUTED)


def draw_phase(axis, rows, topology, coupling, title, opposite_pair):
    samples = [row for row in rows if row["topology"] == topology and float(row["J"]) == coupling]
    time = [float(row["t"]) for row in samples]
    style_axis(axis, title)
    axis.axhline(180, color="#bac5d3", linestyle="--", linewidth=0.9)
    for pair, color, qualifier in ((opposite_pair, OPPOSITE_COLOR, " (opposite phase)"), ("01", OTHER_COLOR, "")):
        gap = [float(row[f"angle_{pair}"]) if row[f"angle_{pair}"] else np.nan for row in samples]
        axis.plot(time, gap, label=f"{pair[0]}–{pair[1]}{qualifier}", color=color, linewidth=2.1)
    axis.set_xlim(0, 8 * np.pi)
    axis.set_ylim(-7, 193)
    axis.set_xticks(np.arange(5) * 2 * np.pi, ["0", r"$2\pi$", r"$4\pi$", r"$6\pi$", r"$8\pi$"])
    axis.set_yticks([0, 45, 90, 135, 180], ["0°", "45°", "90°", "135°", "180°"])
    axis.set_xlabel(r"Time $t$", fontsize=10, color=MUTED)
    axis.set_ylabel("Absolute phase difference", fontsize=10, color=MUTED)
    axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.23), ncols=2,
                fontsize=9, frameon=False, labelcolor=MUTED)


def main():
    short = [row for row in read_csv("summary.csv") if row["window"] == "short"]
    waves = read_csv("waveforms.csv")
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.7))
    fig.subplots_adjust(left=0.087, right=0.973, top=0.75, bottom=0.17, wspace=0.24, hspace=0.52)
    fig.text(0.087, 0.948, "Smaller occupation changes can coexist with opposite phase",
             fontsize=21, weight="bold", color=INK)
    fig.text(0.087, 0.902, r"Four qubits  ·  $\Omega=1$  ·  Initial state $|1000\rangle$  ·  Short window $0\leq t\leq3.688$",
             fontsize=11, color=MUTED)
    fig.text(0.087, 0.863, "Phase differences are shown only where all four analytic-signal envelopes are at least 0.05.",
             fontsize=10, color=MUTED)

    draw_occupation(axes[0, 0], short, "occupation_spread", "A. Short-window occupation spread", "Site spread (time RMS)")
    draw_occupation(axes[0, 1], short, "mean_occupation_amplitude", "B. Short-window occupation amplitude",
                    "Mean half peak-to-peak range")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.080, 0.838), ncols=2,
               frameon=False, fontsize=10, labelcolor=MUTED)

    draw_phase(axes[1, 0], waves, "open", 1.0, "C. Phase differences in the open chain at J = 1.0", "03")
    draw_phase(axes[1, 1], waves, "ring", 0.7, "D. Phase differences in the ring at J = 0.7", "02")
    fig.savefig(DIRECTORY / "phase-alignment.png", dpi=175)


if __name__ == "__main__":
    main()
