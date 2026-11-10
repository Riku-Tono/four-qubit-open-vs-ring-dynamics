"""Overlay the positive and negative local Y-rotation residuals."""

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator


DIRECTORY = Path(__file__).resolve().parent
ANGLES = ("0.025", "0.05", "0.1", "0.2")
COLORS = ("#a45b00", "#00806e", "#1264a0", "#873985")
SIGNS = (("plus", "-", r"$+\theta$"), ("minus", (0, (4, 2)), r"$-\theta$"))
BASELINE = "#657483"
INK, MUTED = "#172637", "#435469"
SHORT, LONG = 3.688, 8 * math.pi
COUPLING = 0.7


def read_waveforms():
    columns = [f"delta_{sign}_{angle}" for angle in ANGLES for sign, _, _ in SIGNS]
    rows = {"open": [], "ring": []}
    with (DIRECTORY / "waveforms.csv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if float(row["J"]) == COUPLING:
                rows[row["topology"]].append(
                    tuple(float(row[key]) for key in ("t", "delta_0", *columns))
                )
    return rows


def draw_window(axis, rows, end, title):
    # A rounded CSV time at the window endpoint may lie just above the bound.
    time, baseline, *curves = zip(*(row for row in rows if row[0] <= end + 1e-7))
    axis.set_title(title, loc="left", color=INK, fontsize=11.5, weight="bold", pad=9)
    axis.set_axisbelow(True)
    axis.grid(color="#e2e9f0", linewidth=0.65)
    axis.plot(time, baseline, color=BASELINE, linestyle=":", linewidth=1.2)
    for index, color in enumerate(COLORS):
        plus, minus = curves[2 * index:2 * index + 2]
        axis.plot(time, plus, color=color, linestyle=SIGNS[0][1], linewidth=1.35)
        axis.plot(time, minus, color=color, linestyle=SIGNS[1][1], linewidth=1.35)

    axis.set_xlim(0, end)
    if end == LONG:
        axis.set_xticks(
            [0, 2 * math.pi, 4 * math.pi, 6 * math.pi, 8 * math.pi],
            ["0", r"$2\pi$", r"$4\pi$", r"$6\pi$", r"$8\pi$"],
        )
    else:
        axis.set_xticks([0, 1, 2, 3, SHORT], ["0", "1", "2", "3", "3.688"])
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5, symmetric=True))
    axis.tick_params(length=0, pad=6, labelsize=9, colors=MUTED)
    for spine in axis.spines.values():
        spine.set_color("#c9d4e0")


def main():
    rows = read_waveforms()
    limit = 1.1 * max(abs(value) for samples in rows.values() for row in samples for value in row[1:])
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 6.85), sharey=True)
    fig.subplots_adjust(left=0.085, right=0.974, top=0.765, bottom=0.105, wspace=0.13, hspace=0.40)
    fig.text(0.085, 0.932, "Signed occupation residuals under opposite Y rotations",
             color=INK, fontsize=20, weight="bold")
    fig.text(0.085, 0.884,
             r"Four qubits   ·   $R_y(\theta)$ on site 0 of $|1000\rangle$   ·   $\Omega=1,\ J=0.7$",
             color=MUTED, fontsize=11)

    angle_handles = [
        Line2D([], [], color=color, linewidth=2.25, label=angle)
        for angle, color in zip(ANGLES, COLORS)
    ]
    fig.legend(handles=angle_handles, ncols=4, loc="upper left",
               bbox_to_anchor=(0.206, 0.862), frameon=False, handlelength=1.75,
               columnspacing=1.25, fontsize=10, labelcolor=MUTED)
    fig.text(0.085, 0.827, r"$|\theta|$ (radians)", fontsize=10, color=MUTED)
    sign_handles = [
        Line2D([], [], color=MUTED, linewidth=1.5, linestyle=style, label=label)
        for _, style, label in SIGNS
    ]
    sign_handles.append(Line2D([], [], color=BASELINE, linestyle=":", linewidth=1.5, label=r"$\theta=0$"))
    fig.legend(handles=sign_handles, ncols=3, loc="upper right", bbox_to_anchor=(0.975, 0.862),
               frameon=False, handlelength=2.0, columnspacing=1.3, fontsize=10, labelcolor=MUTED)

    for index, (topology, label, pair) in enumerate(
        (("open", "Open chain", "3"), ("ring", "Ring", "2"))
    ):
        for column, (end, window) in enumerate(((SHORT, "short window"), (LONG, "long window"))):
            axis = axes[index, column]
            draw_window(axis, rows[topology], end, f"{label} · {window}")
            axis.set_ylim(-limit, limit)
            if index == 1:
                axis.set_xlabel(r"Time $t$", fontsize=10, color=MUTED)
        axes[index, 0].set_ylabel(rf"$R_0+R_{pair}-1$", fontsize=11, color=MUTED)

    fig.savefig(DIRECTORY / "signed-residual-waveforms.png", dpi=180)


if __name__ == "__main__":
    main()
