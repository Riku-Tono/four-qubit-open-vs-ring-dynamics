"""Draw the crossover waveforms and the sampled window-length map."""

import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator


DIRECTORY = Path(__file__).resolve().parent
INK = "#172637"
MUTED = "#435469"
PAIR_STYLES = {
    "03": ("0–3", "#22619c", "#b9d5eb"),
    "12": ("1–2", "#ad4c21", "#efc6b1"),
    "01": ("0–1 / 2–3", "#117366", "#a9dacf"),
    "02": ("0–2 / 1–3", "#6e7184", "#d7d7e3"),
}
CHANNEL_STYLES = {"XX": "#22619c", "YY": "#7a438d", "ZZ": "#117366"}


def read_csv(name):
    with (DIRECTORY / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def values(rows, key):
    return np.array([float(row[key]) for row in rows])


def style_axis(axis):
    axis.tick_params(length=0, pad=6, labelsize=9, colors=MUTED)
    for spine in axis.spines.values():
        spine.set_color("#c9d4e0")


def style_waveform(axis, time, title):
    style_axis(axis)
    axis.set_title(title, loc="left", color=INK, fontsize=12, weight="bold", pad=11)
    axis.set_axisbelow(True)
    axis.axvspan(2, 3, color="#aab6c8", alpha=0.20, linewidth=0)
    axis.grid(color="#e6edf4", linewidth=0.7)
    axis.axhline(0, color="#8294a6", linewidth=0.9)
    axis.set_xlim(time[0], time[-1])
    axis.set_xticks([0, 1, 2, 3])
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
    axis.set_xlabel(r"Time $t$", fontsize=10, color=MUTED)


def draw_waveforms(signed_axis, cumulative_axis, crossing):
    time = values(crossing, "t")
    style_waveform(signed_axis, time, "A. Signed ZZ difference")
    for name, key in (("outer", "01"), ("center", "12")):
        label, color, _ = PAIR_STYLES[key]
        signed_axis.plot(time, values(crossing, f"delta_{name}_ZZ"), color=color, label=label, linewidth=2.2)
    signed_axis.set_ylabel(r"$\Delta C^{ZZ}$", fontsize=11, color=MUTED)
    signed_axis.legend(loc="best", frameon=True, facecolor="white", edgecolor="none",
                       fontsize=9, labelcolor=MUTED, ncols=2)

    style_waveform(cumulative_axis, time, "B. Cumulative outer-pair advantage")
    cumulative = [values(crossing, f"cumulative_{channel}") for channel in CHANNEL_STYLES]
    for channel, series in zip(CHANNEL_STYLES, cumulative):
        cumulative_axis.plot(time, series, color=CHANNEL_STYLES[channel], label=channel, linewidth=2.2)
    cumulative_axis.plot(time, np.sum(cumulative, axis=0), color=MUTED, linestyle="--", label="Sum", linewidth=1.2)
    cumulative_axis.set_ylabel(r"$T^{-1}\int_0^t (\Delta C_{01}^2-\Delta C_{12}^2)\,\mathrm{d}s$",
                              fontsize=10, color=MUTED)
    cumulative_axis.legend(loc="best", frameon=True, facecolor="white", edgecolor="none",
                           fontsize=9, labelcolor=MUTED, ncols=4, columnspacing=0.95, handlelength=1.8)


def draw_window_map(axis, rows, baseline, crossing_j):
    coupling = values(rows, "J")
    windows = values(rows, "T")
    x, x_index = np.unique(coupling, return_inverse=True)
    y, y_index = np.unique(windows, return_inverse=True)
    categories = {name: index for index, name in enumerate(PAIR_STYLES)}
    grid = np.full((len(y), len(x)), np.nan)
    grid[y_index, x_index] = [categories[row["leader"].zfill(2)] for row in rows]
    cmap = ListedColormap([style[2] for style in PAIR_STYLES.values()])
    axis.pcolormesh(x, y, grid, cmap=cmap, shading="nearest", vmin=-0.5, vmax=len(categories) - 0.5)
    axis.axhline(baseline, color=INK, linewidth=1, linestyle=(0, (4, 4)))
    axis.scatter(crossing_j, baseline, color=INK, edgecolor="white", linewidth=0.8, s=42, zorder=4)
    axis.set(xlim=(x[0], x[-1]), ylim=(y[0], y[-1]))
    axis.set_xticks(np.linspace(x[0], x[-1], 6))
    axis.set_yticks(np.linspace(y[0], y[-1], 5), [f"{t:.3f}" for t in np.linspace(y[0], y[-1], 5)])
    style_axis(axis)
    axis.set_title("C. Pair with the largest three-channel RMS for each sampled window", loc="left",
                   color=INK, fontsize=12, weight="bold", pad=11)
    axis.set_xlabel(r"Uniform bond strength $J$", fontsize=11, color=MUTED, labelpad=8)
    axis.set_ylabel(r"Window end $T$", fontsize=11, color=MUTED, labelpad=8)
    axis.text(x[0] + 0.02 * (x[-1] - x[0]), baseline, f"Baseline T = {baseline:.3f}", ha="left", va="center",
              color=INK, fontsize=9, bbox={"facecolor": "white", "edgecolor": "none", "pad": 2})
    present = set(row["leader"].zfill(2) for row in rows)
    legend = [Patch(facecolor=fill, edgecolor="none", label=label)
              for pair, (label, _, fill) in PAIR_STYLES.items() if pair in present]
    axis.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncols=len(legend),
                frameon=False, fontsize=10, labelcolor=MUTED, columnspacing=2)


def main():
    crossing = [row for row in read_csv("waveforms.csv") if row["case"] == "crossing"]
    windows = read_csv("windows.csv")
    crossing_j = float(crossing[0]["J"])
    baseline = float(crossing[-1]["t"])

    fig = plt.figure(figsize=(11.5, 8.4))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.18])
    fig.subplots_adjust(left=0.087, right=0.974, top=0.81, bottom=0.145, hspace=0.55, wspace=0.30)
    signed_axis, cumulative_axis = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
    window_axis = fig.add_subplot(grid[1, :])
    fig.text(0.087, 0.950, "When does the leading pair change?", color=INK, fontsize=23, weight="bold")
    fig.text(0.087, 0.909, r"Ring minus open chain  ·  $\Omega=1$  ·  Initial state $|0000\rangle$",
             color=MUTED, fontsize=11)
    fig.text(0.087, 0.869, f"Panels A and B at the baseline crossover J = {crossing_j:.3f}. Gray shading marks 2 ≤ t ≤ 3.",
             color=MUTED, fontsize=10)
    draw_waveforms(signed_axis, cumulative_axis, crossing)
    draw_window_map(window_axis, windows, baseline, crossing_j)
    fig.savefig(DIRECTORY / "timing-and-window.png", dpi=175)


if __name__ == "__main__":
    main()
