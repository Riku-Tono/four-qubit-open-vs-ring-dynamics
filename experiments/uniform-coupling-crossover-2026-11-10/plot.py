"""Plot the pair RMS crossovers and the outer-pair YY difference."""

import csv
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator


DIRECTORY = Path(__file__).resolve().parent
PAIRS = [
    ("03", "0–3", "#22619c"),
    ("12", "1–2", "#ad4c21"),
    ("01", "0–1 / 2–3", "#117366"),
    ("02", "0–2 / 1–3", "#6e7184"),
]
INK = "#172637"
MUTED = "#435469"
YY_COLOR = "#7a438d"


def read_sweep():
    with (DIRECTORY / "sweep.csv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    coupling = np.array([float(row["J"]) for row in rows])
    pairs = np.array([[float(row[f"pair_{pair}"]) for row in rows] for pair, _, _ in PAIRS])
    yy = np.array([float(row["outer_YY_rms"]) for row in rows])
    return coupling, pairs, yy


def draw_pair_curves(axis, coupling, pairs, crossings):
    intervals = [(coupling[0], crossings[0]["from"])]
    intervals.extend((row["J"], row["to"]) for row in crossings)
    bounds = [position for position, _ in intervals[1:]] + [coupling[-1]]
    styles = {pair: (label, color) for pair, label, color in PAIRS}
    for (start, leader), stop in zip(intervals, bounds):
        label, color = styles[leader]
        axis.axvspan(start, stop, color=color, alpha=0.07, linewidth=0)
        axis.text((start + stop) / 2, 0.95, f"Leader  {label}", transform=axis.get_xaxis_transform(),
                  color=color, fontsize=10, weight="bold", va="top", ha="center")
    for values, (_, label, color) in zip(pairs, PAIRS):
        axis.plot(coupling, values, label=label, color=color, linewidth=2.25)
    for row in crossings:
        position = row["J"]
        height = row["pair_rms"][row["to"]]
        axis.axvline(position, color=MUTED, linewidth=0.9, linestyle=(0, (3, 4)))
        axis.scatter(position, height, color=INK, edgecolor="white", linewidth=0.9, zorder=5, s=47)
        axis.text(position, 0.075, f"{position:.3f}", transform=axis.get_xaxis_transform(),
                  ha="center", va="center", color=INK, fontsize=9,
                  bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5})
    axis.set_ylim(0, 1.22 * np.max(pairs))


def draw_outer_yy(axis, coupling, yy):
    axis.plot(coupling, yy, color=YY_COLOR, linewidth=2.4)
    for index, label, offset, alignment in (
        (np.argmin(yy), "Minimum", (-10, 21), "right"),
        (np.argmax(yy), "Maximum", (-10, -25), "right"),
    ):
        point = coupling[index], yy[index]
        axis.scatter(*point, color=YY_COLOR, edgecolor="white", linewidth=0.9, zorder=5, s=47)
        axis.annotate(f"{label} near J = {point[0]:.2f}", point, xytext=offset,
                      textcoords="offset points", ha=alignment, color=YY_COLOR, fontsize=9,
                      bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.0})
    margin = 0.12 * (np.max(yy) - np.min(yy))
    axis.set_ylim(np.min(yy) - margin, np.max(yy) + margin)


def main():
    coupling, pairs, yy = read_sweep()
    results = json.loads((DIRECTORY / "results.json").read_text(encoding="utf-8"))
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(9.5, 7.1), sharex=True,
                                     gridspec_kw={"height_ratios": [1.55, 1]})
    fig.subplots_adjust(left=0.115, right=0.973, top=0.76, bottom=0.10, hspace=0.34)
    fig.text(0.115, 0.948, "Which pair differs most as the coupling grows?",
             fontsize=22, weight="bold", color=INK)
    fig.text(0.115, 0.900, r"Ring minus open chain  ·  $\Omega=1$  ·  Initial state $|0000\rangle$  ·  $0\leq t\leq3.688$",
             fontsize=11, color=MUTED)

    draw_pair_curves(top, coupling, pairs, results["leader_crossings"])
    draw_outer_yy(bottom, coupling, yy)
    top.set_title("Pair difference across XX, YY and ZZ", loc="left", color=INK, fontsize=12, weight="bold", pad=8)
    bottom.set_title("YY difference on the outer pairs 0–1 and 2–3", loc="left", color=INK,
                     fontsize=12, weight="bold", pad=8)
    handles, labels = top.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.108, 0.872), ncols=4, frameon=False,
               fontsize=10, labelcolor=MUTED, columnspacing=2.2, handlelength=2.2)

    for axis in (top, bottom):
        axis.set_xlim(coupling[0], coupling[-1])
        axis.set_axisbelow(True)
        axis.grid(axis="y", color="#e6edf4", linewidth=0.8)
        axis.tick_params(length=0, pad=6, labelsize=9, colors=MUTED)
        axis.set_ylabel("Time RMS", color=MUTED, fontsize=10, labelpad=10)
        axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
        axis.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        for spine in axis.spines.values():
            spine.set_color("#c9d4e0")
    bottom.set_xticks(np.linspace(1, 2, 6))
    bottom.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    bottom.set_xlabel(r"Uniform bond strength $J$", color=MUTED, fontsize=11, labelpad=8)
    fig.savefig(DIRECTORY / "crossover.png", dpi=175)


if __name__ == "__main__":
    main()
