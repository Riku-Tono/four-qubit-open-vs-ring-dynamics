"""Plot the signed J=2 connected-correlation differences by pair and axis."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DIRECTORY = Path(__file__).resolve().parent
PAIRS = [("01", "0–1 / 2–3"), ("02", "0–2 / 1–3"), ("03", "0–3"), ("12", "1–2")]
CHANNELS = [("XX", "#22619C"), ("YY", "#7A438D"), ("ZZ", "#117366")]


def main():
    with (DIRECTORY / "waveforms.csv").open(encoding="utf-8", newline="") as data:
        rows = [row for row in csv.DictReader(data) if float(row["J"]) == 2.0]
    time = [float(row["t"]) for row in rows]

    fig, axes = plt.subplots(4, 3, figsize=(10, 7.8), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.17, right=0.98, top=0.82, bottom=0.11, wspace=0.14, hspace=0.20)
    fig.text(0.17, 0.955, "Where correlations differ at J = 2", fontsize=22, weight="bold", color="#172637")
    fig.text(0.17, 0.912, r"Ring minus open chain  ·  $\Omega=1$  ·  Initial state $|0000\rangle$",
             fontsize=12, color="#435469")
    fig.text(0.04, 0.465, r"Connected-correlation difference $\Delta C$",
             rotation=90, ha="center", va="center", fontsize=11, color="#435469")

    for row_index, (pair, label) in enumerate(PAIRS):
        for column, (channel, color) in enumerate(CHANNELS):
            axis = axes[row_index, column]
            difference = [float(row[f"ring_{pair}_{channel}"]) - float(row[f"open_{pair}_{channel}"])
                          for row in rows]
            axis.set_axisbelow(True)
            axis.grid(color="#e6edf4", linewidth=0.7)
            axis.axhline(0, color="#8696A8", linewidth=0.85)
            axis.plot(time, difference, color=color, linewidth=2)
            axis.set_xlim(time[0], time[-1])
            axis.set_xticks([0, 1, 2, 3])
            axis.tick_params(length=0, pad=5, labelsize=9, colors="#435469")
            for spine in axis.spines.values():
                spine.set_color("#c9d4e0")
            if row_index == 0:
                axis.set_title(channel, fontsize=12, weight="bold", color=color, pad=9)
            if column == 0:
                axis.set_ylabel(label, rotation=0, ha="right", va="center", labelpad=10,
                                fontsize=10, weight="bold", color="#172637")

    axes[0, 0].set_yticks([-0.4, -0.2, 0, 0.2, 0.4], ["−0.4", "−0.2", "0", "0.2", "0.4"])
    axes[0, 0].set_ylim(-0.48, 0.48)
    fig.supxlabel(r"Time $t$", x=0.575, y=0.035, fontsize=11, color="#435469")
    fig.savefig(DIRECTORY / "j2-waveforms.png", dpi=170)


if __name__ == "__main__":
    main()
