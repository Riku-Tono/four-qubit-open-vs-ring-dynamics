"""Plot the short-window selections and their fixed-condition holdout results."""

import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D


DIRECTORY = Path(__file__).resolve().parent
INK = "#172638"
MUTED = "#4c5c6b"
GRID = "#dce4e8"
BACKGROUND = "#f5f7f8"
PAIR_LABELS = ("01", "02", "03", "12", "13", "23")
CASES = (
    ("unrotated_1000", r"Unrotated $|1000\rangle$", "|1000⟩"),
    ("best_existing", "Best existing condition", "Best existing"),
    ("selected_0.05", "Four local rotations", "Four rotations"),
)
TOPOLOGIES = (("open", "Open chain"), ("ring", "Ring"))
SHORT_END = 3.688
GATE = 0.05
WINDOWS = (("short", "Short", 0.0, SHORT_END), ("holdout", "Holdout", SHORT_END, 8 * np.pi))
SITE_COLORS = ("#14609e", "#c74b22", "#117e65", "#8b49a2")
SITE_STYLES = ("solid", (0, (5, 2.8)), (2, (1.3, 2.0)), (3, (7, 2.5, 1, 2.5)))


def read_csv(name):
    with (DIRECTORY / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def style_wave_axis(axis):
    axis.set_ylim(-0.02, 1.02)
    axis.set_yticks([0, 0.25, 0.5, 0.75, 1])
    axis.set_axisbelow(True)
    axis.grid(color=GRID, linewidth=0.65)
    axis.tick_params(length=0, pad=5, labelsize=9, colors=MUTED)
    for spine in axis.spines.values():
        spine.set_color(GRID)


def plot_waveforms():
    samples = read_csv("waveforms.csv")
    trajectories = {}
    for sample in samples:
        key = sample["topology"], sample["case"]
        point = (float(sample["t"]), *(float(sample[f"R{i}"]) for i in range(4)))
        trajectories.setdefault(key, []).append(point)
    fig, axes = plt.subplots(4, 3, figsize=(12.8, 10), sharey=True)
    fig.subplots_adjust(left=0.125, right=0.975, top=0.83, bottom=0.09, wspace=0.16, hspace=0.40)
    fig.text(0.125, 0.953, "Four-site occupation trajectories", fontsize=23, weight="bold", color=INK)
    fig.text(0.125, 0.917, r"$N=4,\;\Omega=1,\;J=0.7$    ·    Amplitude floor 0.05    ·    Short-window choices are held fixed to $8\pi$.",
             fontsize=11, color=MUTED)
    handles = [Line2D([], [], label=f"Site {site}", color=color, linestyle=style, linewidth=2)
               for site, (color, style) in enumerate(zip(SITE_COLORS, SITE_STYLES))]
    fig.legend(handles=handles, ncols=4, frameon=False, loc="upper left", bbox_to_anchor=(0.118, 0.902),
               labelcolor=MUTED, fontsize=10, columnspacing=1.7)

    for topology_index, (topology, topology_label) in enumerate(TOPOLOGIES):
        for window_index, (window, window_label, start, end) in enumerate(WINDOWS):
            row = 2 * topology_index + window_index
            for column, (case, title, _) in enumerate(CASES):
                axis = axes[row, column]
                t, *rates = np.asarray(sorted(trajectories[(topology, case)])).T
                style_wave_axis(axis)
                for rate, color, style in zip(rates, SITE_COLORS, SITE_STYLES):
                    axis.plot(t, rate, color=color, linestyle=style, linewidth=1.7)
                axis.set_xlim(start, end)
                if window == "short":
                    axis.set_xticks([0, 1, 2, 3, SHORT_END], ["0", "1", "2", "3", "3.688"])
                else:
                    axis.set_xticks(np.arange(1, 5) * 2 * np.pi, [r"$2\pi$", r"$4\pi$", r"$6\pi$", r"$8\pi$"])
                if row == 0:
                    axis.set_title(title, loc="left", weight="bold", fontsize=11, color=INK, pad=9)
                if column == 0:
                    axis.set_ylabel(f"{topology_label}\n{window_label}", rotation=0, ha="right", va="center",
                                    labelpad=13, fontsize=10, weight="bold", color=INK)
    fig.supxlabel(r"Time $t$", x=0.55, y=0.046, fontsize=10, color=MUTED)
    fig.savefig(DIRECTORY / "occupations.png", dpi=175, facecolor="white")
    plt.close(fig)


def phase_label(value):
    if abs(value) >= 179.95:
        return "180°"
    if abs(value) < 0.05:
        return "0°"
    return f"{value:+.1f}°"


def numeric(value):
    return float(value) if value else np.nan


def draw_heatmap(axis, values, colors, minimum, maximum, labels, title, formatter, *, absolute=False):
    palette = LinearSegmentedColormap.from_list(title, colors)
    palette.set_bad("#edf0f3")
    fill = np.abs(values) if absolute else values
    image = axis.imshow(fill, cmap=palette, vmin=minimum, vmax=maximum, aspect="auto", interpolation="nearest")
    for (row, column), value in np.ndenumerate(values):
        rgba = image.cmap(image.norm(fill[row, column]))
        brightness = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
        color = "white" if brightness < 0.52 else INK
        label = formatter(value) if np.isfinite(value) else "—"
        axis.text(column, row, label, ha="center", va="center", color=color, fontsize=8)
    axis.set_xticks(range(len(labels)), labels)
    axis.xaxis.tick_top()
    axis.set_yticks(range(len(values)), [""] * len(values))
    axis.set_title(title, loc="left", color=INK, weight="bold", fontsize=11, pad=30)
    axis.tick_params(axis="both", length=0, pad=6, labelsize=9, colors=MUTED)
    for divider in (2.5, 5.5, 8.5):
        axis.axhline(divider, color="white", linewidth=2)
    for spine in axis.spines.values():
        spine.set_visible(False)


def plot_metrics():
    samples = [row for row in read_csv("comparison.csv") if float(row["gate"]) == GATE]
    records = {(r["topology"], r["window"], r["case"], r["pair"]): r for r in samples}
    order = [(topology, window, case) for topology, _ in TOPOLOGIES for window, *_ in WINDOWS for case, *_ in CASES]
    case_labels = [label for _, _, label in CASES] * 4
    pair_data = [[records[(*key, pair)] for pair in PAIR_LABELS] for key in order]
    solo_data = [pairs[0] for pairs in pair_data]
    fig, axes = plt.subplots(1, 5, figsize=(16.8, 7), gridspec_kw={"width_ratios": (6, 6, 6, 4, 1.55)})
    fig.subplots_adjust(left=0.185, right=0.985, top=0.75, bottom=0.08, wspace=0.105)
    fig.text(0.05, 0.93, "Pairwise phase, locking and occupation amplitude", fontsize=23, weight="bold", color=INK)
    fig.text(0.05, 0.883, "Short and disjoint holdout windows  ·  Main amplitude floor 0.05  ·  Signed lag = phase of first site minus second",
             fontsize=11, color=MUTED)
    fig.text(0.05, 0.847, "Amplitude = half the peak-to-peak range  ·  Coverage = fraction of time all four phase envelopes meet the floor",
             fontsize=11, color=MUTED)
    phase_colors = ("#17664e", "#e6f0df", "#f0e1cb", "#9b423c")
    unit_colors = ("#f4eadf", "#ddece9", "#17664e")
    corr_colors = ("#9b423c", "#f3f0e8", "#17664e")

    for axis, (field, colors, minimum, maximum, title, formatter) in zip(axes[:3], (
        ("lag_deg", phase_colors, 0.0, 180.0, "Mean phase lag", phase_label),
        ("plv", unit_colors, 0.0, 1.0, "Phase-locking value", lambda value: f"{value:.2f}"),
        ("pearson", corr_colors, -1.0, 1.0, "Pearson correlation", lambda value: f"{value:.2f}"),
    )):
        values = np.asarray([[numeric(row[field]) for row in pairs] for pairs in pair_data])
        draw_heatmap(axis, values, colors, minimum, maximum, PAIR_LABELS, title, formatter, absolute=field == "lag_deg")

    amplitudes = np.asarray([[numeric(row[f"amplitude{i}"]) for i in range(4)] for row in solo_data])
    draw_heatmap(axes[3], amplitudes, unit_colors, 0.0, 0.5, [str(i) for i in range(4)], "Site amplitude",
                 lambda value: f"{value:.2f}")
    coverages = np.asarray([[numeric(row["coverage"])] for row in solo_data])
    draw_heatmap(axes[4], coverages, unit_colors, 0.0, 1.0, ["%"], "Coverage",
                 lambda value: f"{100 * value:.0f}")
    axes[0].set_yticklabels(case_labels)
    for group_index, (_, topology_label) in enumerate(TOPOLOGIES):
        for window_index, (_, window_label, *_ends) in enumerate(WINDOWS):
            y = 3 * (2 * group_index + window_index) + 1
            axes[0].text(-0.70, y, f"{topology_label}\n{window_label}", transform=axes[0].get_yaxis_transform(),
                         ha="left", va="center", color=INK, weight="bold", fontsize=10)
    fig.savefig(DIRECTORY / "pair-metrics.png", dpi=175, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    plot_waveforms()
    plot_metrics()
