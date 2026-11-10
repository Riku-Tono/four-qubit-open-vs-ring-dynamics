"""Plot the fixed open-chain states and the verification of their phase persistence."""

import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D


OUTPUT = Path(__file__).resolve().parent
GATE = 0.05
PAIR_LABELS = ("01", "02", "03", "12", "13", "23")
INK = "#172638"
MUTED = "#4c5c6b"
GRID = "#dce4e8"
RED = "#a83e35"
SITE_COLORS = ("#14609e", "#c74b22", "#117e65", "#8b49a2")
SITE_STYLES = ("solid", (0, (5, 2.8)), (2, (1.3, 2.0)), (3, (7, 2.5, 1, 2.5)))
CASES = ("prior_short", "selected_0.05")
CASE_LABELS = ("Previous short-window choice", "New multi-window choice")


def read_csv(name):
    with (OUTPUT / name).open(encoding="utf-8", newline="") as stream:
        yield from csv.DictReader(stream)


def style_axis(axis):
    axis.set_axisbelow(True)
    axis.grid(color=GRID, linewidth=0.7)
    axis.tick_params(length=0, pad=5, labelsize=9, colors=MUTED)
    for spine in axis.spines.values():
        spine.set_color(GRID)


def time_axis(axis, start, end, step):
    first = int(np.ceil(start / (step * np.pi)))
    last = int(np.floor(end / (step * np.pi)))
    multiples = np.arange(first, last + 1) * step
    axis.set_xticks(np.pi * multiples, ["0" if k == 0 else rf"${k}\pi$" for k in multiples])
    axis.set_xlim(start, end)


def signed_angle(value):
    if abs(value) >= 179.995:
        return "180°"
    if abs(value) < 0.005:
        return "0°"
    return f"{value:+.2f}°"


def heatmap(axis, values, labels, title, colors, low, high, formatter, dividers, *, absolute=False):
    palette = LinearSegmentedColormap.from_list(title, colors)
    palette.set_bad("#edf0f3")
    fill = np.abs(values) if absolute else values
    image = axis.imshow(fill, cmap=palette, vmin=low, vmax=high, aspect="auto", interpolation="nearest")
    for (row, column), value in np.ndenumerate(values):
        rgba = image.cmap(image.norm(fill[row, column]))
        brightness = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
        label = formatter(value) if np.isfinite(value) else "—"
        axis.text(column, row, label, ha="center", va="center", fontsize=8.8,
                  color="white" if brightness < 0.52 else INK)
    axis.set_xticks(range(len(labels)), labels)
    axis.xaxis.tick_top()
    axis.set_yticks(range(len(values)), [""] * len(values))
    axis.set_title(title, loc="left", color=INK, weight="bold", fontsize=11, pad=30)
    axis.tick_params(axis="both", length=0, pad=6, labelsize=9, colors=MUTED)
    for row in dividers:
        axis.axhline(row, color="white", linewidth=2)
    for spine in axis.spines.values():
        spine.set_visible(False)


def waveform_data():
    trajectories = {case: [] for case in CASES}
    for row in read_csv("waveforms.csv"):
        case = row["case"]
        if case in trajectories:
            t = float(row["t"])
            occupation = [float(row[f"R{i}"]) for i in range(4)]
            envelope = min(float(row[f"envelope{i}"]) for i in range(4))
            trajectories[case].append((t, *occupation, envelope))
    return {case: np.asarray(sorted(values)).T for case, values in trajectories.items()}


def plot_occupations(windows):
    trajectories = waveform_data()
    search = [row for row in windows if row["role"] == "search"]
    holdout = [row for row in windows if row["role"] == "holdout"]
    time_slices = (search[0], holdout[0])
    final_time = float(holdout[-1]["end"])
    holdout_start = float(holdout[0]["start"])
    split = float(next(row for row in windows if row["window"] == "holdout_B")["start"])
    fig = plt.figure(figsize=(12.8, 9.2))
    layout = fig.add_gridspec(4, 2, left=0.105, right=0.975, bottom=0.085, top=0.80,
                             hspace=0.48, wspace=0.18, height_ratios=(1.65, 1.65, 1, 1))
    fig.text(0.105, 0.953, "Open-chain occupation and amplitude persistence", fontsize=23, weight="bold", color=INK)
    fig.text(0.105, 0.916, r"$N=4,\;\Omega=1,\;J=0.7$     ·     Both initial states are held fixed through $256\pi$.",
             fontsize=11, color=MUTED)
    handles = [Line2D([], [], label=f"Site {site}", color=color, linestyle=style, linewidth=2)
               for site, (color, style) in enumerate(zip(SITE_COLORS, SITE_STYLES))]
    fig.legend(handles=handles, ncols=4, frameon=False, loc="upper left", bbox_to_anchor=(0.10, 0.902),
               labelcolor=MUTED, fontsize=10, columnspacing=1.7)

    for index, (case, label) in enumerate(zip(CASES, CASE_LABELS)):
        t, r0, r1, r2, r3, envelope = trajectories[case]
        for column, interval in enumerate(time_slices):
            axis = fig.add_subplot(layout[index, column])
            style_axis(axis)
            start, end = float(interval["start"]), float(interval["end"])
            keep = (t >= start - 0.021) & (t <= end + 0.021)
            for occupation, color, style in zip((r0, r1, r2, r3), SITE_COLORS, SITE_STYLES):
                axis.plot(t[keep], occupation[keep], color=color, linestyle=style, linewidth=1.5)
            time_axis(axis, start, end, 2)
            axis.set_ylim(-0.025, 1.025)
            axis.set_yticks([0, 0.25, 0.5, 0.75, 1])
            if index == 0:
                axis.set_title(("First selection block", "First untouched block")[column],
                               loc="left", color=INK, weight="bold", fontsize=11, pad=10)
            if column == 0:
                axis.set_ylabel(r"Occupation $R_i(t)$", color=MUTED)
                axis.text(0.02, 0.94, label, transform=axis.transAxes, color=INK, ha="left", va="top",
                          fontsize=9.2, weight="bold", bbox=dict(facecolor="white", edgecolor="none", alpha=0.88))
        axis = fig.add_subplot(layout[index + 2, :])
        style_axis(axis)
        axis.axvspan(holdout_start, split, color="#e9f0f7", zorder=0)
        axis.axvspan(split, final_time, color="#eef5f0", zorder=0)
        axis.plot(t, envelope, color="#1d637c", linewidth=0.65)
        axis.fill_between(t, envelope, GATE, where=envelope < GATE, color=RED, alpha=0.44,
                          interpolate=True, linewidth=0)
        axis.axhline(GATE, color=RED, linewidth=0.9, linestyle="--")
        axis.axvline(holdout_start, color=MUTED, linewidth=0.9)
        axis.axvline(split, color=MUTED, linewidth=0.7, linestyle=":")
        axis.set_ylim(0, max(np.max(trajectories[key][-1]) for key in CASES) * 1.04)
        time_axis(axis, 0, final_time, 32)
        axis.set_ylabel(r"Min $|u_i(t)|$", color=MUTED)
        axis.set_title(f"{label} · minimum phase envelope", loc="left", weight="bold", fontsize=10, color=INK, pad=7)
        if index == 1:
            axis.set_xlabel(r"Time $t$    ·    Red indicates an envelope below 0.05", color=MUTED, labelpad=8)
        else:
            for start, end, stage in ((0, holdout_start, "Selection"), (holdout_start, split, "Holdout A"),
                                      (split, final_time, "Holdout B")):
                axis.text((start + end) / 2, 0.96, stage, transform=axis.get_xaxis_transform(), color=MUTED,
                          ha="center", va="top", fontsize=8.7,
                          bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, pad=1.5))

    fig.savefig(OUTPUT / "occupations-envelope.png", dpi=175, facecolor="white")
    plt.close(fig)


def numeric(value):
    return float(value) if value else np.nan


def plot_metrics(records, windows):
    focus = [row for row in windows if row["role"] in ("search", "holdout_aggregate")]
    heldout = [row for row in windows if row["role"] == "holdout"]
    worst_lag = max(heldout, key=lambda row: max(abs(float(records[CASES[1], row["window"], pair]["lag_deg"] or "inf"))
                                               for pair in PAIR_LABELS))
    worst_plv = min(heldout, key=lambda row: min(float(records[CASES[1], row["window"], pair]["plv"] or 0)
                                                for pair in PAIR_LABELS))
    focus.extend(row for row in heldout if row["window"] in {worst_lag["window"], worst_plv["window"]})
    labels = []
    rows = []
    for row in focus:
        window = row["window"]
        display = window.replace("_", " ").title()
        start, end = (round(float(row[key]) / np.pi) for key in ("start", "end"))
        start_label = "0" if start == 0 else f"{start}π"
        for case, name in zip(CASES, ("prior", "new")):
            rows.append([records[case, window, pair] for pair in PAIR_LABELS])
            labels.append(f"{display}  {start_label}–{end}π  ·  {name}")
    group_dividers = [i - 0.5 for i in range(2, len(rows), 2)]
    summaries = [pairs[0] for pairs in rows]
    fig, axes = plt.subplots(1, 3, figsize=(15.4, 4 + len(rows) * 0.24),
                             gridspec_kw={"width_ratios": (6, 6, 4.3)})
    fig.subplots_adjust(left=0.18, right=0.98, top=0.78, bottom=0.06, wspace=0.09)
    fig.text(0.06, 0.946, "Six-pair phase agreement across selection and holdout", fontsize=23, weight="bold", color=INK)
    fig.text(0.06, 0.909, "Amplitude floor 0.05  ·  Prior = previous short-window choice  ·  New = fixed multi-window choice",
             fontsize=10.5, color=MUTED)
    fig.text(0.06, 0.877, "The 8π holdout rows show the new state's largest absolute mean lag and lowest PLV.",
             fontsize=10.5, color=MUTED)
    fig.text(0.06, 0.846, "Drop = percent of time any envelope is below 0.05  ·  Longest = duration of the longest such dip",
             fontsize=10.5, color=MUTED)

    phase_colors = ("#17664e", "#e2eee4", "#f4ecd5", "#9b423c")
    plv_colors = ("#9b423c", "#f4eadf", "#d6e7e2", "#17664e")
    lag = np.asarray([[numeric(row["lag_deg"]) for row in pairs] for pairs in rows])
    plv = np.asarray([[numeric(row["plv"]) for row in pairs] for pairs in rows])
    heatmap(axes[0], lag, PAIR_LABELS, "Signed mean phase lag", phase_colors, 0, 30, signed_angle,
            group_dividers, absolute=True)
    heatmap(axes[1], plv, PAIR_LABELS, "Phase-locking value", plv_colors, 0.85, 1,
            lambda v: f"{v:.4f}", group_dividers)
    axes[0].set_yticklabels(labels)

    extras = np.asarray([[min(float(r[f"amplitude{i}"]) for i in range(4)),
                          100 * (1 - float(r["coverage"])), float(r["longest_dropout"])] for r in summaries])
    table = axes[2]
    scaled = np.column_stack((np.minimum(extras[:, 0] / 0.15, 1), 1 - np.minimum(extras[:, 1] / 20, 1),
                              1 - np.minimum(extras[:, 2] / 5, 1)))
    heatmap(table, scaled, ("Min amp.", "Drop %", "Longest"), "Amplitude and dropouts", plv_colors, 0, 1,
            lambda _: "", group_dividers)
    for (i, j), value in np.ndenumerate(extras):
        text = f"{value:.3f}" if j == 0 else f"{value:.2f}%" if j == 1 else f"{value:.2f}"
        color = "white" if scaled[i, j] < 0.22 or scaled[i, j] > 0.91 else INK
        table.text(j, i, text, ha="center", va="center", color=color, fontsize=8.8)

    for axis in axes:
        axis.axhline(2 * sum(row["role"] == "search" for row in focus) - 0.5, color=INK, linewidth=1.4)
    fig.savefig(OUTPUT / "pair-metrics.png", dpi=175, facecolor="white")
    plt.close(fig)


def main():
    data = [row for row in read_csv("comparison.csv") if float(row["gate"]) == GATE and row["case"] in CASES]
    records = {(row["case"], row["window"], row["pair"]): row for row in data}
    windows = [row for row in data if row["case"] == CASES[1] and row["pair"] == PAIR_LABELS[0]]
    windows.sort(key=lambda row: float(row["start"]))
    plot_occupations(windows)
    plot_metrics(records, windows)


if __name__ == "__main__":
    main()
