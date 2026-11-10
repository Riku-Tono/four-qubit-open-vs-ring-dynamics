"""Plot the fixed-direction PLV scan and frequency ablation at its onset."""

import argparse
import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator


OUTPUT = Path(__file__).resolve().parent
INK, MUTED, GRID = "#192c3d", "#526675", "#dce4eb"
TEAL, BLUE, RUST, SLATE = "#166e59", "#21679b", "#a44945", "#657785"
MAIN_GATE, PLV_FLOOR = 0.05, 0.95
CARRIER, TARGET = 0.502392496, 0.930636


def read_csv(name):
    with (OUTPUT / name).open(encoding="utf-8", newline="") as file:
        return list(csv.DictReader(file))


def number(row, key):
    value = row[key]
    return float(value) if value not in (None, "") else np.nan


def panel(axis, title, ylabel=None):
    axis.set_axisbelow(True)
    axis.grid(axis="y", color=GRID, linewidth=0.75)
    axis.tick_params(length=0, pad=6, labelsize=9, colors=MUTED)
    axis.spines[["top", "right"]].set_visible(False)
    axis.spines[["bottom", "left"]].set_color(GRID)
    axis.set_title(title, loc="left", color=INK, weight="bold", fontsize=11.5, pad=12)
    if ylabel:
        axis.set_ylabel(ylabel, color=MUTED, fontsize=10)
    axis.xaxis.label.set_color(MUTED)


def header(fig, title, subtitle):
    fig.text(0.095, 0.955, title, color=INK, fontsize=22, weight="bold", va="top")
    fig.text(0.095, 0.897, subtitle, color=MUTED, fontsize=10.5, va="top")


def save(fig, name):
    fig.savefig(OUTPUT / name, dpi=178, facecolor="white")
    plt.close(fig)


def curve(axis, rows, key, color, linestyle="-", label=None, marker="o"):
    ordered = sorted(rows, key=lambda row: float(row["width"]))
    axis.plot([number(row, "width") for row in ordered], [number(row, key) for row in ordered],
              color=color, linestyle=linestyle, marker=marker, markersize=3.1, lw=1.6, label=label)


def grouped_traces(rows):
    by_witness = {}
    for row in rows:
        key = (row["direction"], row["window"], row["pair"])
        by_witness.setdefault(key, []).append(row)
    groups = []
    for (direction, window, pair), trace in sorted(by_witness.items()):
        trace.sort(key=lambda row: float(row["width"]))
        values = np.array([(number(row, "width"), number(row, "plv")) for row in trace])
        for directions, other_window, other_pair, _, other in groups:
            if (window == other_window and pair == other_pair and values.shape == other.shape
                    and np.allclose(values, other, rtol=0, atol=1e-10)):
                directions.append(direction)
                break
        else:
            groups.append(([direction], window, pair, trace, values))
    return groups


def plot_onset(summary, traces):
    by_gate = {gate: sorted([row for row in summary if float(row["gate"]) == gate],
                            key=lambda row: number(row, "width")) for gate in sorted({float(row["gate"]) for row in summary})}
    primary = by_gate[MAIN_GATE]
    first = min(number(row, "width") for row in primary if number(row, "min_plv") < PLV_FLOOR)
    before = max(number(row, "width") for row in primary if number(row, "width") < first)
    by_gate_overlap = all(np.allclose([number(row, "min_plv") for row in series],
                                     [number(row, "min_plv") for row in primary], rtol=0, atol=1e-10)
                          and [row["failing_directions"] for row in series] == [row["failing_directions"] for row in primary]
                          for series in by_gate.values())
    fig, axes = plt.subplot_mosaic([["global", "witness"], ["counts", "counts"]], figsize=(12.4, 7.6),
                                   height_ratios=(1.25, 0.8))
    fig.subplots_adjust(left=0.095, right=0.965, bottom=0.111, top=0.814, wspace=0.26, hspace=0.63)
    caption = "496 fixed directions  ·  Same 24 holdout windows  ·  Six pairs per window"
    if by_gate_overlap:
        caption += "  ·  Identical results at all three floors"
    header(fig, f"The first PLV failure appears at {first:.3f} rad", caption)
    global_axis, witness_axis, counts = (axes[name] for name in ("global", "witness", "counts"))
    panel(global_axis, "A   Lowest PLV over all directions, windows and pairs", "Minimum PLV")
    panel(witness_axis, "B   First failing pairs and windows", "Pair PLV")
    panel(counts, "C   Directions failing at least one window", "Number of directions")
    if by_gate_overlap:
        curve(global_axis, primary, "min_plv", INK)
        curve(counts, primary, "failing_directions", BLUE, marker="s")
    else:
        for gate, color in zip(by_gate, (SLATE, BLUE, TEAL)):
            label = f"Floor {gate:g}"
            curve(global_axis, by_gate[gate], "min_plv", color, label=label)
            curve(counts, by_gate[gate], "failing_directions", color, label=label, marker="s")
        global_axis.legend(frameon=False, fontsize=8)
    width_first = next(row for row in primary if number(row, "width") == first)
    global_axis.scatter(first, number(width_first, "min_plv"), edgecolor=RUST, facecolor="white", s=47, zorder=4)
    groups = grouped_traces([row for row in traces if number(row, "gate") == MAIN_GATE])
    for i, (directions, window, pair, trace, _) in enumerate(groups):
        color = (TEAL, BLUE, SLATE, RUST)[i % 4]
        label = f"{' / '.join(directions)}  ·  {pair}, {window.replace('_', ' ')}"
        curve(witness_axis, trace, "plv", color, linestyle=("-", "--", "-.", ":")[i % 4], label=label)
    witness_axis.legend(loc="upper right", frameon=True, framealpha=1, facecolor="white", edgecolor="white",
                        fontsize=8, labelcolor=MUTED)
    for axis in (global_axis, witness_axis):
        axis.axhline(PLV_FLOOR, color=RUST, linestyle=":", lw=1.3)
        axis.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
        axis.yaxis.set_major_locator(MaxNLocator(5))
    witness_axis.set_xlim(before - 0.0034, first + 0.0034)
    witness_axis.set_xticks(np.arange(round((before - 0.003) * 1000), round((first + 0.003) * 1000) + 1, 2) / 1000)
    witness_axis.relim()
    selected = [value for *_, data in groups for width, value in data if before - 0.0034 <= width <= first + 0.0034]
    if selected:
        lo, hi = min(selected + [PLV_FLOOR]), max(selected + [PLV_FLOOR])
        witness_axis.set_ylim(lo - 0.15 * (hi - lo), hi + 0.25 * (hi - lo))
    counts.yaxis.set_major_locator(MaxNLocator(5, integer=True))
    counts.set_ylim(bottom=-0.3)
    counts.annotate(f"{int(width_first['failing_directions'])} first failures", (first, number(width_first, "failing_directions")),
                    xytext=(-15, 19), textcoords="offset points", ha="right", color=RUST, fontsize=9,
                    arrowprops={"arrowstyle": "-", "color": RUST, "linewidth": 0.8})
    global_axis.set_xlim(number(primary[0], "width") - 0.0007, number(primary[-1], "width") + 0.0007)
    counts.set_xlim(global_axis.get_xlim())
    for axis in axes.values():
        axis.axvspan(before, first, color="#eaf1f4", alpha=0.95, zorder=0)
        axis.set_xlabel("Error half-width (rad)", fontsize=10)
        axis.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    fig.text(0.095, 0.036, f"Shading marks the adjacent tested widths {before:.3f} and {first:.3f} rad. Dotted horizontal lines mark PLV 0.95.",
             color=MUTED, fontsize=9)
    save(fig, "plv-onset.png")


def witness(rows):
    failures = [row for row in rows if row["variant"] == "full" and float(row["gate"]) == MAIN_GATE
                and number(row, "plv") < PLV_FLOOR]
    return min(failures, key=lambda row: (float(row["width"]), round(number(row, "plv"), 12),
                                          row["direction"], row["window"], row["pair"]))


def plot_ablations(rows):
    point = witness(rows)
    keys = ("case", "width", "direction", "window", "pair", "gate")
    group = [row for row in rows if all(row[k] == point[k] for k in keys)]
    main = {row["variant"]: row for row in group if row["variant"] != "omit_line"}
    labels = (("full", "All lines"), ("align_carrier_phase", "Align carrier phases"),
              ("restore_carrier_magnitude", "Restore carrier magnitudes"),
              ("carrier_plus_0930636", "Carrier + 0.930636 only"),
              ("omit_0930636", "Omit 0.930636"),
              ("restore_normalized_0930636", "Reset 0.930636 to nominal ratio"),
              ("carrier_only", "Carrier only"))
    labels = [(key, label) for key, label in labels if key in main]
    values = np.array([number(main[key], "original_mask_plv") for key, _ in labels])
    baseline = number(point, "original_mask_plv")
    width = float(point["width"])
    window = point["window"].replace("_", " ")
    fig, (effects, lines) = plt.subplots(1, 2, figsize=(12.6, 5.6), gridspec_kw={"width_ratios": (1.12, 1)})
    fig.subplots_adjust(left=0.209, right=0.96, bottom=0.19, top=0.74, wspace=0.44)
    header(fig, "Which spectral line drives the PLV loss?",
           f"{point['direction']} at {width:.3f} rad  ·  Pair {point['pair']}  ·  {window}  ·  Original eligible-time mask")

    panel(effects, "A   Targeted spectrum changes")
    effects.grid(axis="y", visible=False)
    effects.grid(axis="x", color=GRID, linewidth=0.75)
    effects.axvline(PLV_FLOOR, color=RUST, linestyle=":", linewidth=1.35)
    for y, ((key, _), value) in enumerate(zip(labels, values)):
        color = RUST if value < PLV_FLOOR else TEAL if key == "omit_0930636" else BLUE
        effects.hlines(y, min(value, PLV_FLOOR), max(value, PLV_FLOOR), color=color, linewidth=2, alpha=0.5)
        effects.scatter(value, y, color=color, s=33, zorder=3)
        effects.annotate(f"{value:.6f}", (value, y), xytext=(7, 0), textcoords="offset points",
                         ha="left", va="center", color=INK, fontsize=8.5)
    effects.set_yticks(range(len(labels)), [label for _, label in labels])
    effects.set_ylim(len(labels) - 0.35, -0.65)
    effects.set_xlim(min(values.min(), PLV_FLOOR) - 0.004, 1.018)
    effects.set_xticks([0.95, 0.96, 0.97, 0.98, 0.99, 1.0])
    effects.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    effects.set_xlabel("Pair PLV")

    modes = [(number(row, "frequency"), number(row, "original_mask_plv") - baseline)
             for row in group if row["variant"] == "omit_line" and abs(number(row, "frequency") - CARRIER) > 1e-6]
    modes = sorted((item for item in modes if abs(item[1]) >= 5e-6), key=lambda item: abs(item[1]), reverse=True)[:6]
    panel(lines, "B   Omit one non-carrier line", "Omitted angular frequency")
    lines.grid(axis="y", visible=False)
    lines.grid(axis="x", color=GRID, linewidth=0.75)
    lines.axvline(0, color=SLATE, lw=0.9)
    for y, (omega, change) in enumerate(modes):
        color = TEAL if abs(omega - TARGET) < 1e-6 else BLUE if change >= 0 else RUST
        lines.barh(y, change, height=0.51, color=color, alpha=0.85)
        lines.annotate(f"{change:+.5f}", (change, y), xytext=(5 if change >= 0 else -5, 0),
                       textcoords="offset points", ha="left" if change >= 0 else "right", va="center",
                       color=INK, fontsize=8.5)
    lines.set_yticks(range(len(modes)), [f"{omega:.6f}" for omega, _ in modes])
    lines.set_ylim(len(modes) - 0.35, -0.65)
    magnitude = max(max((abs(v) for _, v in modes), default=0.001), 0.001)
    lines.set_xlim(min(0, *(v for _, v in modes)) - magnitude * 0.25,
                   max(0, *(v for _, v in modes)) + magnitude * 0.29)
    lines.xaxis.set_major_locator(MaxNLocator(5))
    lines.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    lines.set_xlabel("Change in pair PLV relative to all lines")
    fig.text(0.095, 0.072, "All variants use the original valid times. The dotted line marks PLV 0.95. Panel B shows changes with magnitude ≥ 0.000005.",
             color=MUTED, fontsize=9)
    save(fig, "onset-ablations.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("section", nargs="?", choices=("all", "sweep", "diagnostics"), default="all")
    args = parser.parse_args()
    if args.section in ("all", "sweep"):
        plot_onset(read_csv("width_summary.csv"), read_csv("witness_trace.csv"))
    if args.section in ("all", "diagnostics"):
        plot_ablations(read_csv("mechanism.csv"))


if __name__ == "__main__":
    main()
