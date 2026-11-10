"""Render the independent-angle sweep and the selected open-chain diagnostics."""

import argparse
import csv
from itertools import combinations
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, PercentFormatter


OUTPUT = Path(__file__).resolve().parent
INK, MUTED, GRID = "#192c3d", "#526675", "#dce4eb"
NOMINAL, ALTERED, LIMIT = "#627486", "#126577", "#a44945"
SITE_STYLES = (("#235c9d", "-", "o"), ("#b85128", "--", "s"),
               ("#08755f", "-.", "D"), ("#79568f", ":", "^"))
GATES = (0.025, 0.05, 0.10)
GATE_COLORS = ("#21679b", "#166e59", "#a44945")
PAIRS = tuple(f"{i}{j}" for i, j in combinations(range(4), 2))


def read_csv(name):
    with (OUTPUT / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def column(rows, name):
    return np.array([float(row[name]) if row[name] not in ("", None) else np.nan for row in rows])


def panel(axis, title, ylabel=None):
    axis.set_axisbelow(True)
    axis.grid(axis="y", color=GRID, linewidth=0.75)
    axis.tick_params(length=0, pad=5, labelsize=9, colors=MUTED)
    axis.spines[["top", "right"]].set_visible(False)
    for side in ("bottom", "left"):
        axis.spines[side].set_color(GRID)
    axis.set_title(title, loc="left", color=INK, weight="bold", fontsize=11.5, pad=9)
    if ylabel:
        axis.set_ylabel(ylabel, color=MUTED, fontsize=10)
    axis.xaxis.label.set_color(MUTED)


def header(fig, title, subtitle):
    fig.text(0.085, 0.955, title, color=INK, fontsize=22, weight="bold", va="top")
    fig.text(0.085, 0.908, subtitle, color=MUTED, fontsize=10.5, va="top")


def save(fig, name):
    fig.savefig(OUTPUT / name, dpi=172, facecolor="white")
    plt.close(fig)


def site_legend(fig):
    handles = [Line2D([], [], color=color, linestyle=linestyle, marker=marker, label=f"Site {i}")
               for i, (color, linestyle, marker) in enumerate(SITE_STYLES)]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.078, 0.881), frameon=False,
               ncols=4, fontsize=9.5, labelcolor=MUTED, columnspacing=1.7)


def pretty_case(alias):
    labels = {"first-failure": "First detected failure", "worst-g005": "Worst sampled case",
              "first-coverage": "First detected coverage failure"}
    return labels.get(alias, alias.replace("_", " ").replace("-", " ").capitalize())


def case_rows(rows, name, window=None):
    return [row for row in rows if row["case"] == name and (window is None or row["window"] == window)]


def phase_extent(rows):
    return np.max(np.abs(np.column_stack([column(rows, f"phase_{pair}") for pair in PAIRS])), axis=1)


def plot_waveforms(alias, meta, rows):
    name = pretty_case(alias)
    window = meta["plot_window"]
    gate = float(meta["gate"])
    nominal = case_rows(rows, "nominal", window)
    altered = case_rows(rows, alias, window)
    groups = (nominal, altered)
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.8), sharex=True)
    fig.subplots_adjust(left=0.085, right=0.973, bottom=0.093, top=0.789, wspace=0.23, hspace=0.39)
    width = f"   ·   Error width {float(meta['width']):g} rad" if "width" in meta else ""
    header(fig, f"{name} · four-site occupation and phase",
           f"Fixed initial-angle errors in {window.replace('_', ' ')}   ·   Amplitude floor {gate:g}{width}")
    site_legend(fig)
    for i, (axis, group) in enumerate(zip(axes[0], groups)):
        panel(axis, ("A   Frozen nominal state", "B   Perturbed state")[i], r"Occupation $R_i(t)$" if i == 0 else None)
        time = column(group, "t") / np.pi
        for site, (color, linestyle, _) in enumerate(SITE_STYLES):
            axis.plot(time, column(group, f"R{site}"), color=color, linestyle=linestyle, lw=1.48)
        axis.set_xlim(time[0], time[-1])
        axis.set_ylim(-0.01, 1.01)
        axis.set_yticks(np.arange(5) / 4)

    envelope, phase = axes[1]
    panel(envelope, "C   Weakest envelope across four sites", r"$\min_i|u_i(t)|$")
    panel(phase, "D   Largest instantaneous gap across six pairs", "Degrees")
    for label, group, color, linestyle in (("Nominal", nominal, NOMINAL, "--"),
                                            ("Perturbed", altered, ALTERED, "-")):
        time = column(group, "t") / np.pi
        radius = column(group, "min_envelope")
        envelope.plot(time, radius, color=color, linestyle=linestyle, lw=1.5, label=label)
        phase.plot(time, np.where(radius >= gate, phase_extent(group), np.nan), color=color,
                   linestyle=linestyle, lw=1.4, label=label)
        if label == "Perturbed":
            envelope.fill_between(time, radius, gate, where=radius < gate, interpolate=True,
                                  color=LIMIT, alpha=0.2, linewidth=0)
    envelope.axhline(gate, color=LIMIT, linestyle=":", lw=1.0, label="Envelope floor")
    for axis in (envelope, phase):
        axis.set_ylim(0, axis.get_ylim()[1] * 1.19)
        axis.set_xlabel(r"Time $t/\pi$", fontsize=10)
        axis.xaxis.set_major_locator(MaxNLocator(5, integer=True))
        axis.legend(loc="upper right", ncols=3 if axis is envelope else 2, fontsize=8, labelcolor=MUTED,
                    frameon=True, framealpha=1, facecolor="white", edgecolor="white", columnspacing=1.0)
    phase.text(0.01, 0.97, "Shown when all four envelopes meet the floor", transform=phase.transAxes,
               va="top", color=MUTED, fontsize=8)
    save(fig, f"waveforms-{alias.replace('_', '-')}.png")


def carrier_flags(rows):
    return np.array([row["is_carrier"].lower() in {"true", "1"} for row in rows])


def spectral_lines(axis, rows, residual):
    for site, (color, _, marker) in enumerate(SITE_STYLES):
        selected = [row for row in rows if int(row["site"]) == site]
        omega = column(selected, "omega")
        amplitude = column(selected, "amplitude")
        keep = amplitude > 1e-12
        if residual:
            keep &= ~carrier_flags(selected)
        # The fixed horizontal dodge is visual only. The CSV retains exact frequencies.
        x, y = omega[keep] + (site - 1.5) * 0.0035, amplitude[keep]
        axis.vlines(x, 0, y, color=color, lw=0.95, alpha=0.7)
        axis.scatter(x, y, marker=marker, s=15, color=color, zorder=3)


def carrier_caption(rows):
    carrier = [row for row, selected in zip(rows, carrier_flags(rows)) if selected]
    power = column(carrier, "power_ratio").min()
    phase = column(carrier, "phase_deg")
    differences = (phase[:, None] - phase[None, :] + 180) % 360 - 180
    return f"Carrier phase spread {np.max(np.abs(differences)):.3f}°   ·   Lowest carrier power share {power:.2%}"


def plot_spectrum(alias, meta, rows):
    name = pretty_case(alias)
    groups = (case_rows(rows, "nominal"), case_rows(rows, alias))
    entries = groups[0] + groups[1]
    frequency = column(entries, "omega")
    amplitudes = column(entries, "amplitude")
    ymax = [amplitudes.max(), amplitudes[~carrier_flags(entries)].max()]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.9), sharex=True, sharey="col")
    fig.subplots_adjust(left=0.085, right=0.969, bottom=0.121, top=0.785, wspace=0.25, hspace=0.46)
    header(fig, f"{name} · frequency decomposition",
           r"Positive energy gaps   ·   Common line at $\omega_*=0.502392496$   ·   Site markers are offset slightly for clarity")
    site_legend(fig)
    for i, group in enumerate(groups):
        state = ("Nominal", "Perturbed")[i]
        for j, axis in enumerate(axes[i]):
            suffix = ("full spectrum", "residual lines only")[j]
            panel(axis, f"{chr(65 + 2 * i + j)}   {state} · {suffix}", r"Amplitude $|c_i(\omega)|$")
            spectral_lines(axis, group, residual=bool(j))
            axis.set_xlim(0, frequency.max() + 0.08)
            axis.set_ylim(0, ymax[j] * 1.23)
            axis.yaxis.set_major_locator(MaxNLocator(5))
            if i == 1:
                axis.set_xlabel(r"Angular frequency $\omega$", fontsize=10)
            if j == 0:
                axis.text(0.97, 0.93, carrier_caption(group), transform=axis.transAxes,
                          ha="right", va="top", color=MUTED, fontsize=8.4)
            else:
                residuals = [row for row, carrier in zip(group, carrier_flags(group)) if not carrier]
                peak = max(residuals, key=lambda row: float(row["amplitude"]))
                axis.text(0.97, 0.93, f"Largest residual at ω = {float(peak['omega']):.5f}",
                          transform=axis.transAxes, ha="right", va="top", color=MUTED, fontsize=8.4)
    if "variants" in meta:
        gate = float(meta["gate"])
        variants = meta["variants"]
        variants_to_show = (("full", "Full"), ("carrier_only", "Carrier alone"),
                            ("align_carrier_phase", "Carrier phases aligned"),
                            ("restore_normalized_residual", "Nominal relative residuals"))
        counts = [f"{label} {variants[variant][str(gate)]['pass_windows']}/24" for variant, label in variants_to_show
                  if variant in variants]
        fig.text(0.085, 0.044, f"Window passes at floor {gate:g} with each spectrum's own mask   ·   "
                              + "   ·   ".join(counts), fontsize=8.2, color=MUTED)
    save(fig, f"spectrum-{alias.replace('_', '-')}.png")


def plot_diagnostics():
    data = json.loads((OUTPUT / "diagnostics.json").read_text(encoding="utf-8"))
    waves = read_csv("diagnostic_waveforms.csv")
    spectra = read_csv("diagnostic_spectrum.csv")
    for alias, meta in data["cases"].items():
        if alias != "nominal" and meta.get("plot", True):
            plot_waveforms(alias, meta, waves)
            plot_spectrum(alias, meta, spectra)


def sweep_rows(rows, gate, family):
    selected = [row for row in rows if float(row["gate"]) == gate and row["family"] == family]
    return sorted(selected, key=lambda row: float(row["width"]))


def rate_axis(axis, bottom=0):
    axis.set_ylim(bottom, 1.055)
    axis.yaxis.set_major_formatter(PercentFormatter(1))


def pass_curve(axis, rows, name, color, marker, linestyle="-"):
    width, rate = column(rows, "width"), column(rows, "all_pass_rate")
    axis.plot(width, rate, color=color, marker=marker, linestyle=linestyle, lw=1.55, markersize=3.9, label=name)
    failing = np.flatnonzero(rate < 1)
    if failing.size:
        i = failing[0]
        axis.scatter(width[i], rate[i], facecolor="white", edgecolor=color, marker=marker, s=60, lw=1.6, zorder=4)


def distribution(axis, case_rows, aggregated, metric, higher_is_worse):
    widths = column(aggregated, "width")
    extremes = column(aggregated, metric)
    groups = [column([row for row in case_rows if float(row["width"]) == width], metric) for width in widths]
    groups = [values[np.isfinite(values)] for values in groups]
    low, middle, high = np.array([np.quantile(values, [0.10, 0.50, 0.90]) if values.size
                                 else [np.nan, np.nan, np.nan] for values in groups]).T
    axis.fill_between(widths, low, high, color=ALTERED, alpha=0.13, linewidth=0)
    axis.plot(widths, middle, marker="o", color=ALTERED, markersize=3.7, lw=1.35)
    axis.plot(widths, extremes, marker="s", color=INK, linestyle="--", markersize=3.7, lw=1.2)
    finite = np.r_[low, high, extremes]
    finite = finite[np.isfinite(finite)]
    if finite.size and higher_is_worse:
        axis.set_ylim(bottom=0)


def plot_sweep(summary, cases):
    all_main = sweep_rows(summary, 0.05, "all")
    main_cases = [row for row in cases if float(row["gate"]) == 0.05]
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.4), sharex=True)
    fig.subplots_adjust(left=0.085, right=0.978, bottom=0.108, top=0.79, wspace=0.38, hspace=0.53)
    header(fig, "How far can independent angle errors grow?",
           "Four sites × two angles   ·   Fixed preparation, no retuning   ·   Same 24 held-out windows at every width")
    metrics = (
        ("max_abs_lag", "C   Largest absolute mean lag", "Degrees", True),
        ("min_plv", "D   Lowest phase-locking value", "PLV", False),
        ("min_envelope", "E   Lowest instantaneous envelope", r"$\min_i|u_i(t)|$", False),
        ("min_coverage", "F   Lowest eligible-time fraction", "Coverage", False),
    )
    pass_all, pass_family, *metric_axes = axes.flat
    panel(pass_all, "A   Cases passing all 24 windows", "All-window pass rate")
    panel(pass_family, "B   All-window pass rate by sampling scheme", "Floor 0.10")
    by_gate = [sweep_rows(summary, gate, "all") for gate in GATES]
    rates = [column(rows, "all_pass_rate") for rows in by_gate]
    if all(np.array_equal(rate, rates[0]) for rate in rates[1:]):
        pass_curve(pass_all, by_gate[0], "All floors (0.025, 0.05, 0.10)", GATE_COLORS[1], "o")
    else:
        for gate, rows, color, marker in zip(GATES, by_gate, GATE_COLORS, ("o", "s", "D")):
            pass_curve(pass_all, rows, f"Floor {gate:g}", color, marker)
    rate_axis(pass_all)
    for family, color, marker, line in (("corner", ALTERED, "s", "-"), ("uniform", NOMINAL, "^", "--")):
        pass_curve(pass_family, sweep_rows(summary, 0.10, family), family.capitalize(), color, marker, line)
    rate_axis(pass_family)
    for axis in (pass_all, pass_family):
        axis.legend(loc="lower left", frameon=False, fontsize=8.2, labelcolor=MUTED)

    for axis, (metric, title, ylabel, higher_is_worse) in zip(metric_axes, metrics):
        panel(axis, title, ylabel)
        distribution(axis, main_cases, all_main, metric, higher_is_worse)
    lag, plv, envelope, coverage = metric_axes
    lag.axhline(15, color=LIMIT, lw=1.0, linestyle=":")
    plv.axhline(0.95, color=LIMIT, lw=1.0, linestyle=":")
    envelope.axhline(0.05, color=LIMIT, lw=1.0, linestyle=":")
    envelope.set_ylim(bottom=0)
    coverage.axhline(0.90, color=LIMIT, lw=1.0, linestyle=":")
    coverage.yaxis.set_major_formatter(PercentFormatter(1))
    coverage.set_ylim(top=1.03)
    for axis in axes[-1]:
        axis.set_xlabel("Independent error width (rad)", fontsize=9.5)
    for axis in axes.flat:
        axis.xaxis.set_major_locator(MaxNLocator(7))
    fig.text(0.085, 0.042,
             "Open markers mark the first failure at a tested width.  C–F use floor 0.05.  "
             "Solid = median, shading = 10th–90th percentile, dashed = worst, dotted = criterion or floor.",
             color=MUTED, fontsize=9)
    save(fig, "sweep.png")


def plot_spectral_trends(axes, cases, widths):
    main = [row for row in cases if float(row["gate"]) == 0.05]
    groups = [[row for row in main if float(row["width"]) == width] for width in widths]
    carrier_phase, amplitude, power = axes
    panel(carrier_phase, "D   Carrier phase spread and full-signal lag", "Degrees")
    panel(amplitude, "E   Carrier amplitude and full envelope", "Amplitude")
    panel(power, "F   Weakest site's carrier power share", "Power share")
    comparisons = (
        (carrier_phase, "carrier_max_abs_lag", "max_abs_lag", "Carrier phase spread", "Full mean phase lag", max),
        (amplitude, "carrier_min_amplitude", "min_envelope", "Carrier amplitude", "Full minimum envelope", min),
    )
    for axis, carrier_field, full_field, carrier_label, full_label, extreme in comparisons:
        for field, label, color, marker, linestyle in ((carrier_field, carrier_label, ALTERED, "o", "-"),
                                                       (full_field, full_label, INK, "s", "--")):
            values = [extreme((float(row[field]) for row in group if row[field] not in ("", None)), default=np.nan)
                      for group in groups]
            axis.plot(widths, values, color=color, marker=marker, linestyle=linestyle, markersize=3.4, lw=1.3, label=label)
        axis.legend(loc="best", frameon=False, fontsize=8.0, labelcolor=MUTED)
        axis.set_ylim(bottom=0)
    carrier_phase.axhline(15, color=LIMIT, lw=1, linestyle=":")
    amplitude.axhline(0.05, color=LIMIT, lw=1, linestyle=":")
    per_case = [[min(float(row[f"carrier_power_fraction{i}"]) for i in range(4)) for row in group] for group in groups]
    bottom, median = np.array([(min(values), np.median(values)) for values in per_case]).T
    power.fill_between(widths, bottom, median, color=ALTERED, alpha=0.13, linewidth=0)
    power.plot(widths, median, marker="o", color=ALTERED, markersize=3.5, lw=1.35, label="Median across cases")
    power.plot(widths, bottom, marker="s", color=INK, linestyle="--", markersize=3.5, lw=1.2, label="Worst across cases")
    power.yaxis.set_major_formatter(PercentFormatter(1))
    power.set_ylim(top=1.018)
    power.legend(loc="best", frameon=False, fontsize=8.0, labelcolor=MUTED)


def plot_failure_criteria(summary, cases):
    criteria = (("fail_lag_cases", "Mean phase gap", "#235c9d", "o"),
                ("fail_plv_cases", "PLV", "#166e59", "s"),
                ("fail_coverage_cases", "Eligible time", "#a44945", "^"),
                ("fail_amplitude_cases", "Half-range", "#79568f", "D"))
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.0), sharex=True)
    fig.subplots_adjust(left=0.085, right=0.975, bottom=0.11, top=0.775, wspace=0.38, hspace=0.51)
    header(fig, "Which criteria fail, and what changes in the spectrum?",
           "A–C count cases missing a criterion in at least one window. Criteria can overlap.   ·   D–F use floor 0.05.")
    max_fraction = 0
    for i, (axis, gate) in enumerate(zip(axes[0], GATES)):
        panel(axis, f"{chr(65 + i)}   Amplitude floor {gate:g}", "Cases missing criterion" if i == 0 else None)
        data = sweep_rows(summary, gate, "all")
        width = column(data, "width")
        for field, label, color, marker in criteria:
            fraction = column(data, field) / column(data, "cases")
            max_fraction = max(max_fraction, fraction.max())
            axis.plot(width, fraction, color=color, marker=marker, label=label, lw=1.3, markersize=3.6)
        axis.yaxis.set_major_formatter(PercentFormatter(1))
    for axis in axes[0]:
        axis.set_ylim(-0.02 * max(max_fraction, 0.10), min(1.06, max(0.10, max_fraction * 1.09)))
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="upper left", bbox_to_anchor=(0.077, 0.878),
               frameon=False, ncols=4, fontsize=9, labelcolor=MUTED)
    plot_spectral_trends(axes[1], cases, column(sweep_rows(summary, 0.05, "all"), "width"))
    for axis in axes[1]:
        axis.set_xlabel("Independent error width (rad)", fontsize=9)
        axis.xaxis.set_major_locator(MaxNLocator(7))
    fig.text(0.085, 0.043, "D–E compare separate worst values at each width. The limiting case or site can differ between curves.",
             fontsize=8.8, color=MUTED)
    save(fig, "failure-criteria.png")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", nargs="?", choices=("all", "sweep", "diagnostics"), default="all")
    args = parser.parse_args()
    if args.stage in {"all", "sweep"}:
        summary = read_csv("sweep_summary.csv")
        cases = read_csv("case_summary.csv")
        plot_sweep(summary, cases)
        plot_failure_criteria(summary, cases)
    if args.stage in {"all", "diagnostics"}:
        plot_diagnostics()


if __name__ == "__main__":
    main()
