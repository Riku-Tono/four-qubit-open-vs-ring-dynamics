"""Render the open-chain spectrum, phase ablations and initial-angle error figures."""

import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator


OUTPUT = Path(__file__).resolve().parent
INK, MUTED, GRID = "#172638", "#506070", "#dce4e8"
END, INNER, MISMATCH = "#14609e", "#b7512a", "#17715e"
CONTRAST, SOFT = "#805283", "#a1acb8"
METRICS = (
    ("max_abs_lag", "Largest absolute mean phase gap", "Degrees"),
    ("min_plv", "Lowest phase-locking value", "PLV"),
    ("min_envelope", "Lowest instantaneous envelope", r"$\min_i |u_i|$"),
)
ANGLE_STYLES = (
    ("end_theta", 0, "theta", END, "-", "o", r"End $\theta$"),
    ("end_phi", 0, "phi", END, "--", "s", r"End $\phi$"),
    ("inner_theta", 1, "theta", INNER, "-", "o", r"Interior $\theta$"),
    ("inner_phi", 1, "phi", INNER, "--", "s", r"Interior $\phi$"),
)


def read_csv(name):
    with (OUTPUT / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def as_float(row, field):
    value = row[field]
    return float(value) if value else np.nan


def style_axis(axis, direction="y"):
    axis.set_axisbelow(True)
    axis.grid(axis=direction, color=GRID, linewidth=0.75)
    axis.tick_params(length=0, pad=6, labelsize=9, colors=MUTED)
    axis.xaxis.label.set_color(MUTED)
    axis.yaxis.label.set_color(MUTED)
    axis.spines[["top", "right"]].set_visible(False)
    for spine in ("bottom", "left"):
        axis.spines[spine].set_color(GRID)


def panel_title(axis, title):
    axis.set_title(title, loc="left", color=INK, weight="bold", fontsize=12, pad=12)


def save(fig, filename):
    fig.savefig(OUTPUT / filename, dpi=175, facecolor="white")
    plt.close(fig)


def plot_spectral_lines(axis, rows):
    style_axis(axis)
    omega = np.array([as_float(row, "omega") for row in rows])
    for site, color, marker in ((0, END, "o"), (1, INNER, "s")):
        amplitude = np.array([as_float(row, f"abs{site}") for row in rows])
        visible = amplitude > 1e-12
        x, y = omega[visible], amplitude[visible]
        axis.vlines(x, 0, y, color=color, lw=1.1, alpha=0.7)
        axis.scatter(x, y, color=color, marker=marker, s=14, zorder=3)
    axis.set_xlim(0, max(omega) + 0.16)
    axis.set_ylim(bottom=0)
    axis.set_xlabel(r"Angular frequency $\omega$")
    axis.set_ylabel(r"Occupation coefficient $|c_i(\omega)|$")


def plot_mismatch(axis, selected):
    style_axis(axis, "x")
    panel_title(axis, "C   Residuals that distinguish end and interior")
    omega = [f"{as_float(row, 'omega'):.5f}" for row in selected]
    values = [as_float(row, "normalized_diff_abs") for row in selected]
    bars = axis.barh(omega, values, color=MISMATCH, alpha=0.9, height=0.62)
    axis.invert_yaxis()
    axis.set_xlim(0, max(values) * 1.25)
    axis.bar_label(bars, labels=[f"{v:.4f}" for v in values], padding=4, fontsize=8, color=MUTED)
    axis.set_xlabel(r"Difference in normalized coefficients $|c_0/c_{0*}-c_1/c_{1*}|$")
    axis.set_ylabel(r"Angular frequency $\omega$")


def plot_ablation(axis, rows, leading):
    style_axis(axis, "x")
    panel_title(axis, "D   Effect of retaining or removing a line")
    records = {(row["case"], row["mode"]): row for row in rows}
    labels, values, colors = [], [], []
    for case, label, color in (("full", "Full spectrum", INK),
                               ("dominant_only", "Dominant only", SOFT)):
        row = next(item for item in rows if item["case"] == case)
        labels.append(label)
        values.append(as_float(row, "lag_rms_deg"))
        colors.append(color)
    for case, prefix, color in (("add_one", "Dominant + ", MISMATCH), ("omit_one", "Full − ", CONTRAST)):
        for mode in leading[:3]:
            labels.append(f"{prefix}{as_float(mode, 'omega'):.5f}")
            values.append(as_float(records[case, mode["mode"]], "lag_rms_deg"))
            colors.append(color)
    bars = axis.barh(labels, values, color=colors, alpha=0.9, height=0.62)
    axis.invert_yaxis()
    axis.set_xlim(0, max(values) * 1.25)
    axis.bar_label(bars, labels=[f"{v:.3f}°" for v in values], padding=4, fontsize=8, color=MUTED)
    axis.set_xlabel("RMS fluctuation of the end–interior phase gap (degrees)")


def plot_spectrum():
    rows = read_csv("spectrum.csv")
    dominant = next(row for row in rows if row["is_dominant"].lower() in {"true", "1"})
    residual = [row for row in rows if row is not dominant]
    leading = sorted(residual, key=lambda row: as_float(row, "normalized_diff_abs"), reverse=True)
    frequency = as_float(dominant, "omega")
    phase_gap = (as_float(dominant, "phase0_deg") - as_float(dominant, "phase1_deg") + 180) % 360 - 180

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.0))
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.09, top=0.78, wspace=0.52, hspace=0.66)
    fig.text(0.12, 0.955, "A common line carries the four-site oscillation", fontsize=24, weight="bold", color=INK)
    fig.text(0.12, 0.913,
             r"Frozen mirror product  ·  $R_i(t)=R_i^{\rm DC}+\mathrm{Re}\sum_{\omega>0} c_i(\omega)e^{i\omega t}$",
             fontsize=11, color=MUTED)
    fig.legend(handles=[Line2D([], [], color=END, marker="o", label="End sites 0, 3"),
                        Line2D([], [], color=INNER, marker="s", label="Interior sites 1, 2")],
               loc="upper left", bbox_to_anchor=(0.115, 0.899), frameon=False, ncols=2,
               fontsize=10, labelcolor=MUTED)

    plot_spectral_lines(axes[0, 0], rows)
    panel_title(axes[0, 0], "A   All contributing Bohr frequencies")
    axis = axes[0, 0]
    axis.axvline(frequency, lw=0.9, color=SOFT, ls="--", zorder=0)
    axis.text(0.29, 0.9, rf"$\omega_*={frequency:.5f}$" + "\n"
              + f"End–interior phase gap {phase_gap:+.4f}°", transform=axis.transAxes,
              ha="left", va="top", color=MUTED, fontsize=10)
    axis.yaxis.set_major_locator(MaxNLocator(5))

    plot_spectral_lines(axes[0, 1], residual)
    panel_title(axes[0, 1], "B   Residual lines at enlarged scale")
    axes[0, 1].yaxis.set_major_locator(MaxNLocator(5))

    plot_mismatch(axes[1, 0], leading[:7])
    plot_ablation(axes[1, 1], read_csv("mode_ablation.csv"), leading)
    save(fig, "spectrum.png")


def perturbation_series(rows, baseline, mirror_axis, site, local_axis, single_site):
    if single_site:
        keep = [row for row in rows if row["scope"] == "one_site"
                and row["axis"] == local_axis and int(row["site"]) == site]
    else:
        keep = [row for row in rows if row["scope"] == "one_angle" and row["axis"] == mirror_axis]
    return sorted(keep + [baseline], key=lambda row: as_float(row, "error") if row is not baseline else 0)


def plot_single_errors(axis, rows, baseline, metric, single_site):
    for mirror_axis, site, local_axis, color, line, marker, _ in ANGLE_STYLES:
        series = perturbation_series(rows, baseline, mirror_axis, site, local_axis, single_site)
        x = [as_float(row, "error") if row is not baseline else 0 for row in series]
        y = [as_float(row, metric) for row in series]
        axis.plot(x, y, color=color, linestyle=line, marker=marker, lw=1.35, markersize=4.5)
    errors = [-0.05, -0.01, 0, 0.01, 0.05]
    axis.set_xticks(errors, ["−0.05", "−0.01", "0", "+0.01", "+0.05"])
    axis.set_xlim(-0.056, 0.056)
    axis.set_xlabel("Added angle (rad)")
    axis.axvline(0, color=GRID, lw=1, zorder=0)


def plot_combined_errors(axis, rows, metric):
    values = [[as_float(row, metric) for row in rows if row["scope"] == "combined"
               and np.isclose(as_float(row, "max_abs_error"), error)] for error in (0.01, 0.05)]
    labels = [f"{error:.2f}   ({len(group)} cases)" for error, group in zip((0.01, 0.05), values)]
    axis.boxplot(values, widths=0.49, whis=(0, 100), tick_labels=labels, patch_artist=True,
                 boxprops={"facecolor": "#daece6", "edgecolor": MISMATCH},
                 whiskerprops={"color": MISMATCH}, capprops={"color": MISMATCH},
                 medianprops={"color": INK, "linewidth": 1.2})
    axis.set_xlabel("Maximum absolute error (rad)")


def plot_robustness():
    rows = read_csv("robustness_summary.csv")
    main = [row for row in rows if as_float(row, "gate") == 0.05]
    baseline = next(row for row in main if row["family"] == "nominal")
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 10.4))
    fig.subplots_adjust(left=0.11, right=0.975, bottom=0.125, top=0.795, wspace=0.34, hspace=0.82)
    fig.text(0.11, 0.955, "Sensitivity to preparation-angle errors", fontsize=24, weight="bold", color=INK)
    fig.text(0.11, 0.918,
             r"Initial angles perturbed without reoptimization  ·  Worst values over six pairs and the 24 fixed $8\pi$ windows",
             fontsize=11, color=MUTED)
    fig.legend(handles=[Line2D([], [], color=color, linestyle=line, marker=marker, label=label)
                        for _, _, _, color, line, marker, label in ANGLE_STYLES]
               + [Line2D([], [], color=SOFT, linestyle=":", label="Frozen nominal")],
               loc="upper left", bbox_to_anchor=(0.105, 0.900), ncols=5, frameon=False,
               fontsize=10, labelcolor=MUTED)

    for column, (metric, title, unit) in enumerate(METRICS):
        values = [as_float(row, metric) for row in main]
        low, high = np.nanmin(values), np.nanmax(values)
        span = max(high - low, 0.001)
        panel_title(axes[0, column], title)
        for row, axis in enumerate(axes[:, column]):
            style_axis(axis)
            axis.set_ylabel(unit)
            axis.axhline(as_float(baseline, metric), color=SOFT, lw=1, ls=":", zorder=0)
            axis.yaxis.set_major_locator(MaxNLocator(5))
            if metric == "min_envelope":
                axis.axhline(0.10, color=CONTRAST, lw=1, ls="--", zorder=0)
                axis.set_ylim(min(0.09, low - 0.1 * span), high + 0.17 * span)
            elif metric == "min_plv":
                axis.set_ylim(low - 0.12 * span, min(1.0, high + 0.12 * span))
                axis.ticklabel_format(axis="y", useOffset=False, style="plain")
            else:
                axis.set_ylim(0, high + 0.13 * span)
            if row == 2:
                plot_combined_errors(axis, main, metric)
            else:
                plot_single_errors(axis, main, baseline, metric, single_site=(row == 1))

    row_labels = ("Mirror pair\nOne angle", "Single physical site\nOne angle", "Mirror pairs\nMultiple angles")
    for row, label in enumerate(row_labels):
        pos = axes[row, 0].get_position()
        fig.text(0.020, (pos.y0 + pos.y1) / 2, label, color=INK, weight="bold", fontsize=11,
                 ha="center", va="center", rotation=90)

    counts = [min(int(row["pass_windows"]) for row in rows if as_float(row, "gate") == gate)
              for gate in (0.025, 0.05, 0.1)]
    passed = "    ·    ".join(f"floor {gate}: {count}/24" for gate, count in zip(("0.025", "0.05", "0.10"), counts))
    fig.text(0.11, 0.064, f"Worst case by floor across {len(main)} preparations   ·   {passed}",
             color=INK, fontsize=10.2, weight="bold")
    fig.text(0.11, 0.036, "Main plots use floor 0.05. Purple line marks 0.10. Single-site curves use sites 0 and 1. Whiskers span the tested range.",
             color=MUTED, fontsize=9.5)
    save(fig, "angle-robustness.png")


if __name__ == "__main__":
    plot_spectrum()
    plot_robustness()
