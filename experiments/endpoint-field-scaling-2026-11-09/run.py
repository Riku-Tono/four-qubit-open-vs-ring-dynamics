"""Sweep the endpoint-field magnitude using the previous N=4 model."""

import json
from pathlib import Path
from runpy import run_path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


OUTPUT = Path(__file__).resolve().parent
PREVIOUS = OUTPUT.parent / "endpoint-field-2026-11-09"
MODEL = run_path(str(PREVIOUS / "run.py"))
LINK, END, DT, EPS = (MODEL[key] for key in ("LINK", "END", "DT", "EPS"))
MAGNITUDES = [EPS / 2**i for i in range(3)]
SEAM = MODEL["PAIRS"].index((0, 3))
evolve, connected = MODEL["evolve"], MODEL["connected"]


def seam_shares(time, eps):
    weak_link = connected(evolve(time, LINK, "opposite", eps))
    reference = connected(evolve(time, 0.0, "opposite", eps))
    pair_power = np.sum((weak_link - reference)**2, axis=-2)
    integrals = np.trapezoid(pair_power, x=time, axis=0)
    return integrals[:, SEAM] / integrals.sum(axis=-1)


def measure(dt):
    time = np.arange(round(END / dt) + 1) * dt
    shares = {sign * mag: seam_shares(time, sign * mag)
              for mag in MAGNITUDES for sign in (1, -1)}
    rows = []
    for mag in MAGNITUDES:
        for eps in (mag, -mag):
            qa, qb = map(float, shares[eps])
            delta = qa - qb
            cross_error = qa - float(shares[-eps][1])
            assert abs(cross_error) < 1e-12
            rows.append({"epsilon": eps, "Q_A": qa, "Q_B": qb,
                         "same_sign_delta": delta, "delta_over_epsilon": delta / eps,
                         "cross_sign_error": cross_error})
    return rows


def plot(rows):
    fig, (slope, cross) = plt.subplots(2, 1, figsize=(7.6, 6.3), layout="constrained")
    fig.suptitle(r"Endpoint-field scaling   ($N=4$, $\lambda=0.5$ vs $0$)")
    for sign, color, marker in ((1, "#B44915", "o"), (-1, "#225D97", "x")):
        subset = [row for row in rows if sign * row["epsilon"] > 0]
        x = [abs(row["epsilon"]) for row in subset]
        label = r"$\epsilon>0$" if sign > 0 else r"$\epsilon<0$"
        slope.plot(x, [row["delta_over_epsilon"] for row in subset], marker=marker,
                   linestyle="-" if sign > 0 else "none", color=color, label=label)
        cross.plot(x, [abs(row["cross_sign_error"]) for row in subset], marker=marker,
                   linestyle="none", color=color, label=label)
    slope.set(title="Same-sign initial-state comparison",
              ylabel=r"$[Q_A(\epsilon)-Q_B(\epsilon)]/\epsilon$")
    cross.set(title="Sign reversal symmetry error",
              ylabel=r"$|Q_A(\epsilon)-Q_B(-\epsilon)|$", ylim=(0, None))
    cross.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
    for axis in (slope, cross):
        axis.set_xscale("log", base=2)
        axis.invert_xaxis()
        axis.set_xticks(MAGNITUDES, labels=[f"{x:g}" for x in MAGNITUDES])
        axis.set_xlabel(r"$|\epsilon|$  (halved at each step)")
        axis.grid(alpha=0.22)
        axis.legend(frameon=False)
    fig.savefig(OUTPUT / "q03-scaling.png", dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    fine, coarse = measure(DT / 2), measure(DT)
    baseline = json.loads((PREVIOUS / "results.json").read_text(encoding="utf-8"))
    anchor_error = max(
        abs(row[key] - previous)
        for row in fine[:2]
        for key, previous in zip(("Q_A", "Q_B"),
            baseline["runs"][f"opposite {row['epsilon']:+.2f}"]["q03_dt_half_A_B"])
    )
    assert anchor_error < 1e-12
    grid_error = {key: max(abs(a[key] - b[key]) for a, b in zip(fine, coarse))
                  for key in ("Q_A", "Q_B", "same_sign_delta", "delta_over_epsilon")}
    results = {"lambda": LINK, "reference_lambda": 0.0, "end_time": END,
               "dt": DT / 2, "audit_dt": DT, "fine": fine, "coarse": coarse,
               "max_grid_difference": grid_error, "previous_result_error_at_abs_epsilon_0.05": anchor_error}
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    plot(fine)
    print("Wrote Q03 scaling results and plot to", OUTPUT)
