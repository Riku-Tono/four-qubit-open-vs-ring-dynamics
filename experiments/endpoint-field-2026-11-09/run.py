"""Reproduce the N=4 endpoint-field check with NumPy and Matplotlib."""

import json
from functools import reduce
from itertools import combinations
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


OMEGA, J, LINK = 1.0, 0.7, 0.5
END, DT, EPS = 3.688, 0.002, 0.05
CASES = [(mode, eps) for mode in ("same", "opposite") for eps in (-EPS, EPS)]
PAIRS = list(combinations(range(4), 2))
SITES = range(4)
I = np.eye(2)
PAULI = (
    np.array([[0, 1], [1, 0]]),
    np.array([[0, -1j], [1j, 0]]),
    np.diag([1, -1]),
)
OPS = [[reduce(np.kron, [p if i == j else I for i in SITES]) for j in SITES]
       for p in PAULI]
PAIR_OPS = [[ops[i] @ ops[j] for i, j in PAIRS] for ops in OPS]
OCCUPATION = (np.arange(16)[:, None] >> np.arange(3, -1, -1)) & 1
OUTPUT = Path(__file__).resolve().parent


def evolve(time, link, mode, eps):
    x, _, z = OPS
    h = OMEGA / 2 * sum(x) + J / 2 * (
        z[0] @ z[1] + z[1] @ z[2] + z[2] @ z[3] + link * z[3] @ z[0]
    )
    h += eps * (x[0] + (x[3] if mode == "same" else -x[3]))
    energy, eigenvectors = np.linalg.eigh(h)
    phases = np.exp(-1j * time[:, None, None] * energy)
    # Site 0 is the leftmost bit. The two rows start at |1000> and |0100>.
    return (phases * eigenvectors[[0b1000, 0b0100]]) @ eigenvectors.T


def expect(states, operator):
    return np.sum(states.conj() * (states @ operator.T), axis=-1).real


def connected(states):
    channels = []
    for ops, pair_ops in zip(OPS, PAIR_OPS):
        singles = [expect(states, op) for op in ops]
        channels.append(np.stack([
            expect(states, op) - singles[i] * singles[j]
            for (i, j), op in zip(PAIRS, pair_ops)
        ], axis=-1))
    return np.stack(channels, axis=-2)


def max_abs(values):
    return float(np.max(np.abs(values)))


def run_grid(dt):
    time = np.arange(round(END / dt) + 1) * dt
    runs = {}
    for mode, eps in CASES:
        states = np.array([evolve(time, link, mode, eps) for link in (0.0, LINK)])
        occupation = np.abs(states)**2 @ OCCUPATION
        residual = occupation[:, :, 0, 0] + occupation[:, :, 1, 3] - 1
        delta = connected(states[1]) - connected(states[0])
        pair_power = np.sum(delta**2, axis=-2)
        # The §14 threshold is for w03(t). Integrate the powers over the full grid.
        pair_integral = np.trapezoid(pair_power, x=time, axis=0)
        q = pair_integral[:, PAIRS.index((0, 3))] / np.sum(pair_integral, axis=-1)
        runs[mode, eps] = {
            "states": states, "residual": residual, "q": q,
            "norm_error": max_abs(np.linalg.norm(states, axis=-1) - 1),
        }
    return time, runs


def record_audit(time, fine, coarse):
    short_times = 0.16 / 2.0**np.arange(6)
    indices = np.rint(short_times / (DT / 2)).astype(int)
    assert np.array_equal(short_times, time[indices])
    summary = {"Omega": OMEGA, "J": J, "lambda": LINK, "reference_lambda": 0.0,
               "epsilon_magnitude": EPS, "end_time": END, "dt": DT, "audit_dt": DT / 2,
               "runs": {}, "opposite_comparisons": []}
    for mode, eps in CASES:
        run, base = fine[mode, eps], coarse[mode, eps]
        row = {
            "max_abs_residual_lambda_0": max_abs(run["residual"][0]),
            "max_abs_residual_lambda_0.5": max_abs(run["residual"][1]),
            "max_state_norm_error_dt": base["norm_error"],
            "max_state_norm_error_dt_half": run["norm_error"],
            "max_state_amplitude_change_on_shared_times": max_abs(run["states"][:, ::2] - base["states"]),
            "max_residual_change_on_shared_times": max_abs(run["residual"][:, ::2] - base["residual"]),
            "q03_dt_A_B": base["q"].tolist(),
            "q03_dt_half_A_B": run["q"].tolist(),
            "max_q03_change_when_halving_dt": max_abs(run["q"] - base["q"]),
        }
        assert max(base["norm_error"], run["norm_error"]) < 1e-12
        assert row["max_state_amplitude_change_on_shared_times"] < 1e-12
        assert row["max_residual_change_on_shared_times"] < 1e-12
        assert row["max_q03_change_when_halving_dt"] < 1e-6
        if mode == "same":
            assert max_abs(run["residual"]) < 1e-12
            assert max_abs(base["residual"]) < 1e-12
            assert abs(run["q"][0] - run["q"][1]) < 1e-12
        else:
            ratios = run["residual"][1, indices] / short_times**2
            theory = -2 * OMEGA * eps
            errors = np.abs(ratios - theory)
            assert np.all(np.diff(errors) < 0) and errors[-1] < 1e-6
            row["coefficient_convergence"] = [
                {"t": t, "R_C_over_t_squared": ratio, "theory": theory, "absolute_error": error}
                for t, ratio, error in zip(short_times, ratios, errors)
            ]
            cross = float(run["q"][0] - fine[mode, -eps]["q"][1])
            assert abs(cross) < 1e-12
            assert abs(base["q"][0] - coarse[mode, -eps]["q"][1]) < 1e-12
            summary["opposite_comparisons"].append({
                "epsilon": eps,
                "cross_sign_error_A_eps_minus_B_neg_eps": cross,
                "same_sign_difference_A_eps_minus_B_eps": float(run["q"][0] - run["q"][1]),
            })
        summary["runs"][f"{mode} {eps:+.2f}"] = row
    return summary


def plot_occupation(time, fine):
    fig, (full, zoom) = plt.subplots(2, 1, figsize=(8.0, 6.3), layout="constrained")
    full.plot(time, fine["opposite", EPS]["residual"][1], color="#B44915", label=r"Opposite, $\epsilon=+0.05$")
    full.plot(time, fine["opposite", -EPS]["residual"][1], color="#225D97", label=r"Opposite, $\epsilon=-0.05$")
    full.plot(time, fine["same", EPS]["residual"][1], color="#303B46", linestyle=":", label=r"Same direction, $\epsilon=+0.05$")
    full.set(title=r"Endpoint occupation residual   ($\lambda=0.5$)", xlabel="Time t",
             ylabel=r"$R_C(t)$", xlim=(0, END), ylim=(-0.17, 0.17))
    full.legend(loc="upper left", frameon=False, ncols=2)

    early = time <= 0.4
    zoom.plot(time[early], fine["opposite", EPS]["residual"][1, early], color="#B44915",
              linewidth=2.5, label=r"Opposite, $\epsilon=+0.05$")
    zoom.plot(time[early], -2 * OMEGA * EPS * time[early]**2, color="#303B46", linestyle="--",
              linewidth=1.2, label=r"Short-time prediction $-0.1t^2$")
    zoom.set(title="Short-time detail", xlabel="Time t", ylabel=r"$R_C(t)$", xlim=(0, 0.4))
    zoom.legend(frameon=False, loc="lower left")
    for axis in (full, zoom):
        axis.grid(alpha=0.23)
    fig.savefig(OUTPUT / "occupation.png", dpi=175)
    plt.close(fig)


def plot_q03(comparisons):
    fig, (cross, same) = plt.subplots(1, 2, figsize=(8.8, 3.7), layout="constrained")
    fig.suptitle(r"Integrated seam share $Q_{03}$   ($0\leq t\leq 3.688$, opposite endpoint fields)")
    for axis, key, title, ylabel in (
        (cross, "cross_sign_error_A_eps_minus_B_neg_eps", "Sign reversal symmetry", r"$|Q_A(\epsilon)-Q_B(-\epsilon)|$"),
        (same, "same_sign_difference_A_eps_minus_B_eps", "Same-sign initial-state comparison", r"$Q_A(\epsilon)-Q_B(\epsilon)$"),
    ):
        values = [row[key] for row in comparisons]
        if axis is cross:
            values = [abs(value) for value in values]
        labels = [f"{row['epsilon']:+.2f}" for row in comparisons]
        bars = axis.bar(labels, values, width=0.48, color=["#225D97", "#B44915"])
        axis.axhline(0, color="#303B46", linewidth=0.8)
        axis.set(title=title, xlabel=r"$\epsilon$", ylabel=ylabel)
        axis.grid(axis="y", alpha=0.23)
        axis.set_axisbelow(True)
        if axis is cross:
            axis.set_ylim(0, 1.25 * max(max(values), np.finfo(float).eps))
            axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
            axis.bar_label(bars, fmt="%.2e", padding=3)
        else:
            limit = max(max_abs(values) * 1.30, 1e-3)
            axis.set_ylim(-limit, limit)
            axis.bar_label(bars, fmt="%+.5f", padding=3)
    fig.savefig(OUTPUT / "q03-comparisons.png", dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    coarse_time, coarse = run_grid(DT)
    fine_time, fine = run_grid(DT / 2)
    assert np.array_equal(coarse_time, fine_time[::2])
    audit = record_audit(fine_time, fine, coarse)
    (OUTPUT / "results.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    plot_occupation(fine_time, fine)
    plot_q03(audit["opposite_comparisons"])
    print("Numerical audit passed. Results and plots written to", OUTPUT)
