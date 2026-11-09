"""Compare the N=4 and N=6 endpoint-field seam responses (NumPy >= 2)."""

import json
from functools import reduce
from itertools import combinations
from pathlib import Path

import numpy as np


OMEGA, J, LINK, REFERENCE = 1.0, 0.7, 0.5, 0.0
END, DT = 3.688, 0.001
MAGNITUDES = 0.05 / 2.0**np.arange(5)
FIT_WINDOWS = [(2, 1), (3, 1), (3, 2), (5, 2)]
OUTPUT = Path(__file__).resolve().parent
ARCHIVE = OUTPUT.parent / "endpoint-field-zero-slope-2026-11-09" / "results.json"
I = np.eye(2)
PAULI = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1]))


def max_abs(array):
    return float(np.max(np.abs(array)))


class Chain:
    def __init__(self, size):
        self.size = size
        dim = 2**size
        basis = np.arange(dim)
        self.pairs = list(combinations(range(size), 2))
        self.seam = self.pairs.index((0, size - 1))
        self.left, self.right = np.array(self.pairs).T
        singles = [[reduce(np.kron, [pauli if i == j else I for j in range(size)])
                    for i in range(size)] for pauli in PAULI]
        self.observables = []
        for ops in singles:
            operators = np.array(ops + [ops[i] @ ops[j] for i, j in self.pairs])
            # Each Pauli string acts by one nonzero entry per row.
            columns = np.argmax(np.abs(operators), axis=-1)
            self.observables.append((columns, operators[np.arange(len(operators))[:, None], basis, columns]))

        x, _, z = singles
        self.open_chain = OMEGA / 2 * sum(x) + J / 2 * sum(z[i] @ z[i + 1] for i in range(size - 1))
        self.seam_bond = J / 2 * (z[-1] @ z[0])
        self.perturbation = x[0] - x[-1]
        reflection = np.eye(dim)[:, [int(f"{b:0{size}b}"[::-1], 2) for b in basis]]
        gamma = reduce(np.matmul, z) @ reduce(np.matmul, x[::2])
        self.theta_linear = reflection @ gamma
        initial_a = dim // 2
        mapped = self.theta_linear[:, initial_a]
        initial_b = np.argmax(np.abs(mapped))
        self.initial = [initial_a, int(initial_b)]
        self.initial_labels = [f"{index:0{size}b}" for index in self.initial]
        self.phase = mapped[initial_b]

    def hamiltonian(self, link, epsilon):
        return self.open_chain + link * self.seam_bond + epsilon * self.perturbation

    def evolve(self, time, link, epsilon, differentiate=False):
        energy, vectors = np.linalg.eigh(self.hamiltonian(link, epsilon))
        initial = vectors[self.initial]
        phases = np.exp(-1j * time[:, None] * energy)
        states = (phases[:, None, :] * initial) @ vectors.T
        if not differentiate:
            return states

        # Fréchet derivative of exp(-itH), including the degenerate-energy limit.
        v = vectors.T @ self.perturbation @ vectors
        midpoint = (energy[:, None] + energy) / 2
        frequency = (energy[:, None] - energy) / (2 * np.pi)
        tangents = np.empty_like(states)
        for start in range(0, len(time), 128):
            chunk = slice(start, start + 128)
            t = time[chunk, None, None]
            kernel = -1j * t * np.exp(-1j * t * midpoint) * np.sinc(t * frequency)
            tangents[chunk] = ((kernel * v) @ initial.T).transpose(0, 2, 1) @ vectors.T
        return states, tangents

    def connected(self, states, tangents=None):
        values, slopes = [], []
        conjugates = states.conj()
        tangent_conjugates = None if tangents is None else tangents.conj()
        for columns, coefficients in self.observables:
            means, derivatives = [], []
            for column, coefficient in zip(columns, coefficients):
                applied = states[..., column] * coefficient
                means.append(np.einsum("...i,...i->...", conjugates, applied).real)
                if tangents is not None:
                    derivatives.append(2 * np.einsum("...i,...i->...", tangent_conjugates, applied).real)
            singles, pairs = np.split(np.stack(means, axis=-1), [self.size], axis=-1)
            left, right = singles[..., self.left], singles[..., self.right]
            values.append(pairs - left * right)
            if tangents is not None:
                singles_d, pairs_d = np.split(np.stack(derivatives, axis=-1), [self.size], axis=-1)
                slopes.append(pairs_d - singles_d[..., self.left] * right - left * singles_d[..., self.right])
        values = np.stack(values, axis=-2)
        return values if tangents is None else (values, np.stack(slopes, axis=-2))

    def share(self, time, powers, step):
        integrals = np.trapezoid(powers[::step], x=time[::step], axis=0)
        total, seam = integrals.sum(axis=-1), integrals[:, self.seam]
        return seam / total, total

    def sample(self, time, epsilon):
        linked, opened = (self.evolve(time, link, epsilon) for link in (LINK, REFERENCE))
        delta = self.connected(linked) - self.connected(opened)
        powers = np.sum(delta**2, axis=-2)
        states = np.stack([linked, opened])
        return states, [self.share(time, powers, step)[0] for step in (2, 1)]

    def direct_response(self, time, step, power, power_d):
        shares, total = self.share(time, power, step)
        integral_d = np.trapezoid(power_d[::step], x=time[::step], axis=0)
        numerator_d, denominator_d = integral_d[:, self.seam], integral_d.sum(axis=-1)
        shares_d = (numerator_d - shares * denominator_d) / total
        return {"Q_A_B": shares.tolist(), "numerator_A_B": (shares * total).tolist(),
                "denominator_A_B": total.tolist(), "numerator_derivative_A_B": numerator_d.tolist(),
                "denominator_derivative_A_B": denominator_d.tolist(),
                "dQ_A_B": shares_d.tolist(), "slope": float(shares_d[0] - shares_d[1])}

    def measure(self, time):
        primary, refined = [], []
        norm_error, state_error, operator_error = 0.0, 0.0, 0.0
        for magnitude in MAGNITUDES:
            positive, qp = self.sample(time, magnitude)
            negative, qm = self.sample(time, -magnitude)
            for states, flipped in ((positive, negative), (negative, positive)):
                norm_error = max(norm_error, max_abs(np.linalg.norm(states, axis=-1) - 1))
                mapped = states[:, :, 0].conj() @ self.theta_linear.T
                state_error = max(state_error, max_abs(mapped - self.phase * flipped[:, :, 1]))
            for link in (LINK, REFERENCE):
                u = self.theta_linear
                residual = u @ self.hamiltonian(link, magnitude).conj() @ u.T + self.hamiltonian(link, -magnitude)
                operator_error = max(operator_error, max_abs(residual))
            for rows, plus, minus in zip((primary, refined), qp, qm):
                for epsilon, shares, flipped in ((magnitude, plus, minus), (-magnitude, minus, plus)):
                    qa, qb = map(float, shares)
                    rows.append({"epsilon": float(epsilon), "Q_A": qa, "Q_B": qb,
                                 "same_sign_delta": qa - qb, "delta_over_epsilon": (qa - qb) / epsilon,
                                 "cross_sign_error": qa - float(flipped[1])})

        linked, linked_d = self.connected(*self.evolve(time, LINK, 0, differentiate=True))
        opened, opened_d = self.connected(*self.evolve(time, REFERENCE, 0, differentiate=True))
        delta, delta_d = linked - opened, linked_d - opened_d
        power = np.sum(delta**2, axis=-2)
        power_d = 2 * np.sum(delta * delta_d, axis=-2)
        direct = [self.direct_response(time, step, power, power_d) for step in (2, 1)]
        assert max(norm_error, state_error, operator_error) < 1e-12
        return {"initial_A_B": self.initial_labels, "theta_A_phase": float(self.phase),
                "max_norm_error": norm_error, "max_cross_sign_state_amplitude_error": state_error,
                "max_operator_identity_error": operator_error,
                "primary": analyze(primary, direct[0]), "refined": analyze(refined, direct[1])}


def analyze(rows, direct):
    x = MAGNITUDES**2
    y = np.array([row["delta_over_epsilon"] for row in rows]).reshape(-1, 2).mean(axis=1)
    fits = []
    for points, degree in FIT_WINDOWS:
        coefficients = np.polynomial.polynomial.polyfit(x[-points:], y[-points:], degree)
        fits.append({"points_nearest_zero": points, "degree_in_epsilon_squared": degree,
                     "coefficients": coefficients.tolist(), "intercept": float(coefficients[0])})
    cross_error = max(abs(row["cross_sign_error"]) for row in rows)
    assert cross_error < 1e-12
    assert min(direct["denominator_A_B"]) > 0
    assert abs(fits[2]["intercept"] - direct["slope"]) < 1e-7
    return {"signed_samples": rows, "mean_delta_over_epsilon": y.tolist(), "fits": fits,
            "direct_derivative_at_zero": direct, "max_abs_cross_sign_error": cross_error}


def main():
    fine_dt = DT / 2
    time = np.arange(round(END / fine_dt) + 1) * fine_dt
    n4, n6 = (Chain(size).measure(time) for size in (4, 6))
    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    anchor_error = max(abs(new[key] - old[key]) for grid in ("primary", "refined")
                       for new, old in zip(n4[grid]["signed_samples"], archive[grid], strict=True)
                       for key in ("Q_A", "Q_B"))
    assert n4["initial_A_B"] == ["1000", "0100"] and anchor_error < 1e-12
    assert n6["initial_A_B"] == ["100000", "010100"]
    ratios = {grid: n6[grid]["direct_derivative_at_zero"]["slope"] /
              n4[grid]["direct_derivative_at_zero"]["slope"] for grid in ("primary", "refined")}
    results = {"Omega": OMEGA, "J": J, "lambda": LINK, "reference_lambda": REFERENCE,
               "end_time": END, "dt": DT, "refined_dt": fine_dt, "magnitudes": MAGNITUDES.tolist(),
               "N4": n4, "N6": n6, "N4_archive_max_abs_share_error": anchor_error,
               "slope_ratio_N6_over_N4": ratios}
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print("Primary N4 and N6 slopes:", *(m["primary"]["direct_derivative_at_zero"]["slope"] for m in (n4, n6)))
    print("Primary N6/N4 slope ratio:", ratios["primary"])


if __name__ == "__main__":
    main()
