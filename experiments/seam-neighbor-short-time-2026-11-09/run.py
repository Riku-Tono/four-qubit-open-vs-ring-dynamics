"""Find exact connected-correlation Taylor terms and check finite-time values."""

import csv
import json
from decimal import Decimal, localcontext
from fractions import Fraction
from math import comb, factorial, log2
from pathlib import Path


OUTPUT = Path(__file__).resolve().parent
ARCHIVE = OUTPUT.parent / "seam-neighbor-waveforms-2026-11-09"
PAIRS = {4: ("02", "13"), 6: ("04", "15")}
CHANNELS = ("XX", "YY", "ZZ")
TIMES = tuple(map(Decimal, ("0.0125", "0.025", "0.05", "0.1", "0.2", "0.4")))
H_SCALE, RADIUS, ORDER = 40, 5, 20


class Ising:
    def __init__(self, size, seam, cut_middle=False, *, radius=RADIUS):
        self.dim = 2**size
        self.radius = radius
        self.masks = [1 << (size - 1 - site) for site in range(size)]
        self.z = [[1 - 2 * bool(basis & mask) for basis in range(self.dim)] for mask in self.masks]
        # M = 40H. In this scale Omega/2 = 20, J/2 = 14 and lambda*J/2 = 7.
        bonds = [(i, i + 1, 14) for i in range(size - 1) if not (cut_middle and i == 1)]
        bonds.append((0, size - 1, seam))
        self.diagonal = [sum(k * self.z[i][b] * self.z[j][b] for i, j, k in bonds)
                         for b in range(self.dim)]
        self.initial = [int(b == self.dim // 2) for b in range(self.dim)]
        assert max(map(abs, self.diagonal)) + 20 * size <= H_SCALE * self.radius

    def apply(self, vector):
        return [diagonal * vector[b] + 20 * sum(vector[b ^ mask] for mask in self.masks)
                for b, diagonal in enumerate(self.diagonal)]

    def pauli(self, axis, sites):
        flip = sum(self.masks[i] for i in sites) if axis != "Z" else 0
        signs = [1] * self.dim if axis == "X" else [self.z[sites[0]][b] for b in range(self.dim)]
        if axis != "X" and len(sites) == 2:
            signs = [sign * z for sign, z in zip(signs, self.z[sites[1]])]
        return list(zip((b ^ flip for b in range(self.dim)), signs)), len(sites) if axis == "Y" else 0


def moment_series(powers, observable):
    mapping, y_count = observable
    coefficients = []
    for n in range(len(powers)):
        phase = (1, 0, -1, 0)[(n + y_count) % 4]
        numerator = 0
        if phase:
            for k in range(n + 1):
                bra, ket = powers[n - k], powers[k]
                term = sum(bra[destination] * sign * source for (destination, sign), source in zip(mapping, ket))
                numerator += (-1)**k * comb(n, k) * term
        coefficients.append(Fraction(phase * numerator, H_SCALE**n * factorial(n)))
    return coefficients


def hamiltonian_powers(model, order):
    powers = [model.initial]
    for _ in range(order):
        powers.append(model.apply(powers[-1]))
    return powers


def connected_series(powers, observables):
    combined, left, right = [moment_series(powers, observable) for observable in observables]
    result = [coefficient - sum(left[k] * right[n - k] for k in range(n + 1))
              for n, coefficient in enumerate(combined)]
    assert all(value == 0 for value in result[1::2])
    return result


def leading(series):
    return next((order, value) for order, value in enumerate(series) if value)


def bessel(order, half_argument):
    term = half_argument**order / factorial(order)
    total = term
    k = 1
    while term:
        term *= -(half_argument**2) / (k * (k + order))
        updated = total + term
        if updated == total:
            break
        total = updated
        k += 1
    return total


def evolve(model, t, degree):
    # Chebyshev expansion of exp(-iHt) on [-radius, radius].
    radius = model.radius
    scale = H_SCALE * radius
    real, imag = [Decimal(0)] * model.dim, [Decimal(0)] * model.dim
    older = [Decimal(value) for value in model.initial]
    current = [value / scale for value in model.apply(older)]
    for n in range(degree + 1):
        vector = older if n == 0 else current
        weight = bessel(n, radius * t / 2) * (1 if n == 0 else 2)
        target = real if n % 2 == 0 else imag
        weight *= (1, -1, -1, 1)[n % 4]
        for i, value in enumerate(vector):
            target[i] += weight * value
        if n:
            older, current = current, [2 * value / scale - previous
                                       for value, previous in zip(model.apply(current), older)]
    return real, imag


def moment_value(state, observable):
    real, imag = state
    mapping, y_count = observable
    if y_count == 1:
        return sum(sign * (imag[d] * real[b] - real[d] * imag[b]) for b, (d, sign) in enumerate(mapping))
    phase = 1 if y_count == 0 else -1
    return phase * sum(sign * (real[d] * real[b] + imag[d] * imag[b]) for b, (d, sign) in enumerate(mapping))


def numerics(model, observables, digits, degree):
    with localcontext() as ctx:
        ctx.prec = digits
        states = [evolve(model, t, degree) for t in TIMES]
        norm_error = max(abs(sum(re**2 + im**2 for re, im in zip(*state)) - 1) for state in states)
        results = {}
        for key, operators in observables.items():
            results[key] = []
            for state in states:
                combined, left, right = [moment_value(state, op) for op in operators]
                results[key].append(combined - left * right)
    return results, norm_error


def make_observables(model, pair, channel):
    i, j = map(int, pair)
    return [model.pauli(channel[0], sites) for sites in ((i, j), (i,), (j,))]


def main():
    previous = json.loads((ARCHIVE / "results.json").read_text(encoding="utf-8"))
    assert [previous[key] for key in ("Omega", "J", "lambda", "reference_lambda", "epsilon")] == [1, 0.7, 0.5, 0, 0]
    with (ARCHIVE / "waveforms.csv").open(encoding="utf-8", newline="") as file:
        archived = {Decimal(row["t"]): row for row in csv.DictReader(file) if Decimal(row["t"]) in TIMES}

    cases, series_by_key = [], {}
    max_numeric = max_archive = max_norm = 0.0
    for size, pairs in PAIRS.items():
        for branch, seam in (("linked", 7), ("open", 0)):
            model = Ising(size, seam)
            powers = hamiltonian_powers(model, ORDER)
            observables = {(pair, channel): make_observables(model, pair, channel)
                           for pair in pairs for channel in CHANNELS}
            values, norm = numerics(model, observables, 80, 64)
            refined, refined_norm = numerics(model, observables, 96, 80)
            max_norm = max(max_norm, float(norm), float(refined_norm))
            for (pair, channel), operators in observables.items():
                series = connected_series(powers, operators)
                series_by_key[size, branch, pair, channel] = series
                order, coefficient = leading(series)
                samples = []
                for t, value, check in zip(TIMES, values[pair, channel], refined[pair, channel]):
                    max_numeric = max(max_numeric, float(abs((value - check) / check)))
                    if t in archived:
                        old = float(archived[t][f"N{size}_{pair}_{channel}_{branch}"])
                        max_archive = max(max_archive, abs(float(value) - old))
                    with localcontext() as ctx:
                        ctx.prec = 80
                        predicted = Decimal(coefficient.numerator) / coefficient.denominator * t**order
                        samples.append({"t": float(t), "C": float(value), "ratio": float(value / predicted)})
                assert all(sample["ratio"] > 0 for sample in samples)
                assert abs(samples[0]["ratio"] - 1) < 0.001
                effective_order = log2(abs(samples[1]["C"] / samples[0]["C"]))
                cases.append({"N": size, "pair": pair, "channel": channel, "branch": branch, "order": order,
                              "coefficient": str(coefficient), "coefficient_float": float(coefficient),
                              "series": list(map(str, series)), "samples": samples,
                              "effective_order_smallest_two_times": effective_order})

    differences = []
    for size, pairs in PAIRS.items():
        for pair in pairs:
            for channel in CHANNELS:
                link = series_by_key[size, "linked", pair, channel]
                opened = series_by_key[size, "open", pair, channel]
                order, value = leading([a - b for a, b in zip(link, opened)])
                differences.append({"N": size, "pair": pair, "channel": channel, "order": order, "coefficient": str(value)})

    # Cutting the central bond removes the alternative non-seam route for both N=4 pairs.
    cut = Ising(4, 7, cut_middle=True)
    powers = hamiltonian_powers(cut, ORDER)
    route_check = {}
    for pair in PAIRS[4]:
        order, value = leading(connected_series(powers, make_observables(cut, pair, "XX")))
        opened = series_by_key[4, "open", pair, "XX"][order]
        full = series_by_key[4, "linked", pair, "XX"][order]
        assert order == 6 and value + opened == full
        route_check[pair] = {"open_route": str(opened), "seam_route": str(value), "both_routes": str(full)}

    assert max_numeric < 1e-25 and max_archive < 1e-12 and max_norm < 1e-70
    audit = {"max_relative_change_80_to_96_digits_and_64_to_80_chebyshev_degree": max_numeric,
             "max_prior_waveform_absolute_error": max_archive, "max_norm_squared_error": max_norm,
             "max_ratio_error_by_time": [max(abs(row["samples"][i]["ratio"] - 1) for row in cases)
                                         for i in range(len(TIMES))],
             "max_effective_order_error_smallest_two_times": max(abs(row["effective_order_smallest_two_times"] - row["order"])
                                                                 for row in cases)}
    results = {"Omega": 1, "J": 0.7, "epsilon": 0, "lambda": 0.5, "reference_lambda": 0,
               "initial": {"N4": "1000", "N6": "100000"}, "max_exact_order": ORDER, "cases": cases,
               "leading_branch_differences": differences, "N4_XX_route_check": route_check, "audit": audit}
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
