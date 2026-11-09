"""Check weak-ring Taylor powers, tied routes, and two tuned cancellations exactly."""

import json
import runpy
from fractions import Fraction
from itertools import combinations
from math import factorial, isqrt, sqrt
from pathlib import Path

from coefficients import binomial_row, path_scores


EXPERIMENTS = Path(__file__).resolve().parents[1]
PRIOR = EXPERIMENTS / "seam-neighbor-short-time-2026-11-09"
EIGHT_QUBIT = EXPERIMENTS / "seam-neighbor-eight-qubit-2026-11-09"
CHANNELS = "XYZ"
SEAM = 7  # M = 40 H has ordinary ZZ weight 14 and seam weight 7.
INTERACTION = Fraction(7, 10)


def exact_model(methods, size, order, seam=SEAM, cuts=()):
    model = methods["Ising"](size, seam, radius=size)
    for edge in cuts:
        left, right = model.z[edge], model.z[(edge + 1) % size]
        weight = seam if edge == size - 1 else 14
        model.diagonal = [value - weight * a * b for value, a, b in zip(model.diagonal, left, right)]
    return model, methods["hamiltonian_powers"](model, order)


def series(methods, data, pair, axis, order):
    model, powers = data
    i, j = pair
    operators = [model.pauli(axis, sites) for sites in (pair, (i,), (j,))]
    return methods["connected_series"](powers[:order + 1], operators)


def powers_for_distance(distance):
    if distance == 1:
        return 4, 4, 6
    return 4 * distance - 2, 4 * distance, 4 * distance + 2


def predicted_coefficients(size, i, j, distance, rows):
    m = 2 * distance
    xx_row, yz_row = rows[distance]
    if 2 * distance == size:
        first = path_scores(size, i, distance, xx_row, yz_row)
        second = path_scores(size, j, distance, xx_row, yz_row)
        scores = [a + b for a, b in zip(first, second)]
    else:
        start = i if j - i == distance else j
        scores = path_scores(size, start, distance, xx_row, yz_row)
    factors = (
        INTERACTION**(m - 2) / (4 * m**2 * factorial(m - 1)**2),
        INTERACTION**m / (4 * factorial(m + 2)**2),
        INTERACTION**m / (4 * factorial(m + 2)**2),
    )
    return dict(zip(CHANNELS, (score * factor for score, factor in zip(scores, factors))))


def adjacent_coefficients(size, i, j):
    start = i if j - i == 1 else j
    weights = [-INTERACTION] + [INTERACTION] * (size - 2) + [-INTERACTION / 2]
    central = weights[start]
    exterior = weights[start - 1] + weights[(start + 1) % size]
    coupling = abs(central)
    yy = coupling * (central / 4 + 5 * exterior / 24)
    zz = coupling * (central / 36 + exterior / 24)
    return {"X": -central * yy / coupling, "Y": yy, "Z": zz}


def verify_powers(methods):
    archived = json.loads((EIGHT_QUBIT / "results.json").read_text(encoding="utf-8"))
    previous = {tuple(map(int, row["pair"])): row["leading"]
                for row in archived["pair_cases"] if row["branch"] == "linked"}
    rows = {d: (binomial_row(2*d), binomial_row(2*d+2)) for d in range(2, 7)}
    count = tied_count = formula_count = adjacent_count = 0
    for size in range(3, 13):
        model = exact_model(methods, size, 4 * (size // 2) + 2)
        for i, j in combinations(range(size), 2):
            distance = min(j - i, size - (j - i))
            if distance == 1:
                predicted = adjacent_coefficients(size, i, j)
            elif size != 4:
                predicted = predicted_coefficients(size, i, j, distance, rows)
            else:
                predicted = {}
            for axis, order in zip(CHANNELS, powers_for_distance(distance)):
                coefficients = series(methods, model, (i, j), axis, order)
                assert not any(coefficients[:-1]) and coefficients[-1], (size, i, j, axis, coefficients)
                if axis in predicted:
                    assert coefficients[-1] == predicted[axis], (size, i, j, axis, coefficients[-1], predicted[axis])
                    if distance == 1:
                        adjacent_count += 1
                    else:
                        formula_count += 1
                if size == 8:
                    old = previous[i, j][axis * 2]
                    assert (old["order"], Fraction(old["coefficient"])) == (order, coefficients[-1])
                count += 1
                tied_count += int(2 * distance == size)
        print(f"N={size}: all pair powers match", flush=True)
    print(f"{count} exact channel cases passed at lambda=1/2, including {tied_count} tied-route cases.")
    print(f"{formula_count} leading coefficients match the path-score formulas.")
    print(f"{adjacent_count} leading coefficients match the adjacent-pair formulas.")


def verify_tied_routes(methods):
    for size in (4, 6, 8, 10):
        distance = size // 2
        max_order = 4 * distance + 2
        ring = exact_model(methods, size, max_order)
        for i in range(distance):
            j = i + distance
            plain_path = tuple(range(i, j))
            seam_path = tuple(range(j, size)) + tuple(range(i))
            route_models = []
            for other in (seam_path, plain_path):
                # For d >= 3 an interior cut preserves both exterior bonds of the surviving path.
                # For d = 2 include both single cuts, subtracting their common double-cut term.
                cuts = ((other[1],),) if distance >= 3 else ((other[0],), (other[1],), other)
                route_models.append([exact_model(methods, size, max_order, cuts=cut) for cut in cuts])
            for axis, order in zip(CHANNELS, powers_for_distance(distance)):
                total = series(methods, ring, (i, j), axis, order)[-1]
                routes = []
                for models in route_models:
                    terms = [series(methods, model, (i, j), axis, order)[-1] for model in models]
                    routes.append(terms[0] if distance >= 3 else terms[0] + terms[1] - terms[2])
                plain, seam = routes
                residual = Fraction(0)
                if size == 4 and axis in "YZ":
                    residual = Fraction(2401 * (1 if i else -1), 230400 if axis == "Y" else 5760000)
                assert total == plain + seam + residual, (size, i, j, axis, total, plain, seam)
                if size == 4 or (size == 8 and i == 0):
                    print(f"N={size} pair ({i},{j}) {axis * 2} t^{order}: "
                          f"plain {plain}, seam {seam}, full-cycle {residual}, sum {total}")
    print("Cut decompositions matched all antipodal coefficients at N=4,6,8,10.")


def polynomial_remainder(samples, quadratic):
    # Newton interpolation in k = 14 lambda, with values at consecutive integer k.
    # Keep only the remainder modulo a*lambda**2 + b*lambda + c.
    terms = []
    for order in range(len(samples)):
        terms.append(samples[0] / factorial(order))
        samples = [right - left for left, right in zip(samples, samples[1:])]
    a, b, c = map(Fraction, quadratic)
    constant = linear = Fraction(0)
    for index in reversed(range(len(terms))):
        constant, linear = (
            terms[index] - index * constant - 14 * c * linear / a,
            14 * constant - index * linear - 14 * b * linear / a,
        )
    return constant, linear


def verify_cancellation_roots(methods):
    quadratics = {"Y": (745, -302, -23), "Z": (185, -78, -7)}
    orders = {"Y": 16, "Z": 18}
    samples = {axis: [] for axis in quadratics}
    for k in range(21):
        data = exact_model(methods, 8, 20, seam=k)
        for axis, order in orders.items():
            samples[axis].append(series(methods, data, (0, 4), axis, order + 2))
    # At Taylor order p each coefficient is a polynomial in lambda of degree at most p.
    # These 21 distinct exact samples suffice through order 20.
    for axis, quadratic in quadratics.items():
        order = orders[axis]
        a, b, c = quadratic
        factor = samples[axis][0][order] / c
        assert factor > 0
        for k, coefficients in enumerate(samples[axis]):
            lam = Fraction(k, 14)
            assert not any(coefficients[:order])
            assert coefficients[order] == factor * (a * lam**2 + b * lam + c)
        discriminant = b*b - 4*a*c
        assert isqrt(discriminant)**2 != discriminant
        # An irreducible quadratic cannot share a root with a nonzero linear remainder.
        assert any(polynomial_remainder([row[order + 2] for row in samples[axis]], quadratic))
        root = (-b + sqrt(discriminant)) / (2 * a)
        print(f"N=8 pair (0,4) {axis * 2}: lambda={root:.12f} cancels t^{order}; first power is t^{order + 2}.")


def main():
    methods = runpy.run_path(str(PRIOR / "run.py"))
    verify_powers(methods)
    verify_tied_routes(methods)
    verify_cancellation_roots(methods)


if __name__ == "__main__":
    main()
