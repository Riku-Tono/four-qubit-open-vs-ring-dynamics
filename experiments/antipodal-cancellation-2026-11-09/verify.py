"""Independent rational state-vector checks of antipodal coefficients."""

import runpy
from fractions import Fraction
from math import gcd, lcm, sqrt
from pathlib import Path

from formulas import coefficients


SIZES = (4, 6, 8, 10, 12)
CHANNELS = ("XX", "YY", "ZZ")
LAMBDAS = tuple(Fraction(k, 2) for k in range(4))
PRIOR = Path(__file__).resolve().parents[1] / "seam-neighbor-short-time-2026-11-09" / "run.py"


def interpolate_quadratic(values):
    # The four values are at lambda=0, 1/2, 1, 3/2.
    constant, half, one, three_halves = values
    quadratic = 2 * (one - 2 * half + constant)
    linear = one - constant - quadratic
    assert constant + Fraction(3, 2) * linear + Fraction(9, 4) * quadratic == three_halves
    return constant, linear, quadratic


def dense_coefficients(methods, size):
    distance = size // 2
    powers_by_axis = dict(zip(CHANNELS, (4 * distance - 2, 4 * distance, 4 * distance + 2)))
    samples = {(start, axis): [] for start in range(distance) for axis in CHANNELS}
    for lam in LAMBDAS:
        # The existing exact model scales H by 40. Then the seam is 14*lambda.
        model = methods["Ising"](size, int(14 * lam), radius=size + 1)
        powers = methods["hamiltonian_powers"](model, max(powers_by_axis.values()))
        for (start, axis), values in samples.items():
            end = start + distance
            operators = [model.pauli(axis[0], sites) for sites in ((start, end), (start,), (end,))]
            order = powers_by_axis[axis]
            series = methods["connected_series"](powers[:order + 1], operators)
            assert not any(series[:-1]), (size, start, axis, lam)
            values.append(series[-1])
    # At the candidate orders the path counting bounds the seam degree by two,
    # including the four-site full-cycle term, which is linear in the seam.
    return {key: interpolate_quadratic(values) for key, values in samples.items()}


def primitive(polynomial):
    denominator = lcm(*(value.denominator for value in polynomial))
    integers = tuple(int(value * denominator) for value in polynomial)
    divisor = gcd(*integers)
    return tuple(value // divisor for value in integers)


def roots_for_display(polynomial):
    c, b, a = primitive(polynomial)
    if not a:
        return (-c / b,) if b and 0 < -c / b < 1 else ()
    discriminant = b * b - 4 * a * c
    if discriminant < 0:
        return ()
    roots = {(-b + sign * sqrt(discriminant)) / (2 * a) for sign in (-1, 1)}
    return tuple(sorted(root for root in roots if 0 < root < 1))


def main():
    methods = runpy.run_path(str(PRIOR))
    count = formula_count = 0
    for size in SIZES:
        polynomials = dense_coefficients(methods, size)
        for (start, axis), polynomial in polynomials.items():
            if axis in ("YY", "ZZ"):
                expected = coefficients(size, start, axis)
                assert polynomial == expected, (size, start, axis, polynomial, expected)
                formula_count += 1
            roots = roots_for_display(polynomial)
            assert axis != "XX" or not roots
            if roots:
                roots_text = ", ".join(f"{root:.12f}" for root in roots)
                print(f"N={size}, pair ({start},{start + size // 2}), {axis}: {roots_text}", flush=True)
        count += len(polynomials)
        print(f"N={size}: {len(polynomials)} exact polynomial checks", flush=True)
    print(f"{count} coefficient polynomials checked against rational state-vector Taylor.")
    print(f"{formula_count} YY/ZZ analytic polynomials match the independent Taylor calculation.")


if __name__ == "__main__":
    main()
