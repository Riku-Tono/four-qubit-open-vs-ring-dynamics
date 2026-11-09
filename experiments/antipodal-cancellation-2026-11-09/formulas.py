"""Exact antipodal YY/ZZ coefficient polynomials and their weak-link roots."""

from decimal import Decimal, localcontext
from fractions import Fraction
from math import comb, factorial


AXES = ("YY", "ZZ")
J = Fraction(7, 10)


def pair_weight(length, r, s, axis):
    weight = comb(length, r) * comb(length, s) * (s - r)**2
    if axis == "YY":
        weight *= (length - 1)**2 - (r + s - length)**2
    return weight


def kernel_moments(length):
    return 2**(length - 1), length * 2**(length - 3), (3*length**2 - 2*length) * 2**(length - 5)


def total_weight(length, axis):
    mass, second, fourth = kernel_moments(length)
    total = mass * second
    if axis == "YY":
        total = (length - 1)**2 * total - (mass * fourth - second**2)
    return total


def row_weight(length, r, axis):
    mass, second, fourth = kernel_moments(length)
    offset_squared = (r - length // 2)**2
    row = second + offset_squared * mass
    if axis == "YY":
        row = (length - 1)**2 * row - (fourth - 2*offset_squared*second + offset_squared**2 * mass)
    return comb(length, r) * row


def polynomial(size, start, axis):
    """Return integer (constant, linear, quadratic), with scale Omega^p J^N / ((N+2)!)^2."""
    distance = size // 2
    if size < 4 or size % 2 or not 0 <= start < distance or axis not in AXES:
        raise ValueError("Expected an even ring N >= 4, antipodal start in [0, N/2), and YY or ZZ")
    length = size + 2
    if size == 4:
        if start == 0:
            denominator = 96 if axis == "YY" else 2400
            reduced = (-1, -10, -1)
        elif axis == "YY":
            denominator, reduced = 1440, (-13, 94, -13)
        else:
            denominator, reduced = 21600, (-11, 50, -11)
        return tuple(coefficient * (factorial(length)**2 // denominator) for coefficient in reduced)

    total = total_weight(length, axis)
    row0 = row_weight(length, 0, axis)
    if start == 0:
        row2 = row_weight(length, 2, axis)
        pair02 = pair_weight(length, 0, 2, axis)
        return (-row0 + 2*pair02, total - row0 - 3*row2 + 4*pair02, total - 2*row0 - row2 + 2*pair02)

    h = 2 * (distance - start)
    first = row_weight(length, h, axis)
    second = row_weight(length, h + 2, axis)
    pair = pair_weight(length, h, h + 2, axis)
    quadratic = total - first - 2*second + 2*pair
    linear = -first + 2*pair
    if start == distance - 1:
        return row0, linear + row0 - total, quadratic
    if start == 1:
        return -total + 2*row0, linear, quadratic
    return total, linear, quadratic


def coefficients(size, start, axis):
    """Return the exact Taylor coefficient polynomial at Omega=1 and J=7/10."""
    scale = J**size / factorial(size + 2)**2
    return tuple(value * scale for value in polynomial(size, start, axis))


def weak_root(size, start, axis):
    """Approximate the unique root in (0,1), or return None when the interval is empty."""
    constant, linear, quadratic = polynomial(size, start, axis)
    at_one = constant + linear + quadratic
    if constant * at_one >= 0:
        return None
    with localcontext() as context:
        context.prec = 40
        radical = Decimal(linear**2 - 4*quadratic*constant).sqrt()
        # Choose the crossing direction and rationalize the numerator for small roots.
        if constant < 0:
            return -2*constant / (linear + radical)
        return 2*constant / (-linear + radical)


def evaluate(poly, value):
    constant, linear, quadratic = poly
    return constant + value * (linear + value * quadratic)


def verify_finite_comparisons():
    # Exact separators for YY < ZZ in the near-seam family at small sizes.
    separators = (Fraction(1, 18), Fraction(3, 200), Fraction(1, 200),
                  Fraction(3, 2000), Fraction(9, 20000), Fraction(1, 8000))
    for size, separator in zip(range(6, 18, 2), separators):
        start = size // 2 - 1
        assert evaluate(polynomial(size, start, "YY"), separator) < 0
        assert evaluate(polynomial(size, start, "ZZ"), separator) > 0
    assert evaluate(polynomial(4, 1, "YY"), Fraction(1, 5)) > 0
    assert evaluate(polynomial(4, 1, "ZZ"), Fraction(1, 5)) < 0

    # At these separators the smaller ring's root is higher for both axes.
    for size, denominator in zip(range(4, 18, 2), (10, 30, 100, 300, 1000, 3000, 10000)):
        separator = Fraction(1, denominator)
        for axis in AXES:
            smaller = polynomial(size, size // 2 - 1, axis)
            larger = polynomial(size + 2, size // 2, axis)
            assert evaluate(smaller, separator) * smaller[0] > 0
            assert evaluate(larger, separator) * larger[0] < 0

    # A separator above the near-seam root and below the site-0 root.
    position_separators = {
        "YY": (Fraction(1, 10), Fraction(9, 1000), Fraction(9, 5000), Fraction(21, 50000), Fraction(23, 200000)),
        "ZZ": (Fraction(1, 10), Fraction(9, 1000), Fraction(9, 5000), Fraction(27, 50000), Fraction(147, 1000000)),
    }
    for axis, axis_separators in position_separators.items():
        for size, separator in zip(range(8, 18, 2), axis_separators):
            assert evaluate(polynomial(size, 0, axis), separator) < 0
            assert evaluate(polynomial(size, size // 2 - 1, axis), separator) < 0


def main():
    verify_finite_comparisons()
    print("N, pair, lambda_YY, lambda_ZZ")
    for size in (4, 6, 8, 10, 12, 14, 16, 18, 20, 32, 64):
        distance = size // 2
        for start in (0, distance - 1):
            roots = [weak_root(size, start, axis) for axis in AXES]
            displayed = [format(root, ".12g") if root is not None else "-" for root in roots]
            print(f"{size}, ({start},{start + distance}), {', '.join(displayed)}")


if __name__ == "__main__":
    main()
