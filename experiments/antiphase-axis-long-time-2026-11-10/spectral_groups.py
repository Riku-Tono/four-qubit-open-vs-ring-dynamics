"""Certify equality of cross-parity energy gaps with integer polynomials."""

import math
from fractions import Fraction

import numpy as np


_BRACKET_DENOMINATOR = 10**8
_ENERGY_SCALE = 20


def _characteristic(matrix):
    """Return descending characteristic coefficients by Faddeev–LeVerrier."""
    integer_matrix = np.asarray(matrix, dtype=object)
    size = len(integer_matrix)
    identity = np.eye(size, dtype=object)
    accumulator = identity
    coefficients = [1]
    for order in range(1, size + 1):
        product = integer_matrix @ accumulator
        coefficient, remainder = divmod(-int(np.trace(product)), order)
        if remainder:
            raise ArithmeticError("Noninteger characteristic coefficient")
        coefficients.append(coefficient)
        accumulator = product + coefficient * identity
    return coefficients


def _primitive(polynomial):
    first = next((i for i, value in enumerate(polynomial) if value), len(polynomial))
    coefficients = polynomial[first:]
    if not coefficients:
        return []
    divisor = math.gcd(*coefficients)
    if coefficients[0] < 0:
        divisor = -divisor
    return [value // divisor for value in coefficients]


def _pseudo_remainder(dividend, divisor):
    """Compute a remainder up to a nonzero rational scalar."""
    while len(dividend) >= len(divisor):
        leading = dividend[0]
        remainder = [divisor[0] * value for value in dividend]
        for i, coefficient in enumerate(divisor):
            remainder[i] -= leading * coefficient
        dividend = _primitive(remainder)
    return dividend


def _squarefree_part(polynomial):
    degree = len(polynomial) - 1
    derivative = [(degree - i) * value for i, value in enumerate(polynomial[:-1])]
    gcd, remainder = _primitive(polynomial), _primitive(derivative)
    while remainder:
        gcd, remainder = remainder, _pseudo_remainder(gcd, remainder)

    dividend = polynomial.copy()
    quotient = []
    for i in range(len(dividend) - len(gcd) + 1):
        coefficient, remainder = divmod(dividend[i], gcd[0])
        if remainder:
            raise ArithmeticError("Polynomial gcd did not divide exactly")
        quotient.append(coefficient)
        for j, value in enumerate(gcd):
            dividend[i + j] -= coefficient * value
    if any(dividend):
        raise ArithmeticError("Nonzero remainder in square-free factorization")
    return quotient, gcd


def _scaled_value(polynomial, numerator):
    """Evaluate the sign of a polynomial at numerator / denominator."""
    value = polynomial[0]
    denominator_power = 1
    for coefficient in polynomial[1:]:
        denominator_power *= _BRACKET_DENOMINATOR
        value = value * numerator + coefficient * denominator_power
    return value


def _clusters(intervals):
    groups = []
    by_lower_bound = sorted(enumerate(intervals), key=lambda item: item[1])
    for index, (lower, upper) in by_lower_bound:
        if groups and lower <= groups[-1][1]:
            group = groups[-1]
            group[1] = max(group[1], upper)
            group[2].append(index)
        else:
            groups.append([lower, upper, [index]])
    return groups


def _eigenvalue_brackets(block, eigenvalues):
    if not np.array_equal(block, block.T):
        raise ValueError("Parity block must be symmetric")
    if len(eigenvalues) != len(block):
        raise ValueError("Provide all eigenvalues of each parity block")
    approximate_brackets = []
    for value in eigenvalues:
        center = round(float(value) * _BRACKET_DENOMINATOR)
        approximate_brackets.append((center - 2, center + 2))

    brackets = [None] * len(block)
    groups = _clusters(approximate_brackets)
    polynomial = _characteristic(block)
    squarefree, repeated = _squarefree_part(polynomial)
    if len(groups) != len(squarefree) - 1:
        raise ArithmeticError("Incomplete isolation of parity-block roots")
    multiplicity_factors = [squarefree]
    while len(repeated) > 1:
        squarefree, repeated = _squarefree_part(repeated)
        multiplicity_factors.append(squarefree)

    for lower, upper, members in groups:
        contains_root = [
            _scaled_value(factor, lower) * _scaled_value(factor, upper) < 0
            for factor in multiplicity_factors
        ]
        if not contains_root[0] or sum(contains_root) != len(members):
            raise ArithmeticError("Parity-block root multiplicity not certified")
        for member in members:
            brackets[member] = (lower, upper)
    return brackets, polynomial


def _gap_difference(plus, minus):
    return (
        np.kron(plus, np.eye(len(minus), dtype=int))
        - np.kron(np.eye(len(plus), dtype=int), minus)
    )


def certify_gap_groups(plus_block, minus_block, plus_eigenvalues, minus_eigenvalues):
    """Return signed gap IDs and a certificate for the integer blocks of K=20H.

    The positive IDs increase with the size of the physical gap. Floating-point
    eigenvalues locate narrow rational brackets. Exact signs and degree counts
    certify those brackets and the equality of all gaps assigned the same ID.
    """
    plus_block, minus_block = (np.asarray(block) for block in (plus_block, minus_block))
    if not all(np.issubdtype(block.dtype, np.integer) for block in (plus_block, minus_block)):
        raise TypeError("The parity blocks of K=20H must contain integers")
    plus, plus_poly = _eigenvalue_brackets(plus_block, plus_eigenvalues)
    minus, minus_poly = _eigenvalue_brackets(minus_block, minus_eigenvalues)

    # This Kronecker difference has the multiset of all cross-parity gaps.
    polynomial = _characteristic(_gap_difference(plus_block, minus_block))
    if any(polynomial[1::2]) or not polynomial[-1]:
        raise ArithmeticError("Expected an even gap polynomial with nonzero constant term")
    squared = polynomial[::2]
    squarefree, repeated = _squarefree_part(squared)

    intervals = []
    signs = []
    for plus_lower, plus_upper in plus:
        for minus_lower, minus_upper in minus:
            lower, upper = plus_lower - minus_upper, plus_upper - minus_lower
            if lower > 0:
                signs.append(1)
                intervals.append((lower, upper))
            elif upper < 0:
                signs.append(-1)
                intervals.append((-upper, -lower))
            else:
                raise ArithmeticError("A gap bracket contains zero")

    # Every cluster contains a root, and every positive root is represented.
    # Matching the exact number of distinct squared roots proves that each
    # cluster contains precisely one positive gap.
    groups = _clusters(intervals)
    if len(groups) != len(squarefree) - 1:
        raise ArithmeticError("Incomplete separation of distinct positive gaps")
    ids = np.empty(len(intervals), dtype=int)
    multiplicities = []
    for group_id, (_, _, members) in enumerate(groups, start=1):
        for member in members:
            ids[member] = signs[member] * group_id
        positive = sum(signs[member] == 1 for member in members)
        negative = len(members) - positive
        if positive != negative:
            raise ArithmeticError("An even gap polynomial requires matching ± multiplicities")
        multiplicities.append(positive)

    denominator = _BRACKET_DENOMINATOR * _ENERGY_SCALE
    spacing_numerator = min(right[0] - left[1] for left, right in zip(groups, groups[1:]))
    audit = {
        "plus_characteristic_coefficients": plus_poly,
        "minus_characteristic_coefficients": minus_poly,
        "minus_eigenvalue_enclosure_numerators": [list(bracket) for bracket in minus],
        "parity_energy_enclosure_denominator": _BRACKET_DENOMINATOR,
        "cross_parity_difference_count": len(intervals),
        "squared_gap_polynomial_degree": len(squared) - 1,
        "squared_gap_gcd_degree": len(repeated) - 1,
        "squarefree_squared_gap_coefficients": squarefree,
        "distinct_positive_gap_count": len(groups),
        "distinct_signed_gap_count": 2 * len(groups),
        "positive_gap_multiplicities": multiplicities,
        "physical_gap_enclosure_denominator": denominator,
        "positive_gap_enclosure_numerators": [[lower, upper] for lower, upper, _ in groups],
        "minimum_positive_gap_spacing_numerator": spacing_numerator,
        "minimum_positive_gap_spacing_lower_bound": math.nextafter(
            spacing_numerator / denominator, -math.inf
        ),
    }
    return ids.reshape(len(plus), len(minus)), audit


def _open_candidates(gap_audit):
    # The open-chain S(t) cross element is F02+i F01. The 8x8 Jacobi
    # matrix below is unitarily similar to the scaled Majorana generator i20M.
    edges = [20, 14, 20, 14, 20, 14, 20]
    majorana = np.diag(edges, 1) + np.diag(edges, -1)
    brackets, polynomial = _eigenvalue_brackets(majorana, np.linalg.eigvalsh(majorana))
    active = []
    for lower, upper in brackets:
        if upper < 0:
            continue
        matching = [
            index
            for index, (gap_lower, gap_upper) in enumerate(
                gap_audit["positive_gap_enclosure_numerators"], start=1
            )
            if lower < gap_upper and gap_lower < upper
        ]
        if lower <= 0 or len(matching) != 1:
            raise ArithmeticError("A Majorana frequency did not match exactly one gap bracket")
        active.append(matching[0])
    return sorted(set(active)), {"scaled_majorana_characteristic_coefficients": polynomial}


def _half_turn(block, parity):
    """Represent the ring half-turn and its two integer block restrictions."""
    reps = np.arange(8)
    rotated = ((reps & 3) << 2) | (reps >> 2)
    signs = np.where(rotated >= 8, parity, 1)
    targets = np.minimum(rotated, 15 - rotated)
    transform = np.zeros((8, 8), dtype=int)
    transform[targets, reps] = signs
    if not np.array_equal(transform @ block, block @ transform):
        raise ArithmeticError("Ring block lacks half-turn symmetry")

    polynomials = []
    for sector in (1, -1):
        columns, pivots = [], []
        for site, target in enumerate(targets):
            if site > target or (site == target and signs[site] != sector):
                continue
            column = np.zeros(8, dtype=int)
            column[site] = 1
            if site != target:
                column[target] = sector * signs[site]
            columns.append(column)
            pivots.append(site)
        basis = np.column_stack(columns)
        restricted = (block @ basis)[pivots, :]
        if not np.array_equal(block @ basis, basis @ restricted):
            raise ArithmeticError("Half-turn basis is not invariant")
        polynomials.append(_characteristic(restricted))
    return transform, polynomials


def _ring_candidates(plus, minus, group_ids, gap_audit, diagonal):
    plus_turn, plus_polys = _half_turn(plus, 1)
    minus_turn, minus_polys = _half_turn(minus, -1)
    if (plus_turn[0, 0] != 1 or minus_turn[0, 0] != 1
            or not np.array_equal(plus_turn @ diagonal, diagonal @ minus_turn)):
        raise ArithmeticError("Initial reference or observable lacks half-turn symmetry")
    if minus_polys != [[1, 0, -1584, 0, 160000], [1, 0, -800, 0, 160000]]:
        raise ArithmeticError("Unexpected half-turn spectra in the odd parity block")

    # The T²=-1 sector of K_- has characteristic polynomial (x²-400)².
    scaled_twenty = 20 * _BRACKET_DENOMINATOR
    kept_columns = [
        index
        for index, (lower, upper) in enumerate(gap_audit["minus_eigenvalue_enclosure_numerators"])
        if not (lower < scaled_twenty < upper or lower < -scaled_twenty < upper)
    ]
    if len(kept_columns) != 4:
        raise ArithmeticError("Failed to locate the four T²=+1 roots of K_-")
    active = sorted(set(map(int, np.abs(group_ids[:, kept_columns]).ravel())))
    return active, {"plus_half_turn_characteristic_coefficients": plus_polys,
                    "minus_half_turn_characteristic_coefficients": minus_polys}


def _determinant_modulo(matrix, modulus):
    matrix = [row.copy() for row in matrix]
    determinant = 1
    for column in range(len(matrix)):
        pivot = next((row for row in range(column, len(matrix)) if matrix[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            determinant = -determinant
        diagonal = matrix[column][column]
        determinant = determinant * diagonal % modulus
        inverse = pow(diagonal, -1, modulus)
        for row in range(column + 1, len(matrix)):
            scale = matrix[row][column] * inverse % modulus
            for j in range(column + 1, len(matrix)):
                matrix[row][j] = (matrix[row][j] - scale * matrix[column][j]) % modulus
    return determinant


def _integer_moments(plus, minus, diagonal, powers):
    difference = _gap_difference(plus, minus).astype(object)
    vector = diagonal.astype(object).ravel()
    moments = ([], [])
    # For A=|1000>, B=|0000>, the two orientations are e7⊗e0 and e0⊗e7.
    # Odd moments carry one extra nonzero factor of the frequency.
    for power in range(powers):
        parity = power % 2
        moments[parity].append(int(vector[56] + (2 * parity - 1) * vector[7]))
        vector = difference @ vector
    return moments


def _hankel_determinants(moments, count):
    """Certify a nonzero amplitude at each candidate frequency for both slopes."""
    modulus = 1_000_003
    determinants = [
        _determinant_modulo(
            [[sequence[i + j] % modulus for j in range(count)] for i in range(count)], modulus
        )
        for sequence in reversed(moments)
    ]
    if not all(determinants):
        raise ArithmeticError("Did not certify nonzero amplitudes for all candidate frequencies")
    return modulus, *determinants


def _power_sums(polynomial):
    """Newton sums of the roots of the monic squared-gap polynomial."""
    degree = len(polynomial) - 1
    sums = [degree]
    for power in range(1, 2 * degree):
        value = -sum(polynomial[j] * sums[power - j] for j in range(1, min(power, degree + 1)))
        if power <= degree:
            value -= power * polynomial[power]
        sums.append(value)
    return sums


def _inverse_quadratic(matrix, vector):
    """Evaluate v^T M^-1 v for a positive definite integer matrix."""
    rows = [[Fraction(value) for value in [*row, rhs]] for row, rhs in zip(matrix, vector)]
    size = len(rows)
    for i in range(size):
        for j in range(i + 1, size):
            scale = rows[j][i] / rows[i][i]
            for k in range(i + 1, size + 1):
                rows[j][k] -= scale * rows[i][k]
    solution = [0] * size
    for i in reversed(range(size)):
        remainder = rows[i][-1] - sum(rows[i][j] * solution[j] for j in range(i + 1, size))
        solution[i] = remainder / rows[i][i]
    return sum(rhs * value for rhs, value in zip(vector, solution))


def _exact_gram_diagonals(polynomial, moments):
    # With q the distinct squared gaps and V[i,j]=q[j]**i, the two
    # moment matrices are VV^T and V diag(q) V^T. The integer moments
    # are V times 4 times the Y amplitudes, and V times -4 times the X
    # amplitudes weighted by sqrt(q). A squared sine or cosine averages to 1/2.
    count = len(polynomial) - 1
    sums = _power_sums(polynomial)
    yy, xx = (
        _inverse_quadratic(
            [[sums[i + j + shift] for j in range(count)] for i in range(count)],
            moments[shift][:count],
        ) / 32
        for shift in (0, 1)
    )
    return [xx.numerator, xx.denominator], [yy.numerator, yy.denominator]


def certify_active_support(topology, plus_block, minus_block, group_ids, gap_audit, pair_z):
    """Certify the spectral support of f_X and f_Y for the specified pair.

    The exact selection rules give an upper bound on the support. The Hankel
    determinant of each integer moment sequence certifies that all its
    candidates occur: this matrix is a Vandermonde times the amplitude
    diagonal times the transposed Vandermonde.
    """
    plus, minus = np.asarray(plus_block), np.asarray(minus_block)
    if topology not in ("open", "ring"):
        raise ValueError("Expected open or ring topology")
    other_bit = 1 if topology == "open" else 2
    expected_diagonal = 2 - 2 * ((np.arange(8) & other_bit) != 0)
    values = np.asarray(pair_z)[:8]
    if not np.array_equal(values, expected_diagonal):
        raise ValueError("Expected the prescribed S diagonal on representatives 0 through 7")
    diagonal = np.diag(expected_diagonal)
    if topology == "open":
        active, audit = _open_candidates(gap_audit)
        expected = 4
    else:
        active, audit = _ring_candidates(plus, minus, group_ids, gap_audit, diagonal)
        expected = 10
    if len(active) != expected:
        raise ArithmeticError("Unexpected count of candidate frequencies")
    raw_count = gap_audit["distinct_positive_gap_count"]
    moments = _integer_moments(plus, minus, diagonal, max(4 * expected - 2, 2 * raw_count))
    modulus, determinant_x, determinant_y = _hankel_determinants(moments, expected)
    exact_xx, exact_yy = _exact_gram_diagonals(gap_audit["squarefree_squared_gap_coefficients"], moments)
    inactive = sorted(set(range(1, gap_audit["distinct_positive_gap_count"] + 1)) - set(active))
    return {
        **audit,
        "active_positive_gap_ids": active,
        "inactive_positive_gap_ids": inactive,
        "active_positive_gap_count": expected,
        "hankel_modulus": modulus,
        "hankel_determinant_mod_prime_x": determinant_x,
        "hankel_determinant_mod_prime_y": determinant_y,
        "exact_gram_xx": exact_xx,
        "exact_gram_yy": exact_yy,
    }
