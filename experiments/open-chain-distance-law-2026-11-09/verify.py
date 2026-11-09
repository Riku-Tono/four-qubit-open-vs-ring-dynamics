"""Verify the near-boundary distance law using exact operator commutators."""

import json
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from math import comb, factorial, lcm
from pathlib import Path


CHANNELS = ("XX", "YY", "ZZ")
SIZES = range(4, 13)
OMEGA, J = Fraction(1), Fraction(7, 10)
EXPERIMENTS = Path(__file__).resolve().parents[1]
PRIOR = EXPERIMENTS / "seam-neighbor-short-time-2026-11-09"
EIGHT_QUBIT = EXPERIMENTS / "seam-neighbor-eight-qubit-2026-11-09"


def single_pauli(axis, site):
    # Each operator is a signed i^(k(k-1)/2) times an ordered k-Majorana word.
    if axis == "X":
        return 3 << (2 * site), 1
    prefix = (1 << (2 * site)) - 1
    last = 2 * site + int(axis == "Y")
    return prefix | (1 << last), (-1)**site


def phase(mask):
    degree = mask.bit_count()
    return degree * (degree - 1) // 2


def commuting_product(left, right):
    left_mask, left_sign = left
    right_mask, right_sign = right
    mask = left_mask ^ right_mask
    # Sorting the concatenation anticommutes each right Majorana past the
    # higher-indexed left Majoranas. Repeated Majoranas square to identity.
    swaps = 0
    pending = right_mask
    while pending:
        bit = pending & -pending
        swaps += (left_mask >> bit.bit_length()).bit_count()
        pending ^= bit
    power = (phase(left_mask) + phase(right_mask) - phase(mask) + 2 * swaps) % 4
    assert power in (0, 2)
    return mask, left_sign * right_sign * (1 if power == 0 else -1)


def initial_expectations(size):
    # All Z words are diagonal in |10...0>. All other Pauli words average to zero.
    z_words = [(0, 1)]
    for site in range(size):
        word = single_pauli("Z", site)
        sign = -1 if site == 0 else 1
        products = [commuting_product(previous, word) for previous in z_words]
        z_words.extend((mask, sign * value) for mask, value in products)
    return dict(z_words)


def transpose_hops(mask, weights, bonds):
    # i[M, a_r] = weights[r] a_(r+1) - weights[r-1] a_(r-1).
    # The transpose moves the expectation back to the input of each hop.
    hops = (mask ^ (mask >> 1)) & bonds
    while hops:
        bit = hops & -hops
        weight = weights[bit.bit_length() - 1]
        yield mask ^ bit ^ (bit << 1), -weight if mask & bit else weight
        hops ^= bit


def moment_derivatives(size, operators, max_order, omega, interaction):
    # M = scale * H. Use integer derivatives throughout the commutator recurrence.
    scale = lcm(omega.denominator, interaction.denominator)
    weights = [int(scale * value) for value in (omega, interaction) * (size - 1) + (omega,)]
    bonds = (1 << len(weights)) - 1
    degrees = {mask.bit_count() for mask, _ in operators}
    moments = {operator: [] for operator in operators}
    state = {mask: value for mask, value in initial_expectations(size).items()
             if mask.bit_count() in degrees}
    transitions = {}
    for order in range(max_order + 1):
        for operator in operators:
            mask, sign = operator
            moments[operator].append(sign * state.get(mask, 0))
        if order == max_order:
            break
        updated = {}
        for mask, value in state.items():
            if mask not in transitions:
                transitions[mask] = tuple(transpose_hops(mask, weights, bonds))
            for target, weight in transitions[mask]:
                updated[target] = updated.get(target, 0) + weight * value
        state = {mask: value for mask, value in updated.items() if value}
    return scale, moments


def exact_series(size, max_order, omega, interaction):
    observables = {}
    for start in (0, 1):
        end = start + size - 2
        for channel in CHANNELS:
            left = single_pauli(channel[0], start)
            right = single_pauli(channel[0], end)
            observables[start, channel] = commuting_product(left, right), left, right
    operators = {operator for triple in observables.values() for operator in triple}
    scale, moments = moment_derivatives(size, operators, max_order, omega, interaction)
    series = {}
    for key, (joint, left, right) in observables.items():
        series[key] = [
            Fraction(moments[joint][n] - sum(
                comb(n, k) * moments[left][k] * moments[right][n - k]
                for k in range(n + 1)
            ), scale**n * factorial(n))
            for n in range(max_order + 1)
        ]
    return series


def predicted(distance, start, channel, omega, interaction):
    if channel == "XX":
        q = 2 * distance - 1
        b = 2**(q - 1)
        numerator = b * (b - q*q + q - 2) if start == 0 else b*b
        coefficient = Fraction(numerator, factorial(q)**2 * (q + 1))
        return 4 * distance - 2, coefficient * omega**(2 * distance) * interaction**(2 * distance - 2)

    m = 2 * distance + 1
    denominator = factorial(m + 1) * factorial(m)
    if channel == "YY":
        common = m * 2**(m - 4) * (3*m*m + 7*m - 2)
        correction = 4*m*(2*m - 1) if start == 0 else 2*m*m*(m + 1)
        coefficient = Fraction(correction - common, denominator)
        return 4 * distance, coefficient * omega**(2 * distance) * interaction**(2 * distance)

    common = (m + 2) * 2**(m - 2)
    correction = 4*m if start == 0 else 2*(m + 1)
    coefficient = Fraction(correction - common, denominator)
    return 4 * distance + 2, coefficient * omega**(2 * distance + 2) * interaction**(2 * distance)


def archived_coefficients():
    prior = json.loads((PRIOR / "results.json").read_text(encoding="utf-8"))
    eight = json.loads((EIGHT_QUBIT / "results.json").read_text(encoding="utf-8"))
    archived = {(row["N"], int(row["pair"][0]), row["channel"]):
                (row["order"], Fraction(row["coefficient"]))
                for row in prior["cases"] if row["branch"] == "open"}
    for row in eight["pair_cases"]:
        if row["branch"] == "open" and row["pair"] in ("06", "17"):
            for channel, leading in row["leading"].items():
                archived[8, int(row["pair"][0]), channel] = leading["order"], Fraction(leading["coefficient"])
    assert len(archived) == 18
    return archived


def state_vector_module():
    spec = spec_from_file_location("short_time", PRIOR / "run.py")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def state_vector_series(size, max_order, module):
    model = open_model(size, module)
    powers = module.hamiltonian_powers(model, max_order)
    return {
        (start, channel): module.connected_series(
            powers, module.make_observables(model, f"{start}{start + size - 2}", channel)
        )
        for start in (0, 1)
        for channel in CHANNELS
    }


def open_model(size, module):
    radius = int((size * OMEGA + (size - 1) * J) / 2) + 1
    return module.Ising(size, 0, radius=radius)


def verify(size, omega, interaction, archived, direct=None):
    distance = size - 2
    max_order = 4 * distance + 2
    exact = exact_series(size, max_order, omega, interaction)
    if direct is not None:
        assert exact == state_vector_series(size, max_order, direct), size
    for (start, channel), series in exact.items():
        order, coefficient = predicted(distance, start, channel, omega, interaction)
        assert coefficient != 0
        assert not any(series[:order]), (size, start, channel, series)
        assert series[order] == coefficient, (size, start, channel, series[order], coefficient)
        key = size, start, channel
        if key in archived:
            assert archived[key] == (order, coefficient), (key, archived[key], order, coefficient)


def interval_sum(size, left, right, channel):
    distance = right - left
    m = 2 * distance
    left_offsets = range(0, m + 1, 2) if left else (0,)
    weight = Fraction(0)
    for p in left_offsets:
        right_offsets = range(p, m + 1, 2) if right < size - 1 else (m,)
        for q in right_offsets:
            denominator = factorial(p) * factorial(q + 2) * factorial(m - p + 2) * factorial(m - q)
            term = Fraction((q - p + 2)**2, denominator)
            if channel == "YY":
                term *= (p + q + 1) * (2*m - p - q + 1)
            negative = (left == 0 and q > 0) or (left == 1 and p > 0)
            weight += -term if negative else term
    if channel == "YY":
        return 4 * distance, weight * OMEGA**m * J**m
    return 4 * distance + 2, weight * OMEGA**(m + 2) * J**m


def verify_all_open_pairs(module):
    count = 0
    for size in range(3, 10):
        model = open_model(size, module)
        powers = module.hamiltonian_powers(model, 4 * (size - 1) + 2)
        for i in range(size - 2):
            for j in range(i + 2, size):
                for channel in CHANNELS:
                    series = module.connected_series(powers, module.make_observables(model, f"{i}{j}", channel))
                    if channel == "XX":
                        expected = predicted(j - i, int(i > 0), channel, OMEGA, J)
                    else:
                        expected = interval_sum(size, i, j, channel)
                    assert module.leading(series) == expected, (size, i, j, channel, expected)
                    count += 1
    print(f"{count} all-pair leading terms at d >= 2 match direct state-vector Taylor for N=3 through N=9.")


def verify_lienhard_static_limit(module):
    # Lienhard et al. (2018), published Eqs. (E5)-(E7) with delta = U = 2J.
    interaction = 2 * J
    expected = (
        OMEGA**4 * interaction**2 / 144,
        OMEGA**6 * interaction**4 / 21600,
        OMEGA**8 * interaction**6 / 12700800,
    )
    model = module.Ising(8, 14, radius=7)
    model.initial = [1] + [0] * (model.dim - 1)
    powers = module.hamiltonian_powers(model, 14)
    for distance, coefficient in enumerate(expected, start=1):
        observables = module.make_observables(model, f"0{distance}", "ZZ")
        order, actual = module.leading(module.connected_series(powers, observables))
        assert (order, actual / 4) == (4 * distance + 2, coefficient)
    print("Three Lienhard static-limit density coefficients match direct state-vector Taylor.")


def main():
    archived = archived_coefficients()
    direct = state_vector_module()
    for size in SIZES:
        verify(size, OMEGA, J, archived, direct if size <= 9 else None)
        print(f"N={size}: six exact leading terms verified through t^{4*(size - 2) + 2}", flush=True)
    for size in (4, 5, 8, 9):
        verify(size, Fraction(3, 4), Fraction(-2, 3), {})
    print("All series for N=4 through N=9 match the independent state-vector Taylor method.")
    print("18 archived leading terms match. Four sizes also verified at Omega=3/4, J=-2/3.")
    verify_all_open_pairs(direct)
    verify_lienhard_static_limit(direct)


if __name__ == "__main__":
    main()
