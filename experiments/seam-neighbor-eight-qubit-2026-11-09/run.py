"""Check exact leading powers for all N=8 pairs and propagate the matched pairs."""

import json
import runpy
from decimal import Decimal, localcontext
from itertools import combinations
from math import log2
from pathlib import Path


OUTPUT = Path(__file__).resolve().parent
PRIOR = OUTPUT.parent / "seam-neighbor-short-time-2026-11-09"
SIZE, ORDER, RADIUS = 8, 30, 7
MATCHED_PAIRS = ("06", "17")
PRECISIONS = ((110, 80), (130, 96))


def expected_powers(distance):
    if distance == 1:
        return (4, 4, 6)
    return (4 * distance - 2, 4 * distance, 4 * distance + 2)


def main():
    methods = runpy.run_path(str(PRIOR / "run.py"))
    previous = json.loads((PRIOR / "results.json").read_text(encoding="utf-8"))
    parameters = {key: previous[key] for key in ("Omega", "J", "epsilon", "lambda", "reference_lambda")}
    assert tuple(parameters.values()) == (1, 0.7, 0, 0.5, 0)
    previous_differences = {(row["pair"], row["channel"]): (row["order"], row["coefficient"])
                            for row in previous["leading_branch_differences"] if row["N"] == 6}
    times, channels = methods["TIMES"], methods["CHANNELS"]
    series_by_key, pair_cases, matched_cases, differences = {}, [], [], []
    max_relative_change = max_norm_error = 0.0

    for branch, seam in (("open", 0), ("linked", 7)):
        model = methods["Ising"](SIZE, seam, radius=RADIUS)
        powers = methods["hamiltonian_powers"](model, ORDER)
        matched_observables = {}
        for i, j in combinations(range(SIZE), 2):
            pair = f"{i}{j}"
            separation = j - i
            distance = min(separation, SIZE - separation) if seam else separation
            leading_terms = {}
            for channel in channels:
                operators = methods["make_observables"](model, pair, channel)
                series = methods["connected_series"](powers, operators)
                order, value = methods["leading"](series)
                leading_terms[channel] = {"order": order, "coefficient": str(value)}
                if pair in MATCHED_PAIRS:
                    matched_observables[pair, channel] = operators
                    series_by_key[branch, pair, channel] = series
            assert tuple(term["order"] for term in leading_terms.values()) == expected_powers(distance)
            pair_cases.append({"branch": branch, "pair": pair, "distance": distance, "leading": leading_terms})

        values, norm_error = methods["numerics"](model, matched_observables, *PRECISIONS[0])
        refined, refined_norm_error = methods["numerics"](model, matched_observables, *PRECISIONS[1])
        max_norm_error = max(max_norm_error, float(norm_error), float(refined_norm_error))
        with localcontext() as ctx:
            ctx.prec = PRECISIONS[1][0]
            for pair, channel in matched_observables:
                series = series_by_key[branch, pair, channel]
                order, coefficient = methods["leading"](series)
                samples = []
                for t, value, check in zip(times, values[pair, channel], refined[pair, channel]):
                    max_relative_change = max(max_relative_change, float(abs((value - check) / check)))
                    predicted = Decimal(coefficient.numerator) / coefficient.denominator * t**order
                    samples.append({"t": str(t), "C": str(check), "ratio": float(check / predicted)})
                effective_order = log2(abs(refined[pair, channel][1] / refined[pair, channel][0]))
                matched_cases.append({"branch": branch, "pair": pair, "channel": channel, "order": order,
                                      "coefficient": str(coefficient), "series": list(map(str, series)),
                                      "samples": samples, "effective_order_smallest_two_times": effective_order})

    for pair, previous_pair in zip(MATCHED_PAIRS, ("04", "15")):
        for channel in channels:
            linked = series_by_key["linked", pair, channel]
            opened = series_by_key["open", pair, channel]
            order, value = methods["leading"]([a - b for a, b in zip(linked, opened)])
            assert (order, str(value)) == previous_differences[previous_pair, channel]
            differences.append({"pair": pair, "channel": channel, "order": order, "coefficient": str(value)})

    audit = {
        "max_relative_change_between_precision_settings": max_relative_change,
        "max_norm_squared_error": max_norm_error,
        "max_ratio_error_by_time": {str(t): max(abs(row["samples"][i]["ratio"] - 1) for row in matched_cases)
                                    for i, t in enumerate(times)},
        "max_effective_order_error_smallest_two_times": max(abs(row["effective_order_smallest_two_times"] - row["order"])
                                                             for row in matched_cases),
    }
    assert max_relative_change < 1e-35 and max_norm_error < 1e-100
    assert audit["max_ratio_error_by_time"][str(times[0])] < 0.001
    results = {**parameters, "N": SIZE, "initial": "10000000", "max_exact_order": ORDER,
               "chebyshev_radius": RADIUS, "decimal_digits_and_chebyshev_degree": PRECISIONS,
               "pair_cases": pair_cases, "matched_cases": matched_cases,
               "leading_branch_differences": differences, "audit": audit}
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
