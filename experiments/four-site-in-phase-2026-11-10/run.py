"""Select product states on the short window and audit fixed states on longer windows."""

import argparse
import csv
import json

import numpy as np

from metrics import (GATES, MAIN_GATE, MAX_LAG, MIN_COVERAGE, MIN_PLV, OUTPUT, PAIRS,
                     TOPOLOGIES, WINDOWS, Model, describe, phase_samples, product_state, time_grid)


SEED = 41070
SEARCH_DT, FINAL_DT = 0.02, 0.002
FAMILIES = {
    "independent": ((0, 1, 2, 3), 32768),
    "mirror": ((0, 1, 1, 0), 8192),
    "uniform": ((0, 0, 0, 0), 1024),
}
STARTS, PROPOSALS, STEPS, FINALISTS = 8, 12, 60, 12


def save_json(name, data):
    (OUTPUT / name).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def wrap(angles):
    angles[..., 0] = np.clip(angles[..., 0], 0, np.pi)
    angles[..., 1] = (angles[..., 1] + np.pi) % (2 * np.pi) - np.pi
    return angles


def random_angles(rng, count, sites):
    return np.stack((np.arccos(rng.uniform(-1, 1, (count, sites))),
                     rng.uniform(-np.pi, np.pi, (count, sites))), axis=-1)


def anchors(sites):
    directions = np.array([[0, 0], [np.pi, 0], [np.pi / 2, np.pi / 2], [np.pi / 2, -np.pi / 2]])
    choices = np.indices((len(directions),) * sites).reshape(sites, -1).T
    return directions[choices]


def short_samples(grid, angles, pattern):
    occupation, analytic = grid.sample(product_state(angles[:, pattern, :]))
    envelope, relative = phase_samples(analytic)
    amplitude = np.ptp(occupation, axis=1).min(axis=1) / 2
    weights = np.ones(occupation.shape[1])
    weights[[0, -1]] = 0.5
    low_index = int((1 - MIN_COVERAGE) * len(weights))
    reserve = np.minimum(amplitude, np.partition(envelope, low_index, axis=1)[:, low_index])
    return envelope, relative, amplitude, weights, reserve


def search_ranks(samples, gate):
    envelope, relative, amplitude, weights, reserve = samples
    masked = (envelope >= gate) * weights
    valid = masked.sum(axis=1)
    coverage = valid / weights.sum()
    locking = np.einsum("bt,btp->bp", masked, relative) / np.maximum(valid[:, None], 1e-15)
    feasible = (coverage >= MIN_COVERAGE) & (amplitude >= gate)
    passing = feasible & (np.abs(locking).min(axis=1) >= MIN_PLV) & (
        np.abs(np.angle(locking, deg=True)).max(axis=1) <= MAX_LAG
    )
    # Column order gives the priority of each comparison. Round off numerical phase ties.
    priority = np.where(passing, 2, feasible.astype(int))
    score = np.where(feasible, np.round(locking.real.min(axis=1), 12), reserve - gate)
    return np.column_stack((priority, score, reserve))


def sorted_rows(ranks):
    return np.lexsort(ranks[:, ::-1].T)[::-1]


def refine(rng, grid, pattern, gate, starting):
    current = starting.copy()
    count, sites, _ = current.shape
    for scale in np.geomspace(0.5, 0.003, STEPS):
        noise = rng.normal(size=(count, PROPOSALS, sites, 2))
        noise[:, 0] = 0
        proposals = wrap(current[:, None] + scale * noise).reshape(-1, sites, 2)
        ranks = search_ranks(short_samples(grid, proposals, pattern), gate).reshape(count, PROPOSALS, 3)
        chosen = [sorted_rows(group)[0] for group in ranks]
        current = proposals[np.arange(count) * PROPOSALS + chosen]
    return current


def choose():
    rng = np.random.default_rng(SEED)
    results = {"seed": SEED, "search_dt": SEARCH_DT, "selection_dt": FINAL_DT,
               "starts": STARTS, "proposals_per_step": PROPOSALS, "steps": STEPS,
               "short_window": WINDOWS["short"], "gates": GATES, "topologies": {}}
    for topology in TOPOLOGIES:
        model = Model(topology)
        coarse = model.prepare(time_grid(WINDOWS["short"], SEARCH_DT))
        fine = model.prepare(time_grid(WINDOWS["short"], FINAL_DT))
        candidates = {gate: {} for gate in GATES}
        for family, (pattern, count) in FAMILIES.items():
            pattern = list(pattern)
            sites = max(pattern) + 1
            fixed = anchors(sites)
            angles = np.concatenate((fixed, random_angles(rng, count, sites)))
            ranks = {gate: [] for gate in GATES}
            for group in np.array_split(angles, np.ceil(len(angles) / 256)):
                samples = short_samples(coarse, group, pattern)
                for gate in GATES:
                    ranks[gate].append(search_ranks(samples, gate))
            for gate in GATES:
                rank = np.concatenate(ranks[gate])
                leading = angles[sorted_rows(rank)[:FINALISTS]]
                optimized = refine(rng, coarse, pattern, gate, leading[:STARTS])
                pool = np.concatenate((leading, optimized))
                fine_ranks = search_ranks(short_samples(fine, pool, pattern), gate)
                best = sorted_rows(fine_ranks)[0]
                expanded = pool[best, pattern, :]
                occupation, analytic = fine.sample(product_state(expanded))
                candidates[gate][family] = {
                    "angles_rad": expanded.tolist(), "rank": fine_ranks[best].tolist(),
                    "short": describe(occupation, analytic, gate), "random_draws": count,
                    "anchors": len(fixed), "initial_amplitude_passes": int(np.sum(rank[:, 0] >= 1)),
                    "initial_phase_passes": int(np.sum(rank[:, 0] == 2)),
                }
            print(f"Selected short-window finalists: {topology} {family}", flush=True)
        results["topologies"][topology] = {
            str(gate): {"selected_family": max(families, key=lambda key: families[key]["rank"]),
                        "families": families} for gate, families in candidates.items()
        }
    save_json("selection.json", results)


def rank_result(result):
    return (result["passes"], result["amplitude_pass"],
            round(result["min_in_phase"], 12) if result["min_in_phase"] is not None else -2,
            result["coverage"], min(result["amplitude"]))


def case_states(model, topology, selection):
    from baseline import preparations

    baseline = preparations(topology)
    names = list(baseline)
    grid = model.prepare(time_grid(WINDOWS["short"], FINAL_DT))
    occ, analytic = grid.sample(np.asarray(list(baseline.values())))
    best = {gate: max(range(len(names)), key=lambda i: rank_result(describe(occ[i], analytic[i], gate)))
            for gate in GATES}
    cases = {"unrotated_1000": baseline["1000"], "best_existing": baseline[names[best[MAIN_GATE]]]}
    sources = {"unrotated_1000": "1000", "best_existing": names[best[MAIN_GATE]]}
    for gate in GATES:
        sources[f"best_existing_{gate}"] = names[best[gate]]
        cases[f"best_existing_{gate}"] = baseline[names[best[gate]]]
        chosen = selection[topology][str(gate)]
        cases[f"selected_{gate}"] = product_state(np.asarray(chosen["families"][chosen["selected_family"]]["angles_rad"]))
    return cases, sources


def evaluate():
    selection = json.loads((OUTPUT / "selection.json").read_text(encoding="utf-8"))["topologies"]
    results = {"windows": WINDOWS, "gates": GATES, "dt": FINAL_DT, "topologies": {}}
    comparison = []
    waveforms = []
    for topology in TOPOLOGIES:
        model = Model(topology)
        cases, sources = case_states(model, topology, selection)
        names = list(cases)
        states = np.asarray(list(cases.values()))
        summaries = {case: {} for case in cases}
        audit = {"max_direct_error": 0, "max_coverage_change": 0, "max_amplitude_change": 0,
                 "max_plv_change": 0, "max_lag_change_deg": 0, "max_pearson_change": 0,
                 "pass_changes": [], "max_ring_waveform_spread": 0}
        for window, bounds in WINDOWS.items():
            time = time_grid(bounds, FINAL_DT)
            occupation, analytic = model.prepare(time).sample(states)
            for index, (case, state) in enumerate(cases.items()):
                occ, signal = occupation[index], analytic[index]
                direct = model.direct(state, time[::59])
                audit["max_direct_error"] = max(audit["max_direct_error"], float(np.max(np.abs(occ[::59] - direct))))
                metrics = {gate: describe(occ, signal, gate) for gate in GATES}
                summaries[case][window] = metrics
                if case.startswith("selected_"):
                    gate = float(case.removeprefix("selected_"))
                    finer = model.prepare(time_grid(bounds, FINAL_DT / 2))
                    fine = describe(*finer.sample(state), gate)
                    coarse = metrics[gate]
                    for key, a, b in (("coverage", coarse["coverage"], fine["coverage"]),
                                      ("amplitude", min(coarse["amplitude"]), min(fine["amplitude"])),
                                      ("plv", coarse["min_plv"], fine["min_plv"]),
                                      ("lag_change_deg", coarse["max_abs_lag"], fine["max_abs_lag"]),
                                      ("pearson", coarse["min_pearson"], fine["min_pearson"])):
                        field = "max_lag_change_deg" if key == "lag_change_deg" else f"max_{key}_change"
                        if a is not None and b is not None:
                            audit[field] = max(audit[field], abs(a - b))
                    if coarse["passes"] != fine["passes"]:
                        audit["pass_changes"].append([case, window])
                    if topology == "ring":
                        audit["max_ring_waveform_spread"] = max(audit["max_ring_waveform_spread"],
                                                               float(np.ptp(occ, axis=1).max()))
                for gate, result in metrics.items():
                    for pair, pair_metrics in result["pairs"].items():
                        comparison.append(dict(topology=topology, case=case, gate=gate, window=window, pair=pair,
                                               **pair_metrics, coverage=result["coverage"],
                                               **{f"amplitude{i}": a for i, a in enumerate(result["amplitude"])},
                                               passes=result["passes"]))
            if window == "full":
                plot_indices = np.unique(np.r_[np.arange(0, len(time), 5), len(time) - 1])
                for case in ("unrotated_1000", "best_existing", *(f"selected_{gate}" for gate in GATES)):
                    for t, values in zip(time[plot_indices], occupation[names.index(case), plot_indices]):
                        waveforms.append(dict(topology=topology, case=case, t=t,
                                              **{f"R{i}": value for i, value in enumerate(values)}))
        assert audit["max_direct_error"] < 1e-10
        results["topologies"][topology] = {"baseline_sources": sources, "cases": summaries, "audit": audit}
        print(f"Validated fixed states: {topology}", flush=True)
    save_json("results.json", results)
    for name, rows in (("comparison.csv", comparison), ("waveforms.csv", waveforms)):
        with (OUTPUT / name).open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", nargs="?", choices=("select", "evaluate", "all"), default="all")
    stage = parser.parse_args().stage
    if stage != "evaluate":
        choose()
    if stage != "select":
        evaluate()
