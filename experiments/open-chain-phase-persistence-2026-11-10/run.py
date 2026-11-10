"""Select on eight open-chain windows and evaluate frozen states on new windows."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


OUTPUT = Path(__file__).resolve().parent
PRIOR = OUTPUT.parent / "four-site-in-phase-2026-11-10"
sys.path.insert(0, str(PRIOR))
from metrics import (GATES, MAIN_GATE, MAX_LAG, MIN_COVERAGE, MIN_PLV, Model, describe,
                     phase_samples, product_state, time_grid)


BLOCK = 8 * np.pi
SEARCH_WINDOWS = {f"search_{i + 1:02}": (i * BLOCK, (i + 1) * BLOCK) for i in range(8)}
HOLDOUT_WINDOWS = {f"holdout_{i - 7:02}": (i * BLOCK, (i + 1) * BLOCK) for i in range(8, 32)}
AGGREGATES = {"holdout_A": (8 * BLOCK, 16 * BLOCK), "holdout_B": (16 * BLOCK, 32 * BLOCK)}
REFERENCE_WINDOW = {"prior_short_window": (0, 3.688)}
FAMILIES = {
    "independent": ((0, 1, 2, 3), 16384, 8192),
    "mirror": ((0, 1, 1, 0), 8192, 4096),
    "uniform": ((0, 0, 0, 0), 2048, 1024),
}
SEED = 47051
SEARCH_DT, RERANK_DT, FINAL_DT = 0.12, 0.02, 0.002
STARTS, PROPOSALS, STEPS, FINALISTS = 6, 10, 52, 12
SPECTRAL_SEED = np.array([[1.059, 0], [0.893, 0], [0.893, 0], [1.059, 0]])


def save_json(name, value):
    (OUTPUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prior_selections():
    prior = json.loads((PRIOR / "selection.json").read_text(encoding="utf-8"))["topologies"]["open"]
    return {gate: np.array(entry["families"][entry["selected_family"]]["angles_rad"])
            for gate, entry in prior.items()}


def population(rng, pattern, n_uniform, n_x):
    sites = max(pattern) + 1
    axes = np.array([[0, 0], [np.pi, 0], [np.pi / 2, 0], [np.pi / 2, -np.pi],
                     [np.pi / 2, np.pi / 2], [np.pi / 2, -np.pi / 2]])
    fixed = axes[np.indices((len(axes),) * sites).reshape(sites, -1).T]
    # Project reproducible anchors into each family, including the previous short-window states.
    representatives = [pattern.index(site) for site in range(sites)]
    seeded = np.asarray([SPECTRAL_SEED, *prior_selections().values()])[:, representatives, :]
    uniform = np.stack((np.arccos(rng.uniform(-1, 1, (n_uniform, sites))),
                        rng.uniform(-np.pi, np.pi, (n_uniform, sites))), axis=-1)
    x = rng.choice([-1, 1], size=(n_x, sites)) * (1 - 0.65 * rng.random((n_x, sites))**1.5)
    azimuth = rng.uniform(-np.pi, np.pi, (n_x, sites))
    radius = np.sqrt(1 - x**2)
    biased = np.stack((np.arccos(radius * np.cos(azimuth)), np.arctan2(radius * np.sin(azimuth), x)), axis=-1)
    return np.concatenate((fixed, seeded, uniform, biased)), len(fixed), len(seeded)


class TrainingGrid:
    def __init__(self, model, dt):
        grids = np.array([time_grid(bounds, dt) for bounds in SEARCH_WINDOWS.values()])
        self.windows, self.points = grids.shape
        self.grid = model.prepare(grids.ravel())
        self.weights = np.ones(self.points)
        self.weights[[0, -1]] = 0.5

    def sample(self, angles, pattern):
        occ, signal = self.grid.sample(product_state(angles[:, pattern, :]))
        shape = (len(angles), self.windows, self.points, 4)
        min_env, relative = phase_samples(signal.reshape(shape))
        min_amp = np.ptp(occ.reshape(shape), axis=-2).min(axis=-1) / 2
        low_index = int(0.1 * self.points)
        reserve = np.partition(min_env, low_index, axis=-1)[..., low_index].min(axis=-1)
        return min_env, relative, min_amp, reserve

    def ranks(self, samples, gate):
        envelope, relative, amplitude, reserve = samples
        masked = (envelope >= gate) * self.weights
        valid = masked.sum(axis=-1)
        coverage = valid / self.weights.sum()
        locking = np.einsum("bwt,bwtp->bwp", masked, relative) / np.maximum(valid[..., None], 1e-15)
        worst_lag = np.abs(np.angle(locking, deg=True)).max(axis=-1)
        min_plv = np.abs(locking).min(axis=-1)
        slacks = np.minimum.reduce((amplitude / gate - 1, (coverage - MIN_COVERAGE) / (1 - MIN_COVERAGE),
                                    (min_plv - MIN_PLV) / (1 - MIN_PLV), 1 - worst_lag / MAX_LAG))
        return np.column_stack((np.round(slacks.min(axis=-1), 12),
                                np.round(locking.real.min(axis=(1, 2)), 12), reserve))


def rank_order(ranks):
    margin, phase, reserve = ranks.T
    return np.lexsort((reserve, phase, margin))[::-1]


def refine(rng, grid, pattern, gate, starting):
    current = starting.copy()
    count, sites, _ = current.shape
    for scale in np.geomspace(0.35, 0.0015, STEPS):
        noise = rng.normal(size=(count, PROPOSALS, sites, 2))
        noise[:, 0] = 0
        proposals = (current[:, None] + scale * noise).reshape(-1, sites, 2)
        proposals[..., 0] = np.clip(proposals[..., 0], 0, np.pi)
        proposals[..., 1] = (proposals[..., 1] + np.pi) % (2 * np.pi) - np.pi
        ranks = grid.ranks(grid.sample(proposals, pattern), gate).reshape(count, PROPOSALS, 3)
        chosen = [rank_order(group)[0] for group in ranks]
        current = proposals[np.arange(count) * PROPOSALS + chosen]
    return current


def select():
    rng = np.random.default_rng(SEED)
    model = Model("open")
    coarse, medium, fine = (TrainingGrid(model, dt) for dt in (SEARCH_DT, RERANK_DT, FINAL_DT))
    candidates = {str(gate): {} for gate in GATES}
    for family, (pattern, n_uniform, n_x) in FAMILIES.items():
        angles, n_axes, n_seeded = population(rng, pattern, n_uniform, n_x)
        ranks = {gate: [] for gate in GATES}
        for start in range(0, len(angles), 96):
            samples = coarse.sample(angles[start:start + 96], pattern)
            for gate in GATES:
                ranks[gate].append(coarse.ranks(samples, gate))
        for gate in GATES:
            coarse_ranks = np.concatenate(ranks[gate])
            leading = angles[rank_order(coarse_ranks)[:FINALISTS]]
            improved = refine(rng, coarse, pattern, gate, leading[:STARTS])
            pool = np.concatenate((leading, improved))
            medium_ranks = medium.ranks(medium.sample(pool, pattern), gate)
            pool = pool[rank_order(medium_ranks)[:4]]
            fine_ranks = fine.ranks(fine.sample(pool, pattern), gate)
            best = rank_order(fine_ranks)[0]
            candidates[str(gate)][family] = {
                "angles_rad": pool[best, pattern, :].tolist(), "training_rank": fine_ranks[best].tolist(),
                "uniform_draws": n_uniform, "x_biased_draws": n_x, "axis_anchors": n_axes,
                "additional_anchors": n_seeded, "coarse_pass_all": int(np.count_nonzero(coarse_ranks[:, 0] >= 0)),
            }
        print(f"Training finished: {family}", flush=True)

    frozen = {gate: {"selected_family": max(entries, key=lambda f: entries[f]["training_rank"]),
                     "families": entries} for gate, entries in candidates.items()}
    save_json("selection.json", {
        "seed": SEED, "protocol_sha256": sha256(OUTPUT / "PROTOCOL.md"), "search_windows": SEARCH_WINDOWS,
        "target_steps": [SEARCH_DT, RERANK_DT, FINAL_DT], "starts": STARTS, "proposals_per_step": PROPOSALS,
        "steps": STEPS, "selection": frozen,
    })
    print(f"Frozen selection SHA-256: {sha256(OUTPUT / 'selection.json')}", flush=True)


def longest_dropout(time, envelope, gate):
    below = envelope < gate
    switches = np.flatnonzero(below[:-1] != below[1:])
    # Interpolate threshold crossings between adjacent samples of the joint envelope.
    crossing = time[switches] + ((gate - envelope[switches]) * np.diff(time)[switches]
                                / (envelope[switches + 1] - envelope[switches]))
    start = np.r_[time[0] if below[0] else [], crossing[~below[switches]]]
    end = np.r_[crossing[below[switches]], time[-1] if below[-1] else []]
    if not len(start):
        return 0.0, None
    longest = np.argmax(end - start)
    return float(end[longest] - start[longest]), [float(start[longest]), float(end[longest])]


def describe_window(time, occupation, signal, gate):
    metrics = describe(occupation, signal, gate)
    envelope = np.abs(signal)
    joint = envelope.min(axis=-1)
    weights = np.ones(len(time))
    weights[[0, -1]] = 0.5
    duration, interval = longest_dropout(time, joint, gate)
    return {**metrics, "min_envelope": float(joint.min()), "envelope_p10": float(np.quantile(joint, 0.1)),
            "site_min_envelope": envelope.min(axis=0).tolist(),
            "site_below_gate": np.average(envelope < gate, axis=0, weights=weights).tolist(),
            "longest_dropout": duration, "dropout_interval": interval}


def comparison_rows(case, gate, window, role, bounds, metrics):
    summary = {key: metrics[key] for key in ("coverage", "passes", "min_envelope", "envelope_p10",
                                             "longest_dropout", "joint_within_30_deg")}
    sites = {f"amplitude{i}": value for i, value in enumerate(metrics["amplitude"])}
    sites.update({f"below_gate{i}": value for i, value in enumerate(metrics["site_below_gate"])})
    return [dict(case=case, gate=gate, window=window, role=role, start=bounds[0], end=bounds[1],
                 pair=pair, **pair_metrics, **summary, **sites) for pair, pair_metrics in metrics["pairs"].items()]


def write_waveforms(model, states):
    time = time_grid((0, 32 * BLOCK), 0.02)
    grid = model.prepare(time)
    with (OUTPUT / "waveforms.csv").open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(["t", "case", *(f"R{i}" for i in range(4)), *(f"envelope{i}" for i in range(4))])
        for case, state in states.items():
            occupation, signal = grid.sample(state)
            writer.writerows([f"{t:.9f}", case, *(f"{r:.9f}" for r in occ), *(f"{a:.9f}" for a in env)]
                             for t, occ, env in zip(time, occupation, np.abs(signal)))


def evaluate():
    from certificate import spectral_bound

    frozen_path = OUTPUT / "selection.json"
    frozen = json.loads(frozen_path.read_text(encoding="utf-8"))
    assert frozen["protocol_sha256"] == sha256(OUTPUT / "PROTOCOL.md"), "Protocol changed after selection"
    selections = frozen["selection"]
    angles = {f"selected_{gate}": entry["families"][entry["selected_family"]]["angles_rad"]
              for gate, entry in selections.items()}
    states = {case: product_state(np.asarray(state_angles)) for case, state_angles in
              {"basis_0000": np.zeros((4, 2)), "prior_short": prior_selections()[str(MAIN_GATE)], **angles}.items()}
    model = Model("open")
    windows = {name: {"role": role, "bounds": bounds}
               for role, group in (("search", SEARCH_WINDOWS), ("holdout", HOLDOUT_WINDOWS),
                                   ("holdout_aggregate", AGGREGATES), ("reference", REFERENCE_WINDOW))
               for name, bounds in group.items()}
    cases = {case: {} for case in states}
    audit = {"direct_error": 0.0, "max_lag_change_deg": 0.0, "max_plv_change": 0.0,
             "max_coverage_change": 0.0, "max_amplitude_change": 0.0, "max_dropout_change": 0.0,
             "pass_changes": []}
    all_states = np.asarray(list(states.values()))
    with (OUTPUT / "comparison.csv").open("w", encoding="utf-8", newline="") as output:
        writer = None
        for name, window in windows.items():
            bounds = window["bounds"]
            time = time_grid(bounds, FINAL_DT)
            occ, signal = model.prepare(time).sample(all_states)
            for index, (case, state) in enumerate(states.items()):
                direct = model.direct(state, time[::997])
                audit["direct_error"] = max(audit["direct_error"], float(np.max(np.abs(occ[index, ::997] - direct))))
                results = {str(gate): describe_window(time, occ[index], signal[index], gate) for gate in GATES}
                cases[case][name] = results
                for gate, metrics in results.items():
                    rows = comparison_rows(case, gate, name, window["role"], bounds, metrics)
                    if writer is None:
                        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
                        writer.writeheader()
                    writer.writerows(rows)

            # Re-evaluate each frozen state at its own gate, and the main state at all gates.
            fine_time = time_grid(bounds, FINAL_DT / 2)
            fine_grid = model.prepare(fine_time)
            for gate in GATES:
                case = f"selected_{gate}"
                fine_occ, fine_signal = fine_grid.sample(states[case])
                for check_gate in GATES if gate == MAIN_GATE else (gate,):
                    fine = describe_window(fine_time, fine_occ, fine_signal, check_gate)
                    coarse = cases[case][name][str(check_gate)]
                    for field, key in (("max_lag_change_deg", "max_abs_lag"), ("max_plv_change", "min_plv"),
                                       ("max_coverage_change", "coverage"), ("max_dropout_change", "longest_dropout")):
                        if coarse[key] is not None and fine[key] is not None:
                            audit[field] = max(audit[field], abs(coarse[key] - fine[key]))
                    amp_change = np.max(np.abs(np.asarray(coarse["amplitude"]) - fine["amplitude"]))
                    audit["max_amplitude_change"] = max(audit["max_amplitude_change"], float(amp_change))
                    if coarse["passes"] != fine["passes"]:
                        audit["pass_changes"].append([case, name, check_gate])
            print(f"Evaluated {name}", flush=True)

    assert audit["direct_error"] < 1e-10
    bounds = {case: spectral_bound(model, state) for case, state in states.items()}
    save_json("results.json", {"selection_sha256": sha256(frozen_path), "dt": FINAL_DT, "windows": windows,
                              "cases": cases, "audit": audit, "spectral_bounds": bounds})
    write_waveforms(model, states)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("select", "evaluate"))
    stage = parser.parse_args().stage
    if stage == "select":
        select()
    else:
        evaluate()
