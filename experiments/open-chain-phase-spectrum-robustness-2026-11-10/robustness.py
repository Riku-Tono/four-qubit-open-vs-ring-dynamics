"""Evaluate fixed angle errors on the previous open-chain holdout windows."""

import csv
import gzip
import hashlib
from itertools import product
import json
from pathlib import Path
from runpy import run_path
import sys

import numpy as np


OUTPUT = Path(__file__).resolve().parent
PRIOR = OUTPUT.parent / "open-chain-phase-persistence-2026-11-10"
sys.path.insert(0, str(OUTPUT.parent / "four-site-in-phase-2026-11-10"))
from metrics import (GATES, LEFT, MAIN_GATE, MAX_LAG, MIN_COVERAGE, MIN_PLV, PAIRS, RIGHT,
                     Model, describe, phase_samples, product_state, time_grid)


AXES = ("end_theta", "end_phi", "inner_theta", "inner_phi")
ERRORS = (-0.05, -0.01, 0, 0.01, 0.05)
GROUP_FIELDS = tuple(f"delta_{axis}" for axis in AXES)
SUMMARY_METRICS = ("max_abs_lag", "min_plv", "min_envelope", "min_coverage", "min_half_amplitude",
                   "worst_instantaneous_lag")
UPPER = {"max_abs_lag", "worst_instantaneous_lag"}
BATCH = 24


def sources():
    selection_path = PRIOR / "selection.json"
    checksum = hashlib.sha256(selection_path.read_bytes()).hexdigest()
    selection = json.loads(selection_path.read_text(encoding="utf-8"))["selection"][str(MAIN_GATE)]
    nominal = np.asarray(selection["families"][selection["selected_family"]]["angles_rad"])
    results = json.loads((PRIOR / "results.json").read_text(encoding="utf-8"))
    assert checksum == results["selection_sha256"]
    windows = {name: entry["bounds"] for name, entry in results["windows"].items() if entry["role"] == "holdout"}
    assert len(windows) == 24
    return nominal, windows, results["dt"], results["cases"][f"selected_{MAIN_GATE}"], checksum


def conditions():
    metadata, deltas = [], []

    def add(case, family, scope, site, axis, error, groups, changes):
        metadata.append(dict(case=case, family=family, scope=scope, site=site, axis=axis, error=error,
                             **dict(zip(GROUP_FIELDS, groups)), max_abs_error=np.abs(changes).max(),
                             n_changed=int(np.count_nonzero(groups if family == "mirror_grid" else changes))))
        deltas.append(changes)

    add("nominal", "nominal", "nominal", None, None, None, (0,) * 4, np.zeros((4, 2)))
    for index, offsets in enumerate(product(ERRORS, repeat=4)):
        changed = [i for i, offset in enumerate(offsets) if offset]
        if not changed:
            continue
        single = len(changed) == 1
        expanded = np.asarray(offsets).reshape(2, 2)[[0, 1, 1, 0]]
        add(f"m{index:03}", "mirror_grid", "one_angle" if single else "combined", None,
            AXES[changed[0]] if single else None, offsets[changed[0]] if single else None, offsets, expanded)
    for site, axis, error in product(range(4), range(2), (value for value in ERRORS if value)):
        expanded = np.zeros((4, 2))
        expanded[site, axis] = error
        add(f"s{len(metadata) - 624:02}", "single_site", "one_site", site, ("theta", "phi")[axis], error,
            (None,) * 4, expanded)
    assert len(metadata) == 657
    return metadata, np.asarray(deltas)


def measure(occupation, analytic):
    """Batch version of the shared phase and amplitude evaluator, reusing the phase samples across gates."""
    envelope = np.abs(analytic)
    min_env, relative = phase_samples(analytic)
    worst_real = relative.real.min(axis=-1)
    weights = np.ones(occupation.shape[1])
    weights[[0, -1]] = 0.5
    total = weights.sum()
    amplitude = np.ptp(occupation, axis=1) / 2
    centered = occupation - np.average(occupation, weights=weights, axis=1)[:, None, :]
    variance = np.average(centered**2, weights=weights, axis=1)
    covariance = np.average(centered[..., LEFT] * centered[..., RIGHT], weights=weights, axis=1)
    normalization = np.sqrt(variance[:, LEFT] * variance[:, RIGHT])
    pearson = np.divide(covariance, normalization, out=np.full_like(covariance, np.nan),
                        where=(variance[:, LEFT] > 1e-26) & (variance[:, RIGHT] > 1e-26))
    unmasked = np.einsum("t,btp->bp", weights, relative) / total
    common = {"amplitude": amplitude, "site_min_envelope": envelope.min(axis=1),
              "min_envelope": min_env.min(axis=1), "pearson": pearson}
    for gate in GATES:
        valid = min_env >= gate
        partial = ~np.all(valid, axis=1)
        locking = unmasked.copy()
        coverage = np.ones(len(occupation))
        if partial.any():
            masked = weights * valid[partial]
            duration = masked.sum(axis=1)
            numerator = np.einsum("bt,btp->bp", masked, relative[partial])
            locking[partial] = np.divide(numerator, duration[:, None], out=np.full_like(numerator, np.nan),
                                         where=duration[:, None] > 0)
            coverage[partial] = duration / total
        plv = np.minimum(1, np.abs(locking))
        lag = np.where(plv > 1e-12, np.angle(locking, deg=True), np.nan)
        instantaneous = np.rad2deg(np.arccos(np.clip(np.min(np.where(valid, worst_real, 1), axis=1), -1, 1)))
        instantaneous[coverage == 0] = np.nan
        yield gate, {**common, "lag_deg": lag, "plv": plv, "in_phase": locking.real, "coverage": coverage,
                     "max_abs_lag": np.max(np.abs(lag), axis=1), "min_plv": plv.min(axis=1),
                     "worst_instantaneous_lag": instantaneous,
                     "passes": ((amplitude.min(axis=1) >= gate) & (coverage >= MIN_COVERAGE)
                                & (np.abs(lag).max(axis=1) <= MAX_LAG) & (plv.min(axis=1) >= MIN_PLV))}


def scalar(value):
    return value.item() if np.isfinite(value) else None


def window_row(meta, gate, window, bounds, data, i):
    return dict(
        case=meta["case"], gate=gate, window=window, start=bounds[0], end=bounds[1],
        **{key: scalar(data[key][i]) for key in ("passes", "coverage", "min_envelope", "max_abs_lag",
                                                 "min_plv", "worst_instantaneous_lag")},
        **{f"{metric}{j}": scalar(value) for metric in ("amplitude", "site_min_envelope")
           for j, value in enumerate(data[metric][i])},
        **{f"{metric}_{pair}": scalar(value) for metric in ("lag_deg", "plv", "pearson", "in_phase")
           for pair, value in zip(PAIRS, data[metric][i])},
    )


def summarize(meta, gate, rows, bound):
    metrics = {"max_abs_lag": "max_abs_lag", "min_plv": "min_plv", "min_envelope": "min_envelope",
               "min_coverage": "coverage", "worst_instantaneous_lag": "worst_instantaneous_lag"}
    values = {}
    for key, source in metrics.items():
        defined = [row for row in rows if row[source] is not None]
        selected = ((max if key in UPPER else min)(defined, key=lambda row: row[source]) if defined else None)
        values[key] = selected[source] if selected else None
        values[f"{key}_window"] = selected["window"] if selected else None
    return {**meta, "gate": gate, "pass_windows": sum(row["passes"] for row in rows), **values,
            "min_half_amplitude": min(row[f"amplitude{i}"] for row in rows for i in range(4)),
            "spectral_min_envelope_lower": bound.get("min_envelope_lower"),
            "spectral_max_instantaneous_lag_upper": bound.get("max_abs_lag_upper_deg"),
            "spectral_min_plv_lower": bound.get("min_plv_lower")}


def group_summary(summaries):
    groups = {"mirror_all": lambda row: row["family"] == "mirror_grid",
              "single_all": lambda row: row["family"] == "single_site"}
    for scope, size in product(("one_angle", "combined", "one_site"), (0.01, 0.05)):
        groups[f"{scope}_{size}"] = lambda row, s=scope, e=size: row["scope"] == s and row["max_abs_error"] == e
    worst = {}
    for key, keep in groups.items():
        worst[key] = {}
        for gate in GATES:
            subset = [row for row in summaries if row["gate"] == gate and keep(row)]
            selected = {}
            for metric in (*SUMMARY_METRICS, "spectral_min_envelope_lower", "spectral_max_instantaneous_lag_upper"):
                defined = [row for row in subset if row[metric] is not None]
                highest = metric in UPPER or metric == "spectral_max_instantaneous_lag_upper"
                if defined:
                    selected[metric] = (max if highest else min)(defined, key=lambda row: row[metric])
            worst[key][str(gate)] = {"conditions": len(subset),
                                     "conditions_all_pass": sum(row["pass_windows"] == 24 for row in subset),
                                     "min_pass_windows": min(row["pass_windows"] for row in subset),
                                     "worst": selected}
    return worst


def update_audit(audit, row, reference):
    for key in ("max_abs_lag", "min_plv", "min_envelope", "coverage", "worst_instantaneous_lag",
                *(f"amplitude{i}" for i in range(4)),
                *(f"{m}_{p}" for m in ("lag_deg", "plv") for p in PAIRS)):
        if key in reference and row[key] is not None and reference[key] is not None:
            delta = abs(row[key] - reference[key])
            audit[key] = max(audit.get(key, 0), delta)
    if row["passes"] != reference["passes"]:
        audit.setdefault("pass_changes", []).append((row["case"], row["gate"], row["window"]))


def check_fine(model, states, metadata, windows, dt, window_results, groups):
    cases = {"nominal"}
    for gates in groups.values():
        for entry in gates.values():
            cases.update(row["case"] for row in entry["worst"].values())
    indices = [i for i, row in enumerate(metadata) if row["case"] in cases]
    audit = {"cases": [metadata[i]["case"] for i in indices], "pass_changes": []}
    for name, bounds in windows.items():
        grid = model.prepare(time_grid(bounds, dt / 2))
        for start in range(0, len(indices), BATCH // 2):
            selected = indices[start:start + BATCH // 2]
            for gate, data in measure(*grid.sample(states[selected])):
                for j, index in enumerate(selected):
                    meta = metadata[index]
                    row = window_row(meta, gate, name, bounds, data, j)
                    update_audit(audit, row, window_results[meta["case"], gate][name])
    return audit


def main():
    nominal, windows, dt, prior, checksum = sources()
    metadata, deltas = conditions()
    states = product_state(nominal + deltas)
    model = Model("open")
    certificate = run_path(str(PRIOR / "certificate.py"))["spectral_bound"]
    bounds_by_case = {meta["case"]: certificate(model, state) for meta, state in zip(metadata, states)}
    window_results = {(meta["case"], gate): {} for meta in metadata for gate in GATES}
    prior_audit = {"pass_changes": []}
    engine_audit = {"pass_changes": []}
    with gzip.open(OUTPUT / "robustness_windows.csv.gz", "wt", encoding="utf-8", newline="", compresslevel=5) as file:
        writer = None
        for name, bounds in windows.items():
            grid = model.prepare(time_grid(bounds, dt))
            for start in range(0, len(states), BATCH):
                occ, analytic = grid.sample(states[start:start + BATCH])
                for gate, data in measure(occ, analytic):
                    for j, meta in enumerate(metadata[start:start + BATCH]):
                        row = window_row(meta, gate, name, bounds, data, j)
                        if start == 0 and j == 0:
                            reference = prior[name][str(gate)]
                            flat = {**reference, **{f"amplitude{i}": a for i, a in enumerate(reference["amplitude"])},
                                    **{f"{m}_{p}": d[m] for p, d in reference["pairs"].items() for m in ("lag_deg", "plv")}}
                            update_audit(prior_audit, row, flat)
                        if name == next(iter(windows)) and start == 0 and j < 2:
                            reference = describe(occ[j], analytic[j], gate)
                            flat = {**reference, **{f"amplitude{i}": a for i, a in enumerate(reference["amplitude"])},
                                    **{f"{m}_{p}": d[m] for p, d in reference["pairs"].items() for m in ("lag_deg", "plv")}}
                            update_audit(engine_audit, row, flat)
                        if writer is None:
                            writer = csv.DictWriter(file, fieldnames=list(row))
                            writer.writeheader()
                        writer.writerow(row)
                        window_results[meta["case"], gate][name] = row
            print(f"Evaluated {name}", flush=True)

    summaries = [summarize(meta, gate, list(window_results[meta["case"], gate].values()), bounds_by_case[meta["case"]])
                 for meta in metadata for gate in GATES]
    baseline = {row["gate"]: row for row in summaries if row["case"] == "nominal"}
    for row in summaries:
        for metric in ("pass_windows", *SUMMARY_METRICS):
            a, b = row[metric], baseline[row["gate"]][metric]
            row[f"change_{metric}"] = a - b if a is not None and b is not None else None
    with (OUTPUT / "robustness_summary.csv").open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    worst = group_summary(summaries)
    fine = check_fine(model, states, metadata, windows, dt, window_results, worst)
    for audit in (prior_audit, engine_audit):
        assert not audit["pass_changes"] and max(value for key, value in audit.items() if key != "pass_changes") < 1e-9
    result = {"nominal_angles_rad": nominal.tolist(), "prior_selection_sha256": checksum, "windows": windows,
              "dt": dt, "gates": GATES, "condition_count": len(metadata), "nominal": baseline,
              "nominal_spectral_bound": bounds_by_case["nominal"], "groups": worst,
              "audits": {"previous_nominal": prior_audit, "batch_against_shared": engine_audit, "fine_grid": fine}}
    (OUTPUT / "robustness.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
