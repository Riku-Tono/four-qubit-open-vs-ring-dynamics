"""Scan fixed, asymmetric errors of eight initial-state angles without optimization."""

import csv
import gzip
from itertools import product
import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
REFERENCE = run_path(str(OUTPUT.parent / "open-chain-phase-spectrum-robustness-2026-11-10" / "robustness.py"))
GATES, MAIN_GATE, MIN_COVERAGE, MAX_LAG, MIN_PLV = (
    REFERENCE[k] for k in ("GATES", "MAIN_GATE", "MIN_COVERAGE", "MAX_LAG", "MIN_PLV")
)
LEFT, RIGHT, PAIRS = (REFERENCE[k] for k in ("LEFT", "RIGHT", "PAIRS"))
Model, product_state, time_grid = (REFERENCE[k] for k in ("Model", "product_state", "time_grid"))
WIDTHS = (0.01, 0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.25, 0.30, 0.40, 0.50)
SEED, UNIFORM_SIZE = 41107, 256
DT, TRACKED_LINE = 0.01, 0.502392496
AXES = tuple(f"delta_{axis}{site}" for site in range(4) for axis in ("theta", "phi"))
FAILURES = ("lag", "plv", "coverage", "amplitude", "undefined")


def samples():
    signs = np.array(list(product((-1, 1), repeat=8))).reshape(-1, 4, 2)
    asymmetric = ~np.all(signs == signs[:, ::-1], axis=(1, 2))
    random = np.random.default_rng(SEED).uniform(-1, 1, (UNIFORM_SIZE, 4, 2))
    assert asymmetric.sum() == 240 and np.all(np.any(random != random[:, ::-1], axis=(1, 2)))
    labels = [(f"c{i:03}", "corner") for i in np.flatnonzero(asymmetric)]
    labels += [(f"u{i:03}", "uniform") for i in range(len(random))]
    return labels, np.concatenate((signs[asymmetric], random))


def spectrum(model, states):
    initial = states @ model.vectors
    transition = initial[:, model.a].conj() * initial[:, model.b]
    return np.add.reduceat(transition[:, None, :] * model.observed, model.group_starts, axis=-1)


def spectral_records(coefficients, carrier):
    line = coefficients[..., carrier]
    magnitude = np.abs(coefficients)
    power = (magnitude**2).sum(axis=-1)
    line_amp = magnitude[..., carrier]
    records = {"carrier_amplitude": line_amp, "carrier_phase_deg": np.angle(line, deg=True),
               "carrier_power_fraction": np.divide(line_amp**2, power, out=np.full_like(power, np.nan), where=power > 0),
               "residual_l1": magnitude.sum(axis=-1) - line_amp,
               "carrier_largest": line_amp >= magnitude.max(axis=-1)}
    for i in range(len(coefficients)):
        yield {**{f"{key}{site}": value.item() for key, block in records.items()
                  for site, value in enumerate(block[i])},
               "carrier_min_amplitude": float(line_amp[i].min()),
               "carrier_max_abs_lag": float(np.max(np.abs(np.angle(line[i, LEFT] * line[i, RIGHT].conj(), deg=True))))}


def measure(signal):
    envelope = np.abs(signal)
    minimum = envelope.min(axis=-1)
    unit = np.divide(signal, envelope, out=np.zeros_like(signal), where=envelope > 0)
    relative = unit[..., LEFT] * unit[..., RIGHT].conj()
    weights = np.ones(signal.shape[1])
    weights[[0, -1]] = 0.5
    total = weights.sum()
    amplitudes = np.ptp(signal.real, axis=1) / 2
    unmasked = np.einsum("t,btp->bp", weights, relative) / total
    common = {"amplitudes": amplitudes, "min_half_amplitude": amplitudes.min(axis=1),
              "site_min_envelope": envelope.min(axis=1), "min_envelope": minimum.min(axis=1)}
    records = {}
    for gate in GATES:
        valid = minimum >= gate
        partial = ~np.all(valid, axis=1)
        locking = unmasked.copy()
        coverage = np.ones(len(signal))
        site_coverage = np.ones((len(signal), 4))
        if partial.any():
            masked = valid[partial] * weights
            duration = masked.sum(axis=1)
            numerator = np.einsum("bt,btp->bp", masked, relative[partial])
            locking[partial] = np.divide(numerator, duration[:, None], out=np.full_like(numerator, np.nan),
                                         where=duration[:, None] > 0)
            coverage[partial] = duration / total
            site_coverage[partial] = np.einsum("t,bti->bi", weights, envelope[partial] >= gate) / total
        plv = np.minimum(1, np.abs(locking))
        lag = np.where(plv > 1e-12, np.angle(locking, deg=True), np.nan)
        max_lag, min_plv = np.max(np.abs(lag), axis=1), plv.min(axis=1)
        failures = dict(fail_lag=np.any(np.abs(lag) > MAX_LAG, axis=1), fail_plv=np.any(plv < MIN_PLV, axis=1),
                        fail_coverage=coverage < MIN_COVERAGE, fail_amplitude=common["min_half_amplitude"] < gate,
                        fail_undefined=~np.all(np.isfinite(lag), axis=1))
        margin = np.minimum.reduce((common["min_half_amplitude"] / gate - 1, (coverage - MIN_COVERAGE) / (1 - MIN_COVERAGE),
                                    (min_plv - MIN_PLV) / (1 - MIN_PLV), 1 - max_lag / MAX_LAG))
        records[gate] = {**common, **failures, "lag_deg": lag, "plv": plv, "max_abs_lag": max_lag,
                         "min_plv": min_plv, "coverage": coverage, "site_coverage": site_coverage,
                         "margin": np.where(np.isfinite(margin), margin, -100.0),
                         "passes": ~np.logical_or.reduce(list(failures.values()))}
    return records


def evaluate_window(model, coefficients, bounds, dt, batch):
    kernel = model.prepare(time_grid(bounds, dt)).kernel
    chunks = [measure((coefficients[start:start + batch] @ kernel).swapaxes(-1, -2))
              for start in range(0, len(coefficients), batch)]
    return {gate: {key: np.concatenate([chunk[gate][key] for chunk in chunks]) for key in chunks[0][gate]}
            for gate in GATES}


def threshold_neighborhood(data):
    near = np.zeros(len(data[MAIN_GATE]["passes"]), dtype=bool)
    for gate, row in data.items():
        near |= ((np.abs(row["max_abs_lag"] - MAX_LAG) <= 0.15) | (np.abs(row["min_plv"] - MIN_PLV) <= 0.002)
                 | (np.abs(row["coverage"] - MIN_COVERAGE) <= 0.003)
                 | (np.abs(row["min_half_amplitude"] - gate) <= 0.0001))
    return near


def finite(value):
    return value.item() if np.isfinite(value) else None


def window_rows(name, bounds, cases, data, steps):
    for gate, row in data.items():
        for i, case in enumerate(cases):
            yield dict(case=case["case"], width=case["width"], family=case["family"], window=name, gate=gate,
                       start=bounds[0], end=bounds[1], target_dt=float(steps[i]),
                       **{key: finite(row[key][i]) for key in ("passes", "coverage", "max_abs_lag", "min_plv",
                                                               "min_half_amplitude", "min_envelope", "margin",
                                                               *(f"fail_{k}" for k in FAILURES))},
                       **{f"{k}{j}": finite(v) for k in ("amplitudes", "site_min_envelope", "site_coverage")
                          for j, v in enumerate(row[k][i])},
                       **{f"{k}_{pair}": finite(v) for k in ("lag_deg", "plv") for pair, v in zip(PAIRS, row[k][i])})


def summarize(cases, windows):
    summaries = []
    names = list(windows)
    for gate in GATES:
        by_window = [window[gate] for window in windows.values()]
        def collect(key):
            return np.stack([window[key] for window in by_window], axis=1)

        lag, plv, margin = (collect(key) for key in ("max_abs_lag", "min_plv", "margin"))
        passes = collect("passes").sum(axis=1)
        failures = {key: collect(f"fail_{key}") for key in FAILURES}
        minima = {key: collect(key).min(axis=1) for key in ("min_envelope", "min_half_amplitude", "site_min_envelope", "site_coverage")}
        minima["min_coverage"] = collect("coverage").min(axis=1)
        for i, case in enumerate(cases):
            worst = int(np.argmin(margin[i]))
            summaries.append(dict(**case, gate=gate, pass_windows=int(passes[i]), pass_rate=float(passes[i] / len(names)),
                                  all_pass=bool(passes[i] == len(names)), worst_window=names[worst],
                                  worst_margin=float(margin[i, worst]),
                                  max_abs_lag=finite(np.max(lag[i], where=np.isfinite(lag[i]), initial=-np.inf)),
                                  min_plv=finite(np.min(plv[i], where=np.isfinite(plv[i]), initial=np.inf)),
                                  **{f"fail_{key}_windows": int(values[i].sum()) for key, values in failures.items()},
                                  **{f"first_fail_{key}_window": names[np.argmax(values[i])] if values[i].any() else None
                                     for key, values in failures.items()},
                                  **{key: float(minima[key][i]) for key in ("min_envelope", "min_half_amplitude", "min_coverage")},
                                  **{f"{key}{j}": float(v) for key in ("site_min_envelope", "site_coverage")
                                     for j, v in enumerate(minima[key][i])}))
    return summaries


def aggregate(rows, width):
    for gate, family in product(GATES, ("all", "corner", "uniform") if width else ("all",)):
        subset = [row for row in rows if row["gate"] == gate and (family == "all" or row["family"] == family)]
        pass_count = sum(row["all_pass"] for row in subset)
        yield dict(width=width, gate=gate, family=family, cases=len(subset), all_pass_cases=pass_count,
                   all_pass_rate=pass_count / len(subset), min_pass_windows=min(row["pass_windows"] for row in subset),
                   max_abs_lag=max((row["max_abs_lag"] for row in subset if row["max_abs_lag"] is not None), default=None),
                   min_plv=min((row["min_plv"] for row in subset if row["min_plv"] is not None), default=None),
                   **{key: min(row[key] for row in subset) for key in ("min_envelope", "min_coverage", "min_half_amplitude")},
                   **{f"fail_{key}_cases": sum(row[f"fail_{key}_windows"] > 0 for row in subset) for key in FAILURES})


def write_csv(name, rows):
    with (OUTPUT / name).open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def diagnose(rows, widths):
    def worst(cases):
        return min(cases, key=lambda row: (row["pass_windows"], row["worst_margin"], row["case"]))

    first, onset, manifest = {}, {}, {}
    for gate in GATES:
        gate_rows = [row for row in rows if row["width"] > 0 and row["gate"] == gate]
        failing = [row for row in gate_rows if not row["all_pass"]]
        width = min((row["width"] for row in failing), default=None)
        first[str(gate)] = {"width": width, "previous_tested_width": max((w for w in widths if width and w < width), default=None),
                            "representative": worst([row for row in failing if row["width"] == width]) if width else None,
                            "failing_cases": sum(row["width"] == width for row in failing)}
        onset[str(gate)] = {key: min((row["width"] for row in gate_rows if row[f"fail_{key}_windows"] > 0), default=None)
                            for key in FAILURES}

    def add(alias, row, role, plot_window=None):
        if row is not None and not any(x["case"] == row["case"] and x["gate"] == row["gate"] for x in manifest.values()):
            manifest[alias] = {"offsets_rad": np.array([row[key] for key in AXES]).reshape(4, 2).tolist(),
                               "gate": row["gate"], "plot_window": plot_window or row["worst_window"], "case": row["case"],
                               "width": row["width"], "direction": row["direction"], "role": role}

    add("first-failure", first[str(MAIN_GATE)]["representative"], "first failure at main floor")
    for gate, alias in ((0.025, "first-g0025"), (0.10, "first-g010")):
        candidate = first[str(gate)]["representative"]
        main = first[str(MAIN_GATE)]["representative"]
        if candidate and (not main or candidate["case"] != main["case"]):
            add(alias, candidate, "first failure at alternate floor")
    overall = worst([row for row in rows if row["width"] > 0 and row["gate"] == MAIN_GATE])
    add("worst-g005", overall, "lowest pass count at main floor")
    coverage = [row for row in rows if row["width"] == onset["0.1"]["coverage"] and row["gate"] == 0.10
                and row["fail_coverage_windows"] > 0]
    if coverage:
        weak = min(coverage, key=lambda row: row["min_coverage"])
        add("first-coverage", weak, "lowest coverage at first coverage-failure width, floor 0.10",
            weak["first_fail_coverage_window"])
    return first, onset, overall, manifest


def main():
    nominal, windows, fine_dt, reference, checksum = REFERENCE["sources"]()
    model = Model("open")
    line = int(np.argmin(np.abs(model.frequencies - TRACKED_LINE)))
    assert abs(model.frequencies[line] - TRACKED_LINE) < 1e-9
    labels, directions = samples()
    all_rows, width_rows, tested = [], [], []
    audit = {"refined_case_windows": 0, "refinement_pass_flips": 0, "nominal_max_error": 0.0,
             "direct_max_error": 0.0}
    with gzip.open(OUTPUT / "window_metrics.csv.gz", "wt", encoding="utf-8", newline="", compresslevel=4) as file:
        writer = None
        for width in (0, *WIDTHS):
            if width > 0.30 and all(any(row["gate"] == gate and not row["all_pass"] for row in all_rows) for gate in GATES):
                break
            if width:
                tested.append(width)
            deviations = width * directions if width else np.zeros((1, 4, 2))
            states = product_state(nominal + deviations)
            coefficients = spectrum(model, states)
            direct_time = np.array([0, windows["holdout_01"][0], windows["holdout_24"][1]])
            sample = 0.5 + (coefficients @ model.prepare(direct_time).kernel).swapaxes(-1, -2).real
            direct = np.stack([model.direct(state, direct_time) for state in states])
            audit["direct_max_error"] = max(audit["direct_max_error"], float(np.max(np.abs(direct - sample))))
            meta = []
            for (label, family), offset, spec in zip(labels if width else [("nominal", "nominal")], deviations,
                                                     spectral_records(coefficients, line)):
                case = f"w{width:.3f}-{label}" if width else label
                meta.append(dict(case=case, width=width, family=family, direction=label,
                                 max_abs_error=float(np.abs(offset).max()), **dict(zip(AXES, offset.ravel())), **spec))
            by_window = {}
            for name, bounds in windows.items():
                steps = np.full(len(meta), DT if width else fine_dt)
                data = evaluate_window(model, coefficients, bounds, steps[0], 48)
                if width:
                    near = threshold_neighborhood(data)
                    if near.any():
                        refined = evaluate_window(model, coefficients[near], bounds, fine_dt, 18)
                        audit["refined_case_windows"] += int(near.sum())
                        for gate in GATES:
                            audit["refinement_pass_flips"] += int(np.count_nonzero(data[gate]["passes"][near] != refined[gate]["passes"]))
                            # Common arrays are shared across the gate dictionaries.
                            for key, values in refined[gate].items():
                                data[gate][key][near] = values
                        steps[near] = fine_dt
                else:
                    for gate, values in data.items():
                        prior = reference[name][str(gate)]
                        assert values["passes"][0] == prior["passes"]
                        audit["nominal_max_error"] = max(audit["nominal_max_error"],
                            *(abs(values[key][0] - prior[key]) for key in ("coverage", "min_envelope", "max_abs_lag", "min_plv")))
                by_window[name] = data
                output = window_rows(name, bounds, meta, data, steps)
                if writer is None:
                    row = next(output)
                    writer = csv.DictWriter(file, fieldnames=list(row))
                    writer.writeheader()
                    writer.writerow(row)
                writer.writerows(output)
            rows = summarize(meta, by_window)
            all_rows.extend(rows)
            width_rows.extend(aggregate(rows, width))
            print(f"Finished radius {width:.3f}", flush=True)

    assert audit["direct_max_error"] < 1e-10 and audit["nominal_max_error"] < 1e-9
    write_csv("case_summary.csv", all_rows)
    write_csv("sweep_summary.csv", width_rows)
    first, onsets, worst, manifest = diagnose(all_rows, tested)
    provenance = {"selection_sha256": checksum, "nominal_angles_rad": nominal.tolist(), "windows": windows,
                  "target_dt": DT, "refined_dt": fine_dt, "carrier_frequency": float(model.frequencies[line])}
    (OUTPUT / "diagnostic_cases.json").write_text(json.dumps({**provenance, "cases": manifest}, indent=2, allow_nan=False) + "\n",
                                                  encoding="utf-8")
    (OUTPUT / "sweep.json").write_text(json.dumps({**provenance, "widths": tested, "seed": SEED,
                                                  "directions_per_width": len(directions), "audit": audit,
                                                  "first_failure": first, "onset_by_criterion": onsets,
                                                  "worst_main_gate": worst}, indent=2, allow_nan=False) + "\n",
                                       encoding="utf-8")


if __name__ == "__main__":
    main()
