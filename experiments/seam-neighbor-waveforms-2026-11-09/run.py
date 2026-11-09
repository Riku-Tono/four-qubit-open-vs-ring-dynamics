"""Resolve seam-neighbor correlation differences into branch and cross terms."""

import json
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
MODEL = run_path(str(OUTPUT.parent / "endpoint-field-six-qubit-2026-11-09" / "run.py"))
ARCHIVE = OUTPUT.parent / "seam-pair-time-window-2026-11-09" / "results.json"
END, DT = MODEL["END"], MODEL["DT"]
LINK, REFERENCE = MODEL["LINK"], MODEL["REFERENCE"]
PAIRS = {4: ("02", "13"), 6: ("04", "15")}
CHANNELS = ("XX", "YY", "ZZ")
BRANCHES = ("linked", "open")
TERMS = ("linked_squared", "open_squared", "cross", "difference_squared")


def waveforms(time, size):
    chain = MODEL["Chain"](size)
    indices = [chain.pairs.index(tuple(map(int, pair))) for pair in PAIRS[size]]
    waves = []
    for link in (LINK, REFERENCE):
        states = chain.evolve(time, link, 0)[:, :1]
        assert np.max(np.abs(np.linalg.norm(states, axis=-1) - 1)) < 1e-12
        waves.append(chain.connected(states)[:, 0].transpose(0, 2, 1)[:, indices])
    return waves


def integrate(time, waves):
    linked, opened = waves
    # Last axis stores the two squared terms, signed cross term, and direct integral.
    powers = np.stack((linked**2, opened**2, -2 * linked * opened, (linked - opened)**2), axis=-1)
    return np.trapezoid(powers, x=time, axis=0)


def pack_integrals(values):
    row = dict(zip(TERMS, map(float, values)))
    linked, opened, cross, _ = values
    row["overlap_cosine"] = float(-cross / (2 * np.sqrt(linked * opened)))
    return row


def main():
    fine_dt = DT / 2
    fine_time = np.arange(round(END / fine_dt) + 1) * fine_dt
    time = fine_time[::2]
    # Keep all samples in the CSV, and a thinner grid plus the exact endpoint in JSON for display.
    preview = np.unique(np.r_[np.arange(0, len(time), 16), len(time) - 1])
    display = {"time": time[preview].round(3).tolist()}
    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    baseline = next(row for row in archive["checkpoints"] if row["T"] == END)
    columns, samples = ["t"], [time]
    integrals, refined, totals = {}, {}, {}
    closure_error = archive_error = grid_error = 0.0

    for size, pairs in PAIRS.items():
        key = f"N{size}"
        fine_waves = waveforms(fine_time, size)
        waves = [wave[::2] for wave in fine_waves]
        primary, fine = integrate(time, waves), integrate(fine_time, fine_waves)
        for terms in (primary, fine):
            closure_error = max(closure_error, float(np.max(np.abs(terms[..., :3].sum(axis=-1) - terms[..., 3]))))
        grid_error = max(grid_error, float(np.max(np.abs(primary - fine))))
        archive_error = max(archive_error, max(float(np.max(np.abs(primary[i, :, 3] - baseline[f"{key}_pair_channels"][pair])))
                                               for i, pair in enumerate(pairs)))
        integrals[key], refined[key], display[key] = {}, {}, {}
        totals[key] = {"both_pairs": pack_integrals(primary.sum(axis=(0, 1)))}

        for i, pair in enumerate(pairs):
            integrals[key][pair], refined[key][pair], display[key][pair] = {}, {}, {}
            totals[key][pair] = pack_integrals(primary[i].sum(axis=0))
            for j, channel in enumerate(CHANNELS):
                integrals[key][pair][channel] = pack_integrals(primary[i, j])
                refined[key][pair][channel] = pack_integrals(fine[i, j])
                display[key][pair][channel] = {}
                for branch, wave in zip(BRANCHES, waves):
                    column = wave[:, i, j]
                    columns.append(f"{key}_{pair}_{channel}_{branch}")
                    samples.append(column)
                    display[key][pair][channel][branch] = column[preview].tolist()

    assert max(archive_error, closure_error) < 1e-12
    assert grid_error < 1e-6
    audit = {"max_decomposition_closure_error": closure_error, "max_archived_pair_channel_error": archive_error,
             "max_integral_term_change_dt_halving": grid_error}
    results = {"Omega": MODEL["OMEGA"], "J": MODEL["J"], "lambda": LINK, "reference_lambda": REFERENCE,
               "epsilon": 0, "T": END, "dt": DT, "audit_dt": fine_dt,
               "initial": {"N4": "1000", "N6": "100000"}, "pair_mapping": list(zip(PAIRS[4], PAIRS[6])),
               "integrals": integrals, "refined_integrals": refined, "totals": totals, "audit": audit,
               "waveforms": display}
    np.savetxt(OUTPUT / "waveforms.csv", np.column_stack(samples), delimiter=",",
               fmt=["%.3f"] + ["%.16g"] * (len(columns) - 1), header=",".join(columns), comments="")
    (OUTPUT / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"totals": totals, "audit": audit}, indent=2))


if __name__ == "__main__":
    main()
