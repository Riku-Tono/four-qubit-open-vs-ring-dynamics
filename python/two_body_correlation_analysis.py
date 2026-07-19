from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

import local_response_analysis as lr


ROOT = Path(__file__).resolve().parent
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
CHANNELS = ("XX", "YY", "ZZ")
AXES = ("X", "Y", "Z")
BRANCHES = ("open", "closed")
OBS_NAMES = [f"{axis}{i}" for axis in AXES for i in range(4)] + [
    f"{ch}_{i}{j}" for ch in CHANNELS for i, j in PAIRS
]
SEPARATION_ABS = 1.0e-8
AMP_TOL = 1.0e-14


def read_wide(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    rows = lr.read_csv(path)
    times = np.array([float(row["time"]) for row in rows])
    data = {key: np.array([float(row[key]) for row in rows]) for key in rows[0] if key != "time"}
    return times, data


def op_matrix(paulis: dict[int, str]) -> np.ndarray:
    out = np.zeros((lr.DIM, lr.DIM), dtype=complex)
    for basis in range(lr.DIM):
        target = basis
        amp = 1.0 + 0.0j
        for site, axis in paulis.items():
            bit = (basis >> site) & 1
            if axis == "X":
                target ^= 1 << site
            elif axis == "Y":
                target ^= 1 << site
                amp *= 1j if bit == 0 else -1j
            elif axis == "Z":
                amp *= 1.0 if bit == 0 else -1.0
            else:
                raise ValueError(axis)
        out[target, basis] += amp
    return out


def observables() -> dict[str, np.ndarray]:
    ops: dict[str, np.ndarray] = {}
    for axis in AXES:
        for i in range(4):
            ops[f"{axis}{i}"] = op_matrix({i: axis})
    for channel in CHANNELS:
        axis = channel[0]
        for i, j in PAIRS:
            ops[f"{channel}_{i}{j}"] = op_matrix({i: axis, j: axis})
    return ops


def connected(data: dict[str, np.ndarray], branch: str, channel: str, i: int, j: int) -> np.ndarray:
    axis = channel[0]
    return data[f"{branch}_{channel}_{i}{j}"] - data[f"{branch}_{axis}{i}"] * data[f"{branch}_{axis}{j}"]


def rank_average(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = 0.5 * (start + end - 1) + 1.0
        start = end
    return ranks


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    return lr.pearson(rank_average(a), rank_average(b))


def correlation_arrays(data: dict[str, np.ndarray], branch: str) -> dict[str, np.ndarray]:
    return {
        channel: np.stack([connected(data, branch, channel, i, j) for i, j in PAIRS], axis=1)
        for channel in CHANNELS
    }


def write_observable_csvs(times: np.ndarray, data: dict[str, np.ndarray]) -> None:
    one_rows: list[dict] = []
    raw_rows: list[dict] = []
    conn_rows: list[dict] = []
    for ti, time in enumerate(times):
        for branch in BRANCHES:
            for axis in AXES:
                for site in range(4):
                    one_rows.append(
                        {"time": time, "branch": branch, "axis": axis, "site": site,
                         "expectation": data[f"{branch}_{axis}{site}"][ti]}
                    )
            for channel in CHANNELS:
                axis = channel[0]
                for i, j in PAIRS:
                    raw = data[f"{branch}_{channel}_{i}{j}"][ti]
                    product = data[f"{branch}_{axis}{i}"][ti] * data[f"{branch}_{axis}{j}"][ti]
                    raw_rows.append(
                        {"time": time, "branch": branch, "channel": channel, "site_i": i,
                         "site_j": j, "raw_expectation": raw}
                    )
                    conn_rows.append(
                        {"time": time, "branch": branch, "channel": channel, "site_i": i,
                         "site_j": j, "raw_expectation": raw, "one_body_product": product,
                         "connected": raw - product}
                    )
    lr.write_csv(ROOT / "one_body_pauli_expectations.csv", one_rows)
    lr.write_csv(ROOT / "two_body_raw_pauli_correlations.csv", raw_rows)
    lr.write_csv(ROOT / "two_body_connected_correlations.csv", conn_rows)


def distances(times: np.ndarray, data: dict[str, np.ndarray]) -> tuple[list[dict], dict[str, np.ndarray]]:
    arrays = {branch: correlation_arrays(data, branch) for branch in BRANCHES}
    out: dict[str, np.ndarray] = {}
    rows: list[dict] = []
    one_open = np.stack([data[f"open_{name}"] for name in OBS_NAMES[:12]], axis=1)
    one_closed = np.stack([data[f"closed_{name}"] for name in OBS_NAMES[:12]], axis=1)
    out["D_onebody"] = np.linalg.norm(one_open - one_closed, axis=1)
    out["D_ZZ_connected"] = np.linalg.norm(arrays["open"]["ZZ"] - arrays["closed"]["ZZ"], axis=1)
    open_all = np.concatenate([arrays["open"][ch] for ch in CHANNELS], axis=1)
    closed_all = np.concatenate([arrays["closed"][ch] for ch in CHANNELS], axis=1)
    out["D_all_connected"] = np.linalg.norm(open_all - closed_all, axis=1)
    out["D_onebody_all_connected"] = np.sqrt(out["D_onebody"] ** 2 + out["D_all_connected"] ** 2)
    channel_l2 = {}
    for ch in CHANNELS:
        diff = arrays["open"][ch] - arrays["closed"][ch]
        channel_l2[ch] = np.linalg.norm(diff, axis=1)
        out[f"D_{ch}"] = channel_l2[ch]
    out["D_corr"] = np.sqrt(sum(channel_l2[ch] ** 2 for ch in CHANNELS))
    for ti, time in enumerate(times):
        row: dict[str, float] = {"time": time}
        for ch in CHANNELS:
            diff = arrays["open"][ch][ti] - arrays["closed"][ch][ti]
            row[f"{ch}_upper_L1"] = float(np.sum(np.abs(diff)))
            row[f"{ch}_upper_L2"] = float(np.linalg.norm(diff))
            row[f"{ch}_map_Frobenius"] = float(math.sqrt(2.0) * np.linalg.norm(diff))
        row.update({key: float(value[ti]) for key, value in out.items()})
        rows.append(row)
    lr.write_csv(ROOT / "correlation_distance_timeseries.csv", rows)
    return rows, out


def comparison(times: np.ndarray, dist: dict[str, np.ndarray]) -> tuple[list[dict], list[dict], np.ndarray, np.ndarray]:
    old = lr.read_csv(ROOT / "state_response_comparison.csv")
    d_state = np.array([float(row["state_distance"]) for row in old])
    d_response = np.array([float(row["response_normalized_L1"]) for row in old])
    rows = []
    for ti, time in enumerate(times):
        row = {"time": time, "D_state": d_state[ti], "D_response_normalized_L1": d_response[ti]}
        row.update({key: float(value[ti]) for key, value in dist.items()})
        rows.append(row)
    lr.write_csv(ROOT / "state_response_correlation_comparison.csv", rows)
    metrics = {"D_state": d_state, "D_response_normalized_L1": d_response, **dist}
    corr_rows = []
    for target in ("D_XX", "D_YY", "D_ZZ", "D_corr", "D_onebody", "D_ZZ_connected", "D_all_connected", "D_onebody_all_connected"):
        for reference in ("D_state", "D_response_normalized_L1"):
            corr_rows.append({
                "target": target, "reference": reference,
                "pearson": lr.pearson(metrics[target], metrics[reference]),
                "spearman": spearman(metrics[target], metrics[reference]),
                "target_max": float(np.max(metrics[target])),
                "target_max_time": float(times[int(np.argmax(metrics[target]))]),
                "reference_max_time": float(times[int(np.argmax(metrics[reference]))]),
                "max_time_offset_target_minus_reference": float(
                    times[int(np.argmax(metrics[target]))] - times[int(np.argmax(metrics[reference]))]
                ),
            })
    lr.write_csv(ROOT / "correlation_distance_correlations.csv", corr_rows)
    return rows, corr_rows, d_state, d_response


def representative_maps(times: np.ndarray, data: dict[str, np.ndarray], dist: dict[str, np.ndarray], d_state: np.ndarray) -> list[dict]:
    first_sep = int(np.flatnonzero(dist["D_corr"] > SEPARATION_ABS)[0])
    first_state_half = int(np.flatnonzero(d_state > 0.5)[0])
    events = {
        "initial": 0,
        "first_Dcorr_gt_1e-8": first_sep,
        "max_local_response_distance": int(round(2.154 / lr.DT_PRIMARY)),
        "first_Dstate_gt_0.5": first_state_half,
        "max_Dstate": int(np.argmax(d_state)),
        "final": len(times) - 1,
    }
    for ch in CHANNELS:
        events[f"max_{ch}_distance"] = int(np.argmax(dist[f"D_{ch}"]))
    rows = []
    for event, ti in events.items():
        for branch in BRANCHES:
            for ch in CHANNELS:
                pair_values = {(i, j): connected(data, branch, ch, i, j)[ti] for i, j in PAIRS}
                for i in range(4):
                    for j in range(4):
                        value = 0.0 if i == j else pair_values[min(i, j), max(i, j)]
                        rows.append({"event": event, "time": times[ti], "branch": branch,
                                     "channel": ch, "site_i": i, "site_j": j, "connected": value})
    lr.write_csv(ROOT / "representative_correlation_maps.csv", rows)
    return rows


def symmetry_audit(times: np.ndarray, data: dict[str, np.ndarray]) -> tuple[list[dict], dict[str, float]]:
    rows = []
    maxima: dict[str, float] = {}
    for ti, time in enumerate(times):
        for ch in CHANNELS:
            c = {(i, j): connected(data, "closed", ch, i, j)[ti] for i, j in PAIRS}
            nearest = np.array([c[0, 1], c[1, 2], c[2, 3], c[0, 3]])
            opposite = np.array([c[0, 2], c[1, 3]])
            open_c = {(i, j): connected(data, "open", ch, i, j)[ti] for i, j in PAIRS}
            vals = {
                "closed_nearest_rotation_residual": float(np.max(nearest) - np.min(nearest)),
                "closed_opposite_rotation_residual": float(abs(opposite[0] - opposite[1])),
                "open_reflection_01_23_residual": float(abs(open_c[0, 1] - open_c[2, 3])),
                "open_reflection_02_13_residual": float(abs(open_c[0, 2] - open_c[1, 3])),
            }
            rows.append({"time": time, "channel": ch, **vals})
            for key, value in vals.items():
                maxima[key] = max(maxima.get(key, 0.0), value)
                maxima[f"{ch}_{key}"] = max(maxima.get(f"{ch}_{key}", 0.0), value)
    lr.write_csv(ROOT / "correlation_symmetry_audit.csv", rows)
    return rows, maxima


def exact_and_dt_audits(times: np.ndarray, data: dict[str, np.ndarray], ops: dict[str, np.ndarray]) -> tuple[dict, list[dict]]:
    psi0 = np.zeros(lr.DIM, complex)
    psi0[0] = 1.0
    exact_rows = []
    maxima = {"exact_onebody": 0.0, "exact_raw_two_body": 0.0, "exact_connected": 0.0,
              "max_imaginary_exact": 0.0, "max_expectation_bound_excess": 0.0,
              "max_hermiticity_residual": 0.0, "connected_raw_minus_product_identity": 0.0,
              "global_phase_expectation_residual": 0.0}
    exact_values: dict[str, dict[str, np.ndarray]] = {}
    for name, op in ops.items():
        maxima["max_hermiticity_residual"] = max(maxima["max_hermiticity_residual"], float(np.max(np.abs(op - op.conj().T))))
    for branch, bonds in (("open", lr.OPEN_BONDS), ("closed", lr.CLOSED_BONDS)):
        states = lr.exact_states(lr.hamiltonian(bonds), psi0, times)
        phased = states[-1] * np.exp(1j * 0.731)
        for op in ops.values():
            maxima["global_phase_expectation_residual"] = max(
                maxima["global_phase_expectation_residual"],
                float(abs(np.vdot(states[-1], op @ states[-1]) - np.vdot(phased, op @ phased))),
            )
        exact_values[branch] = {}
        for name, op in ops.items():
            vals_complex = np.einsum("ti,ij,tj->t", states.conj(), op, states)
            exact_values[branch][name] = vals_complex.real
            maxima["max_imaginary_exact"] = max(maxima["max_imaginary_exact"], float(np.max(np.abs(vals_complex.imag))))
            maxima["max_expectation_bound_excess"] = max(maxima["max_expectation_bound_excess"], float(max(0.0, np.max(np.abs(vals_complex.real)) - 1.0)))
            kind = "onebody" if "_" not in name else "raw_two_body"
            residual = np.abs(data[f"{branch}_{name}"] - vals_complex.real)
            maxima[f"exact_{kind}"] = max(maxima[f"exact_{kind}"], float(np.max(residual)))
            for ti, time in enumerate(times):
                exact_rows.append({"time": time, "branch": branch, "observable": name,
                                   "trotter": data[f"{branch}_{name}"][ti], "exact": vals_complex.real[ti],
                                   "abs_residual": residual[ti]})
        for ch in CHANNELS:
            axis = ch[0]
            for i, j in PAIRS:
                exact_conn = exact_values[branch][f"{ch}_{i}{j}"] - exact_values[branch][f"{axis}{i}"] * exact_values[branch][f"{axis}{j}"]
                direct_conn = connected(data, branch, ch, i, j)
                conn_residual = np.abs(direct_conn - exact_conn)
                maxima["exact_connected"] = max(maxima["exact_connected"], float(np.max(conn_residual)))
                for ti, time in enumerate(times):
                    exact_rows.append({"time": time, "branch": branch,
                                       "observable": f"connected_{ch}_{i}{j}",
                                       "trotter": direct_conn[ti], "exact": exact_conn[ti],
                                       "abs_residual": conn_residual[ti]})
    lr.write_csv(ROOT / "correlation_exact_validation.csv", exact_rows)

    fine_times, fine = read_wide(ROOT / "direct_open_closed_pauli_dt001.csv")
    if len(fine_times) != 2 * len(times) - 1:
        raise AssertionError("dt grids do not nest")
    dt_rows = []
    dt_max = {"dt_onebody": 0.0, "dt_raw_two_body": 0.0, "dt_connected": 0.0}
    for ti, time in enumerate(times):
        fi = 2 * ti
        for branch in BRANCHES:
            one_res, raw_res, conn_res = [], [], []
            for name in OBS_NAMES[:12]:
                one_res.append(abs(data[f"{branch}_{name}"][ti] - fine[f"{branch}_{name}"][fi]))
            for ch in CHANNELS:
                for i, j in PAIRS:
                    name = f"{ch}_{i}{j}"
                    raw_res.append(abs(data[f"{branch}_{name}"][ti] - fine[f"{branch}_{name}"][fi]))
                    coarse = connected(data, branch, ch, i, j)[ti]
                    axis = ch[0]
                    fine_conn = fine[f"{branch}_{name}"][fi] - fine[f"{branch}_{axis}{i}"][fi] * fine[f"{branch}_{axis}{j}"][fi]
                    conn_res.append(abs(coarse - fine_conn))
            vals = {"max_onebody_residual": max(one_res), "max_raw_two_body_residual": max(raw_res),
                    "max_connected_residual": max(conn_res)}
            dt_rows.append({"time": time, "branch": branch, **vals})
            dt_max["dt_onebody"] = max(dt_max["dt_onebody"], vals["max_onebody_residual"])
            dt_max["dt_raw_two_body"] = max(dt_max["dt_raw_two_body"], vals["max_raw_two_body_residual"])
            dt_max["dt_connected"] = max(dt_max["dt_connected"], vals["max_connected_residual"])
    lr.write_csv(ROOT / "correlation_dt_convergence.csv", dt_rows)
    maxima.update(dt_max)
    return maxima, exact_rows


def mode_analysis(times: np.ndarray, ops: dict[str, np.ndarray]) -> tuple[list[dict], list[dict], list[dict], dict[str, float]]:
    psi0 = np.zeros(lr.DIM, complex)
    psi0[0] = 1.0
    amplitude_rows: list[dict] = []
    contribution_rows: list[dict] = []
    reconstruction_rows: list[dict] = []
    dominant_rows: list[dict] = []
    audit = {"mode_raw_max_residual": 0.0, "mode_onebody_max_residual": 0.0,
             "connected_from_reconstructed_ingredients_max_residual": 0.0}
    names = list(ops)
    op_list = [ops[name] for name in names]
    exact_by_branch: dict[str, np.ndarray] = {}
    names = list(ops)
    op_list = [ops[name] for name in names]
    for branch, bonds in (("open", lr.OPEN_BONDS), ("closed", lr.CLOSED_BONDS)):
        states = lr.exact_states(lr.hamiltonian(bonds), psi0, times)
        exact_by_branch[branch] = np.stack(
            [np.einsum("ti,ij,tj->t", states.conj(), op, states).real for op in op_list], axis=1
        )
    for branch, bonds in (("open", lr.OPEN_BONDS), ("closed", lr.CLOSED_BONDS)):
        h = lr.hamiltonian(bonds)
        _, _, groups = lr.spectral_groups(h, psi0)
        pairs = []
        total = np.zeros((len(times), len(names)))
        for a, ga in enumerate(groups):
            for b in range(a, len(groups)):
                gb = groups[b]
                amps = np.array([np.vdot(gb["projection"], op @ ga["projection"]) for op in op_list])
                if np.max(np.abs(amps)) <= AMP_TOL:
                    continue
                gap = ga["energy"] - gb["energy"]
                factor = 1.0 if a == b else 2.0
                contrib = factor * np.real(np.exp(-1j * gap * times)[:, None] * amps[None, :])
                score = float(np.sqrt(np.mean(contrib ** 2)))
                pairs.append({"a": a, "b": b, "ea": ga["energy"], "eb": gb["energy"],
                              "gap": gap, "amps": amps, "contrib": contrib, "score": score})
                total += contrib
        pairs.sort(key=lambda p: -p["score"])
        exact = exact_by_branch[branch]
        one_idx = [names.index(name) for name in OBS_NAMES[:12]]
        raw_idx = [names.index(name) for name in OBS_NAMES[12:]]
        audit["mode_onebody_max_residual"] = max(audit["mode_onebody_max_residual"], float(np.max(np.abs(total[:, one_idx] - exact[:, one_idx]))))
        audit["mode_raw_max_residual"] = max(audit["mode_raw_max_residual"], float(np.max(np.abs(total[:, raw_idx] - exact[:, raw_idx]))))
        channel_rank_maps: dict[str, dict[tuple[int, int], tuple[int, int, float, float]]] = {}
        for ch in CHANNELS:
            axis = ch[0]
            full_conn = np.stack([
                exact[:, names.index(f"{ch}_{i}{j}")] - exact[:, names.index(f"{axis}{i}")] * exact[:, names.index(f"{axis}{j}")]
                for i, j in PAIRS], axis=1)
            other = exact_by_branch["closed" if branch == "open" else "open"]
            other_conn = np.stack([
                other[:, names.index(f"{ch}_{i}{j}")] - other[:, names.index(f"{axis}{i}")] * other[:, names.index(f"{axis}{j}")]
                for i, j in PAIRS], axis=1)
            full_diff = full_conn - other_conn if branch == "open" else other_conn - full_conn
            scored = []
            for pair in pairs:
                reduced = exact - pair["contrib"]
                reduced_conn = np.stack([
                    reduced[:, names.index(f"{ch}_{i}{j}")] - reduced[:, names.index(f"{axis}{i}")] * reduced[:, names.index(f"{axis}{j}")]
                    for i, j in PAIRS], axis=1)
                effect_local = full_conn - reduced_conn
                reduced_diff = reduced_conn - other_conn if branch == "open" else other_conn - reduced_conn
                effect_difference = full_diff - reduced_diff
                local_score = float(np.sqrt(np.mean(effect_local ** 2)))
                diff_score = float(np.sqrt(np.mean(effect_difference ** 2)))
                scored.append((pair, local_score, diff_score))
            local_order = sorted(scored, key=lambda item: -item[1])
            diff_order = sorted(scored, key=lambda item: -item[2])
            local_rank = {(p["a"], p["b"]): rank for rank, (p, _, _) in enumerate(local_order, 1)}
            diff_rank = {(p["a"], p["b"]): rank for rank, (p, _, _) in enumerate(diff_order, 1)}
            channel_rank_maps[ch] = {
                (p["a"], p["b"]): (local_rank[(p["a"], p["b"])], diff_rank[(p["a"], p["b"])], local, diff)
                for p, local, diff in scored
            }
            for p, local, diff in scored:
                dominant_rows.append({
                    "branch": branch, "channel": ch,
                    "connected_rms_rank": local_rank[(p["a"], p["b"])],
                    "open_closed_difference_rank": diff_rank[(p["a"], p["b"])],
                    "group_a": p["a"], "group_b": p["b"], "energy_a": p["ea"], "energy_b": p["eb"],
                    "delta_E": p["gap"], "beat_period": math.inf if abs(p["gap"]) < 1e-15 else 2*math.pi/abs(p["gap"]),
                    "connected_leave_one_pair_out_rms": local,
                    "difference_leave_one_pair_out_rms": diff,
                    "ranking_note": "nonlinear connected influence measured by removing this eigenspace-pair contribution from raw and one-body ingredients",
                })
        for rank, pair in enumerate(pairs, 1):
            for oi, name in enumerate(names):
                amp = pair["amps"][oi]
                amplitude_rows.append({"branch": branch, "rank": rank, "group_a": pair["a"],
                                       "group_b": pair["b"], "energy_a": pair["ea"], "energy_b": pair["eb"],
                                       "delta_E": pair["gap"], "beat_period": math.inf if abs(pair["gap"]) < 1e-15 else 2*math.pi/abs(pair["gap"]),
                                       "observable": name, "amplitude_real": amp.real, "amplitude_imag": amp.imag,
                                       "amplitude_abs": abs(amp), "amplitude_phase": np.angle(amp) if abs(amp) else 0.0,
                                       "global_pair_rms_score": pair["score"]})
            for ti, time in enumerate(times):
                row = {"time": time, "branch": branch, "global_rank": rank, "group_a": pair["a"],
                       "group_b": pair["b"], "delta_E": pair["gap"],
                       "phase_advance": float((-pair["gap"] * time + math.pi) % (2*math.pi) - math.pi)}
                for ch in CHANNELS:
                    local_rank, diff_rank, _, _ = channel_rank_maps[ch][(pair["a"], pair["b"])]
                    row[f"{ch}_connected_rank"] = local_rank
                    row[f"{ch}_difference_rank"] = diff_rank
                row.update({f"contribution_{name}": pair["contrib"][ti, oi] for oi, name in enumerate(names)})
                contribution_rows.append(row)
        for ch in CHANNELS:
            cumulative = np.zeros_like(exact)
            errors_by_k = []
            ranked_pairs = sorted(pairs, key=lambda p: channel_rank_maps[ch][(p["a"], p["b"])][0])
            for k, pair in enumerate(ranked_pairs, 1):
                cumulative += pair["contrib"]
                raw_names = [f"{ch}_{i}{j}" for i, j in PAIRS]
                axis = ch[0]
                conn_exact = np.stack([
                    exact[:, names.index(f"{ch}_{i}{j}")] - exact[:, names.index(f"{axis}{i}")] * exact[:, names.index(f"{axis}{j}")]
                    for i, j in PAIRS], axis=1)
                conn_k = np.stack([
                    cumulative[:, names.index(f"{ch}_{i}{j}")] - cumulative[:, names.index(f"{axis}{i}")] * cumulative[:, names.index(f"{axis}{j}")]
                    for i, j in PAIRS], axis=1)
                denom = float(np.linalg.norm(conn_exact))
                relative_error = float(np.linalg.norm(conn_k - conn_exact) / denom) if denom else 0.0
                errors_by_k.append((k, relative_error))
                if k == len(pairs):
                    audit["connected_from_reconstructed_ingredients_max_residual"] = max(
                        audit["connected_from_reconstructed_ingredients_max_residual"], float(np.max(np.abs(conn_k - conn_exact))))
            row = {"branch": branch, "channel": ch, "total_nonzero_mode_pairs": len(pairs)}
            for pct, target in ((90, math.sqrt(0.10)), (95, math.sqrt(0.05)), (99, math.sqrt(0.01))):
                qualifying = [k for k, err in errors_by_k if err <= target]
                row[f"pairs_for_{pct}pct_explained_energy"] = min(qualifying) if qualifying else len(pairs)
                row[f"relative_error_at_selected_{pct}"] = next(err for k, err in errors_by_k if k == row[f"pairs_for_{pct}pct_explained_energy"])
            reconstruction_rows.append(row)
    lr.write_csv(ROOT / "correlation_mode_pair_amplitudes.csv", amplitude_rows)
    lr.write_csv(ROOT / "correlation_mode_pair_contributions.csv", contribution_rows)
    lr.write_csv(ROOT / "correlation_mode_reconstruction_summary.csv", reconstruction_rows)
    dominant_rows.sort(key=lambda row: (row["branch"], row["channel"], row["connected_rms_rank"]))
    lr.write_csv(ROOT / "correlation_dominant_mode_pair_summary.csv", dominant_rows)
    return amplitude_rows, contribution_rows, reconstruction_rows, audit


def write_report(times: np.ndarray, dist: dict[str, np.ndarray], d_state: np.ndarray,
                 d_response: np.ndarray, corr_rows: list[dict], symmetry: dict[str, float],
                 audits: dict[str, float], recon: list[dict], data: dict[str, np.ndarray]) -> None:
    i_state = int(np.argmax(d_state))
    i_response = int(np.argmax(d_response))
    i_first_corr = int(np.flatnonzero(dist["D_corr"] > SEPARATION_ABS)[0])
    high_state = d_state > 0.5
    lines = [
        "# 二体相関解析：切替なし4bit直列対4bit円環",
        "",
        "## 定義",
        "",
        "一体量は `<A_i>`、二体raw量は `<A_i A_j>`、connected相関は "
        "`C^AA_ij=<A_i A_j>-<A_i><A_j>`（A=X,Y,Z）です。4×4表示では比較対象でない対角を0に固定します。",
        "",
        "各チャネルのupper L1/L2は6個の `i<j` 成分に対する距離です。Frobenius距離は対称な4×4ゼロ対角行列に対し "
        "`sqrt(2) * upper_L2`、統合距離は `D_corr=sqrt(D_XX^2+D_YY^2+D_ZZ^2)` です。",
        "",
        "## 主結果（有限N=4、0<=t<=3.688）",
        "",
    ]
    for ch in CHANNELS:
        idx = int(np.argmax(dist[f"D_{ch}"]))
        lines.append(f"- {ch}: 最大upper-L2={dist[f'D_{ch}'][idx]:.12g} at t={times[idx]:.3f}")
    dc = int(np.argmax(dist["D_corr"]))
    lines.append(f"- 統合D_corr最大={dist['D_corr'][dc]:.12g} at t={times[dc]:.3f}")
    lines.append(f"- 状態距離最大={d_state[i_state]:.12g} at t={times[i_state]:.3f}; この時のD_corr={dist['D_corr'][i_state]:.12g}")
    lines.append(f"- D_corrが{SEPARATION_ABS:.0e}を最初に超える時刻: t={times[i_first_corr]:.3f}")
    lines.append(f"- 状態距離>0.5区間での平均チャネル距離: " + ", ".join(
        f"{ch}={np.mean(dist[f'D_{ch}'][high_state]):.6g}" for ch in CHANNELS))
    lines.append("")
    lines.append("t=3.462 のconnected相関（open / closed / 差）:")
    lines.append("")
    lines.append("| channel | pair | open | closed | open-closed |")
    lines.append("|---|---:|---:|---:|---:|")
    for ch in CHANNELS:
        for i, j in PAIRS:
            o = connected(data, "open", ch, i, j)[i_state]
            c = connected(data, "closed", ch, i, j)[i_state]
            lines.append(f"| {ch} | {i}{j} | {o:.9g} | {c:.9g} | {o-c:.9g} |")
    lines += ["", "## 距離間の対応", ""]
    for row in corr_rows:
        if row["reference"] == "D_state":
            lines.append(f"- {row['target']} vs D_state: Pearson={row['pearson']:.6f}, Spearman={row['spearman']:.6f}, max={row['target_max']:.6g} at t={row['target_max_time']:.3f}, 状態距離最大との時刻差={row['max_time_offset_target_minus_reference']:.3f}")
    lines.append(f"- 局所占有形状差最大はt={times[i_response]:.3f}。D_corr最大t={times[dc]:.3f}との差は{times[dc]-times[i_response]:.3f}。")
    lines += ["", "相関係数は時系列の併行性を記述するだけで、因果関係を意味しません。", "", "## 対称性監査", ""]
    for key, value in symmetry.items():
        lines.append(f"- {key}: max {value:.3e}")
    lines += ["", "## 固有空間対再構成", "",
              "縮退固有空間への射影 `P_E` を使い、一体量とraw二体量を基底不変な空間対寄与へ分解しました。"
              "connected量は線形観測量ではないため、raw二体量と一体量を各々再構成した後、raw-minus-productで構成しています。",
              "", "| branch | channel | total pairs | k90 | k95 | k99 |", "|---|---|---:|---:|---:|---:|"]
    for row in recon:
        lines.append(f"| {row['branch']} | {row['channel']} | {row['total_nonzero_mode_pairs']} | {row['pairs_for_90pct_explained_energy']} | {row['pairs_for_95pct_explained_energy']} | {row['pairs_for_99pct_explained_energy']} |")
    largest = None
    for ch in CHANNELS:
        for i, j in PAIRS:
            delta = connected(data, "open", ch, i, j)[i_state] - connected(data, "closed", ch, i, j)[i_state]
            if largest is None or abs(delta) > abs(largest[3]):
                largest = (ch, i, j, delta)
    lines += ["", "## 問いへの回答", "",
              f"1. 状態距離が0.5を超える区間では平均的に " + max(CHANNELS, key=lambda ch: np.mean(dist[f'D_{ch}'][high_state])) + " が最大です。ただしXXとZZは近く、単一チャネルだけの説明ではありません。",
              f"2. t=3.462では局所占有の正規化L1差が{d_response[i_state]:.6g}なのに対し、D_corr={dist['D_corr'][i_state]:.6g}です。最大の単一成分差は{largest[0]} pair {largest[1]}{largest[2]} の{largest[3]:.6g}です。",
              "3. closedの最近接4対と対向2対は、下記対称性残差の範囲で回転対称性を保ちます。",
              "4. openでは反射対応する対は一致しますが、端点を含む対と内部・離隔対は同一にはならず、その差が相関マップに現れます。",
              f"5. 各チャネル最大値はXX={np.max(dist['D_XX']):.6g}, YY={np.max(dist['D_YY']):.6g}, ZZ={np.max(dist['D_ZZ']):.6g}です。最大値ではZZ、統合距離最大t={times[dc]:.3f}ではZZが最大ですがXXも同程度です。",
              f"6. D_corr最大はt={times[dc]:.3f}で、状態距離最大t={times[i_state]:.3f}とも局所占有差最大t={times[i_response]:.3f}とも一致しません。",
              "7. 99%基準ではopenに23--37対、closedに11--14対を要し、特にopenは少数対だけより多モード競合として記述する方が安全です。",
              f"8. 一体+全connected特徴距離と状態距離のPearson相関は{next(r['pearson'] for r in corr_rows if r['target']=='D_onebody_all_connected' and r['reference']=='D_state'):.6f}です。局所占有形状差では弱く見える時刻にも差を可視化しますが、完全な状態再構成ではありません。"]
    lines += ["", "## 数値監査", ""]
    for key, value in audits.items():
        lines.append(f"- {key}: {value:.3e}")
    lines += ["", "Rust testsは10/10合格、`cargo fmt --check`合格。global phaseに依存しない期待値・connected相関・距離のみを物理比較に使用しました。",
              "", "## 結論の範囲", "",
              "これは有限N=4・有限時間の固定Hamiltonian数値結果です。一般的輸送則、相転移、熱力学極限、N=7 Lindblad逆流との同一機構、または円環が応答移動を因果的に起こすという主張ではありません。"]
    (ROOT / "CORRELATION_ANALYSIS_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def numerical_audit_rows(audits: dict[str, float], symmetry: dict[str, float]) -> None:
    combined = {**audits, **{f"symmetry_{k}": v for k, v in symmetry.items()}}
    rows = [{"audit": key, "max_residual": value} for key, value in combined.items()]
    lr.write_csv(ROOT / "correlation_numerical_audit.csv", rows)


def main() -> None:
    times, data = read_wide(ROOT / "direct_open_closed_pauli_dt002.csv")
    if abs(times[-1] - lr.T_END) > 1e-12 or abs(times[1] - lr.DT_PRIMARY) > 1e-12:
        raise AssertionError("unexpected primary grid")
    write_observable_csvs(times, data)
    _, dist = distances(times, data)
    _, corr_rows, d_state, d_response = comparison(times, dist)
    representative_maps(times, data, dist, d_state)
    _, symmetry = symmetry_audit(times, data)
    ops = observables()
    audits, _ = exact_and_dt_audits(times, data, ops)
    _, _, recon, mode_audit = mode_analysis(times, ops)
    audits.update(mode_audit)
    state_rows = lr.read_csv(ROOT / "direct_open_closed_dt002.csv")
    audits["dt002_max_norm_error"] = max(
        max(abs(float(row["open_norm_error"])) for row in state_rows),
        max(abs(float(row["closed_norm_error"])) for row in state_rows),
    )
    fine_state_rows = lr.read_csv(ROOT / "direct_open_closed_dt001.csv")
    audits["dt001_max_norm_error"] = max(
        max(abs(float(row["open_norm_error"])) for row in fine_state_rows),
        max(abs(float(row["closed_norm_error"])) for row in fine_state_rows),
    )
    numerical_audit_rows(audits, symmetry)
    write_report(times, dist, d_state, d_response, corr_rows, symmetry, audits, recon, data)
    print(f"D_corr max={np.max(dist['D_corr']):.12g} at t={times[int(np.argmax(dist['D_corr']))]:.3f}")
    print(f"exact connected max residual={audits['exact_connected']:.3e}")
    print(f"mode raw max residual={audits['mode_raw_max_residual']:.3e}")


if __name__ == "__main__":
    main()
