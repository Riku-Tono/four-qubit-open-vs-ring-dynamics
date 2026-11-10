"""Compare finite-window FFT phases with the energy-gap reference."""

import json
from pathlib import Path

import numpy as np

from run import ANGLES, COUPLINGS, DT, GATE, LEFT, PAIRS, RIGHT, TOPOLOGIES, WINDOWS, prepare, sample, time_grid


OUTPUT = Path(__file__).resolve().parent
EDGE_FRACTION = 0.1


def finite_fft_signal(centered):
    """Keep the known-DC-centered real part and FFT its Hilbert quadrature."""
    count = len(centered)
    multiplier = -1j * np.sign(np.fft.fftfreq(count))
    if count % 2 == 0:
        multiplier[count // 2] = 0  # A real Nyquist bin has zero quadrature.
    quadrature = np.fft.ifft(np.fft.fft(centered, axis=0) * multiplier[:, None], axis=0).real
    return centered + 1j * quadrature


def compare(reference, fft_signal, selection):
    reference, fft_signal = reference[selection], fft_signal[selection]
    reference_valid = np.all(np.abs(reference) >= GATE, axis=1)
    fft_valid = np.all(np.abs(fft_signal) >= GATE, axis=1)
    common = reference_valid & fft_valid
    result = {
        "samples": int(len(reference)),
        "reference_coverage": float(reference_valid.mean()),
        "fft_coverage": float(fft_valid.mean()),
        "common_coverage": float(common.mean()),
        "gate_disagreement_fraction": float(np.mean(reference_valid != fft_valid)),
        "max_analytic_signal_error": float(np.max(np.abs(fft_signal - reference))),
        "pairs": {},
    }
    if not common.any():
        return result

    ref_unit, fft_unit = (np.exp(1j * np.angle(values[common])) for values in (reference, fft_signal))
    ref_pair = ref_unit[:, LEFT] * ref_unit[:, RIGHT].conj()
    fft_pair = fft_unit[:, LEFT] * fft_unit[:, RIGHT].conj()
    error = np.abs(np.angle(fft_pair * ref_pair.conj(), deg=True))
    ref_mean, fft_mean = ref_pair.mean(axis=0), fft_pair.mean(axis=0)
    lag_diff = np.angle(fft_mean * ref_mean.conj(), deg=True)
    for index, pair in enumerate(PAIRS):
        result["pairs"][pair] = {
            "pointwise_lag_mae_deg": float(error[:, index].mean()),
            "pointwise_lag_max_error_deg": float(error[:, index].max()),
            "reference_lag_deg": float(np.angle(ref_mean[index], deg=True)),
            "fft_lag_deg": float(np.angle(fft_mean[index], deg=True)),
            "fft_minus_reference_lag_deg": float(lag_diff[index]),
            "reference_plv": float(abs(ref_mean[index])),
            "fft_plv": float(abs(fft_mean[index])),
            "fft_minus_reference_plv": float(abs(fft_mean[index]) - abs(ref_mean[index])),
        }
    return result


def main():
    results = {
        "phase_gate": GATE,
        "dt_target": DT,
        "edge_fraction_per_end": EDGE_FRACTION,
        "sampling": "left-inclusive uniform grid, right endpoint excluded",
        "reference": "positive Hamiltonian energy gaps",
        "fft": "periodic Hilbert quadrature with no padding and known spectral DC",
        "phase_comparison": "uniformly weighted common all-site envelope gate",
        "cases": {},
    }
    for topology, link, _ in TOPOLOGIES:
        by_j = {}
        for coupling in COUPLINGS:
            prepared = prepare(coupling, link)
            by_window = {}
            for name, bounds in WINDOWS.items():
                time = time_grid(bounds, DT)[:-1]
                position = (time - bounds[0]) / (bounds[1] - bounds[0])
                center = (position >= EDGE_FRACTION) & (position < 1 - EDGE_FRACTION)
                regions = {"whole": np.ones(len(time), dtype=bool), "edges": ~center, "center": center}
                _, oscillation, _ = sample(prepared, time, angles=ANGLES)
                by_angle = {}
                for angle, reference in zip(ANGLES, oscillation):
                    fft_signal = finite_fft_signal(reference.real)
                    by_angle[str(angle)] = {
                        region: compare(reference, fft_signal, mask) for region, mask in regions.items()
                    }
                by_window[name] = {"window": bounds, "dt": float(time[1] - time[0]), "by_angle": by_angle}
            by_j[str(coupling)] = by_window
        results["cases"][topology] = by_j

    # The uncoupled long window is an integer number of periods for the FFT grid.
    for topology, _, _ in TOPOLOGIES:
        periodic = results["cases"][topology]["0.0"]["long"]["by_angle"]
        assert max(case["whole"]["max_analytic_signal_error"] for case in periodic.values()) < 1e-10
    (OUTPUT / "endpoint-audit.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
