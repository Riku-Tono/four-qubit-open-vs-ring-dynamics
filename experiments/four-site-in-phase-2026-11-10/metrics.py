"""Occupation waveforms and common phase metrics for the fixed N=4 model."""

from itertools import combinations
from pathlib import Path
from runpy import run_path

import numpy as np


OUTPUT = Path(__file__).resolve().parent
REFERENCE = run_path(str(OUTPUT.parent / "endpoint-field-six-qubit-2026-11-09" / "run.py"))
TOPOLOGIES = {"open": 0, "ring": 1}
WINDOWS = {"short": (0, 3.688), "full": (0, 8 * np.pi), "holdout": (3.688, 8 * np.pi),
           "late": (4 * np.pi, 8 * np.pi), "extended": (8 * np.pi, 64 * np.pi)}
GATES = (0.025, 0.05, 0.10)
MAIN_GATE = 0.05
MIN_COVERAGE, MAX_LAG, MIN_PLV = 0.90, 15, 0.95
LEFT, RIGHT = np.asarray(list(combinations(range(4), 2))).T
PAIRS = tuple(f"{i}{j}" for i, j in zip(LEFT, RIGHT))
OCCUPATION = (np.arange(16)[:, None] >> np.arange(3, -1, -1)) & 1
ZERO_GAP = 1e-10


class Model:
    def __init__(self, topology):
        chain = REFERENCE["Chain"](4)
        self.hamiltonian = chain.hamiltonian(TOPOLOGIES[topology], 0)
        self.energy, self.vectors = np.linalg.eigh(self.hamiltonian)
        gaps = self.energy[:, None] - self.energy[None, :]
        observed = np.einsum("sa,si,sb->abi", self.vectors, OCCUPATION, self.vectors)
        self.dc = observed * (np.abs(gaps) <= ZERO_GAP)[..., None]
        a, b = np.where((gaps > ZERO_GAP) & (np.max(np.abs(observed), axis=-1) > 1e-13))
        order = np.argsort(gaps[a, b])
        self.a, self.b = a[order], b[order]
        frequencies = gaps[self.a, self.b]
        self.group_starts = np.r_[0, np.flatnonzero(np.diff(frequencies) > ZERO_GAP) + 1]
        self.frequencies = frequencies[self.group_starts]
        self.observed = 2 * observed[self.a, self.b].T

    def prepare(self, time):
        return Grid(self, time)

    def direct(self, state, time):
        evolved = (np.exp(-1j * time[:, None] * self.energy) * (state @ self.vectors)) @ self.vectors.T
        return np.abs(evolved)**2 @ OCCUPATION


class Grid:
    def __init__(self, model, time):
        self.model = model
        self.kernel = np.exp(1j * model.frequencies[:, None] * time)

    def sample(self, states):
        model = self.model
        initial = states @ model.vectors
        products = initial[..., model.a].conj() * initial[..., model.b]
        contributions = products[..., None, :] * model.observed
        grouped = np.add.reduceat(contributions, model.group_starts, axis=-1)
        analytic = (grouped @ self.kernel).swapaxes(-1, -2)
        dc = np.einsum("...a,abi,...b->...i", initial.conj(), model.dc, initial).real
        return dc[..., None, :] + analytic.real, analytic


def time_grid(bounds, dt=0.002):
    start, end = bounds
    return np.linspace(start, end, round((end - start) / dt) + 1)


def product_state(angles):
    """Map (...,4,2) Bloch polar and azimuthal angles to (...,16)."""
    theta, phi = np.moveaxis(angles, -1, 0)
    spinor = np.stack((np.cos(theta / 2), np.exp(1j * phi) * np.sin(theta / 2)), axis=-1)
    return np.einsum("...a,...b,...c,...d->...abcd", *np.moveaxis(spinor, -2, 0)).reshape(
        (*angles.shape[:-2], 16)
    )


def phase_samples(analytic):
    envelope = np.abs(analytic)
    unit = np.divide(analytic, envelope, out=np.zeros_like(analytic), where=envelope > 0)
    relative = unit[..., LEFT] * unit[..., RIGHT].conj()
    return envelope.min(axis=-1), relative


def describe(occupation, analytic, gate):
    """Measure six pairs with a common amplitude mask and whole-window Pearson r."""
    min_envelope, relative = phase_samples(analytic)
    valid = min_envelope >= gate
    weights = np.ones(len(occupation))
    weights[[0, -1]] = 0.5
    coverage = float(np.average(valid, weights=weights))
    amplitude = np.ptp(occupation, axis=0) / 2

    centered = occupation - np.average(occupation, axis=0, weights=weights)
    variance = np.average(centered**2, axis=0, weights=weights)
    covariance = np.average(centered[:, LEFT] * centered[:, RIGHT], axis=0, weights=weights)
    normalization = np.sqrt(variance[LEFT] * variance[RIGHT])
    pearson = np.divide(covariance, normalization, out=np.full(6, np.nan),
                        where=(variance[LEFT] > 1e-26) & (variance[RIGHT] > 1e-26))

    locking = (np.average(relative[valid], axis=0, weights=weights[valid]) if np.any(valid)
               else np.full(6, np.nan + 1j * np.nan))
    plv = np.minimum(1, np.abs(locking))
    lag = np.where(plv > 1e-12, np.angle(locking, deg=True), np.nan)
    near_zero = np.all(np.abs(np.angle(relative, deg=True)) <= 30, axis=-1)

    def finite(value):
        return float(value) if np.isfinite(value) else None

    amplitude_pass = np.min(amplitude) >= gate and coverage >= MIN_COVERAGE
    return {
        "amplitude": amplitude.tolist(), "coverage": coverage,
        "min_in_phase": finite(np.min(locking.real)), "max_abs_lag": finite(np.max(np.abs(lag))),
        "min_plv": finite(np.min(plv)), "min_pearson": finite(np.min(pearson)),
        "joint_within_30_deg": float(np.average(valid & near_zero, weights=weights)),
        "amplitude_pass": bool(amplitude_pass),
        "passes": bool(amplitude_pass and np.all(np.abs(lag) <= MAX_LAG) and np.all(plv >= MIN_PLV)),
        "pairs": {pair: {"lag_deg": finite(angle), "plv": finite(p), "pearson": finite(r),
                         "in_phase": finite(re)}
                  for pair, angle, p, r, re in zip(PAIRS, lag, plv, pearson, locking.real)},
    }
