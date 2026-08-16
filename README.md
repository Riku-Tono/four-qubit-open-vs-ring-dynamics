# Four-qubit open chain vs. closed ring, no switch, no quench

**One-paragraph summary.** This package records a direct numerical comparison of two fixed four-qubit topologies — an open chain and a closed ring — evolved from the same initial state `|0000>` under time-independent Hamiltonians, with no bond switching and no quench at any time. Three families of observables are compared: the global state distance, the site-resolved local occupation response `R_i(t)=<n_i>`, and one-body plus two-body Pauli correlations. The central observation is that **at times where the normalized local occupation distributions of the two topologies are close, large differences still remain in the one-body expectations and in the two-body connected correlations.** However, we do **not** claim that this finite set of observables accounts for the full state difference: the one-body-plus-connected feature distance is a distance between finitely many measured quantities, not a reconstruction of the state distance.

---

## 1. What is being compared

Two branches, both starting at `t=0` from `|0000>`, both evolved under their own fixed Hamiltonian for the whole observation window:

| Branch | Bonds |
|---|---|
| open (chain) | `(0,1), (1,2), (2,3)` |
| closed (ring) | open bonds plus `(3,0)` |

There is no bond switching, no quench, and no time-dependent term in either branch. Site labels `0,1,2,3` mean the same physical sites in both branches; the ring differs only by the extra periodic bond `(3,0)`.

This is a **topology comparison**. It is a separate experiment from the pre-existing local-quench study with `t_switch=1.688`. Anything drawn from that older study appears only in the clearly marked supplementary section (§12) and uses the post-switch elapsed time `tau`, which must never be identified with the time origin `t` used everywhere else in this README.

## 2. Hamiltonian and fixed conditions

- `N = 4`
- `Omega = 1.0`
- `J = 0.7`
- initial state `|0000>`
- time window `0 <= t <= 3.688`
- primary step `dt = 0.002`; audit step `dt = 0.001`
- symmetric second-order Trotter propagation; exact diagonalization used for validation and for the modal analysis

## 3. Numerical method and definitions

**State distance.** The global-phase-invariant distance `sqrt(2 - 2 |<psi_open|psi_closed>|)` is used throughout.

**Local response.** The primary site-resolved observable is

```
O_i = n_i = |1_i><1_i| = (I - Z_i)/2
R_i(t) = <psi(t)| n_i |psi(t)> = (1 - <Z_i>)/2
```

`R_i` is the probability that site `i` is in state 1, so `R_i` lies in `[0,1]` and is global-phase invariant. This is a definition introduced here as a common quantity for both the fixed and the quench analyses; it is not a claim that a site-resolved response existed in the older artifacts.

For nonzero total response, `p_i = R_i / sum_j R_j`. Raw L1/L2 and normalized L1/L2 distances are all stored; the primary distribution-shape distance quoted in this README is the **normalized L1**, and total variation is one half of it. Raw and normalized distances are reported separately so that a difference in distribution *shape* is never confused with a difference in *total* response. At `t=0` all responses vanish, so `p_i` is stored as zero and PR/IPR are recorded as NaN.

**One-body and two-body Pauli quantities.** For every `i<j` the package stores the one-body expectations `<X_i>, <Y_i>, <Z_i>` and the raw two-body expectations `<X_iX_j>, <Y_iY_j>, <Z_iZ_j>`. The connected correlations are defined as

```
C^AA_ij = <A_i A_j> - <A_i><A_j>,   A = X, Y, Z
```

Raw and connected quantities are stored separately and must not be conflated. In the displayed 4x4 connected-correlation maps **the diagonal is set to zero by convention**; the diagonal is not an on-site variance and is not interpreted.

**Correlation distances.** For each channel, `upper_L1` and `upper_L2` are taken over the six `i<j` entries of the open-minus-closed difference. For a symmetric zero-diagonal matrix the Frobenius distance is `sqrt(2)*upper_L2`. The combined correlation distance is

```
D_corr = sqrt(D_XX_upper_L2^2 + D_YY_upper_L2^2 + D_ZZ_upper_L2^2)
```

**Feature distances.** Ordinary unnormalized Euclidean distances over, respectively: the 12 one-body values, the six ZZ-connected values, all 18 connected values, and the concatenation of 12 one-body plus 18 connected values. **These are finite-observable feature distances. They are not reconstructions of the full state distance.** Connected correlations here are correlation functions of specified Pauli observables and are not used as entanglement measures.

## 4. State distance — what it shows

| Quantity | Value |
|---|---:|
| maximum state distance in window | `0.546934944036` |
| time of maximum | `t = 3.462` |
| state distance at `t = 3.688` | `0.544161993982` |
| fidelity at `t = 3.688` | `0.725808344260` |

The two branches separate progressively from `t=0` without any switching event, reach maximum separation at `t=3.462`, and remain near that level to the end of the window.

## 5. Local response distribution — what it shows

Site-resolved responses at the time of maximum state distance, `t = 3.462`:

| site | open | closed | open − closed |
|---:|---:|---:|---:|
| 0 | 0.704891575 | 0.558997688 | +0.145893887 |
| 1 | 0.650208418 | 0.558997688 | +0.091210730 |
| 2 | 0.650208418 | 0.558997688 | +0.091210730 |
| 3 | 0.704891575 | 0.558997688 | +0.145893887 |

Closed is uniform across all four sites. Open is higher than closed at every site, with the end sites above the interior sites.

Key times along the window:

| event | t | state distance | normalized L1 | raw L1 |
|---|---:|---:|---:|---:|
| initial | 0.000 | 0 | 0 | 0 |
| max sets first differ | 0.008 | 1.58e-5 | 4.41e-6 | 2.82e-10 |
| state distance > 0.1 | 0.674 | 0.100172 | 0.025660 | 0.010756 |
| state distance > 0.25 | 1.190 | 0.250207 | 0.066673 | 0.075171 |
| response L1 maximum | 2.154 | 0.461112 | **0.111999** | 0.297007 |
| state distance > 0.5 | 2.494 | 0.500048 | 0.106085 | 0.348019 |
| state distance maximum | 3.462 | **0.546935** | 0.040354 | 0.474209 |
| final | 3.688 | 0.544162 | 0.020227 | **0.529384** |

At `t=3.462` the normalized L1 distance is `0.040354` (total variation `0.020177`) while the state distance is `0.546935`. The *shape* of the normalized occupation distribution is comparatively close there, even though the raw L1 difference `0.474209` shows a substantial difference in total response.

Maximum-position behaviour: fixed closed keeps all four sites equal at every saved time (fourfold tie throughout, zero position changes). Fixed open has `R_0=R_3` and `R_1=R_2` by reflection symmetry; after leaving the initial near-zero fourfold tie at `t=0.008`, the end pair `0|3` is maximal for the rest of the window, with no migration, reflection, or round trip. Ties are resolved with tolerance `max(1e-10, 1e-8*|R_max|)`; the persistent ties here are symmetry-induced, not numerical noise.

Correspondence with the state distance: `Pearson(D_state, normalized L1) = 0.688010`, `Spearman = 0.439491`. Correlation coefficients here describe co-variation of time series only and carry no causal content.

## 6. One-body and two-body Pauli correlations — what they show

Channel maxima over the window:

| channel | max upper-L2 | t |
|---|---:|---:|
| XX | 0.604701641381 | 3.340 |
| YY | 0.469912653201 | 3.190 |
| ZZ | 0.665558112764 | 3.688 |
| combined `D_corr` | **0.988619981363** | **3.200** |

`D_corr` first exceeds `1e-8` at `t = 0.012`.

At the time of maximum state distance, `t = 3.462`, the combined correlation distance is `D_corr = 0.974637690628`, while the normalized local occupation L1 distance is only `0.040354`. Connected correlations at that time:

| channel | pair | open | closed | open − closed |
|---|---:|---:|---:|---:|
| XX | 01 | -0.196609683 | 0.00676489627 | -0.203374579 |
| XX | 02 | -0.0235895718 | -0.0297264872 | 0.00613691543 |
| XX | 03 | 0.478921356 | 0.00676489627 | **0.47215646** |
| XX | 12 | -0.231523258 | 0.00676489627 | -0.238288154 |
| XX | 13 | -0.0235895718 | -0.0297264872 | 0.00613691543 |
| XX | 23 | -0.196609683 | 0.00676489627 | -0.203374579 |
| YY | 01 | 0.37779505 | 0.316782279 | 0.0610127711 |
| YY | 02 | 0.213084588 | 0.368699724 | -0.155615136 |
| YY | 03 | 0.0486186388 | 0.316782279 | -0.26816364 |
| YY | 12 | 0.56456502 | 0.316782279 | 0.247782741 |
| YY | 13 | 0.213084588 | 0.368699724 | -0.155615136 |
| YY | 23 | 0.37779505 | 0.316782279 | 0.0610127711 |
| ZZ | 01 | 0.416403264 | 0.464558696 | -0.0481554322 |
| ZZ | 02 | 0.145375256 | 0.456845621 | -0.311470366 |
| ZZ | 03 | 0.0238689717 | 0.464558696 | -0.440689725 |
| ZZ | 12 | 0.391416953 | 0.464558696 | -0.0731417431 |
| ZZ | 13 | 0.145375256 | 0.456845621 | -0.311470366 |
| ZZ | 23 | 0.416403264 | 0.464558696 | -0.0481554322 |

The single largest connected-correlation difference at `t=3.462` is the XX channel, pair `(0,3)`, at `0.472156` — that is, the pair joined by the ring's extra bond.

Averaged over the interval where the state distance exceeds 0.5, the mean channel distances are `XX = 0.564273`, `YY = 0.40003`, `ZZ = 0.63102`. ZZ is the largest on average, but XX is of comparable size, so the separation in this window is not described by any single channel alone.

## 7. Correspondence between state distance, local response, and correlation distances

| feature distance | Pearson vs `D_state` | Spearman vs `D_state` | max | t of max | offset from `t=3.462` |
|---|---:|---:|---:|---:|---:|
| `D_XX` | 0.904029 | 0.997396 | 0.604702 | 3.340 | -0.122 |
| `D_YY` | 0.941816 | 0.958606 | 0.469913 | 3.190 | -0.272 |
| `D_ZZ` | 0.952728 | 0.976435 | 0.665558 | 3.688 | +0.226 |
| `D_corr` | **0.953794** | **0.990560** | 0.98862 | 3.200 | -0.262 |
| `D_onebody` | 0.969556 | 0.998150 | 0.667066 | 3.688 | +0.226 |
| `D_ZZ_connected` | 0.952728 | 0.976435 | 0.665558 | 3.688 | +0.226 |
| `D_all_connected` | 0.953794 | 0.990560 | 0.98862 | 3.200 | -0.262 |
| `D_onebody_all_connected` | 0.981143 | 0.998150 | 1.18283 | 3.688 | +0.226 |

The three maxima do not coincide: `D_corr` peaks at `t=3.200`, the state distance at `t=3.462`, and the normalized local occupation difference at `t=2.154` (offset of 1.046 from the `D_corr` peak). These quantities track each other closely in rank order but are distinct functions of time with distinct extrema.

All coefficients above are descriptive measures of co-variation between time series. They do not license any causal statement, in either direction, about the ring bond and the observed differences.

## 8. Symmetries of the two branches

The closed branch retains the rotation symmetry of the ring together with the rotation-invariant initial state `|0000>`; the open branch retains reflection symmetry. Both are verified numerically as residuals:

| symmetry residual | max |
|---|---:|
| `closed_nearest_rotation_residual` | 7.438e-15 |
| `closed_opposite_rotation_residual` | 5.773e-15 |
| `open_reflection_01_23_residual` | 7.161e-15 |
| `open_reflection_02_13_residual` | 6.078e-15 |
| `XX_closed_nearest_rotation_residual` | 4.996e-15 |
| `XX_closed_opposite_rotation_residual` | 5.176e-15 |
| `XX_open_reflection_01_23_residual` | 6.217e-15 |
| `XX_open_reflection_02_13_residual` | 3.220e-15 |
| `YY_closed_nearest_rotation_residual` | 6.634e-15 |
| `YY_closed_opposite_rotation_residual` | 5.773e-15 |
| `YY_open_reflection_01_23_residual` | 7.161e-15 |
| `YY_open_reflection_02_13_residual` | 3.761e-15 |
| `ZZ_closed_nearest_rotation_residual` | 7.438e-15 |
| `ZZ_closed_opposite_rotation_residual` | 4.774e-15 |
| `ZZ_open_reflection_01_23_residual` | 6.051e-15 |
| `ZZ_open_reflection_02_13_residual` | 6.078e-15 |

All symmetry residuals are at or below roughly `1e-14`. In the closed branch the four nearest-neighbour pairs agree with each other and the two opposite pairs agree with each other. In the open branch reflection-related pairs agree, but pairs containing an end site differ from interior and separated pairs, and that difference is what appears in the correlation maps.

## 9. Reconstruction from degenerate eigenspace pairs

Degenerate energies are grouped into eigenspaces and the analysis uses eigenspace projectors `P_E`, with contributions `A_EF,i = Tr(P_E rho_0 P_F O_i)`. For `E != F` the real contribution of an unordered pair is `2 Re[A_EF,i exp(-i(E-F)t)]`; the `E=F` diagonal term is counted once. **Because degenerate levels are grouped before ranking, the resulting weights do not depend on the arbitrary choice of eigenvectors inside a degenerate space.**

Connected correlations are not linear in the state, so they cannot be decomposed additively over eigenspace pairs. The procedure is: decompose the raw two-body terms and the one-body terms linearly over eigenspace pairs, reconstruct each, and only then form the connected quantity as **raw minus the product of the reconstructed one-body terms**. Because that product is nonlinear, the connected ranking score is a leave-one-eigenspace-pair-out RMS influence and is not claimed to be a unique additive decomposition.

Eigenspace pairs required to reach a given fraction of explained energy:

| branch | channel | total pairs | k90 | k95 | k99 |
|---|---|---:|---:|---:|---:|
| open | XX | 55 | 11 | 14 | 23 |
| open | YY | 55 | 19 | 22 | 37 |
| open | ZZ | 55 | 15 | 25 | 31 |
| closed | XX | 21 | 7 | 7 | 11 |
| closed | YY | 21 | 11 | 12 | 14 |
| closed | ZZ | 21 | 8 | 8 | 14 |

At the 99% level the open branch requires 23 to 37 pairs and the closed branch 11 to 14. The open branch in particular is more safely described as **multi-mode competition** than as a few-mode structure.

For the local response alone the counts are smaller: open reaches 95% of the total response with 6 eigenspace pairs and 99% with 7, while the spatially non-uniform component reaches 95% with 4 pairs but needs 15 for 99%. The closed response is spatially uniform, so it has no pairs carrying a maximum-position shift.

## 10. Numerical audit

All audits below are reported as observed residuals; none of them is an interpretation.

**Correlation analysis audit**

| check | residual / result |
|---|---:|
| `exact_onebody` | 7.140e-07 |
| `exact_raw_two_body` | 6.221e-07 |
| `exact_connected` | 5.473e-07 |
| `max_imaginary_exact` | 1.405e-16 |
| `max_expectation_bound_excess` | 8.882e-16 |
| `max_hermiticity_residual` | 0.000e+00 |
| `connected_raw_minus_product_identity` | 0.000e+00 |
| `global_phase_expectation_residual` | 2.225e-16 |
| `dt_onebody` | 5.355e-07 |
| `dt_raw_two_body` | 4.666e-07 |
| `dt_connected` | 4.104e-07 |
| `mode_raw_max_residual` | 6.661e-16 |
| `mode_onebody_max_residual` | 9.437e-16 |
| `connected_from_reconstructed_ingredients_max_residual` | 1.443e-15 |
| `dt002_max_norm_error` | 6.402e-13 |
| `dt001_max_norm_error` | 7.814e-13 |

**State and local-response audit**

| check | residual / result |
|---|---:|
| `cargo fmt --all -- --check` | PASS |
| Current Rust test suite | 10 passed / 0 failed |
| exact vs. Trotter local response | 2.82e-07 |
| exact vs. Trotter normalized distribution L1 | 4.90e-07 |
| `dt=0.002` vs `dt=0.001` local response | 2.12e-07 |
| `dt=0.002` vs `dt=0.001` normalized distribution L1 | 3.68e-07 |
| exact vs. Trotter maximum-set mismatches | 0 |
| step-size maximum-set mismatches | 0 |
| exact state-distance check | 4.59e-08 |
| mode sum, open | 9.99e-16 |
| mode sum, closed | 5.55e-16 |
| `dt=0.002` vs `dt=0.001` state distance, open | 4.8514666e-07 |
| `dt=0.002` vs `dt=0.001` state distance, closed | 7.2695678e-07 |
| global-phase invariance | Rust test PASS |

Exact-diagonalization agreement, step-size convergence, norm conservation, Hermiticity, reality of expectation values, Pauli expectation bounds, global-phase invariance, and mode-sum completeness all pass. Only global-phase-invariant quantities — expectation values, connected correlations, and the distances built from them — are used for physical comparison.

## 11. How to reproduce and audit

Tested with Python 3.14.5, NumPy 2.5.0, and Rust/Cargo 1.97.0. NumPy is pinned in `requirements.txt`; Cargo crates are included in `vendor/`, so the Rust stages do not require network access.

```text
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python verify_reproduction.py
```

On macOS/Linux, use `.venv/bin/python` instead. The verifier creates a clean temporary working directory, executes formatting, all 10 Rust tests, both Trotter runs, and all three Python analyses. It then compares all 40 regenerated CSV files against the packaged references with absolute tolerance `2e-12` and relative tolerance `1e-10`, and writes `reproduction_audit.json` with the environment, command results, comparison maxima, and final verdict.

The Rust binary produces the time series and state distances at both step sizes; the Python scripts produce the modal, local-response, and correlation analyses. The modal audit uses exact diagonalization with degenerate eigenspaces grouped before ranking. Cargo/test PASS results are recorded by the verifier and are not hard-coded into the numerical CSV audit.

Everything needed to re-run is in the archive:

- `four_qubit_direct_open_closed_two_body_correlation_reproducible_2026-07-20.zip`

## 12. Key output files

**Reports**

- `four_qubit_direct_open_closed_RUN_RESULTS.md` — run conditions, state distance, fidelity, step-size and exact audits
- `LOCAL_RESPONSE_ANALYSIS_REPORT.md` / `four_qubit_local_response_analysis_REPORT.md` — local response `R_i(t)`, maximum-position statistics, modal decomposition
- `CORRELATION_ANALYSIS_REPORT.md` — one-body and two-body Pauli correlations, channel distances, symmetry audit, eigenspace-pair reconstruction

**State and mode CSVs**

- `direct_open_closed_dt002.csv`, `direct_open_closed_dt001.csv`, `direct_open_closed_dt_comparison.csv`
- `direct_open_closed_summary.csv`
- `direct_open_closed_mode_weights.csv`, `direct_open_closed_top6_pair_phases.csv`
- `direct_open_closed_exact_validation.csv`
- `direct_open_closed_pauli_dt002.csv`, `direct_open_closed_pauli_dt001.csv`

**Local-response CSVs**

- `direct_open_closed_local_response_dt002.csv`, `direct_open_closed_local_response_dt001.csv`
- `local_response_timeseries.csv`, `local_response_max_position_timeseries.csv`
- `local_response_distribution_distance.csv`, `state_response_comparison.csv`, `local_response_key_times.csv`
- `local_response_mode_pair_contributions.csv`, `dominant_mode_pair_summary.csv`, `local_response_mode_reconstruction_summary.csv`
- `local_response_exact_validation.csv`, `local_response_dt_convergence.csv`, `local_response_numerical_audit.csv`
- `state_response_correlations.csv`

**Correlation CSVs**

- `one_body_pauli_expectations.csv`, `two_body_raw_pauli_correlations.csv`, `two_body_connected_correlations.csv`
- `correlation_distance_timeseries.csv`, `state_response_correlation_comparison.csv`, `correlation_distance_correlations.csv`
- `representative_correlation_maps.csv`, `correlation_symmetry_audit.csv`
- `correlation_mode_pair_amplitudes.csv`, `correlation_mode_pair_contributions.csv`, `correlation_dominant_mode_pair_summary.csv`, `correlation_mode_reconstruction_summary.csv`
- `correlation_exact_validation.csv`, `correlation_dt_convergence.csv`, `correlation_numerical_audit.csv`

**Supplementary quench artifacts (separate experiment).** The package also contains a re-analysis of the pre-existing `t_switch=1.688` local-quench run using the same `R_i=<n_i>` definition, in `quench_local_response_timeseries.csv` and `fixed_vs_quench_response_comparison.csv`. That run is parameterized by post-switch elapsed time `tau`, which is a different time origin from `t`; the two must not be placed on a common axis. In that re-analysis both quench branches keep the maximum set `0|3` over `0 <= tau <= 2` with zero position changes, and the maximum normalized response difference is `0.01747`, smaller than the fixed-topology value `0.1120`. The older artifacts contained no site-resolved response definition, so earlier descriptions of response evolution there cannot be treated as a reproduction of this quantity.

## 13. What can and cannot be claimed

**Supported by the numbers in this package**

- Two fixed topologies starting from the same state separate without any switching event, reaching state distance `0.546935` at `t=3.462`.
- At `t=3.462` the normalized local occupation distributions are comparatively close (normalized L1 `0.040354`) while the combined connected-correlation distance is `0.974638`, and the largest single connected difference is XX pair `(0,3)` at `0.472156`.
- `D_corr` and the state distance co-vary strongly (Pearson `0.953794`, Spearman `0.990560`) but peak at different times (`3.200` vs `3.462`), and neither coincides with the local-occupation peak at `2.154`.
- Ring rotation symmetry and chain reflection symmetry hold to residuals of order `1e-14` or below.
- 99% reconstruction requires 23–37 eigenspace pairs on the open side and 11–14 on the closed side.
- All listed exactness, convergence, norm, Hermiticity, reality, bound, phase-invariance, and mode-sum audits pass.

**Explicitly not claimed**

- That the ring bond *causes* any transfer of correlation. Correlation coefficients here are descriptive only.
- That two-body correlations account for the full state difference. The feature distances are distances between finitely many observables.
- That the ZZ correlation difference is the cause of the state distance. On average ZZ is the largest channel above state distance 0.5, but XX is comparable.
- That entanglement has been detected. Connected correlations are correlation functions of chosen Pauli observables, not entanglement measures.
- Any general boundary mechanism, transport law, phase transition, or thermodynamic-limit statement.
- Any identity of mechanism with the `N=7` Lindblad backflow study.

**Scope.** This is a finite `N=4`, finite-time (`0 <= t <= 3.688`), fixed-Hamiltonian numerical record at `Omega=1.0`, `J=0.7`, from the initial state `|0000>`. Every statement above is scoped to those conditions and to the observables defined in §3.
---

## 14. Subsequent weak-link and symmetry extensions (2026-08-17)

The original fixed-endpoint record in §§1–13 is preserved unchanged above. The analyses in §§14–18 are later extensions built from that sealed baseline. They do not replace the direct open-versus-closed comparison.

The main extension uses

```text
H(lambda) = Omega/2 sum_i X_i
          + J/2 [Z0Z1 + Z1Z2 + Z2Z3 + lambda Z3Z0]
```

with the same `N=4`, `Omega=1`, `J=0.7`, time window, and site ordering. The sampled grid is `lambda=0,0.05,...,1`; `lambda=0` is the open-chain reference and `lambda=1` is the closed-ring endpoint. This is a closing-bond-strength sweep, not a claim of continuous topology.

For connected correlations, the signed branch difference is

```text
Delta C_ij^AA(t,lambda) = C_ij^AA(t,lambda) - C_ij^AA(t,0).
```

The pointwise squared-difference share of pair `ij` is

```text
w_ij = sum_A (Delta C_ij^AA)^2
       / sum_{k<l,A} (Delta C_kl^AA)^2,
```

and `Q_ij` denotes the corresponding time-integrated numerator divided by the time-integrated denominator. Pointwise ratios are undefined where the all-channel difference norm is at or below `1e-8`. Squaring discards sign, so every concentration result must be read together with the signed-channel results.

### 14.1 Weak-link sweep

The endpoint checks against the sealed open/closed baseline all pass. At `lambda=1`, the maximum state distance is `0.546934958021` at `t=3.462`, and its time RMS is `0.393869262797`.

At `lambda=1`, the time-integrated seam shares for pair `(0,3)` are:

| channel | integrated seam share |
|---|---:|
| XX | 0.412739 |
| YY | 0.282331 |
| ZZ | 0.548788 |
| all channels | 0.451198 |

At `t=3.462`, the all-channel pointwise seam share is `0.514835`. The signed seam differences are `-0.472156818` in XX, `+0.268163791` in YY, and `+0.440689566` in ZZ. The seam is therefore prominent in squared-difference space, but the three channels do not move in one common signed direction.

Using the predeclared rule `Q03 > 1/3` and requiring it to hold for at least 75% of nonzero `lambda` values and at least 50% of the valid nonzero-`lambda` time grid, the observed fractions are `0.550` and `0.761`. The mechanical classification is therefore **limited**, not broadly robust.

### 14.2 Pair-space redistribution

The integrated top pair changes from the reflection-tied pair set `02|13` for `lambda=0.05–0.15` to the seam pair `03` for `lambda=0.20–1.00`. Linear interpolation between sampled values places the seam-versus-`02/13` crossings near `lambda=0.166`. The separate fixed-threshold event `Q03=1/3` occurs later, near `lambda=0.452`.

At `lambda=1`, the integrated six-pair decomposition is:

| pair | integrated share |
|---|---:|
| 01 | 0.063122 |
| 02 | 0.132206 |
| 03 | 0.451198 |
| 12 | 0.158144 |
| 13 | 0.132206 |
| 23 | 0.063122 |

When `w03<1/3`, the energy-weighted non-seam remainder is split `0.170` to the two seam-adjacent pairs `01+23` and `0.830` to `12+02+13`. Thus the data do not show a simple seam-to-adjacent-to-remaining sequence. These are concurrent partitions of a squared observable-difference denominator, not evidence of propagation.

Reflection-related equalities are numerically retained: `01=23` and `02=13`, with integrated-share residuals at approximately `1e-15`.

## 15. Signed versus squared connected-correlation differences

The signed audit separates the integral

```text
I = integral Delta C dt
```

from the absolute integral

```text
M = integral |Delta C| dt
```

and reports the signed balance `P=I/M` where defined. At `lambda=1`, the seam results are:

| channel | Q03 | I | M | P | sign changes |
|---|---:|---:|---:|---:|---:|
| XX | 0.412739 | -0.609628 | 0.609628 | -1.000000 | 0 |
| YY | 0.282331 | +0.355200 | 0.385902 | +0.920442 | 2 |
| ZZ | 0.548788 | +0.913784 | 0.913784 | +1.000000 | 0 |

The large seam squared contribution is therefore compatible with persistent negative XX and predominantly positive YY/ZZ. `Q` measures magnitude concentration; `P` measures signed consistency and cancellation.

Sign-canceling components do occur. The strongest fractional cancellation is at `lambda=0.30`, pair `23`, channel YY, with `M=0.0496291` and `P=-0.00503727`, but its magnitude ranks only 309th of 360 components. Among components with channel share at least `0.1`, the smallest observed `|P|` is about `0.662`. The data therefore do not support the stronger statement that the dominant squared contributions are strongly canceled.

## 16. End–interior occupation contrast versus seam concentration

For the original `|0000>` initial state, define the signed one-body contrast

```text
E(t,lambda) = [R0 + R3 - R1 - R2] / 2.
```

The time RMS of `E` is nonincreasing across the sampled `lambda` grid, falling from `0.0848382` at `lambda=0` to a numerical residual `6.922e-17` at `lambda=1`. In contrast, integrated `Q03` is nondecreasing over its 20 defined nonzero-`lambda` values, from `0.235105` at `lambda=0.05` to `0.451198` at `lambda=1`.

Across the 20 summaries, `Spearman(E_RMS,Q03)=-1.000000` and `Pearson=-0.980965`. Across all valid `(lambda,t)` points, `Spearman(|E|,w03)=-0.840902` and `Pearson=-0.814780`. These are descriptive relationships on a finite, serially related grid; they are not independent-sample inference or causal evidence.

The predeclared quartile check—`|E|` at or below its 25th percentile and `w03` at or above its 75th percentile—selects 5,387 points, or `14.6657%` of valid points, spanning `lambda=0.20–1.00` and `t=0.012–1.276`. Dividing `|E|` by total occupation preserves the suppression trend where defined, so the raw trend is not explained only by a uniform change in total occupation.

The `|E|`, `w03`, and state-distance peaks generally occur at different times. None of these three quantities is used as a proxy for either of the others.

## 17. Initial-state sensitivity

Only two additional initial states were introduced, with all Hamiltonian, grid, operator, floor, and per-state `lambda=0` reference definitions held fixed.

| initial state | E_RMS at lambda=0 | E_RMS at lambda=1 | Q03 at lambda=0.05 | Q03 at lambda=1 | state-distance RMS at lambda=1 |
|---|---:|---:|---:|---:|---:|
| `|0000>` | 0.084838 | 0.000000 | 0.235105 | 0.451198 | 0.393869 |
| `|1000>` | 0.247139 | 0.249492 | 0.020359 | 0.109157 | 0.457350 |
| `|0100>` | 0.275810 | 0.249492 | 0.020359 | 0.109157 | 0.457350 |

For `|1000>`, `E_RMS` rises slightly rather than being suppressed, while `Q03` is nonmonotonic and remains below the `|0000>` curve at every defined `lambda`. For `|0100>`, `E_RMS` decreases from a larger open-endpoint value. Both asymmetric single-site states reach the same three endpoint values at `lambda=1`, but equality of these selected metrics is not equality of the full states.

The full `|1000>` and `|0100>` curves agree to `4.382e-15` for `Q03` and `1.221e-15` for time-RMS state distance. Their `E_RMS` curves differ for `lambda<1`. Thus the suppression of `E` together with growth of `Q03`, seen for `|0000>`, is initial-state sensitive rather than a universal Hamiltonian-only pattern.

### 17.1 Exact symmetry explanation

The `Q03` and state-distance equalities for `|1000>` and `|0100>` are enforced by an antiunitary symmetry.

Let `R` implement the reflection `r=(0 3)(1 2)`, and define

```text
Gamma = Z0 Z1 Z2 Z3 X0 X2.
```

For every sampled `lambda`—indeed for the complete Hamiltonian family—

```text
[R,H(lambda)] = 0,
Gamma H(lambda) Gamma^-1 = -H(lambda),
R Gamma |1000> = -|0100>.
```

Because `H(lambda)` is real in the computational basis, complex conjugation `K` gives the antiunitary map `Theta=R Gamma K`, which commutes with the same-time evolution:

```text
Theta exp[-i H(lambda)t] = exp[-i H(lambda)t] Theta.
```

The same `Theta` maps both the `lambda=0` and `lambda>0` evolutions, so it preserves the overlap magnitude used in the state distance. For connected correlations it only permutes pairs and changes channel/site signs; squaring removes those signs, the all-pair denominator is invariant, and the seam pair `03` is fixed by reflection. Therefore pointwise `w03` and integrated `Q03` are symmetry invariants.

The occupation contrast is not invariant because some `Z_i` map to `-Z_j`, which sends `R_i=(1-Z_i)/2` to `1-R_j`. At `lambda=1`, the saved arrays satisfy `E_0100(t)=-E_1000(t)` to `8.33e-16`, explaining the equal endpoint RMS values without implying pointwise equality of the occupations.

The operator identities have zero matrix residual. Stored-array residuals are `4.21e-15` for the all-channel difference norm, `1.00e-15` for reconstructed seam squared energy, and `4.40e-15` for integrated `Q03`. A larger pointwise ratio discrepancy near the `1e-8` floor is roundoff amplification, not symmetry breaking.

## 18. Extension artifacts, audit status, and updated boundary

All extension stages report **PASS** for their numerical audits and independent regeneration checks. The sealed baseline and previously completed stage artifacts were hash-checked before and after the relevant work and were not modified.

Primary local reports:

- [Weak-link sweep](C:/Users/yauki/Documents/Codex/2026-08-17/pe/outputs/four_qubit_weak_link_sweep_2026-08-17/weak_link_sweep/REPORT.md)
- [Pair-space redistribution](C:/Users/yauki/Documents/Codex/2026-08-17/pe/outputs/four_qubit_weak_link_sweep_2026-08-17/pair_space_redistribution/REPORT.md)
- [Signed versus squared differences](C:/Users/yauki/Documents/Codex/2026-08-17/vs/outputs/signed_vs_squared/REPORT.md)
- [End–interior contrast versus seam concentration](C:/Users/yauki/Documents/Codex/2026-08-17/vs/outputs/end_interior_vs_seam/REPORT.md)
- [`|1000>` comparison](C:/Users/yauki/Documents/Codex/2026-08-17/vs/outputs/initial_1000_three_metric_comparison/REPORT.md)
- [`|0100>` three-state comparison](C:/Users/yauki/Documents/Codex/2026-08-17/vs/outputs/initial_0100_three_state_comparison/REPORT.md)
- [Antiunitary symmetry explanation](C:/Users/yauki/Documents/Codex/2026-08-17/vs/outputs/symmetry_explanation_1000_0100/REPORT.md)

The added results support only finite-`N=4`, finite-time, sampled-`lambda`, selected-observable statements under the stated Hamiltonian and initial states. They do not establish a continuous topology change, a phase transition, transport, information flow, correlation transfer, causation, entanglement, a general initial-state law, or thermodynamic-limit behavior. In particular, a larger squared contribution is not restated as a larger signed correlation, and suppression of `E` is not treated as identification of the full state with a ring state.


