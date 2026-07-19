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
| Rust tests (state/response stage) | 8 passed / 0 failed |
| Rust tests (correlation stage) | 10 passed / 0 failed |
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

```text
cargo fmt --all -- --check
cargo test --release --offline
cargo run --release --offline
python analyze_modes.py
python local_response_analysis.py
python two_body_correlation_analysis.py
```

The Rust binary produces the time series and state distances at both step sizes; the Python scripts produce the modal, local-response, and correlation analyses. The modal audit uses exact diagonalization with degenerate eigenspaces grouped before ranking.

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
