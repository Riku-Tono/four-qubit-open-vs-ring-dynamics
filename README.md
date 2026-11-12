# Four-qubit open chain vs. closed ring, no switch, no quench

**One-paragraph summary.** This package records a direct numerical comparison of two fixed four-qubit topologies — an open chain and a closed ring — evolved from the same initial state `|0000>` under time-independent Hamiltonians, with no bond switching and no quench at any time. Three families of observables are compared: the global state distance, the site-resolved local occupation response `R_i(t)=<n_i>`, and one-body plus two-body Pauli correlations. The central observation is that **at times where the normalized local occupation distributions of the two topologies are close, large differences still remain in the one-body expectations and in the two-body connected correlations.** However, we do **not** claim that this finite set of observables accounts for the full state difference: the one-body-plus-connected feature distance is a distance between finitely many measured quantities, not a reconstruction of the state distance.

---

**Short-time extension.** Sections 21–24 use the initial state `|10...0>` to study distance-dependent connected correlations on open chains and weak rings. Section 24 classifies the weak-link strengths that cancel the leading YY or ZZ term for every antipodal pair of an even ring. [Distance laws and cancellation roots](#21-eight-qubit-check-of-the-distance-rule-2026-11-09)

**Rotation-axis extension.** Sections 25–27 track the occupation-sum residual of an antiphase pair after rotating site 0 of the initial state $\lvert1000\rangle$. In the later-window tests, the short-window axis beats both pure axes only in the first open-chain window. The long-time optimum is Y for the open chain and X for the full ring. [Rotation-axis results](#25-complementary-occupations-and-signed-site-0-rotations-2026-11-10)

**In-phase occupations and preparation errors.** Sections 29–31 follow a four-site product state whose open-chain occupations stay in phase through the 24 held-out windows ending at $256\pi$. The first PLV failure on the tested 0.001 rad grid occurs at 0.091 rad, on pair 03 in four of 496 fixed error directions. [In-phase results](#29-four-site-in-phase-occupations-and-long-window-validation-2026-11-10) · [PLV boundary](#31-independent-angle-errors-and-the-first-plv-failure-2026-11-10)

**Finite-shot graph identification.** The November 13 study extends the comparison to all 64 labeled four-qubit graphs. With independent 5% readout flips and a fixed measurement schedule, the certificates bound the optimum for 95% success on each graph between 53 and 99 shots. [Section 32](#32-finite-shot-identification-of-all-64-four-qubit-graphs-2026-11-13)

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



## 19. Endpoint transverse-field perturbations (2026-11-09)

This extension tests the symmetry in §17.1 with `N=4`, `Omega=1`, `J=0.7`, and `0 <= t <= 3.688`. The weak-link branch has `lambda=0.5` and the open reference has `lambda=0`. Each comparison applies the same perturbation to both branches and uses initial states `A=|1000>` and `B=|0100>`.

### 19.1 Symmetry and the short-time occupation residual

Write `H_same = H + epsilon*(X0 + X3)` and `H_opp = H + epsilon*(X0 - X3)`. Since `Theta` sends `X0` to `-X3` and `X3` to `-X0`, the Hamiltonians obey

```text
Theta H_same(lambda,epsilon) Theta^-1 = -H_same(lambda,epsilon)
Theta H_opp(lambda,epsilon) Theta^-1 = -H_opp(lambda,-epsilon)
```

Define the occupation residual using the §3 occupations as `R_C(t) = R0^A(t) + R3^B(t) - 1`. For the same-direction fields the symmetry gives `R_C(t)=0` exactly. The maximum observed `|R_C|` at `lambda=0.5` and `epsilon=+/-0.05` was `2.78e-15`.

For opposite-direction fields, the short-time expansion is `R_C(t) = -2*Omega*epsilon*t^2 + O(t^4)`. At `epsilon=+0.05`, successive halvings of the observation time reduced the error in `R_C/t^2` by about a factor of four, reaching `-0.099999350` at `t=0.005` against the predicted `-0.1`. The negative perturbation gave `+0.099999349` against `+0.1`.

### 19.2 Seam share and linear response

`Q_A(epsilon)` and `Q_B(epsilon)` are the integrated `Q03` of §14 for the two initial states. Under opposite-direction fields, the sign-reversing symmetry gives `Q_A(epsilon)=Q_B(-epsilon)` exactly where the share is defined. The largest absolute discrepancy across the ten signed samples below was `1.08e-15`.

The same-sign difference `Delta Q(epsilon) = Q_A(epsilon) - Q_B(epsilon)` is therefore odd in `epsilon`. At `epsilon=+0.05`, the two shares are `0.075032997` and `0.017521337`. The table gives the difference for positive `epsilon` and the response averaged over both signs. Negative `epsilon` reverses the sign of `Delta Q`.

| `abs(epsilon)` | `epsilon^2` | `Delta Q` at positive `epsilon` | `Delta Q/epsilon` |
|---:|---:|---:|---:|
| 0.05 | 2.5e-3 | 0.057511660 | 1.150233194 |
| 0.025 | 6.25e-4 | 0.030345067 | 1.213802682 |
| 0.0125 | 1.5625e-4 | 0.015374866 | 1.229989276 |
| 0.00625 | 3.90625e-5 | 0.007712839 | 1.234054240 |
| 0.003125 | 9.765625e-6 | 0.003859599 | 1.235071625 |

The intercept of `Delta Q/epsilon` against `epsilon^2` estimates **`Delta Q'(0) = 1.23541`** for this Hamiltonian and observation window. A line through the two smallest magnitudes gives `1.235410753`. A quadratic in `epsilon^2` through the three smallest gives `1.235410855`, in agreement with a separate calculation that differentiated the propagator at `epsilon=0` on the same time grid.

### 19.3 Numerical audit and reproduction

The scripts evaluate each state directly by diagonalizing the 16-by-16 Hamiltonian. The initial perturbation check gives a maximum state-norm error of `1.55e-15`. Halving the trapezoidal integration interval from `0.002` to `0.001` changes its `Q03` values by at most `3.90e-8`. In the zero-field extrapolation, halving the interval from `0.001` to `0.0005` shifts the fitted intercepts by about `1.12e-7`.

Run the three `run.py` scripts in order from their experiment folders with NumPy 2.0 or later and Matplotlib. Each folder contains its numerical output and figures.

- Endpoint-field check. [Report](experiments/endpoint-field-2026-11-09/REPORT.md) · [Occupation residual](experiments/endpoint-field-2026-11-09/occupation.png) · [Seam-share comparisons](experiments/endpoint-field-2026-11-09/q03-comparisons.png)
- Perturbation halving. [Report](experiments/endpoint-field-scaling-2026-11-09/REPORT.md) · [Scaling and symmetry error](experiments/endpoint-field-scaling-2026-11-09/q03-scaling.png)
- Zero-field extrapolation. [Report](experiments/endpoint-field-zero-slope-2026-11-09/REPORT.md) · [Response versus epsilon squared](experiments/endpoint-field-zero-slope-2026-11-09/q03-zero-slope.png)

## 20. Six-qubit extension and short-time seam-neighbor correlations (2026-11-09)

The extension uses `Omega=1`, `J=0.7`, and the same weak-link strength `lambda=0.5` against the open reference `lambda=0`. Sites are numbered from the left. The initial state A has only the left endpoint excited in both sizes, `|1000>` and `|100000>`. Unless stated otherwise, integrated values use `T=3.688` and trapezoidal spacing `dt=0.001`.

### 20.1 Antiunitary pairing and the zero-field slope

For even `N`, let `R` be site reflection, `K` complex conjugation in the computational basis, and `Gamma_N = (product of all Z_i)(product of X_i on even sites)`. With `Theta_N = R Gamma_N K`, the perturbation `epsilon*(X_0 - X_(N-1))` gives

```text
Theta_N H(lambda,epsilon) Theta_N^-1 = -H(lambda,-epsilon)
Theta_N exp[-it H(lambda,epsilon)] Theta_N^-1 = exp[-it H(lambda,-epsilon)]
Theta_4 |1000> = -|0100>
Theta_6 |100000> = |010100>
```

Conjugation by `Gamma_N` flips the transverse fields and the ZZ bonds, since every bond including the seam joins opposite site parities. Reflection swaps the endpoint fields. Antilinearity changes the sign of `i`, so the transformed evolution is at the same time. Reflection fixes the seam pair and preserves the sum of squared correlation differences. Consequently `Q_A(epsilon) = Q_B(-epsilon)` for the integrated seam share at both sizes. The maximum six-qubit sign-reversal discrepancy over the two time grids was `6.6e-15` or less.

The same-sign difference is odd in `epsilon`. Define `I_ij = integral sum_a (C_ij^aa(lambda=0.5) - C_ij^aa(lambda=0))^2 dt`, `S = sum_(i<j) I_ij`, and `Q_A = I_seam/S` for initial state A. Its slope at zero field is

```text
chi_N = d[Q_A(epsilon) - Q_B(epsilon)]/d epsilon at epsilon=0
      = 2 I_seam'/S - 2 I_seam S'/S^2
```

| Size | `2 I_seam'/S` | `-2 I_seam S'/S^2` | `chi_N` |
|---:|---:|---:|---:|
| 4 | 1.183264934 | +0.052145921 | 1.235410855 |
| 6 | 2.502318543 | -0.259012235 | 2.243306307 |

The propagator derivative and the extrapolation in `epsilon^2` agree within `8.3e-10` for six qubits on the primary grid. The normalized slope grows by a factor of `1.81584`. The unnormalized seam derivative falls from `0.21316` to `0.19190`, while `S(0)` falls from `0.360287` to `0.153382`. The comparison is for these initial states, the given observation window, and the all-pair normalization at each size.

### 20.2 Which pairs account for the change in the denominator?

At `epsilon=0`, match the two sites nearest each end by `0,1,2,3 -> 0,1,4,5`. Count the nine six-qubit pairs involving a central site separately. The two pairs with a two-edge path through the seam account for `89.20%` of the net reduction in `S` at `T=3.688`.

| Pair in N=4 | Matched pair in N=6 | `I_4` | `I_6` | Reduction `I_4 - I_6` |
|:---|:---|---:|---:|---:|
| 02 | 04 | 0.111023 | 0.021517 | 0.089506 |
| 13 | 15 | 0.121227 | 0.026164 | 0.095063 |
| 03 (seam) | 05 (seam) | 0.014125 | 0.014155 | -0.000030 |
| All pairs | All pairs | 0.360287 | 0.153382 | 0.206906 |

On all 2,001 sampled values of `T` from `2.688` to `4.688`, these two matched pairs have the two largest decreases. Their combined decrease accounts for `88.68–101.38%` of the net change. Above 100%, the other pairs partially offset the decrease. Grouping all pairs at shortest distance two on the linked ring gives `75.36–88.71%` of the net change over the same window, without requiring a one-to-one match.

### 20.3 Waveforms and the signed cross term

For each matched pair and Pauli channel, write `I = U + V + W` with `U = integral C_0.5^2 dt`, `V = integral C_0^2 dt`, and `W = -2 integral C_0.5 C_0 dt`. Summed over both matched pairs and three channels, the six-qubit minus four-qubit changes are

| Weak-link square `Delta U` | Open-chain square `Delta V` | Cross term `Delta W` | Squared difference `Delta I` |
|---:|---:|---:|---:|
| -0.848591 | -0.317507 | +0.981529 | -0.184569 |

Both branch squares shrink by about 92%. The cross term becomes less negative overall and partially offsets their decrease. For `02 -> 04` in XX, `W` changes from `-0.54018` to `+0.00244`. For `13 -> 15`, YY and ZZ together account for about 85% of the decrease in `I`. The waveform report gives all 12 signed decompositions and the complete 24-trace data.

### 20.4 Short-time series and the longer open-chain path

At `epsilon=0`, the open-chain path for each matched pair changes from two edges to four. The linked branch retains a two-edge path through the seam. Exact rational arithmetic, with the product of the one-site expectations subtracted as a series, gives the first nonzero powers of the connected correlations.

| Size and branch | XX | YY | ZZ |
|:---|---:|---:|---:|
| N=4, linked or open | 6 | 8 | 10 |
| N=6, linked | 6 | 8 | 10 |
| N=6, open | 14 | 16 | 18 |

For the XX branch difference `C_0.5 - C_0`, the leading `t^6` coefficients agree across sizes, `-49/1800` for `02 -> 04` and `+49/1800` for `13 -> 15`. Cutting the central bond in the four-qubit linked graph verifies the separate two-edge route contributions to this coefficient. The finite-window XX integral difference therefore depends on higher-order waveform behavior. The YY and ZZ branch-difference coefficients are smaller in magnitude for the six-qubit matched pairs.

The exact leading terms were checked against a separate Pauli-commutator calculation and finite-time Chebyshev propagation. At `t=0.0125`, the largest relative discrepancy `|C(t)/(a t^p) - 1|` was `1.69e-4` across all 24 cases. The report lists the exact coefficient of every leading term.

### 20.5 Data and reproduction

After running the three §19 stages in order, run `python -B run.py` from each of these experiment folders in the order listed. The first three require NumPy 2.0 or later. The short-time script uses the Python standard library.

- Six-qubit symmetry and response. [Report](experiments/endpoint-field-six-qubit-2026-11-09/REPORT.md) · [Results](experiments/endpoint-field-six-qubit-2026-11-09/results.json)
- Pair weights and time-window scan. [Report](experiments/seam-pair-time-window-2026-11-09/REPORT.md) · [Results](experiments/seam-pair-time-window-2026-11-09/results.json)
- Branch waveforms and cross terms. [Report](experiments/seam-neighbor-waveforms-2026-11-09/REPORT.md) · [Results](experiments/seam-neighbor-waveforms-2026-11-09/results.json) · [Waveforms CSV](experiments/seam-neighbor-waveforms-2026-11-09/waveforms.csv)
- Exact short-time series. [Report](experiments/seam-neighbor-short-time-2026-11-09/REPORT.md) · [Results](experiments/seam-neighbor-short-time-2026-11-09/results.json)

## 21. Eight-qubit check of the distance rule (2026-11-09)

At the same `Omega=1`, `J=0.7`, and `epsilon=0`, the eight-qubit initial state is `|10000000>`. The matched pairs extend from `02 -> 04 -> 06` and `13 -> 15 -> 17`. The first nonzero powers of the connected correlations are the same for both pairs in each row.

| Size and branch | Shortest path in edges | XX | YY | ZZ |
|:---|---:|---:|---:|---:|
| N=4, open | 2 | 6 | 8 | 10 |
| N=6, open | 4 | 14 | 16 | 18 |
| N=8, open | 6 | 22 | 24 | 26 |
| N=4, 6, 8, linked | 2 | 6 | 8 | 10 |

For the branch difference `C_0.5 - C_0`, the `t^6` coefficients in XX remain `-49/1800` and `+49/1800` across all three sizes. Between `N=6` and `N=8`, the leading coefficients of the branch difference also agree in YY and ZZ.

An exact scan through order 30 covered all 28 pairs in each eight-qubit branch. At shortest distance `d >= 2`, the leading orders are `(4d-2, 4d, 4d+2)` in every case. Nearest-neighbor pairs have orders `(4,4,6)`. For the two matched pairs on the open chain, the report also derives closed-form XX coefficients and proves the leading order `4d-2` for every even `N >= 4` with the same one-excitation initial state. The all-axis distance rule is verified here at finite size.

Run `python -B run.py` in the eight-qubit experiment folder. The script uses the Python standard library and the preceding short-time experiment's code and results. [Report and derivation](experiments/seam-neighbor-eight-qubit-2026-11-09/REPORT.md) · [Exact coefficients and numerical checks](experiments/seam-neighbor-eight-qubit-2026-11-09/results.json)

## 22. Open-chain distance law, independent check, and prior work (2026-11-09)

For the open chain with uniform nonzero transverse field and bond strength and the initial state `|10...0>`, every pair at distance `d >= 2` has first nonzero connected-correlation orders `(4d-2, 4d, 4d+2)` in XX, YY, and ZZ, respectively. The proof applies to both even and odd chain sizes. It extends the XX calculation in §21 by counting the minimum Majorana hops that produce the same nonempty spin-flip interval on both sides of the connected correlator. The signed contributions at that order sum to a nonzero coefficient. The report gives a finite sum for any such pair and closed forms for `(0,N-2)` and `(1,N-1)`. The linked branch is covered by the finite-size results in §21.

Exact Majorana commutators verify the matched-pair formulas at `N=4,...,12`. For `N=4,...,9`, every computed matched-pair coefficient through `4d+2` also agrees with the direct state-vector expansion in the computational basis. The direct expansion independently confirms the leading terms for all 84 pairs at distance at least two across open chains of sizes `N=3,...,9`.

The longitudinal exponent `4d+2` has direct precedent in the short-time Rydberg density correlations of Lienhard et al. Their first three leading coefficients agree with a direct series calculation after specializing their equations to a uniform static ring with an all-aligned initial state. The literature comparison separates that precedent and the known Wick formulas from the coefficients for the single flipped endpoint in the open chain. Priority for the three-axis formulas remains undetermined. [Lienhard et al. (2018), §V.B and Appendix E](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.8.021070)

Run `python -B verify.py` in the distance-law folder to repeat the exact checks. [Proof and literature comparison](experiments/open-chain-distance-law-2026-11-09/REPORT.md) · [Verification script](experiments/open-chain-distance-law-2026-11-09/verify.py)

## 23. Weak-ring distance law and antipodal cancellation (2026-11-09)

For the ring with `lambda=1/2`, uniform nonzero `Omega` and `J`, and initial state `|10...0>`, the leading connected-correlation orders are `(4d-2, 4d, 4d+2)` for all `N >= 3` and ring distances `d >= 2`. Adjacent pairs have orders `(4,4,6)`. The proof uses an analytic bound at `d >= 400`, exhaustive integer checks of the shorter-distance configurations, and exact small-ring cases. A separate state-vector calculation verifies 855 channel cases through `N=12`.

At antipodal separation, the leading coefficient is the sum of the two shortest-route contributions, including the immediate exterior bonds for YY and ZZ. The four-site antipodal YY and ZZ coefficients also receive a full-cycle term. At `N=8`, pair `(0,4)` has opposite ZZ route contributions at `lambda=1/2`, with ratio `-85/86` for the seam route relative to the ordinary route.

At `Omega=1, J=0.7`, tuning the seam to `(151+32 sqrt(39))/745` raises the leading YY order for `(0,4)` from 16 to 18. Tuning it to `(39+16 sqrt(11))/185` raises the ZZ order from 18 to 20. Both cancellations were checked using exact coefficient polynomials. The antipodal XX coefficient is nonzero throughout `0 < lambda <= 1`.

Run `python -B coefficients.py` for the finite part of the proof and `python -B verify.py` for the independent Taylor checks. [Proof and cancellation examples](experiments/weak-ring-distance-law-2026-11-09/REPORT.md) · [Exact Taylor verification](experiments/weak-ring-distance-law-2026-11-09/verify.py)

## 24. Antipodal cancellation roots for all even sizes (2026-11-09)

For even `N=2d >= 4` and the same one-flipped-site initial state, the candidate YY and ZZ coefficients at orders `2N` and `2N+2` are quadratic in the seam strength. The coefficient formulas give the following classification in `0 < lambda < 1`.

| Antipodal pair | Sizes with one root per axis |
|:---|:---|
| `(0,d)` | Even `N >= 8` |
| `(d-1,N-1)` | Even `N >= 4` |
| Other antipodal pairs | None |

Each cancelling pair has distinct simple roots with `lambda_YY < lambda_ZZ`. Each series strictly decreases as `N` increases. For either position the asymptotic roots are `(N+2)/2^(N+1)` for ZZ and `(3/4)(N+2)/2^(N+1)` for YY. At the same size and axis, the site-0 pair has the larger root whenever both exist. The antipodal XX coefficient is nonzero throughout `0 < lambda <= 1`. Its leading contribution uses only the bonds on the two routes, while YY and ZZ also depend on the immediately exterior bonds.

| `N` | Pair | `lambda_YY` | `lambda_ZZ` |
|---:|:---|---:|---:|
| 6 | `(2,5)` | 0.050191310 | 0.062533422 |
| 8 | `(0,4)` | 0.470926089 | 0.497654036 |
| 8 | `(3,7)` | 0.013492741 | 0.016852776 |
| 10 | `(0,5)` | 0.012253477 | 0.016005605 |

The general formulas match all 40 YY and ZZ antipodal polynomials from an exact state-vector calculation at `N=4,6,8,10,12`. Run `python -B formulas.py` for root values and exact finite comparisons, and `python -B verify.py` for the state-vector calculation. [General formulas and proof](experiments/antipodal-cancellation-2026-11-09/REPORT.md) · [Exact verification](experiments/antipodal-cancellation-2026-11-09/verify.py)

## 25. Complementary occupations and signed site-0 rotations (2026-11-10)

This extension starts from $\lvert1000\rangle$ with $N=4$ and $\Omega=1$. The open chain has bonds $01,12,23$. The full ring adds $30$ with the same $J/2$ coefficient as the other ZZ bonds. A local rotation changes the preparation of site 0 before the time-independent evolution. The residual is $\delta=R_0+R_k-1$, with $k=3$ for the open chain and $k=2$ for the ring.

Without rotation, the pair sum is exactly one at all times. The open-chain proof uses the linear Majorana evolution at either endpoint. The ring proof combines a half-turn, the product of all Z operators, an even-site X flip, and complex conjugation. Removing the spectral DC term gives opposite analytic signals for the two sites, so their phases differ by $\pi$ wherever the signals are nonzero. [Occupation and phase measurements](experiments/phase-alignment-2026-11-10/REPORT.md) · [Exact complementarity proofs](experiments/antiphase-constraints-2026-11-10/REPORT.md)

Write $A=\lvert1000\rangle$, $B=\lvert0000\rangle$, $S(t)=e^{itH}(Z_0+Z_k)e^{-itH}$, and $C(t)=\langle A|S(t)|B\rangle$. In the preparation $\cos(\theta/2)A+e^{i\phi}\sin(\theta/2)B$, X rotation has $\phi=-\pi/2$ and Y rotation has $\phi=\pi$. Complementarity gives $\langle A|S(t)|A\rangle=0$, so the signed-angle decomposition is exact.

~~~text
delta_even = -sin(theta/2)^2 <B|S(t)|B>/2
delta_odd  = sin(theta) L_phi(t)
L_phi(t)  = -(Re C(t) cos(phi) - Im C(t) sin(phi))/2
~~~

X selects the imaginary part with coefficient $-1/2$ and Y selects the real part with coefficient $+1/2$. The leading short-time terms in $L_X,L_Y$ are $-\Omega t/2$ and $d_0\Omega Jt^2/4$, respectively, where $d_0=1$ for the chain and $2$ for the ring. At $J=0$, the Y odd component vanishes exactly while the X odd component persists. At $J=0.7$, $\lvert\theta\rvert=0.20$, and $0\le t\le3.688$, the sampled maximum absolute odd residual for X and Y is $0.082660$ and $0.053852$ on the chain, and $0.060063$ and $0.076201$ on the ring.

The rotation comparison reports phase concentration separately from the occupation-sum residual. On the chain at $\theta=+0.20$ in the short window, X gives a higher phase-locking value (0.999997 versus 0.996181) while Y gives a mean absolute phase difference closer to $180^\circ$ ($171.98^\circ$ versus $169.05^\circ$). Both axes were evaluated on the same time samples where all four sites passed the envelope threshold. [Y rotation](experiments/antiphase-y-rotation-2026-11-10/REPORT.md) · [Signed Y rotation](experiments/antiphase-y-signed-2026-11-10/REPORT.md) · [X/Y comparison and derivation](experiments/antiphase-xy-rotation-2026-11-10/REPORT.md) · [Signed residuals](experiments/antiphase-xy-rotation-2026-11-10/residuals.csv) · [Phase metrics](experiments/antiphase-xy-rotation-2026-11-10/phases.csv)

## 26. Finite-window rotation axes and later-window validation (2026-11-10)

For a rotation axis at angle $\beta$ from X toward Y in the XY plane, $\phi=-\pi/2-\beta$. Let $f=(L_X,L_Y)^{\mathsf T}$, $v=(\cos\beta,\sin\beta)^{\mathsf T}$, and $G_W=\langle ff^{\mathsf T}\rangle_W$ on a time window $W$. The time RMS of the angular slope is $\sqrt{v^{\mathsf T}G_Wv}$. The eigenvector of the smaller eigenvalue gives the optimal axis modulo $180^\circ$. The same axis minimizes the odd residual RMS at the tested nonzero rotation angles. It also minimizes the total RMS pooled over the two signs because pooling the mean squares cancels the even–odd cross term. A single signed angle retains that term.

At $J=0.7$, the axes selected on $[0,3.688]$ were $58.584^\circ$ for the chain and $33.015^\circ$ for the ring. The later-window test uses the saved full-precision angles throughout. The baseline in each window is the pure X or Y axis with the lower slope RMS.

| Topology | Evaluation window | Fixed-axis slope RMS | Better pure axis | Reduction versus that axis |
|:---|:---|---:|:---|---:|
| Open | $[0,3.688]$ | 0.068484 | Y | +59.19% |
| Open | $(3.688,8\pi]$ | 0.149640 | Y | +3.08% |
| Open | $(8\pi,16\pi]$ | 0.181878 | Y | −4.36% |
| Ring | $[0,3.688]$ | 0.097438 | X | +41.89% |
| Ring | $(3.688,8\pi]$ | 0.105333 | X | −9.06% |
| Ring | $(8\pi,16\pi]$ | 0.133724 | X | −16.69% |

**The short-window advantage did not persist.** The fixed axis beat both pure axes only in the first later window on the chain. The same win/loss pattern held for the total residual at each tested sign and magnitude ($\pm0.10,\pm0.20$). The reports also separate occupation amplitude from cancellation between the two site traces. In the first later chain window at $\theta=-0.20$, the fixed axis widened the gap between the two site standard deviations relative to X. Stronger anticorrelation contributed to a small reduction in total RMS, from $0.030858$ to $0.030664$. The separately fitted later-window axes serve as diagnostics for how far the fixed axes are from the window minima. [Axis sweep and optimization](experiments/antiphase-axis-optimization-2026-11-10/REPORT.md) · [Later-window report](experiments/antiphase-axis-holdout-2026-11-10/REPORT.md) · [RMS comparison](experiments/antiphase-axis-holdout-2026-11-10/rms-comparison.png) · [Fixed and diagnostic axes](experiments/antiphase-axis-holdout-2026-11-10/regret.csv)

## 27. Long-time limit of the occupation-residual Gram matrix (2026-11-10)

At the same $J=0.7,\Omega=1$, the spectral calculation adds the coefficients of exactly equal signed energy gaps before squaring. Integer characteristic polynomials and rational isolating intervals distinguish equal gaps from close but distinct ones. Symmetry leaves four positive frequencies in the open-chain slope signals and ten in the ring slope signals. $L_X$ is a sum of sines and $L_Y$ a sum of cosines, so their product has zero infinite-time mean, including the terms with equal frequency. The two diagonal entries were independently evaluated in rational arithmetic.

| Topology | Exact $G_{\infty,XX}$ | Exact $G_{\infty,YY}$ | Optimal axis | Minimum slope RMS | Reduction versus the better pure axis |
|:---|---:|---:|:---|---:|---:|
| Open | $318475/9741742$ | $758225/29225226$ | Y | 0.161072022 | 0% |
| Ring | $3290125/236511872$ | $5858025/236511872$ | X | 0.117945052 | 0% |

The limiting eigenvalue gaps are $0.006748$ on the chain and $0.010857$ on the ring. The limiting axes are unique modulo $180^\circ$. The chain has the smaller eigenvalue separation and its optimal direction is more sensitive to a perturbation of comparable size in the Gram matrix.

| Cumulative endpoint $T$ | Open axis / reduction | Ring axis / reduction |
|:---|:---|:---|
| $8\pi$ | $56.35^\circ$ / 10.12% | $20.44^\circ$ / 8.67% |
| $16\pi$ | $68.74^\circ$ / 2.38% | $11.06^\circ$ / 2.53% |
| $32\pi$ | $79.16^\circ$ / 0.52% | $12.76^\circ$ / 2.97% |
| $64\pi$ | $83.84^\circ$ / 0.14% | $12.43^\circ$ / 2.28% |

The ring's finite-window axis and reduction vary nonmonotonically. The difference between the distinct ring frequencies $0.34698748$ and $0.36127247$ gives a beat period of about $439.85$, and this difference-frequency term accounts for about 87% of $G_{XY}(64\pi)$. On the prespecified later window $(32\pi,64\pi]$, the axes chosen in the infinite-time limit were held fixed. The open-chain Y RMS was $0.160511075$ versus X at $0.176888200$. The ring X RMS was $0.119572545$ versus Y at $0.155916812$. In both cases the fixed axis remained the better pure axis on that window. [Long-time derivation and audit](experiments/antiphase-axis-long-time-2026-11-10/REPORT.md) · [Convergence figure](experiments/antiphase-axis-long-time-2026-11-10/long-time-convergence.png) · [Finite-window results](experiments/antiphase-axis-long-time-2026-11-10/convergence.csv) · [Fixed-axis results](experiments/antiphase-axis-long-time-2026-11-10/late_holdout.csv)

### Reproduction

With NumPy 2 or later, run `python -B run.py` in each of the four folders below in order. The later-window validation reads the full-precision short-window axes from the optimization folder's `optima.csv`. Run `python -B plot.py` in each folder with Matplotlib to regenerate its figures.

1. [Signed X/Y rotation](experiments/antiphase-xy-rotation-2026-11-10/run.py)
2. [Finite-window axis optimization](experiments/antiphase-axis-optimization-2026-11-10/run.py)
3. [Fixed-axis later-window validation](experiments/antiphase-axis-holdout-2026-11-10/run.py)
4. [Long-time limit and cumulative windows](experiments/antiphase-axis-long-time-2026-11-10/run.py) · [Exact gap and moment calculations](experiments/antiphase-axis-long-time-2026-11-10/spectral_groups.py)

## 28. Uniform-coupling crossover in connected correlations (2026-11-10)

With $N=4$, $\Omega=1$, and $\lvert0000\rangle$, both branches use the same coefficient $J/2$ on every existing ZZ bond. On $[0,3.688]$, rank each pair by the time RMS of its three-axis connected-correlation difference. For $1\le J\le2$, the largest pair changes from the closing bond $03$ to the central bond $12$ near $J=1.26765$, then to the tied outer bonds $01,23$ near $J=1.57203$. At the second crossover, the outer bonds' larger ZZ contribution offsets their smaller XX and YY contributions. The crossing shifts with the observation window. [Uniform-coupling sweep](experiments/uniform-coupling-sweep-2026-11-10/REPORT.md) · [Crossover scan](experiments/uniform-coupling-crossover-2026-11-10/REPORT.md) · [Time-window decomposition](experiments/crossover-time-window-2026-11-10/REPORT.md)

## 29. Four-site in-phase occupations and long-window validation (2026-11-10)

For the time-independent model at $N=4$, $\Omega=1$, and $J=0.7$, the search allows independent local Bloch angles for all four sites. For each window, all sites must have an occupation half-range of at least $g$, and all four analytic-signal envelopes must reach $g$ simultaneously for at least 90% of the window. On those times, all six pairs must have a circular mean phase difference within $\pm15^\circ$ and a phase-locking value (PLV) of at least 0.95. The tested floors are $g=0.025,0.05,0.10$. [Four-site search and baselines](experiments/four-site-in-phase-2026-11-10/REPORT.md)

On the ring, identical local preparations give identical occupation traces by translation symmetry. The first open-chain search chose a state on $[0,3.688]$ whose PLV and envelope coverage fell below threshold on the disjoint window $[3.688,8\pi]$. A follow-up search selected a mirror-symmetric product on eight $8\pi$ windows over $[0,64\pi]$, then held it fixed for 24 new windows over $[64\pi,256\pi]$. [Long-window report](experiments/open-chain-phase-persistence-2026-11-10/REPORT.md) · [Frozen initial angles](experiments/open-chain-phase-persistence-2026-11-10/selection.json) · [Occupation and envelope figure](experiments/open-chain-phase-persistence-2026-11-10/occupations-envelope.png)

| Fixed open-chain preparation, $g=0.05$ | Held-out windows passed | Largest absolute mean lag | Lowest PLV | Lowest window coverage |
|:---|---:|---:|---:|---:|
| Selected on the short window | 9/24 | $7.602^\circ$ | 0.939657 | 68.14% |
| Selected across eight windows | 24/24 | $0.201^\circ$ | 0.999365 | 100% |

The eight-window state also passed all 24 windows at the other two floors. The minimum envelope on the held-out interval was 0.241902. Its Bloch angles $(\theta,\phi)$ in radians are $(1.028786,-0.012331)$ at sites 0 and 3, and $(0.848414,-0.010693)$ at sites 1 and 2. The selection file stores the full precision.

## 30. Spectral mechanism and angle-error controls (2026-11-10)

Grouping equal energy gaps yields 20 positive frequencies across the four occupation observables. For the frozen open-chain state, the shared line $\omega_*=0.502392496$ carries 99.879% of the oscillatory power at the ends and 99.826% at the interior sites. Its amplitudes are 0.255065 and 0.336237, with an end–interior phase difference of $-0.00177^\circ$. The sum of the other amplitudes gives a numerical envelope lower bound of 0.241788 across sites. A phase bound using those residual magnitudes independently gives $8.081^\circ$, so the shared line is strong enough to account for persistent phase agreement. The weaker lines set the instantaneous fluctuations. The interior-only line at $1.761139$ accounts for about 59.4% of their first-order phase variance in the frozen state. [Spectral derivation](experiments/open-chain-phase-spectrum-robustness-2026-11-10/REPORT.md) · [Frequency coefficients](experiments/open-chain-phase-spectrum-robustness-2026-11-10/spectrum.csv) · [Spectrum and line removal](experiments/open-chain-phase-spectrum-robustness-2026-11-10/spectrum.png)

The fixed-state error controls include 625 preparations on a mirror-symmetric angle grid and 32 single-site perturbations, with angle offsets of up to 0.05 rad. All 657 preparations pass all 24 windows at all three floors. The worst mean phase difference is $0.919^\circ$, the lowest PLV is 0.993658, and the minimum envelope is 0.215882. [Angle-error results](experiments/open-chain-phase-spectrum-robustness-2026-11-10/robustness_summary.csv)

## 31. Independent angle errors and the first PLV failure (2026-11-10)

The independent-error sweep applies fixed directions to all eight Bloch angles without retuning the preparation. Each width uses the same 240 reflection-breaking sign corners and 256 independent uniform draws. The coarser sweep passed through 0.075 rad and found six failing directions at 0.100 rad. The refined sweep evaluates all 496 directions at each of 26 widths from 0.075 to 0.100 rad in 0.001 rad steps, using the same 24 windows and three floors. [Independent-error sweep](experiments/open-chain-independent-angle-errors-2026-11-10/REPORT.md) · [Refined PLV scan](experiments/open-chain-plv-onset-2026-11-10/REPORT.md) · [All 26 widths](experiments/open-chain-plv-onset-2026-11-10/width_summary.csv)

Every direction passes through 0.090 rad. At **0.091 rad**, four directions each fail one window, always on endpoint pair $03$. All three floors give the same result and the common envelope coverage is 100%.

| Direction IDs | Failing window | PLV at 0.090 rad | PLV at 0.091 rad |
|:---|:---|---:|---:|
| `c090`, `c165` | $[104\pi,112\pi]$ | 0.950412 | 0.949324 |
| `c015`, `c240` | $[152\pi,160\pi]$ | 0.950633 | 0.949551 |

Removing the $0.930636$ line produces the largest PLV recovery, to 0.997769 and 0.997730 in the two failing windows, and restores all 24 windows for the four directions. Keeping only the carrier and $0.930636$ also passes, with PLVs 0.954651 and 0.955004 in those windows. The weaker endpoint lines at $1.332896$ and $1.604652$ contribute to the threshold crossing. The boundary describes the tested directions and width grid. [First failing conditions](experiments/open-chain-plv-onset-2026-11-10/first_witnesses.csv) · [PLV onset figure](experiments/open-chain-plv-onset-2026-11-10/plv-onset.png) · [Line-removal figure](experiments/open-chain-plv-onset-2026-11-10/onset-ablations.png)

### Reproducing the in-phase results

Use NumPy 2 or later and Matplotlib. Each command runs from its experiment folder. The frozen selections and compact outputs are included with the code. The three larger arrays listed below can be regenerated by the corresponding command.

| Experiment | Commands in order | Large generated output |
|:---|:---|:---|
| [Four-site search](experiments/four-site-in-phase-2026-11-10/REPORT.md) | `python -B baseline.py`, `python -B run.py`, `python -B plot.py` | |
| [Long-window validation](experiments/open-chain-phase-persistence-2026-11-10/REPORT.md) | `python -B run.py evaluate`, `python -B plot.py` | `waveforms.csv` from `evaluate` |
| [Spectral mechanism and angle controls](experiments/open-chain-phase-spectrum-robustness-2026-11-10/REPORT.md) | `python -B spectrum.py`, `python -B majorana.py`, `python -B robustness.py`, `python -B plot.py` | |
| [Independent errors](experiments/open-chain-independent-angle-errors-2026-11-10/REPORT.md) | `python -B sweep.py`, `python -B diagnostics.py`, `python -B plot.py` | `window_metrics.csv.gz` from `sweep.py` |
| [Refined PLV onset](experiments/open-chain-plv-onset-2026-11-10/REPORT.md) | `python -B scan.py`, `python -B diagnostics.py`, `python -B plot.py` | `all_windows.npz` from `scan.py` |

To repeat the eight-window selection, run `python -B run.py select` in the long-window validation folder before `evaluate`. The published `selection.json` is the reference for the downstream perturbation experiments.

## 32. Finite-shot identification of all 64 four-qubit graphs (2026-11-13)

The comparison extends to all $2^6=64$ labeled graphs on four sites. Each independent shot starts in $\lvert0000\rangle$ and evolves under $H_G=\tfrac12\sum_i X_i+\tfrac7{20}\sum_{ij\in G} Z_i Z_j$ for a time in $[0,3.688]$. The shot measures all four sites in the same X, Y, or Z basis. Each readout bit flips independently with known probability $p$. The schedule fixes the settings and shot counts in advance. Success requires correctly recovering all six edge indicators with probability at least 95% for each true graph.

| Readout flip probability | Lower bound on required shots | Certified sufficient shots |
|:---|---:|---:|
| 0% | 17 | 97 |
| 2% | 22 | 130 |
| 5% | **53** | **99** |

At $p=5\%$, the 99-shot design measures Z at $t=2.122$ for 88 shots and at $t=3.531$ for 11. The decoder sums a fixed bias and the integer scores for each observed bit string, then returns the unique highest-scoring graph. Counting ties as failures, the worst-case error probability is at most **4.999064%**. The certificate uses rational probability bounds and integer convolution to account for overlaps between competing graphs' error events. [99-shot design, decoder, and certificates](experiments/Fable-Codex-2026-11-13/Codex/four_qubit_p5_99_audit.zip)

The lower bound covers the full continuous time window. With 52 shots, if the true graph lacks only edge $23$, the probabilities of outputting six particular rival graphs must sum to at least **5.2736919284%** to meet the 95% success requirement on those rivals. Those six outputs are mutually exclusive and already exceed the 5% error allowance on the true graph. The same bound rules out fewer shots because a decoder could ignore any extra observations. [Derivation and replay commands (Japanese)](experiments/Fable-Codex-2026-11-13/REPORT.md)

The experiment folder also contains the original graph-enumeration and finite-shot search scripts, Monte Carlo counts, and earlier certificates. [Experiment files](experiments/Fable-Codex-2026-11-13) · [Monte Carlo counts by graph](experiments/Fable-Codex-2026-11-13/Codex/graph_success_counts.csv)
