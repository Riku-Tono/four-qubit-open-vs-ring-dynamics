# Four-qubit direct open-vs-closed comparison

This package compares two fixed topologies from the same initial state, with no switch or quench.

It also contains a local-response audit and a same-observable reanalysis of the existing t_switch=1.688 quench artifacts.

## Fixed conditions

- N = 4
- Omega = 1.0
- J = 0.7
- initial state = `|0000>`
- time range = `0 <= t <= 3.688`
- open bonds = `(0,1), (1,2), (2,3)`
- closed bonds = open bonds plus `(3,0)`
- symmetric second-order Trotter steps `dt=0.002` and `dt=0.001`

The open and closed states both start at time zero. This is a topology comparison, separate from the existing `t_switch=1.688` local-quench experiment.

## Run

```text
cargo fmt --all -- --check
cargo test --release --offline
cargo run --release --offline
python analyze_modes.py
python local_response_analysis.py
python two_body_correlation_analysis.py
```

The Python modal audit uses exact diagonalization. Degenerate eigenspaces are grouped before ranking, so weights do not depend on an arbitrary eigenvector basis inside a degenerate space.

## Local-response definition

The older quench source did not contain a site-resolved R_i; it stored average magnetization, bond correlation, and energy quantities. This package therefore introduces one common definition for both fixed and quench analyses:

O_i = n_i = |1_i><1_i| = (I-Z_i)/2

R_i(t) = <n_i>

R_i is the probability that site i is in state 1, so it lies in [0,1] and is global-phase invariant. This is a newly fixed common definition, not a claim that a prior site-response implementation existed.

For nonzero total response, p_i=R_i/sum_j R_j. Saved distance columns include raw L1/L2 and normalized L1/L2; the primary distribution-shape distance in the report is normalized L1. Total variation is one half of normalized L1.

- IPR = sum_i R_i^2 / (sum_i R_i)^2
- PR = 1/IPR

At t=0, where all responses vanish, normalized p_i is stored as zero and PR/IPR are undefined (NaN).

Maximum-site ties use max(1e-10, 1e-8 * abs(R_max)) as the absolute/relative tolerance. Every site within that tolerance of the maximum is recorded. The same threshold is used for the near_tie classification.

Site labels are 0,1,2,3 in both topologies. Closed has periodic boundary bond (3,0).

## Two-body correlation definition

For every `i<j`, the package saves all one-body expectations `<X_i>`, `<Y_i>`, `<Z_i>` and raw two-body expectations `<X_iX_j>`, `<Y_iY_j>`, `<Z_iZ_j>`. The connected quantities are

`C^AA_ij = <A_i A_j> - <A_i><A_j>`, for `A=X,Y,Z`.

Raw and connected values are stored separately. The displayed 4 by 4 connected-correlation maps use zero on the diagonal; the diagonal is not interpreted as the on-site variance. For each channel, `upper_L1` is the sum of absolute open/closed differences over the six `i<j` entries and `upper_L2` is their Euclidean norm. The symmetric zero-diagonal matrix Frobenius distance is `sqrt(2)*upper_L2`. The combined distance is

`D_corr = sqrt(D_XX_upper_L2^2 + D_YY_upper_L2^2 + D_ZZ_upper_L2^2)`.

Feature distances are ordinary unnormalized Euclidean distances over, respectively: 12 one-body values, six ZZ-connected values, all 18 connected values, and the concatenation of the 12 one-body plus 18 connected values. These are finite-observable feature distances, not reconstructions of the full state distance.

The first nonzero map event uses `D_corr > 1e-8`. Mode analysis groups degenerate energies and uses eigenspace projections. Raw two-body and one-body terms are decomposed linearly; connected correlations are then reconstructed as raw minus the product of the reconstructed one-body terms. Because that product is nonlinear, the connected-ranking score is a leave-one-eigenspace-pair-out RMS influence, and is not claimed to be a unique additive decomposition.

## Outputs

- `direct_open_closed_dt002.csv`
- `direct_open_closed_dt001.csv`
- `direct_open_closed_dt_comparison.csv`
- `direct_open_closed_summary.csv`
- `direct_open_closed_mode_weights.csv`
- `direct_open_closed_top6_pair_phases.csv`
- `direct_open_closed_exact_validation.csv`
- `direct_open_closed_pauli_dt002.csv`
- `direct_open_closed_pauli_dt001.csv`

Additional local-response outputs:

- direct_open_closed_local_response_dt002.csv
- direct_open_closed_local_response_dt001.csv
- local_response_timeseries.csv
- local_response_max_position_timeseries.csv
- local_response_distribution_distance.csv
- state_response_comparison.csv
- local_response_key_times.csv
- local_response_mode_pair_contributions.csv
- dominant_mode_pair_summary.csv
- local_response_mode_reconstruction_summary.csv
- local_response_exact_validation.csv
- local_response_dt_convergence.csv
- quench_local_response_timeseries.csv
- fixed_vs_quench_response_comparison.csv
- state_response_correlations.csv
- local_response_numerical_audit.csv
- LOCAL_RESPONSE_ANALYSIS_REPORT.md

Two-body correlation outputs:

- one_body_pauli_expectations.csv
- two_body_raw_pauli_correlations.csv
- two_body_connected_correlations.csv
- correlation_distance_timeseries.csv
- state_response_correlation_comparison.csv
- correlation_distance_correlations.csv
- representative_correlation_maps.csv
- correlation_symmetry_audit.csv
- correlation_mode_pair_amplitudes.csv
- correlation_mode_pair_contributions.csv
- correlation_dominant_mode_pair_summary.csv
- correlation_mode_reconstruction_summary.csv
- correlation_exact_validation.csv
- correlation_dt_convergence.csv
- correlation_numerical_audit.csv
- CORRELATION_ANALYSIS_REPORT.md

## Scope

This is a finite N=4, finite-time numerical comparison. It does not establish a general boundary mechanism, thermodynamic-limit result, phase transition, transport law, or correspondence with the N=7 Lindblad backflow.
