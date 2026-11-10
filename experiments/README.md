# Additional experiments

[日本語](README.md)

This index groups the experiments from November 9 and 10, 2026 by topic. The reports describe the calculation settings, results, validation, and saved data. The four-qubit baseline and the August extensions each have a reproducibility ZIP on their release page. [July baseline](https://github.com/Riku-Tono/four-qubit-open-vs-ring-dynamics/releases/tag/v0.1.0-n4-fixed) · [August extensions](https://github.com/Riku-Tono/four-qubit-open-vs-ring-dynamics/releases/tag/v0%2C2)

- [Response to transverse fields at the endpoints](#response-to-transverse-fields-at-the-endpoints)
- [Waveforms and short-time expansions near the seam](#waveforms-and-short-time-expansions-near-the-seam)
- [Distance laws and cancellation in connected correlations](#distance-laws-and-cancellation-in-connected-correlations)
- [Uniform coupling strength and correlation differences](#uniform-coupling-strength-and-correlation-differences)
- [Antiphase pairs and local rotations](#antiphase-pairs-and-local-rotations)
- [In-phase oscillations and preparation-angle errors on four sites](#in-phase-oscillations-and-preparation-angle-errors-on-four-sites)

## How to run

Download or clone the entire repository. Some experiments use code or saved results from neighboring directories. Run each listed script in order from the corresponding experiment directory, with `python -B` before the script name and any arguments. The scripts regenerate result files, so use a working copy if you want to refer back to the archived results.

To verify the open-chain distance law, start from the repository root and run the following commands.

```sh
cd experiments/open-chain-distance-law-2026-11-09
python -B verify.py
```

Most numerical experiments use NumPy 2.x. The plotting scripts and the first three endpoint-field experiments also require Matplotlib. The tables identify the experiments that use the Python standard library. See each report for the calculation settings and required inputs.

## Response to transverse fields at the endpoints

| Report | What it examines | Scripts in order |
| --- | --- | --- |
| [Four-qubit endpoint fields and antiunitary symmetry](endpoint-field-2026-11-09/REPORT.md) | How the signs of the transverse fields at the endpoints affect correlation differences | `run.py` |
| [Opposite endpoint fields and scaling of the seam pair's share](endpoint-field-scaling-2026-11-09/REPORT.md) | How the seam pair's contribution varies with field strength | `run.py` |
| [Linear response of the seam pair's share for four qubits](endpoint-field-zero-slope-2026-11-09/REPORT.md) | Extrapolating the first-order response to a weak field | `run.py` |
| [Linear response of the seam pair's share for six qubits](endpoint-field-six-qubit-2026-11-09/REPORT.md) | Comparison with the four-qubit case and a check at two integration step sizes | `run.py` |

## Waveforms and short-time expansions near the seam

| Report | What it examines | Scripts in order |
| --- | --- | --- |
| [All-pair integrals for four and six qubits across time windows](seam-pair-time-window-2026-11-09/REPORT.md) | Which pairs contribute to the overall difference and how their contributions vary with the observation window | `run.py` |
| [Correlation waveforms for two pairs near the seam](seam-neighbor-waveforms-2026-11-09/REPORT.md) | Open-chain and ring waveforms and the contributions to the integrated squared difference | `run.py` |
| [Short-time expansions for two pairs near the seam](seam-neighbor-short-time-2026-11-09/REPORT.md) | How path length relates to the leading terms of the exact series | `run.py` (standard library) |

## Distance laws and cancellation in connected correlations

The scripts in this section use the Python standard library and saved experiment files.

| Report | What it examines | Scripts in order |
| --- | --- | --- |
| [Path length and leading orders of connected correlations for eight qubits](seam-neighbor-eight-qubit-2026-11-09/REPORT.md) | Leading orders for an eight-qubit open chain and ring | `run.py` |
| [Short-time distance law for open Ising chains](open-chain-distance-law-2026-11-09/REPORT.md) | How distance determines the leading orders of XX, YY, and ZZ | `verify.py` |
| [Short-time distance law for rings with a weak link](weak-ring-distance-law-2026-11-09/REPORT.md) | Coefficients from the two paths and the resulting distance law | `coefficients.py` → `verify.py` |
| [Weak-link strengths that cancel leading terms for antipodal pairs](antipodal-cancellation-2026-11-09/REPORT.md) | Conditions under which the leading YY or ZZ term vanishes in even-sized rings | `formulas.py` → `verify.py` |

## Uniform coupling strength and correlation differences

| Report | What it examines | Scripts in order |
| --- | --- | --- |
| [Uniform coupling strength and connected correlation differences between the open chain and ring](uniform-coupling-sweep-2026-11-10/REPORT.md) | How the differences vary with coupling strength for each pair and axis | `run.py` → `plot.py` |
| [Coupling strength at which the pair with the largest correlation difference changes](uniform-coupling-crossover-2026-11-10/REPORT.md) | Crossover between the outer pairs and the central pair | `run.py` → `plot.py` |
| [Time-window analysis of the crossover between outer and central pairs](crossover-time-window-2026-11-10/REPORT.md) | How the ranking of correlation differences depends on the observation end time | `run.py` → `plot.py` |

## Antiphase pairs and local rotations

| Report | What it examines | Scripts in order |
| --- | --- | --- |
| [Occupations, amplitudes, and phases from the initial state &#124;1000⟩](phase-alignment-2026-11-10/REPORT.md) | Phase differences and complementary occupations in the open chain and ring | `run.py` → `plot.py` |
| [Antiphase constraints for four qubits](antiphase-constraints-2026-11-10/REPORT.md) | The constraint that each relevant pair's occupations sum to one and the numerical residuals | `run.py` |
| [Y-axis rotation at site 0 and the antiphase constraint](antiphase-y-rotation-2026-11-10/REPORT.md) | Residuals and phases after rotating the initial state | `run.py` → `endpoint_audit.py` |
| [Antiphase pair residuals under positive and negative Y rotations](antiphase-y-signed-2026-11-10/REPORT.md) | Components that are even and odd under reversal of the rotation angle | `run.py` → `plot.py` |
| [X- and Y-axis rotations of qubit 0](antiphase-xy-rotation-2026-11-10/REPORT.md) | How the residual of the occupation sum depends on the rotation axis | `run.py` → `plot.py` |
| [Rotation axis that minimizes antiphase pair residuals](antiphase-axis-optimization-2026-11-10/REPORT.md) | The optimal axis over a finite observation window | `run.py` → `plot.py` |
| [Testing the short-window rotation axis in later windows](antiphase-axis-holdout-2026-11-10/REPORT.md) | How the axis selected in a short observation window performs in subsequent windows | `run.py` → `plot.py` |
| [Long-time limit of the optimal rotation axis for antiphase pairs](antiphase-axis-long-time-2026-11-10/REPORT.md) | How the optimal axis changes as the observation window grows and what it approaches | `run.py` → `plot.py` |

## In-phase oscillations and preparation-angle errors on four sites

| Report | What it examines | Scripts in order |
| --- | --- | --- |
| [Searching for in-phase occupations using local rotations on four sites](four-site-in-phase-2026-11-10/REPORT.md) | A search for initial states whose site occupations oscillate in phase | `baseline.py` → `run.py select` → `run.py evaluate` → `plot.py` |
| [How long in-phase occupations persist on the open four-site chain](open-chain-phase-persistence-2026-11-10/REPORT.md) | Validation of the selected state in later observation windows | `run.py evaluate` → `plot.py` |
| [Frequencies of in-phase oscillations and robustness to angle errors](open-chain-phase-spectrum-robustness-2026-11-10/REPORT.md) | Frequency components and sensitivity to changes in the preparation angles | `spectrum.py` → `majorana.py` → `robustness.py` → `plot.py` |
| [In-phase oscillations under independent preparation-angle errors](open-chain-independent-angle-errors-2026-11-10/REPORT.md) | Phase alignment when the angle error differs from site to site | `sweep.py` → `diagnostics.py` → `plot.py` |
| [First tested error magnitude at which PLV falls below 0.95](open-chain-plv-onset-2026-11-10/REPORT.md) | Where the phase-locking value first falls below the threshold on the tested error grid | `scan.py` → `diagnostics.py` → `plot.py` |

The long-time in-phase validation uses the published `selection.json`. To repeat the search for the initial state, run `python -B run.py select` first in a separate working copy.
