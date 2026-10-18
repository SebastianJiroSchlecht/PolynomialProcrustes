# Reproducibility audit (18 October 2026)

Baseline: `c0dcaadfb8635a00d2136c67aa973f71a8f468d9`, with PEVD Toolbox pinned at `a8f4704beb362a1c6a0d2d82a8f80c60fbbe7cb3`. Paper: [Weiss, Schlecht and Moonen, IEEE TSP 74, 904–919 (2026)](https://doi.org/10.1109/TSP.2026.3663449); the [accepted manuscript](https://strathprints.strath.ac.uk/95550/) was inspected for figure coverage and experimental parameters.

## Confirmed source defects addressed

- `ProcrustesMetrics` subtracted `eye(3)` regardless of matrix size. This prevented the two-channel Figure 12 cache from being computed and broke the six-channel ensemble metrics.
- `ProcrustesRandomMatrixSingleTest` assigned `PMetrics` but read `PMetric`, causing an undefined-variable error. It also put the solver's final FFT length and convergence mismatch into the opposite columns from its own documented output contract.
- The advertised one-argument solver call referenced undefined `NOrd`. It now defaults to allpass order 32. The two-input dimension check compared A's row count to B's column count, incorrectly rejecting compatible rectangular matrices.
- The solver removed both spectral-zero bins and repeated-singular-value bins, but chose the partial-IDFT reconstruction only when the first kind was present. Reconstruction now checks the actual retained bin set.
- The large-matrix single test defined allpass order 64 but passed literal 16; it now passes its configured order. The large ensemble harness omitted dimension 96 and only requested five seeds per dimension; it now requests seven dimensions and 100 seeds each.
- The random ensemble harness requested 1,000 runs rather than the paper's 10,000. It now requests seeds 0–9999. A sixth optional solver argument exposes the convergence tolerance; the random single test uses the paper's 1e-4 threshold, while the solver default remains 1e-8.
- Recomputations previously appended to the shipped result files. They now write separate `*-recomputed.txt` files. Archived data are preserved, including their historical column ordering.
- Figure 3 wrote an unnumbered filename; Figure 20 supplied five tick labels for three tick positions. Both are corrected.

## Evidence gaps and limits

The shipped first ensemble contains 1,020 rows, 1,001 unique seeds, and 19 repeated seeds, rather than evidence for 10,000 independent instances. Computing its unmodified rows gives mean PU error 2.1908756e-10 (sample standard deviation 8.8158880e-11) and mean normalised LS excess 4.9761787e-4 (sample standard deviation 1.3136302e-4). These are close to the paper's rounded statistics but do not verify its claimed sample count. No duplicate rows have been silently removed, and no missing trials have been invented.

The second ensemble contains all 700 distinct (seed, dimension) pairs expected by the plotting script. Its archived results have not been regenerated during this audit. Historical computation times cannot be reproduced exactly on different hardware.

MAT-file caches in the author's local checkout are ignored by Git. The fresh-checkout workflow exercises the uncached paths. Figures 8 and 11 are explanatory diagrams, not supplied MATLAB plots. Successful script execution by itself does not establish agreement with every curve, theorem, or numerical claim in the paper. In particular, fixed index windows in Figures 12–16 and 19–20 deserve comparison with the paper if solver output lengths change.

The email mentioning `FileList.txt` and three third-party comparison-method folders matches a local polynomial Gram–Schmidt/QR project, not this repository. Those files were not imported into Procrustes. PEVD Toolbox has its own MIT license; redistribution rights for any other third-party implementations must be established from their actual licenses or permission text rather than inferred from permission to use them.

## Validation

Numerical regressions and all supplied plotting scripts are run by `.github/workflows/reproducibility.yml` in MATLAB R2024b with Signal Processing Toolbox. The workflow retains generated EPS/PNG plots and a reproduction log. Local MATLAB execution required MathWorks sign-in and did not reach script execution. Full ensemble regeneration is excluded from the automated run; workflow results are recorded in the pull request.
