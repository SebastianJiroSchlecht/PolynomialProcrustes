# Reproducibility audit (18 October 2026)

Baseline: `c0dcaadfb8635a00d2136c67aa973f71a8f468d9`, with PEVD Toolbox pinned at `a8f4704beb362a1c6a0d2d82a8f80c60fbbe7cb3`. Paper: [Weiss, Schlecht and Moonen, IEEE TSP 74, 904–919 (2026)](https://doi.org/10.1109/TSP.2026.3663449); the [accepted manuscript](https://strathprints.strath.ac.uk/95550/) was inspected for figure coverage and experimental parameters.

## Confirmed source defects addressed

- `ProcrustesMetrics` subtracted `eye(3)` regardless of matrix size. This prevented the two-channel Figure 12 cache from being computed and broke the six-channel ensemble metrics.
- `ProcrustesRandomMatrixSingleTest` assigned `PMetrics` but read `PMetric`, causing an undefined-variable error. It also put the solver's final FFT length and convergence mismatch into the opposite columns from its own documented output contract.
- The advertised one-argument solver call referenced undefined `NOrd`. It now defaults to allpass order 32. The two-input dimension check compared A's row count to B's column count, incorrectly rejecting compatible rectangular matrices.
- Allpass design formed a dense 4097-by-4097 diagonal weight matrix for every filter. Direct row weighting now applies the same equations without that allocation or dense matrix multiplication; the numerical regressions cover the weighted allpass design.
- The fresh 96-channel example failed on a zero crossing between the final and first frequency bins. Crossing detection now includes the circular boundary and refines each crossing within its actual adjacent retained-bin interval. The random-matrix generator likewise records its boundary crossing once, rather than once per interior bin.
- The solver removed both spectral-zero bins and repeated-singular-value bins, but chose the partial-IDFT reconstruction only when the first kind was present. Reconstruction now checks the actual retained bin set.
- The large-matrix single test defined allpass order 64 but passed literal 16; it now passes its configured order. The large ensemble harness omitted dimension 96 and only requested five seeds per dimension; it now requests seven dimensions and 100 seeds each.
- The random ensemble harness requested 1,000 runs rather than the paper's 10,000. It now requests seeds 0–9999. A sixth optional solver argument exposes the convergence tolerance; the random single test uses the paper's 1e-4 threshold, while the solver default remains 1e-8.
- Recomputations previously appended to the shipped result files. They now write separate `*-recomputed.txt` files. Archived data are preserved, including their historical column ordering.
- Zero refinement now stops at a numerically exact root, avoiding changes in which side of a fixed allpass-design grid point is selected by roundoff. The Figure 16 FFT buffer accommodates the returned support while retaining its original plotting grid. Auto-updating legends in Figures 7 and 16 are disabled so later annotations do not appear as extra data series.
- Plotting scripts assumed fixed solver output lengths and offsets, including an older causal Figure 12 cache. Windows now derive zero lag from the returned support instead of hard-coded FFT lengths.
- Figures 19–20 selected different matrix entries from the paper. The plotted rows and columns now match the paper (rows 4–6, columns 29–32). Figure 12 now explicitly uses the example's 1e-10 convergence threshold.
- Figure 3 wrote an unnumbered filename; Figure 20 supplied five tick labels for three tick positions. Both are corrected.

## Evidence gaps and limits

The shipped first ensemble contains 1,020 rows, 1,001 unique seeds, and 19 repeated seeds, rather than evidence for 10,000 independent instances. Computing its unmodified rows gives mean PU error 2.1908756e-10 (sample standard deviation 8.8158880e-11) and mean normalised LS excess 4.9761787e-4 (sample standard deviation 1.3136302e-4). These are close to the paper's rounded statistics but do not verify its claimed sample count. No duplicate rows have been silently removed, and no missing trials have been invented.

The second ensemble contains all 700 distinct (seed, dimension) pairs expected by the plotting script. Its archived results have not been regenerated during this audit. Historical computation times cannot be reproduced exactly on different hardware.

MAT-file caches in the author's local checkout are ignored by Git. The fresh-checkout workflow exercises the uncached paths. Figures 8 and 11 are explanatory diagrams, not supplied MATLAB plots. Successful script execution by itself does not establish agreement with every curve, theorem, or numerical claim in the paper. The figure windows now derive the zero-lag offset from the solver output length; the physical plotting ranges are retained.

Visual comparison of the recomputed Figures 19–20 with the accepted manuscript still shows different coefficient patterns and amplitudes, even after selecting the published rows 4–6 and columns 29–32. Their execution is now verified, but agreement with the published random realization is not. The original realization, seed/scaling settings, or final figure-generation files need author confirmation; this PR does not infer replacement data from the plot.

The email mentioning `FileList.txt` and three third-party comparison-method folders matches a local polynomial Gram–Schmidt/QR project, not this repository. Those files were not imported into Procrustes. PEVD Toolbox has its own MIT license; redistribution rights for any other third-party implementations must be established from their actual licenses or permission text rather than inferred from permission to use them.

## Validation

Numerical regressions and all supplied plotting scripts are run by `.github/workflows/reproducibility.yml` in MATLAB R2024b with Signal Processing Toolbox. The workflow retains generated EPS/PNG plots and a reproduction log. Local MATLAB execution required MathWorks sign-in and did not reach script execution. Full ensemble regeneration is excluded from the automated run; workflow results are recorded in the pull request.

The final code run passed on GitHub: [MATLAB workflow 39154534321](https://github.com/SebastianJiroSchlecht/PolynomialProcrustes/actions/runs/39154534321), code commit `bac8304`. All numerical regressions and all ten plotting/evaluation scripts completed from a fresh recursive checkout with no MAT-file caches; 20 EPS files were generated. Rendered plots were visually spot-checked against the manuscript, including Figures 12, 16 and 19–20. The three-channel example returned LS mismatch 1.09395447151, agreeing with the reported 1.0940 to rounding, and PU error 4.73e-23. The full Monte Carlo ensembles were not rerun. Later audit-text and whitespace-only changes do not alter the tested algorithms.
