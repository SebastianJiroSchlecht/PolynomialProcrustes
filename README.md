# <img src="Media/Logo_PP.png" alt="Logo" width="50" height="50"> Polynomial Procrustes <img src="Media/Logo_PP.png" alt="Logo" width="50" height="50">


This repository contains the accompanying files for the paper

S. Weiss, S.J. Schlecht, and M. Moonen: "Best Least Squares Paraunitary Approximation: Analytic Procrustes Problem," IEEE Trans. Signal Process., vol. 74, pp. 904–919, 2026. [doi:10.1109/TSP.2026.3663449](https://doi.org/10.1109/TSP.2026.3663449)

This repository supplies MATLAB scripts for the numerical figures and archived ensemble results. The inventory and limitations below distinguish plotting archived results from recomputing the experiments.

The code relies on the PEVD Toolbox, which can be found here: 
[http://pevd-toolbox.eee.strath.ac.uk/](http://pevd-toolbox.eee.strath.ac.uk/). The toolbox is also included as a submodule.

## Setup and checks

Clone recursively so the pinned PEVD Toolbox is available:

```sh
git clone --recurse-submodules https://github.com/SebastianJiroSchlecht/PolynomialProcrustes.git
cd PolynomialProcrustes
```

For an existing clone, run `git submodule update --init --recursive`. In MATLAB, start in the repository root and run:

```matlab
addpath(fullfile(pwd,'PEVDToolbox','General'));
addpath(fullfile(pwd,'PEVDToolbox','Decompositions'));
addpath('tests');
test_reproducibility;
audit_figures;
```

The checks target MATLAB R2024b with Signal Processing Toolbox (`impz`, `dftmtx`). `audit_figures` executes every supplied paper plotting script in a separate workspace and writes `reproduction-log.txt`. It does not run the full Monte Carlo ensembles. The GitHub Actions workflow runs these checks from a fresh recursive checkout and saves the EPS plots as an artifact.

## Paper reproduction inventory

| Paper output | Script | Notes |
| --- | --- | --- |
| Figures 1–2 | `Generate_Figs12.m` | Computes the solution. |
| Figure 3 | `Generate_Fig3.m` | Analytic singular values. |
| Figures 4–7 | `Generate_Figs4to7.m` | Truncated Fourier switching approximation. |
| Figure 8 | No script supplied | Schematic target phase/switching response in the paper. |
| Figures 9–10 | `Generate_Figs910.m` | Allpass design and convergence. |
| Figure 11 | No script supplied | Algorithm flowchart. |
| Figure 12 | `Generate_Fig12.m` | Reuses `Generate_Fig12.mat` if present. |
| Figures 13–16 | `Generate_Figs13to16.m` | Reuses `Qhat1024.mat` if present. |
| Figure 17 | `Generate_Fig17.m` | Example singular values for M=32. |
| Figure 18 | `Generate_Fig18.m` | Reads the archived large-matrix ensemble. |
| Figures 19–20 | `Generate_Figs1920.m` | Reuses `LargeMatrixM96Example.mat` if present. |
| First ensemble statistics, Section VI-C | `ProcrustesRandomMatrixEnsembleEval.m` | Reads the archived random-matrix ensemble. |

MAT-file caches are ignored by Git and are not shipped. To test computation rather than reuse a previous result, run in a fresh checkout or move the caches elsewhere first. Plot scripts write EPS files in the working directory. Figures 12–16 and 19–20 include fixed sample-index windows from the original experiments; these must be checked if solver convergence or parameters change.

## Ensemble limitations

The archived `LargeMatrixEnsembleResults.txt` contains 700 distinct (seed, dimension) pairs: 100 seeds for each M in {2,4,8,16,32,64,96}, matching the paper's sample counts. Its runtime column is historical and will depend on the machine and MATLAB release.

The archived `EnsembleResults.txt` contains 1,020 rows but only 1,001 distinct seeds (0–1000), including 19 repeated seeds. This does **not** provide the 10,000 distinct instances described in Section VI-C. Its existing columns are seed, theoretical LS bound, final FFT length, convergence error, paraunitarity error, and achieved LS error. Evaluation of this archive approximately recovers the reported means (PU: 2.19e-10; relative LS excess: 4.98e-4), but cannot independently establish the claimed sample count. The archived file is preserved unchanged.

`ProcrustesRandomMatrixEnsembleTest` now requests 10,000 seeds and the paper's 1e-4 convergence threshold; `ProcrustesLargeMatrixEnsembleTest` requests all seven dimensions and 100 seeds each. They write separate `*-recomputed.txt` files so the archived results remain intact. They append to those files: move a previous recomputation away before starting a fresh full run. To evaluate newly computed results, explicitly change the input filename in the relevant evaluation/plot script. New single-run outputs use the documented column order: seed, bound, convergence error, FFT length, PU error, LS error. Do not concatenate them with the older archive, whose middle two columns are reversed.

These full ensembles can take substantial time and memory. The automated checks exercise one random instance and all plotting scripts, not all 10,700 simulations. See `REPRODUCIBILITY.md` for the audit findings and remaining evidence gaps.

## Dependencies and third-party code

The pinned PEVD Toolbox is a submodule with its own MIT license and copyright notice; retain that notice. The root MIT license applies to this repository's code. No `FileList.txt` or three separately supplied comparison-method directories are present in this repository. Permission to use code does not document permission to redistribute it under this repository's license. Before adding external implementations, retain their original license and attribution and record the upstream source and version; where redistribution terms are not documented, link to the upstream implementation instead.
