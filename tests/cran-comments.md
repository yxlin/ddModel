## Test environments
- local macOS (Apple clang), R 4.4.1
- Ubuntu 22.04 (GCC), R 4.4.1 (GitHub Actions)
- Windows Server 2022 (Rtools), R 4.4.1 (GitHub Actions)
- win-builder (devel, release)

## R CMD check results
0 errors | 0 warnings | 1 note

* Note: New submission.

## Summary
`ddModel` provides fast R/C++ implementations for the Diffusion Decision Model (DDM):
- density, distribution, random sampling, and likelihood evaluation
- vectorised functions for large datasets
- optional integration with the `ggdmc` ecosystem for hierarchical inference

## Dependencies & System requirements
- Imports: Rcpp (>= 1.0.7), RcppArmadillo (>= 0.10.7.5.0)
- SystemRequirements: C++17
- Suggests: ggdmc, ggdmcModel, ggdmcPrior (examples that use these are wrapped in
  `if (requireNamespace("ggdmcModel", quietly = TRUE)) { ... }` etc.)

## Examples & Run time
- Examples are short and self-contained (no I/O, no network, no parallelism).
- Typical runtime < 5 seconds on win-builder.
- All randomness uses local `set.seed()` and does not modify global RNGKind.

## URLs / DOIs
- All URLs verified with `R CMD check --as-cran`.
- DOIs formatted as `https://doi.org/...`.

## Data / Size
- No large data files in the tarball.
- README images are under `man/figures/` and small (<1MB).

## Native code
- C++ via Rcpp/RcppArmadillo only.
- Compiles with CRAN’s default toolchains on all platforms.
- No OpenMP and no non-portable flags.

## Notes for reviewers
- This is a first submission. There are no reverse dependencies yet.
- The package does not write to user directories and uses `tempdir()` where needed.
- No internet access at build or run time.



- First submission of ddModel — R + C++ tools for the Diffusion Decision Model (DDM): density/distribution/simulation and likelihood evaluation.

- Tested on macOS, clang20, Windows (Github Action via rhub).

- Uses Rcpp/RcppArmadillo; SystemRequirements: C++17.

- R CMD check: 0 errors | 0 warnings | 0 note (new submission).