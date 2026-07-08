# CRAN comments for xplaineff

This is the first CRAN submission of `xplaineff`.

## Test results

Local checks were run with
`env LC_ALL=C OMP_NUM_THREADS=1 R CMD check --as-cran xplaineff_0.1.0.tar.gz`.

- macOS Ventura 13.1, R 4.5.0: `0 errors | 0 warnings | 4 notes`.
- Windows Server 2022, R-release 4.6.1 (win-builder): `0 errors | 0 warnings | 1 note`.
- Windows Server 2022, R-devel r90199 (win-builder): `0 errors | 0 warnings | 1 note`.
- `checking CRAN incoming feasibility ... NOTE`: new submission.
- `checking for future file timestamps ... NOTE`: unable to verify current time in the local check environment.
- `checking top-level files ... NOTE`: pandoc is not installed in the local check environment.
- `checking HTML version of manual ... NOTE`: local HTML Tidy is not recent enough, and package `V8` is unavailable.

## Package contents

The package tarball excludes non-package directories such as `simulation/`, `scripts/`, `paper/`, and `figures/`.

## Dependencies

`xplaineff` depends on R (>= 4.3.0) and uses `Rcpp`/`RcppArmadillo` for compiled code.
