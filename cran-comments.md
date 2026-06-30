# CRAN comments for gadget

This is the first CRAN submission of `gadget`.

## Test results

Local checks were run with `LC_ALL=C R CMD check --no-install --no-manual --no-vignettes`.

- macOS Ventura 13.1, R 4.5.0: `0 errors | 0 warnings | 1 note`.
- The remaining note is `checking for future file timestamps ... NOTE unable to verify current time`.

## Package contents

The package tarball excludes non-package directories such as `simulation/`, `scripts/`, `paper/`, and `figures/`.

## Dependencies

`gadget` depends on R (>= 4.3.0) and uses `Rcpp`/`RcppArmadillo` for compiled code.
