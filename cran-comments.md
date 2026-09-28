# CRAN comments for xplaineff 0.1.1

This is a maintenance release of `xplaineff` (previous version 0.1.0).

## Changes

- Adds the introductory vignette `vignette("xplaineff")` with a full workflow on synthetic and
  bike-sharing data. It was left out of the 0.1.0 submission and is now complete; the code chunks that
  need `mlr3`, `mlr3learners`, `ranger`, or `ISLR2` are only evaluated when these suggested packages are
  installed.
- Fixes two memory issues in the compiled split search: effect matrices passed from R were aliased and
  modified in place when they contained `NaN`, and coerced non-double inputs could be freed while still in
  use. Both now work on a copy.
- Excludes observations with a missing split-feature value consistently from both child nodes and from
  the parent totals, guards missing interval indices in the ALE sweep, and skips non-finite local effects.
- Validates `n_quantiles >= 1` up front, makes the "model is required" error of `AleStrategy$fit()`
  reachable, and places observation points of categorical features on the plotted levels in regional PD
  plots.
- Enables roxygen markdown, restores the `GadgetTree` help topic to the index, aligns the documentation
  with the code, and removes unused internal helpers and imports.
- Adds three co-authors to `Authors@R`.

See `NEWS.md` for the complete list.

## Test environments

- Local: macOS Ventura 13.1, R 4.5.0, `R CMD check --as-cran --no-manual` on the source tarball built
  with vignettes.
- GitHub Actions (`R-CMD-check.yaml`): macOS (release), Windows (release), Ubuntu (R 4.3, release, devel).

## R CMD check results

Local: 0 errors | 0 warnings | 3 notes.

- `checking for future file timestamps ... NOTE: unable to verify current time`: the local check
  environment has no network access.
- `checking top-level files ... NOTE`: `pandoc` is not on the `PATH` of the local check environment, so
  `README.md` and `NEWS.md` could not be checked.
- `checking compiled code ... NOTE`: the local sandbox prevents `nm` from writing its cache file; the
  message is unrelated to the package.

The suggested package `xgboost` was not installed locally; its three tests are skipped in that case.

## Dependencies

`xplaineff` depends on R (>= 4.3.0) and uses `Rcpp`/`RcppArmadillo` for compiled code. There are no reverse
dependencies on CRAN.
