# CRAN comments for xplaineff

This is the first CRAN submission of `xplaineff`.

## Test results

Local checks were run with:

```sh
env LC_ALL=C OMP_NUM_THREADS=1 \
  PATH=/Applications/quarto/bin/tools:/opt/homebrew/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin \
  R_LIBS=/Users/zzz/Downloads/LMU_Work/GADGET/.r-lib \
  R CMD check --as-cran xplaineff_0.1.0.tar.gz
```

- macOS Ventura 13.1, R 4.5.0: `0 errors | 0 warnings | 2 notes`.
- Windows Server 2022, R-release 4.6.1 (win-builder): `0 errors | 0 warnings | 1 note`.
- Windows Server 2022, R-devel r90242 (win-builder): `0 errors | 0 warnings | 1 note`.
- `checking CRAN incoming feasibility ... NOTE`: new submission.
- `checking HTML version of manual ... NOTE`: package `V8` is unavailable in the local check environment.

## Package contents

The package tarball excludes non-package directories such as `simulation/`, `scripts/`, `paper/`, `figures/`,
and `.r-lib/`.

## Dependencies

`xplaineff` depends on R (>= 4.3.0) and uses `Rcpp`/`RcppArmadillo` for compiled code.
