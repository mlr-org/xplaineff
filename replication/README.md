# xplaineff replication bundle

This directory contains the independent entry point for regenerating the simulation outputs used in the paper.
It is excluded from the CRAN source package by `.Rbuildignore`; it is intended for reviewers and researchers
who want to rerun the computational experiments from a clean checkout.

## Quick smoke run

From the repository root:

```sh
bash replication/run_replication.sh smoke
```

The default smoke run executes the small xplaineff global runtime workflow with the toy model and current `r` and
`cpp` implementations. It writes its output under `replication/output/smoke/`. The run records `sessionInfo.txt` and
`run-config.txt` before starting.

Additional smoke stages can be requested explicitly:

```sh
STAGES=global,categorical bash replication/run_replication.sh smoke
```

The regional and structural stages require the external Python/reticulate dependencies described below.

## Full run

```sh
bash replication/run_replication.sh full
```

The full run executes the global runtime, regional runtime, categorical recovery, and structural recovery workflows.
It uses the paper benchmark grids and 20 repetitions for runtime benchmarks, and 30 seeds for recovery benchmarks.
The output is written to `replication/output/full/` and includes raw CSV files, summaries, diagnostics, figures,
the run configuration, and the R session information.

Stages can be run separately when a machine does not have every optional dependency:

```sh
STAGES=global bash replication/run_replication.sh full
STAGES=categorical,structural bash replication/run_replication.sh full
```

Set `OUTPUT_DIR` to write to a different location. Set `XPLAINEFF_BENCH_LOAD_ALL=false` to benchmark an installed
version of `xplaineff` instead of the checkout. By default, the wrapper uses the checkout when R build tools are
available and falls back to an installed `xplaineff` package otherwise.

## Fixed settings

The runtime data generator uses seed `21` by default. The runtime benchmarks use the same seed for the random forest
backends. The categorical and structural recovery scripts use their deterministic seed sequences beginning at
`1001`, as specified in their benchmark implementations. Override the runtime data seed with `SEED=<integer>` when
performing a deliberate sensitivity analysis.

The full runtime grids are:

- sample size: `N = 1000, 5000, 10000, 20000`, with `D = 20` fixed;
- feature dimension: `D = 10, 20, 50, 100`, with `N = 10000` fixed;
- resolution: `10, 20, 50`, with `N = 10000`, `D = 20` fixed;
- regional split depth: `2, 5, 8, 10`, with `N = 10000`, `D = 20` fixed.

## Dependencies

The R workflows require R >= 4.3 and the package dependencies in `DESCRIPTION`, plus `data.table`, `ggplot2`, and
`scales` for summaries. Running from the checkout uses `pkgload`; running against an installed package can disable
that path with `XPLAINEFF_BENCH_LOAD_ALL=false`.

The regional and structural comparisons additionally require Python 3, `numpy`, `pandas`, `scikit-learn`, and
`effector`, as well as the R package `reticulate` for the reticulate-backed xplaineff regional benchmark. The paper's
recorded comparison environment used Python 3.13.13, scikit-learn 1.9.0, effector 0.5.0, and reticulate 1.46.0.
The exact versions detected on the current machine are written to `sessionInfo.txt` for every run.

## Output map

- `data/`: generated input data used by the benchmark stages;
- `global_r_runtime/`: raw global benchmark CSV files, summary, and completeness diagnostics;
- `regional_runtime/`: raw regional benchmark CSV files and summary;
- `categorical_recovery/`: raw categorical recovery CSV files and summaries;
- `structural_recovery/`: raw structural recovery CSV files and summaries;
- `figures/`: figures generated from the corresponding summaries;
- `sessionInfo.txt`: R version, platform, and selected package versions;
- `run-config.txt`: mode, stages, and runtime seed used for the run.

The output directory is ignored by git. The source scripts remain under `simulation/`, so the bundle does not duplicate
the implementation or depend on a hidden local path.
