# Simulation Directory

This directory contains the scripts and generated artifacts for the paper simulations.
The generated data and CSV results are regenerable and are ignored by git; publication figures are synced to
`paper/figures/` and kept in git.

## Main Workflows

- `run_runtime_benchmark.sh` runs the coordinated global and regional runtime workflow.
  It generates shared sweep-cell data, runs the global R package benchmark, and runs the regional
  xplaineff/effector benchmark.
  Each run writes to `results/runtime_runs/<RUN_ID>/`, so incomplete or repeated runs do not overwrite earlier output.
  Use `bash simulation/run_runtime_benchmark.sh smoke` for a small wiring check.
- `run_global_r_runtime.sh` runs the global R feature-effect benchmark.
  Publication mode uses 20 timed repetitions after one untimed warm-up.
  It runs three sweeps: sample size `N`, feature dimension `D`, and PD/ALE resolution.
  Runtime trend plots show median timings with interquartile-range ribbons.
  Use `GLOBAL_PACKAGES=...`, `GLOBAL_IMPLS=...`, or `GLOBAL_METHODS=...` to restrict the benchmark grid.
  For example, the current xplaineff-only automatic-path run is:

  ```sh
  GLOBAL_PACKAGES=xplaineff GLOBAL_IMPLS=auto bash simulation/run_global_r_runtime.sh
  ```

- `run_regional_runtime.sh` runs the regional runtime benchmark.
  It generates the shared sweep data, runs xplaineff with the reticulate-backed sklearn path, runs effector, and writes
  precompute, split-search, and total-runtime summaries.
  The historical compact ALE path is explicit through `ALE_COMPACT=true`, which is the default for this wrapper.
- `benchmark_global_r_runtime.R` and `summarize_global_r_runtime.R` are the R entry points behind the global runtime
  wrapper.
- `benchmark_regional_runtime_xplaineff.R`, `benchmark_regional_runtime_xplaineff_reticulate.R`, and
  `benchmark_regional_runtime_effector.py` are the regional runtime benchmark entry points.
- `plot_regional_runtime_paperformat.R` redraws the paper-format regional runtime comparison figures from a completed
  regional summary CSV.
- `run_categorical_recovery.sh` and `run_structural_recovery.sh` run the appendix diagnostics.
  They write raw CSV files under `results/categorical_recovery/` and `results/structural_recovery/`, and sync
  publication figures to `paper/figures/`.

## Paper Runtime Inputs

The current paper figures use integrated summary files rather than a single raw run directory.
The source files are kept under `results/runtime_runs/20260616_214256/`.

- Global runtime figure:
  `global_r_runtime/summary_paper_with_xpl_auto_20260719.csv`.
  This combines the original R-package global benchmark with the newer xplaineff automatic-path run.
- Regional split-search figure:
  `regional_runtime/summary_paper_with_xpl_pruning_effector050_20260719.csv`.
  This combines xplaineff reticulate/sklearn rows with effector 0.5.0 rows.
- Regional total-runtime figure:
  `regional_runtime/summary_paper_with_xpl_pruning_effector050_20260719.csv`.
  This uses the same regional summary and the `total_*` columns.

To redraw the global paper figure from the integrated summary:

```sh
GLOBAL_RUN=simulation/results/runtime_runs/20260616_214256
Rscript simulation/summarize_global_r_runtime.R \
  --summary "${GLOBAL_RUN}/global_r_runtime/summary_paper_with_xpl_auto_20260719.csv" \
  --figdir "${GLOBAL_RUN}/paper_figures" \
  --paper-figdir paper/figures
```

To redraw the regional paper figures from the integrated summary:

```sh
REGIONAL_RUN=simulation/results/runtime_runs/20260616_214256
Rscript simulation/plot_regional_runtime_paperformat.R \
  --summary "${REGIONAL_RUN}/regional_runtime/summary_paper_with_xpl_pruning_effector050_20260719.csv" \
  --figdir "${REGIONAL_RUN}/paper_figures" \
  --paper-figdir paper/figures \
  --tag local_compare
```

This command writes `regional_precompute_runtime_local_compare.png`,
`regional_split_runtime_local_compare.png`, and `regional_total_runtime_local_compare.png`.
The paper currently includes the split-search and total-runtime figures.

## Current Runtime Artifacts

- `results/runtime_runs/20260616_214256/` is the consolidated paper runtime directory.
  Its `global_r_runtime/` folder keeps the global raw CSVs, diagnostics, and integrated paper summary.
  Its `regional_runtime/` folder keeps the regional paper summary.
  Its `paper_figures/` folder keeps paper-format copies generated from the integrated summaries.
- `results/runtime_runs/server_reticulate_sklearn_compact_reps20_20260714_0100/` is the xplaineff regional raw source.
  The `README_pruning_20260719.md` file documents how the active-effect pruning replacement rows were fused for the
  paper summary.
- `results/runtime_runs/server_effector050_py313_reps20_20260715_1441/` is the effector 0.5.0 regional raw source.
- `results/runtime_runs/bikeshare_exhaustive_categorical_20260722/` stores the one-off bikeshare check for the
  categorical exhaustive split option.
  It is useful for auditing the examples around `categorical_split = "exhaustive"`, but it is not used to draw the
  runtime figures.
- `results/runtime_runs/archive_scripts_20260719/` stores one-off profiling scripts kept for auditability.
  They are not part of the publication workflow.
- `results/runtime_runs/current_precompute_path_logic_20260717.md` is a working note about path-selection logic.
  It is not an input to any figure.

Temporary smoke, probe, copy, or one-off diagnostic outputs should not be kept under `results/runtime_runs/` after use.
Use clearly named scratch directories while debugging, then remove them before committing.
If a one-off script is needed for auditability, keep it inside the ignored result directory that it produced.

## Regenerable Data

- `data/global_r_runtime/` contains generated input data shared by the current global and regional runtime benchmarks.
- `data/structural_recovery/` contains generated input data for the structural-recovery benchmark.

## Diagnostic Helpers

- `diagnose_global_r_runtime.R` checks whether a global runtime run has all expected cells and repetitions.
- `probe_effector_wrapper_official.py` checks the official effector wrapper behavior and is not part of the publication
  runtime workflow.
- `probe_bikeshare_exhaustive_categorical.R` reruns the bikeshare categorical exhaustive split check and writes to
  `results/runtime_runs/bikeshare_exhaustive_categorical_20260722/`.
