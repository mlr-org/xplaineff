#!/bin/bash
# Run the publication efficiency benchmark.
# Usage:
#   bash simulation/run_global_r_runtime.sh          # publication global-effect grid
#   bash simulation/run_global_r_runtime.sh smoke    # tiny dependency / wiring check

set -e
cd "$(dirname "$0")/.."

MODE="${1:-publication}"

# Keep every benchmark single-threaded and avoid Intel/OpenMP shared-memory failures in sandboxed shells.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OMP_THREAD_LIMIT="${OMP_THREAD_LIMIT:-1}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-FALSE}"
export KMP_INIT_AT_FORK="${KMP_INIT_AT_FORK:-FALSE}"
export KMP_AFFINITY="${KMP_AFFINITY:-disabled}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export DATATABLE_NUM_THREADS="${DATATABLE_NUM_THREADS:-1}"
export RCPP_PARALLEL_NUM_THREADS="${RCPP_PARALLEL_NUM_THREADS:-1}"
export GADGET_BENCH_LOAD_ALL="${GADGET_BENCH_LOAD_ALL:-true}"
CORES="${CORES:-1}"
PARALLEL_SUB="${PARALLEL_SUB:-false}"

if [ "$MODE" = "smoke" ]; then
  N_VEC="500"
  D_VEC="10"
  FIXED_N="500"
  FIXED_D="10"
  N_GRID_VEC="20"
  N_INTERVALS_VEC="20"
  REPS=1
  PREDICT_REPS=2
  OUTDIR="simulation/results/global_r_runtime_smoke"
  FIGDIR="simulation/results/paper_figures_smoke"
  MODELS="${BENCHMARK_MODELS:-toy}"
else
  N_VEC="1000,2500,5000,10000"
  D_VEC="5,10,20"
  FIXED_N="1000"
  FIXED_D="10"
  N_GRID_VEC="10,20,50"
  N_INTERVALS_VEC="10,20,50"
  REPS=20
  PREDICT_REPS=20
  OUTDIR="simulation/results/global_r_runtime"
  FIGDIR="simulation/results/paper_figures"
  MODELS="${BENCHMARK_MODELS:-rf,toy}"
fi

N_GRID=20
N_INTERVALS=20
DATADIR="simulation/data/global_r_runtime"
PAPER_FIGDIR="paper/figures"
INCLUDE_MLR3=false
DATA_N_VEC="${N_VEC},${FIXED_N}"
DATA_D_VEC="${D_VEC},${FIXED_D}"
case ",${MODELS}," in
  *",mlr3_rf,"*) INCLUDE_MLR3=true ;;
esac

echo "=== Benchmark mode: ${MODE} ==="
echo "    N_VEC=${N_VEC}  D_VEC=${D_VEC}  fixed_N=${FIXED_N}  fixed_D=${FIXED_D}"
echo "    n_grid=${N_GRID}  n_intervals=${N_INTERVALS}"
echo "    n_grid_vec=${N_GRID_VEC}  n_intervals_vec=${N_INTERVALS_VEC}"
echo "    global_models=${MODELS}"
echo "    threads: OMP_NUM_THREADS=${OMP_NUM_THREADS}  OMP_THREAD_LIMIT=${OMP_THREAD_LIMIT}"
echo "    openmp: KMP_INIT_AT_FORK=${KMP_INIT_AT_FORK}  KMP_AFFINITY=${KMP_AFFINITY}"
echo "    R: DATATABLE_NUM_THREADS=${DATATABLE_NUM_THREADS}  GADGET_BENCH_LOAD_ALL=${GADGET_BENCH_LOAD_ALL}"
echo "    DATADIR=${DATADIR}  OUTDIR=${OUTDIR}  FIGDIR=${FIGDIR}"

echo "1. Generating benchmark data..."
Rscript simulation/generate_runtime_data.R \
  --outdir "${DATADIR}" \
  --N-vec "${DATA_N_VEC}" \
  --D-vec "${DATA_D_VEC}"

echo "2. Running global R package benchmark..."
if [ "$PARALLEL_SUB" = "true" ]; then
  echo "    parallel mode: 3 sub_experiments concurrently, K=1 each"
  PIDS=""
  for SUB in vs_N vs_D vs_res; do
    SKIP_BL="true"
    if [ "$SUB" = "vs_N" ]; then SKIP_BL="false"; fi
    LOG="${OUTDIR}/_run_${SUB}.log"
    mkdir -p "${OUTDIR}"
    Rscript simulation/benchmark_global_r_runtime.R \
      --datadir "${DATADIR}" \
      --outdir "${OUTDIR}" \
      --reps "${REPS}" \
      --predict-reps "${PREDICT_REPS}" \
      --N-vec "${N_VEC}" \
      --D-vec "${D_VEC}" \
      --fixed-N "${FIXED_N}" \
      --fixed-D "${FIXED_D}" \
      --n-grid "${N_GRID}" \
      --n-intervals "${N_INTERVALS}" \
      --n-grid-vec "${N_GRID_VEC}" \
      --n-intervals-vec "${N_INTERVALS_VEC}" \
      --models "${MODELS}" \
      --cores 1 \
      --sub-experiments "${SUB}" \
      --skip-predict-baseline "${SKIP_BL}" \
      --output-suffix "${SUB}" \
      > "${LOG}" 2>&1 &
    PIDS="$PIDS $!"
  done
  wait $PIDS
  echo "    merging per-sub_experiment CSVs..."
  Rscript -e '
    args <- commandArgs(trailingOnly = TRUE)
    outdir <- args[1]
    models <- strsplit(args[2], ",", fixed = TRUE)[[1]]
    subs <- c("vs_N", "vs_D", "vs_res")
    for (m in models) {
      parts <- list()
      for (s in subs) {
        p <- file.path(outdir, sprintf("global_r_runtime_%s_%s.csv", m, s))
        if (file.exists(p)) parts[[s]] <- read.csv(p)
      }
      if (length(parts)) {
        merged <- do.call(rbind, parts)
        write.csv(merged, file.path(outdir, sprintf("global_r_runtime_%s.csv", m)), row.names = FALSE)
        for (s in names(parts)) file.remove(file.path(outdir, sprintf("global_r_runtime_%s_%s.csv", m, s)))
        cat(sprintf("Merged: global_r_runtime_%s.csv (%d rows)\n", m, nrow(merged)))
      }
    }
  ' "${OUTDIR}" "${MODELS}"
else
  Rscript simulation/benchmark_global_r_runtime.R \
    --datadir "${DATADIR}" \
    --outdir "${OUTDIR}" \
    --reps "${REPS}" \
    --predict-reps "${PREDICT_REPS}" \
    --N-vec "${N_VEC}" \
    --D-vec "${D_VEC}" \
    --fixed-N "${FIXED_N}" \
    --fixed-D "${FIXED_D}" \
    --n-grid "${N_GRID}" \
    --n-intervals "${N_INTERVALS}" \
    --n-grid-vec "${N_GRID_VEC}" \
    --n-intervals-vec "${N_INTERVALS_VEC}" \
    --models "${MODELS}" \
    --cores "${CORES}"
fi

echo "3. Summarizing and plotting..."
Rscript simulation/summarize_global_r_runtime.R \
  --indir "${OUTDIR}" \
  --figdir "${FIGDIR}" \
  --fixed-D "${FIXED_D}" \
  --models "${MODELS}" \
  --include-mlr3 "${INCLUDE_MLR3}"

if [ "$MODE" = "publication" ] && [ -f "${FIGDIR}/global_r_methods.png" ]; then
  mkdir -p "${PAPER_FIGDIR}"
  cp "${FIGDIR}/global_r_methods.png" "${PAPER_FIGDIR}/global_r_methods.png"
  echo "Synced global runtime figure to ${PAPER_FIGDIR}/global_r_methods.png"
fi

echo "Done. Global raw CSVs and summary.csv in ${OUTDIR}/; global figure in ${FIGDIR}/"
echo "Use simulation/summarize_split_search_runtime_large.R to regenerate the split-search benchmark figure."
