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
PARALLEL_SUB="${PARALLEL_SUB:-true}"
GENERATE_DATA="${GENERATE_DATA:-true}"
GLOBAL_SUB_JOBS="${GLOBAL_SUB_JOBS:-3}"
case "${GLOBAL_SUB_JOBS}" in
  ''|*[!0-9]*) GLOBAL_SUB_JOBS=3 ;;
esac
if [ "${GLOBAL_SUB_JOBS}" -lt 1 ]; then GLOBAL_SUB_JOBS=1; fi
if [ "${GLOBAL_SUB_JOBS}" -gt 3 ]; then GLOBAL_SUB_JOBS=3; fi

if [ "$MODE" = "smoke" ]; then
  N_VEC="500"
  D_VEC="10"
  FIXED_N="500"
  FIXED_D="10"
  N_GRID=10
  N_INTERVALS=10
  N_GRID_VEC="10"
  N_INTERVALS_VEC="10"
  REPS=1
  OUTDIR="simulation/results/global_r_runtime_smoke"
  FIGDIR="simulation/results/paper_figures_smoke"
  MODELS="${BENCHMARK_MODELS:-toy}"
else
  N_VEC="5000,10000,25000,50000"
  D_VEC="10,20,50,100"
  FIXED_N="10000"
  FIXED_D="20"
  N_GRID=20
  N_INTERVALS=20
  N_GRID_VEC="10,20,50"
  N_INTERVALS_VEC="10,20,50"
  REPS=30
  OUTDIR="simulation/results/global_r_runtime"
  FIGDIR="simulation/results/paper_figures"
  MODELS="${BENCHMARK_MODELS:-rf,toy}"
fi

DATADIR="simulation/data/global_r_runtime"
PAPER_FIGDIR="paper/figures"
INCLUDE_MLR3=false
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
echo "    PARALLEL_SUB=${PARALLEL_SUB}  GLOBAL_SUB_JOBS=${GLOBAL_SUB_JOBS}  GENERATE_DATA=${GENERATE_DATA}"

if [ "$GENERATE_DATA" = "true" ]; then
  echo "1. Generating benchmark data..."
  Rscript simulation/generate_runtime_data.R \
    --outdir "${DATADIR}" \
    --N-vec "${N_VEC}" \
    --D-vec "${D_VEC}" \
    --fixed-N "${FIXED_N}" \
    --fixed-D "${FIXED_D}"
else
  echo "1. Skipping benchmark data generation; using ${DATADIR}"
fi

echo "2. Running global R package benchmark..."
if [ "$PARALLEL_SUB" = "true" ]; then
  echo "    parallel mode: up to ${GLOBAL_SUB_JOBS} sub_experiment(s) concurrently, K=1 each"
  PIDS=()
  SUB_NAMES=()
  wait_for_sub_batch() {
    local failed=0
    local idx
    for idx in "${!PIDS[@]}"; do
      if ! wait "${PIDS[$idx]}"; then
        echo "    ERROR: sub_experiment ${SUB_NAMES[$idx]} failed; see ${OUTDIR}/_run_${SUB_NAMES[$idx]}.log" >&2
        failed=1
      fi
    done
    PIDS=()
    SUB_NAMES=()
    if [ "${failed}" -ne 0 ]; then
      exit 1
    fi
  }
  for SUB in vs_N vs_D vs_res; do
    LOG="${OUTDIR}/_run_${SUB}.log"
    mkdir -p "${OUTDIR}"
    Rscript simulation/benchmark_global_r_runtime.R \
      --datadir "${DATADIR}" \
      --outdir "${OUTDIR}" \
      --reps "${REPS}" \
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
      --output-suffix "${SUB}" \
      > "${LOG}" 2>&1 &
    PIDS+=("$!")
    SUB_NAMES+=("${SUB}")
    if [ "${#PIDS[@]}" -ge "${GLOBAL_SUB_JOBS}" ]; then
      wait_for_sub_batch
    fi
  done
  wait_for_sub_batch
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
echo "Use simulation/run_runtime_benchmark.sh for the coordinated global, regional, and diagnostic workflow."
