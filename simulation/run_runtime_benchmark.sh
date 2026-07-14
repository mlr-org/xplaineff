#!/bin/bash
# Run the coordinated global, regional, and diagnostic runtime benchmarks.
# Usage:
#   bash simulation/run_runtime_benchmark.sh
#   bash simulation/run_runtime_benchmark.sh smoke

set -e
cd "$(dirname "$0")/.."

MODE="${1:-publication}"

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
export PYTHONHASHSEED="${PYTHONHASHSEED:-21}"
export XPLAINEFF_BENCH_LOAD_ALL="${XPLAINEFF_BENCH_LOAD_ALL:-true}"
export PARALLEL_SUB="${PARALLEL_SUB:-true}"
export GLOBAL_SUB_JOBS="${GLOBAL_SUB_JOBS:-3}"
export PARALLEL_PHASES="${PARALLEL_PHASES:-true}"
export PARALLEL_REGIONAL_PACKAGES="${PARALLEL_REGIONAL_PACKAGES:-true}"
export SYNC_PAPER_FIGURES="${SYNC_PAPER_FIGURES:-true}"

DATADIR="simulation/data/global_r_runtime"
PAPER_FIGDIR="paper/figures"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
RUN_ROOT="${RUN_ROOT:-simulation/results/runtime_runs/${RUN_ID}}"

if [ "$MODE" = "smoke" ]; then
  N_VEC="500"
  D_VEC="10"
  FIXED_N="500"
  FIXED_D="10"
  RES_VEC="10"
  RESOLUTION="10"
  N_SPLIT_VEC="2"
  N_SPLIT="2"
  N_QUANTILES="19"
  EFFECTOR_NUMERICAL_GRID_SIZE="20"
  REPS="1"
  MODELS="${BENCHMARK_MODELS:-toy}"
  GLOBAL_OUTDIR="${RUN_ROOT}/global_r_runtime_smoke"
  REGIONAL_OUTDIR="${RUN_ROOT}/regional_runtime_smoke"
  FIGDIR="${RUN_ROOT}/paper_figures_smoke"
else
  N_VEC="1000,5000,10000,20000"
  D_VEC="10,20,50,100"
  FIXED_N="10000"
  FIXED_D="20"
  RES_VEC="10,20,50"
  RESOLUTION="20"
  N_SPLIT_VEC="2,5,8,10"
  N_SPLIT="2"
  N_QUANTILES="19"
  EFFECTOR_NUMERICAL_GRID_SIZE="20"
  REPS="20"
  MODELS="${BENCHMARK_MODELS:-rf,toy}"
  GLOBAL_OUTDIR="${RUN_ROOT}/global_r_runtime"
  REGIONAL_OUTDIR="${RUN_ROOT}/regional_runtime"
  FIGDIR="${RUN_ROOT}/paper_figures"
fi

echo "=== Runtime benchmark mode: ${MODE} ==="
echo "    RUN_ID=${RUN_ID}  RUN_ROOT=${RUN_ROOT}"
echo "    N_VEC=${N_VEC}  D_VEC=${D_VEC}  fixed_N=${FIXED_N}  fixed_D=${FIXED_D}"
echo "    RES_VEC=${RES_VEC}  resolution=${RESOLUTION}  n_split_vec=${N_SPLIT_VEC}"
echo "    n_quantiles=${N_QUANTILES}  effector_numerical_grid_size=${EFFECTOR_NUMERICAL_GRID_SIZE}"
echo "    reps=${REPS}  models=${MODELS}  PARALLEL_SUB=${PARALLEL_SUB}  GLOBAL_SUB_JOBS=${GLOBAL_SUB_JOBS}"
echo "    PARALLEL_PHASES=${PARALLEL_PHASES}  PARALLEL_REGIONAL_PACKAGES=${PARALLEL_REGIONAL_PACKAGES}"
echo "    model backend: R xplaineff uses default_predict_fun dispatch; sklearn n_jobs=1"

run_global_benchmark() {
  echo "2a. Running global benchmark..."
  if [ "$MODE" = "smoke" ]; then
    GENERATE_DATA=false BENCHMARK_MODELS="${MODELS}" GLOBAL_OUTDIR="${GLOBAL_OUTDIR}" GLOBAL_FIGDIR="${FIGDIR}" \
      SYNC_PAPER_FIGURES=false bash simulation/run_global_r_runtime.sh smoke
  else
    GENERATE_DATA=false BENCHMARK_MODELS="${MODELS}" GLOBAL_OUTDIR="${GLOBAL_OUTDIR}" GLOBAL_FIGDIR="${FIGDIR}" \
      SYNC_PAPER_FIGURES=false bash simulation/run_global_r_runtime.sh publication
  fi
}

run_regional_xplaineff() {
  echo "2b. Running regional benchmark: xplaineff..."
  Rscript simulation/benchmark_regional_runtime_xplaineff.R \
    --datadir "${DATADIR}" \
    --outdir "${REGIONAL_OUTDIR}" \
    --reps "${REPS}" \
    --N-vec "${N_VEC}" \
    --D-vec "${D_VEC}" \
    --fixed-N "${FIXED_N}" \
    --fixed-D "${FIXED_D}" \
    --resolution "${RESOLUTION}" \
    --resolution-vec "${RES_VEC}" \
    --n-split "${N_SPLIT}" \
    --n-split-vec "${N_SPLIT_VEC}" \
    --n-quantiles "${N_QUANTILES}" \
    --models "${MODELS}"

}

run_regional_effector() {
  echo "2c. Running regional benchmark: effector..."
  python3 simulation/benchmark_regional_runtime_effector.py \
    --datadir "${DATADIR}" \
    --outdir "${REGIONAL_OUTDIR}" \
    --reps "${REPS}" \
    --N-vec "${N_VEC}" \
    --D-vec "${D_VEC}" \
    --fixed-N "${FIXED_N}" \
    --fixed-D "${FIXED_D}" \
    --resolution "${RESOLUTION}" \
    --resolution-vec "${RES_VEC}" \
    --n-split "${N_SPLIT}" \
    --n-split-vec "${N_SPLIT_VEC}" \
    --numerical-features-grid-size "${EFFECTOR_NUMERICAL_GRID_SIZE}" \
    --models "${MODELS}"
}

run_regional_benchmark() {
  if [ "$PARALLEL_REGIONAL_PACKAGES" = "true" ]; then
    mkdir -p "${REGIONAL_OUTDIR}"
    run_regional_xplaineff > "${REGIONAL_OUTDIR}/_run_xplaineff.log" 2>&1 &
    XPLAINEFF_PID=$!
    run_regional_effector > "${REGIONAL_OUTDIR}/_run_effector.log" 2>&1 &
    EFFECTOR_PID=$!

    FAILED=0
    if ! wait "${XPLAINEFF_PID}"; then
      echo "ERROR: regional xplaineff benchmark failed; see ${REGIONAL_OUTDIR}/_run_xplaineff.log" >&2
      FAILED=1
    fi
    if ! wait "${EFFECTOR_PID}"; then
      echo "ERROR: regional effector benchmark failed; see ${REGIONAL_OUTDIR}/_run_effector.log" >&2
      FAILED=1
    fi
    if [ "${FAILED}" -ne 0 ]; then
      exit 1
    fi
  else
    run_regional_xplaineff
    run_regional_effector
  fi
  echo "2d. Summarizing regional benchmark..."
  Rscript simulation/summarize_regional_runtime.R \
    --indir "${REGIONAL_OUTDIR}" \
    --figdir "${FIGDIR}"
}

echo "1. Generating shared benchmark data..."
Rscript simulation/generate_runtime_data.R \
  --outdir "${DATADIR}" \
  --N-vec "${N_VEC}" \
  --D-vec "${D_VEC}" \
  --fixed-N "${FIXED_N}" \
  --fixed-D "${FIXED_D}"

if [ "$PARALLEL_PHASES" = "true" ]; then
  echo "2. Running global and regional benchmarks concurrently..."
  run_global_benchmark &
  GLOBAL_PID=$!
  run_regional_benchmark &
  REGIONAL_PID=$!

  FAILED=0
  if ! wait "${GLOBAL_PID}"; then
    echo "ERROR: global benchmark failed." >&2
    FAILED=1
  fi
  if ! wait "${REGIONAL_PID}"; then
    echo "ERROR: regional benchmark failed." >&2
    FAILED=1
  fi
  if [ "${FAILED}" -ne 0 ]; then
    exit 1
  fi
else
  run_global_benchmark
  run_regional_benchmark
fi

if [ "$MODE" = "publication" ] && [ "$SYNC_PAPER_FIGURES" = "true" ]; then
  mkdir -p "${PAPER_FIGDIR}"
  if [ -f "${FIGDIR}/global_r_methods.png" ]; then
    cp "${FIGDIR}/global_r_methods.png" "${PAPER_FIGDIR}/global_r_methods.png"
    echo "Synced global_r_methods.png to ${PAPER_FIGDIR}/"
  fi
  if [ -f "${FIGDIR}/regional_split_methods_linear.png" ]; then
    cp "${FIGDIR}/regional_split_methods_linear.png" "${PAPER_FIGDIR}/regional_split_methods_linear.png"
    echo "Synced regional_split_methods_linear.png to ${PAPER_FIGDIR}/"
  fi
fi

echo "Done."
