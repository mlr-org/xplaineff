#!/bin/bash
# Run the current regional runtime benchmark.
# Usage:
#   bash simulation/run_regional_runtime.sh
#   bash simulation/run_regional_runtime.sh smoke

set -e
cd "$(dirname "$0")/.."

MODE="${1:-publication}"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
RUN_ROOT="${RUN_ROOT:-simulation/results/runtime_runs/${RUN_ID}}"

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

if [ "$MODE" = "smoke" ]; then
  N_VEC="${N_VEC:-500}"
  D_VEC="${D_VEC:-10}"
  FIXED_N="${FIXED_N:-500}"
  FIXED_D="${FIXED_D:-10}"
  RESOLUTION="${RESOLUTION:-10}"
  RES_VEC="${RES_VEC:-10}"
  N_SPLIT="${N_SPLIT:-2}"
  N_SPLIT_VEC="${N_SPLIT_VEC:-2}"
  REPS="${REPS:-1}"
  MODELS="${BENCHMARK_MODELS:-toy}"
  OUTDIR="${REGIONAL_OUTDIR:-${RUN_ROOT}/regional_runtime_smoke}"
  FIGDIR="${REGIONAL_FIGDIR:-${RUN_ROOT}/figures_smoke}"
else
  N_VEC="${N_VEC:-1000,5000,10000,20000}"
  D_VEC="${D_VEC:-10,20,50,100}"
  FIXED_N="${FIXED_N:-10000}"
  FIXED_D="${FIXED_D:-20}"
  RESOLUTION="${RESOLUTION:-20}"
  RES_VEC="${RES_VEC:-10,20,50}"
  N_SPLIT="${N_SPLIT:-2}"
  N_SPLIT_VEC="${N_SPLIT_VEC:-2,5,8,10}"
  REPS="${REPS:-20}"
  MODELS="${BENCHMARK_MODELS:-rf,toy}"
  OUTDIR="${REGIONAL_OUTDIR:-${RUN_ROOT}/regional_runtime}"
  FIGDIR="${REGIONAL_FIGDIR:-${RUN_ROOT}/figures}"
fi

DATADIR="${REGIONAL_DATADIR:-simulation/data/global_r_runtime}"
PAPER_FIGDIR="${PAPER_FIGDIR:-paper/figures}"
SYNC_PAPER_FIGURES="${SYNC_PAPER_FIGURES:-true}"
PARALLEL_REGIONAL_PACKAGES="${PARALLEL_REGIONAL_PACKAGES:-true}"
GENERATE_DATA="${GENERATE_DATA:-true}"
N_QUANTILES="${N_QUANTILES:-19}"
MIN_NODE_SIZE="${MIN_NODE_SIZE:-50}"
EFFECTOR_NUMERICAL_GRID_SIZE="${EFFECTOR_NUMERICAL_GRID_SIZE:-20}"
RF_N_JOBS="${RF_N_JOBS:-1}"
ALE_COMPACT="${ALE_COMPACT:-true}"

echo "=== Regional runtime benchmark mode: ${MODE} ==="
echo "    RUN_ID=${RUN_ID}  RUN_ROOT=${RUN_ROOT}"
echo "    N_VEC=${N_VEC}  D_VEC=${D_VEC}  fixed_N=${FIXED_N}  fixed_D=${FIXED_D}"
echo "    resolution=${RESOLUTION}  resolution_vec=${RES_VEC}  n_split_vec=${N_SPLIT_VEC}"
echo "    reps=${REPS}  models=${MODELS}  ale_compact=${ALE_COMPACT}  rf_n_jobs=${RF_N_JOBS}"
echo "    OUTDIR=${OUTDIR}  FIGDIR=${FIGDIR}"

if [ "$GENERATE_DATA" = "true" ]; then
  echo "1. Generating shared regional benchmark data..."
  Rscript simulation/generate_runtime_data.R \
    --outdir "${DATADIR}" \
    --N-vec "${N_VEC}" \
    --D-vec "${D_VEC}" \
    --fixed-N "${FIXED_N}" \
    --fixed-D "${FIXED_D}"
else
  echo "1. Skipping data generation; using ${DATADIR}"
fi

run_xplaineff() {
  echo "2a. Running regional benchmark: xplaineff reticulate compact..."
  Rscript simulation/benchmark_regional_runtime_xplaineff_reticulate.R \
    --datadir "${DATADIR}" \
    --outdir "${OUTDIR}" \
    --reps "${REPS}" \
    --N-vec "${N_VEC}" \
    --D-vec "${D_VEC}" \
    --fixed-N "${FIXED_N}" \
    --fixed-D "${FIXED_D}" \
    --resolution "${RESOLUTION}" \
    --resolution-vec "${RES_VEC}" \
    --n-split "${N_SPLIT}" \
    --n-split-vec "${N_SPLIT_VEC}" \
    --min-node-size "${MIN_NODE_SIZE}" \
    --n-quantiles "${N_QUANTILES}" \
    --models "${MODELS}" \
    --rf-n-jobs "${RF_N_JOBS}" \
    --ale-compact "${ALE_COMPACT}" \
    --output-suffix "reticulate_sklearn_compact_reps${REPS}"
}

run_effector() {
  echo "2b. Running regional benchmark: effector..."
  python3 simulation/benchmark_regional_runtime_effector.py \
    --datadir "${DATADIR}" \
    --outdir "${OUTDIR}" \
    --reps "${REPS}" \
    --N-vec "${N_VEC}" \
    --D-vec "${D_VEC}" \
    --fixed-N "${FIXED_N}" \
    --fixed-D "${FIXED_D}" \
    --resolution "${RESOLUTION}" \
    --resolution-vec "${RES_VEC}" \
    --n-split "${N_SPLIT}" \
    --n-split-vec "${N_SPLIT_VEC}" \
    --min-node-size "${MIN_NODE_SIZE}" \
    --numerical-features-grid-size "${EFFECTOR_NUMERICAL_GRID_SIZE}" \
    --models "${MODELS}" \
    --output-suffix "effector050_reps${REPS}"
}

mkdir -p "${OUTDIR}"
if [ "$PARALLEL_REGIONAL_PACKAGES" = "true" ]; then
  run_xplaineff > "${OUTDIR}/_run_xplaineff_reticulate_compact.log" 2>&1 &
  XPLAINEFF_PID=$!
  run_effector > "${OUTDIR}/_run_effector.log" 2>&1 &
  EFFECTOR_PID=$!

  FAILED=0
  if ! wait "${XPLAINEFF_PID}"; then
    echo "ERROR: xplaineff regional benchmark failed; see ${OUTDIR}/_run_xplaineff_reticulate_compact.log" >&2
    FAILED=1
  fi
  if ! wait "${EFFECTOR_PID}"; then
    echo "ERROR: effector regional benchmark failed; see ${OUTDIR}/_run_effector.log" >&2
    FAILED=1
  fi
  if [ "${FAILED}" -ne 0 ]; then
    exit 1
  fi
else
  run_xplaineff
  run_effector
fi

echo "3. Summarizing regional benchmark..."
Rscript simulation/summarize_regional_runtime.R \
  --indir "${OUTDIR}" \
  --figdir "${FIGDIR}" \
  --paper-figdir ""

echo "4. Plotting regional paper-format figures..."
PAPER_SYNC_ARGS=()
if [ "$SYNC_PAPER_FIGURES" != "true" ]; then
  PAPER_SYNC_ARGS=(--no-paper-sync)
fi
Rscript simulation/plot_regional_runtime_paperformat.R \
  --summary "${OUTDIR}/summary.csv" \
  --figdir "${FIGDIR}" \
  --paper-figdir "${PAPER_FIGDIR}" \
  --tag "regional_runtime" \
  "${PAPER_SYNC_ARGS[@]}"

echo "Done. Regional raw CSVs and summary.csv in ${OUTDIR}/; figures in ${FIGDIR}/"
