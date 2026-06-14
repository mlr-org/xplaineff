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
export GADGET_BENCH_LOAD_ALL="${GADGET_BENCH_LOAD_ALL:-true}"
export PARALLEL_SUB="${PARALLEL_SUB:-true}"
export GLOBAL_SUB_JOBS="${GLOBAL_SUB_JOBS:-3}"

DATADIR="simulation/data/global_r_runtime"
PAPER_FIGDIR="paper/figures"

if [ "$MODE" = "smoke" ]; then
  N_VEC="500"
  D_VEC="10"
  FIXED_N="500"
  FIXED_D="10"
  RES_VEC="10"
  RESOLUTION="10"
  N_SPLIT_VEC="2"
  N_SPLIT="2"
  REPS="1"
  MODELS="${BENCHMARK_MODELS:-toy}"
  GLOBAL_OUTDIR="simulation/results/global_r_runtime_smoke"
  REGIONAL_OUTDIR="simulation/results/regional_runtime_smoke"
  DIAGNOSTIC_OUTDIR="simulation/results/ranger_layout_sensitivity_smoke"
  FIGDIR="simulation/results/paper_figures_smoke"
  RUN_DIAGNOSTIC="${RUN_DIAGNOSTIC:-false}"
else
  N_VEC="5000,10000,25000,50000"
  D_VEC="10,20,50,100"
  FIXED_N="10000"
  FIXED_D="20"
  RES_VEC="10,20,50"
  RESOLUTION="20"
  N_SPLIT_VEC="2,5,8,10"
  N_SPLIT="2"
  REPS="30"
  MODELS="${BENCHMARK_MODELS:-rf,toy}"
  GLOBAL_OUTDIR="simulation/results/global_r_runtime"
  REGIONAL_OUTDIR="simulation/results/regional_runtime"
  DIAGNOSTIC_OUTDIR="simulation/results/ranger_layout_sensitivity"
  FIGDIR="simulation/results/paper_figures"
  RUN_DIAGNOSTIC="${RUN_DIAGNOSTIC:-true}"
fi

echo "=== Runtime benchmark mode: ${MODE} ==="
echo "    N_VEC=${N_VEC}  D_VEC=${D_VEC}  fixed_N=${FIXED_N}  fixed_D=${FIXED_D}"
echo "    RES_VEC=${RES_VEC}  resolution=${RESOLUTION}  n_split_vec=${N_SPLIT_VEC}"
echo "    reps=${REPS}  models=${MODELS}  PARALLEL_SUB=${PARALLEL_SUB}  GLOBAL_SUB_JOBS=${GLOBAL_SUB_JOBS}"

echo "1. Generating shared benchmark data..."
Rscript simulation/generate_runtime_data.R \
  --outdir "${DATADIR}" \
  --N-vec "${N_VEC}" \
  --D-vec "${D_VEC}" \
  --fixed-N "${FIXED_N}" \
  --fixed-D "${FIXED_D}"

echo "2. Running global benchmark..."
if [ "$MODE" = "smoke" ]; then
  GENERATE_DATA=false BENCHMARK_MODELS="${MODELS}" bash simulation/run_global_r_runtime.sh smoke
else
  GENERATE_DATA=false BENCHMARK_MODELS="${MODELS}" bash simulation/run_global_r_runtime.sh publication
fi

echo "3. Running regional benchmark: gadget..."
Rscript simulation/benchmark_regional_runtime_gadget.R \
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
  --models "${MODELS}"

echo "4. Running regional benchmark: effector..."
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
  --models "${MODELS}"

echo "5. Summarizing regional benchmark..."
Rscript simulation/summarize_regional_runtime.R \
  --indir "${REGIONAL_OUTDIR}" \
  --figdir "${FIGDIR}"

if [ "$RUN_DIAGNOSTIC" = "true" ]; then
  echo "6. Running ranger layout-sensitivity diagnostic..."
  Rscript simulation/benchmark_ranger_layout_sensitivity.R \
    --datadir "${DATADIR}" \
    --outdir "${DIAGNOSTIC_OUTDIR}" \
    --reps "${REPS}" \
    --warmup 2
  Rscript simulation/summarize_ranger_layout_sensitivity.R \
    --indir "${DIAGNOSTIC_OUTDIR}"
fi

if [ "$MODE" = "publication" ]; then
  mkdir -p "${PAPER_FIGDIR}"
  for FIG in global_r_methods.png regional_split_methods.png regional_total_methods.png; do
    if [ -f "${FIGDIR}/${FIG}" ]; then
      cp "${FIGDIR}/${FIG}" "${PAPER_FIGDIR}/${FIG}"
      echo "Synced ${FIG} to ${PAPER_FIGDIR}/"
    fi
  done
fi

echo "Done."
