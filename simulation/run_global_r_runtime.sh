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

if [ "$MODE" = "smoke" ]; then
  N_VEC="500"
  REPS=1
  PREDICT_REPS=2
  OUTDIR="simulation/results/global_r_runtime_smoke"
  FIGDIR="simulation/results/paper_figures_smoke"
  MODELS="${BENCHMARK_MODELS:-toy}"
else
  N_VEC="500,1000,2500,5000"
  REPS=20
  PREDICT_REPS=20
  OUTDIR="simulation/results/global_r_runtime"
  FIGDIR="simulation/results/paper_figures"
  MODELS="${BENCHMARK_MODELS:-rf,toy}"
fi

D=10
N_GRID=20
N_INTERVALS=20
DATADIR="simulation/data/global_r_runtime"
PAPER_FIGDIR="paper/figures"
INCLUDE_MLR3=false
case ",${MODELS}," in
  *",mlr3_rf,"*) INCLUDE_MLR3=true ;;
esac

echo "=== Benchmark mode: ${MODE} ==="
echo "    N_VEC=${N_VEC}  D=${D}  n_grid=${N_GRID}  n_intervals=${N_INTERVALS}"
echo "    global_models=${MODELS}"
echo "    threads: OMP_NUM_THREADS=${OMP_NUM_THREADS}  OMP_THREAD_LIMIT=${OMP_THREAD_LIMIT}"
echo "    openmp: KMP_INIT_AT_FORK=${KMP_INIT_AT_FORK}  KMP_AFFINITY=${KMP_AFFINITY}"
echo "    R: DATATABLE_NUM_THREADS=${DATATABLE_NUM_THREADS}  GADGET_BENCH_LOAD_ALL=${GADGET_BENCH_LOAD_ALL}"
echo "    DATADIR=${DATADIR}  OUTDIR=${OUTDIR}  FIGDIR=${FIGDIR}"

echo "1. Generating benchmark data..."
Rscript simulation/generate_runtime_data.R \
  --outdir "${DATADIR}" \
  --N-vec "${N_VEC}" \
  --D-vec "${D}"

echo "2. Running global R package benchmark..."
Rscript simulation/benchmark_global_r_runtime.R \
  --datadir "${DATADIR}" \
  --outdir "${OUTDIR}" \
  --reps "${REPS}" \
  --predict-reps "${PREDICT_REPS}" \
  --N-vec "${N_VEC}" \
  --D "${D}" \
  --n-grid "${N_GRID}" \
  --n-intervals "${N_INTERVALS}" \
  --models "${MODELS}"

echo "3. Summarizing and plotting..."
Rscript simulation/summarize_global_r_runtime.R \
  --indir "${OUTDIR}" \
  --figdir "${FIGDIR}" \
  --fixed-D "${D}" \
  --models "${MODELS}" \
  --include-mlr3 "${INCLUDE_MLR3}"

if [ "$MODE" = "publication" ] && [ -f "${FIGDIR}/global_r_methods.png" ]; then
  mkdir -p "${PAPER_FIGDIR}"
  cp "${FIGDIR}/global_r_methods.png" "${PAPER_FIGDIR}/global_r_methods.png"
  echo "Synced global runtime figure to ${PAPER_FIGDIR}/global_r_methods.png"
fi

echo "Done. Global raw CSVs and summary.csv in ${OUTDIR}/; global figure in ${FIGDIR}/"
echo "Use simulation/results/split_search_figures/regional_methods.png for the split-search benchmark figure."
