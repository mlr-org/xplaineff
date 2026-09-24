#!/usr/bin/env bash
# Reproduce the simulation outputs used by the xplaineff paper.
# Usage:
#   bash replication/run_replication.sh smoke
#   bash replication/run_replication.sh full
#   STAGES=global,regional bash replication/run_replication.sh full

set -euo pipefail
cd "$(dirname "$0")/.."

MODE="${1:-smoke}"
case "${MODE}" in
  smoke|full) ;;
  *)
    echo "Usage: bash replication/run_replication.sh [smoke|full]" >&2
    exit 2
    ;;
esac

if [[ -n "${STAGES+x}" ]]; then
  STAGES="${STAGES}"
elif [[ "${MODE}" == "smoke" ]]; then
  STAGES="global"
else
  STAGES="global,regional,categorical,structural"
fi

OUTPUT_DIR="${OUTPUT_DIR:-replication/output/${MODE}}"
SEED="${SEED:-21}"
CORES="${CORES:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OMP_THREAD_LIMIT="${OMP_THREAD_LIMIT:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export DATATABLE_NUM_THREADS="${DATATABLE_NUM_THREADS:-1}"
export RCPP_PARALLEL_NUM_THREADS="${RCPP_PARALLEL_NUM_THREADS:-1}"
export PYTHONHASHSEED="${PYTHONHASHSEED:-21}"
if [[ -z "${XPLAINEFF_BENCH_LOAD_ALL+x}" ]]; then
  if Rscript -e 'if (requireNamespace("pkgbuild", quietly = TRUE) && pkgbuild::has_build_tools()) quit(status = 0) else quit(status = 1)' >/dev/null 2>&1; then
    export XPLAINEFF_BENCH_LOAD_ALL=true
  else
    export XPLAINEFF_BENCH_LOAD_ALL=false
  fi
else
  export XPLAINEFF_BENCH_LOAD_ALL
fi

mkdir -p "${OUTPUT_DIR}"
Rscript replication/session_info.R "${OUTPUT_DIR}/sessionInfo.txt"
printf 'mode=%s\nstages=%s\nseed=%s\nxplaineff_bench_load_all=%s\n' \
  "${MODE}" "${STAGES}" "${SEED}" "${XPLAINEFF_BENCH_LOAD_ALL}" > "${OUTPUT_DIR}/run-config.txt"

has_stage() {
  case ",${STAGES}," in
    *,"$1",*) return 0 ;;
    *) return 1 ;;
  esac
}

generate_runtime_data() {
  local n_vec d_vec fixed_n fixed_d runtime_datadir
  if [[ "${MODE}" == "smoke" ]]; then
    n_vec="500"
    d_vec="10"
    fixed_n="500"
    fixed_d="10"
  else
    n_vec="1000,5000,10000,20000"
    d_vec="10,20,50,100"
    fixed_n="10000"
    fixed_d="20"
  fi
  runtime_datadir="${OUTPUT_DIR}/data/global_r_runtime"
  Rscript simulation/generate_runtime_data.R \
    --outdir "${runtime_datadir}" \
    --seed "${SEED}" \
    --N-vec "${n_vec}" \
    --D-vec "${d_vec}" \
    --fixed-N "${fixed_n}" \
    --fixed-D "${fixed_d}"
}

run_global() {
  echo "[replication] global runtime (${MODE})"
  generate_runtime_data
  local global_packages="${GLOBAL_PACKAGES:-}"
  local global_impls="${GLOBAL_IMPLS:-}"
  local benchmark_models="${BENCHMARK_MODELS:-rf,toy}"
  if [[ "${MODE}" == "smoke" ]]; then
    global_packages="${GLOBAL_PACKAGES:-xplaineff}"
    global_impls="${GLOBAL_IMPLS:-r,cpp}"
    benchmark_models="${BENCHMARK_MODELS:-toy}"
  fi
  RUN_ROOT="${OUTPUT_DIR}/runtime" \
    GLOBAL_DATADIR="${OUTPUT_DIR}/data/global_r_runtime" \
    GLOBAL_OUTDIR="${OUTPUT_DIR}/global_r_runtime" \
    GLOBAL_FIGDIR="${OUTPUT_DIR}/figures/global" \
    GENERATE_DATA=false \
    SYNC_PAPER_FIGURES=false \
    BENCHMARK_MODELS="${benchmark_models}" \
    GLOBAL_PACKAGES="${global_packages}" \
    GLOBAL_IMPLS="${global_impls}" \
    bash simulation/run_global_r_runtime.sh "${MODE}"
}

run_regional() {
  echo "[replication] regional runtime (${MODE})"
  generate_runtime_data
  RUN_ROOT="${OUTPUT_DIR}/runtime" \
    REGIONAL_DATADIR="${OUTPUT_DIR}/data/global_r_runtime" \
    REGIONAL_OUTDIR="${OUTPUT_DIR}/regional_runtime" \
    REGIONAL_FIGDIR="${OUTPUT_DIR}/figures/regional" \
    GENERATE_DATA=false \
    SYNC_PAPER_FIGURES=false \
    BENCHMARK_MODELS="${BENCHMARK_MODELS:-rf,toy}" \
    bash simulation/run_regional_runtime.sh "${MODE}"
}

run_categorical() {
  local outdir="${OUTPUT_DIR}/categorical_recovery"
  local figdir="${OUTPUT_DIR}/figures/categorical"
  local n_seeds="30"
  local sweeps="leakage K N D slope_mag"
  local extra=()
  if [[ "${MODE}" == "smoke" ]]; then
    n_seeds="1"
    sweeps="leakage"
    extra=(
      --leakage-vec "0,0.1,1"
      --K-vec "6"
      --N-vec "500"
      --D-vec "3"
      --group-frac-vec "0.5"
      --slope-mag-vec "1"
    )
  fi
  echo "[replication] categorical recovery (${MODE})"
  mkdir -p "${outdir}" "${figdir}"
  IFS=' ' read -r -a sweep_array <<< "${sweeps}"
  for sweep in "${sweep_array[@]}"; do
    Rscript simulation/benchmark_categorical_recovery.R \
      --sweep "${sweep}" \
      --outdir "${outdir}" \
      --n-seeds "${n_seeds}" \
      --cores "${CORES}" \
      --dgp-type digit \
      "${extra[@]}"
  done
  Rscript simulation/summarize_categorical_recovery.R \
    --indir "${outdir}" \
    --figdir "${figdir}" \
    --sweeps "$(printf '%s' "${sweeps}" | tr ' ' ',')"
}

run_structural() {
  local outdir="${OUTPUT_DIR}/structural_recovery"
  local figdir="${OUTPUT_DIR}/figures/structural"
  local datadir="${OUTPUT_DIR}/data/structural_recovery"
  local n_seeds="30"
  local n_vec="200,500,1000,5000"
  local d_vec="5,10,20"
  local variants="num_0,num_04,cat"
  if [[ "${MODE}" == "smoke" ]]; then
    n_seeds="1"
    n_vec="200"
    d_vec="5"
    variants="num_0"
  fi
  echo "[replication] structural recovery (${MODE})"
  bash simulation/run_structural_recovery.sh \
    --datadir "${datadir}" \
    --outdir "${outdir}" \
    --figdir "${figdir}" \
    --paper-figdir "" \
    --n-seeds "${n_seeds}" \
    --N-vec "${n_vec}" \
    --D-vec "${d_vec}" \
    --variants "${variants}"
}

if has_stage global; then run_global; fi
if has_stage regional; then run_regional; fi
if has_stage categorical; then run_categorical; fi
if has_stage structural; then run_structural; fi

echo "[replication] completed: ${OUTPUT_DIR}"
