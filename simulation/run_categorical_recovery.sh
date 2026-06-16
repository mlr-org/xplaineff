#!/bin/bash
# Run the categorical split-feature recovery benchmark (GADGET-only).
# Usage:
#   bash simulation/run_categorical_recovery.sh           # publication: all sweeps, 30 seeds, digit DGP
#   bash simulation/run_categorical_recovery.sh smoke     # tiny wiring check
#   CORES=8 bash simulation/run_categorical_recovery.sh   # parallel over seeds
#   SWEEPS="leakage K N D group_frac slope_mag" bash simulation/run_categorical_recovery.sh
#   DGP_TYPE=group SWEEPS=leakage bash simulation/run_categorical_recovery.sh   # group DGP, leakage robustness

set -e
cd "$(dirname "$0")/.."

MODE="${1:-publication}"
CORES="${CORES:-1}"
SWEEPS="${SWEEPS:-leakage K N D slope_mag}"
GROUP_FRAC="${GROUP_FRAC:-0.3333333333333333}"
GROUP_FRAC_VEC="${GROUP_FRAC_VEC:-0.1666666666666667,0.3333333333333333,0.5}"
DGP_TYPE="${DGP_TYPE:-digit}"
MAX_PARTITION_K="${MAX_PARTITION_K:-12}"
SUMMARY_SWEEPS="$(printf "%s" "${SWEEPS}" | tr ' ' ',')"

# Keep model-side libraries single-threaded; parallelism is over seeds only.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export DATATABLE_NUM_THREADS="${DATATABLE_NUM_THREADS:-1}"

if [ "$MODE" = "smoke" ]; then
  N_SEEDS=3
  OUTDIR="simulation/results/categorical_recovery_smoke/${DGP_TYPE}"
  FIGDIR="simulation/results/paper_figures_smoke"
  EXTRA="--leakage-vec 0,0.05,0.1,0.25,1 --K-vec 6,8 --N-vec 500,2000 --D-vec 3,10 --group-frac-vec 0.1666666666666667,0.5 --slope-mag-vec 1,4"
else
  N_SEEDS=30
  OUTDIR="simulation/results/categorical_recovery/${DGP_TYPE}"
  FIGDIR="simulation/results/paper_figures"
  EXTRA=""
fi

echo "=== Categorical recovery benchmark: mode=${MODE} dgp=${DGP_TYPE} seeds=${N_SEEDS} cores=${CORES} ==="
echo "    OUTDIR=${OUTDIR}"
echo "    SWEEPS=${SWEEPS}"
echo "    GROUP_FRAC=${GROUP_FRAC}"
echo "    GROUP_FRAC_VEC=${GROUP_FRAC_VEC}"
echo "    MAX_PARTITION_K=${MAX_PARTITION_K}"

for SWEEP in ${SWEEPS}; do
  echo "--- sweep: ${SWEEP} ---"
  Rscript simulation/benchmark_categorical_recovery.R \
    --sweep "${SWEEP}" \
    --outdir "${OUTDIR}" \
    --n-seeds "${N_SEEDS}" \
    --cores "${CORES}" \
    --dgp-type "${DGP_TYPE}" \
    --max-partition-K "${MAX_PARTITION_K}" \
    --group-frac "${GROUP_FRAC}" \
    --group-frac-vec "${GROUP_FRAC_VEC}" \
    ${EXTRA}
done

echo "--- summarize and plot ---"
Rscript simulation/summarize_categorical_recovery.R \
  --indir "${OUTDIR}" \
  --figdir "${FIGDIR}" \
  --sweeps "${SUMMARY_SWEEPS}"

echo "Done. Raw CSVs in ${OUTDIR}/; figures in ${FIGDIR}/"
