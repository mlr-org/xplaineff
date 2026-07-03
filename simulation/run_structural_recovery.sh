#!/bin/bash
# Run the full structural recovery benchmark.

set -euo pipefail
cd "$(dirname "$0")/.."

all_args=("$@")
datadir="simulation/data/structural_recovery"
outdir="simulation/results/structural_recovery"
figdir="simulation/results/paper_figures"
paper_figdir="paper/figures"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --datadir)
      datadir="$2"
      shift 2
      ;;
    --outdir)
      outdir="$2"
      shift 2
      ;;
    --figdir)
      figdir="$2"
      shift 2
      ;;
    --paper-figdir)
      paper_figdir="$2"
      shift 2
      ;;
    *)
      shift
      ;;
  esac
done

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-xplaineff"
mkdir -p "$MPLCONFIGDIR"

echo "1. Generating shared accuracy datasets..."
if [[ ${#all_args[@]} -gt 0 ]]; then
  Rscript simulation/generate_structural_recovery_data.R "${all_args[@]}"
else
  Rscript simulation/generate_structural_recovery_data.R
fi

echo "2. xplaineff accuracy..."
if [[ ${#all_args[@]} -gt 0 ]]; then
  Rscript simulation/benchmark_structural_recovery_xplaineff.R "${all_args[@]}"
else
  Rscript simulation/benchmark_structural_recovery_xplaineff.R
fi

echo "3. effector accuracy..."
if [[ ${#all_args[@]} -gt 0 ]]; then
  python3 -u simulation/benchmark_structural_recovery_effector.py "${all_args[@]}"
else
  python3 -u simulation/benchmark_structural_recovery_effector.py
fi

echo "4. Summarize and refresh paper figures..."
Rscript simulation/summarize_structural_recovery.R \
  --indir "$outdir" \
  --figdir "$figdir" \
  --paper-figdir "$paper_figdir"

echo "Done. Results in $outdir and figures in $figdir."
