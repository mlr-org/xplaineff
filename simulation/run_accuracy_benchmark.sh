#!/bin/bash
# Run the full structural recovery benchmark.

set -euo pipefail
cd "$(dirname "$0")/.."

all_args=("$@")
datadir="simulation/data/accuracy"
outdir="simulation/results/accuracy"
figdir="simulation/results/figures"
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

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-gadget"
mkdir -p "$MPLCONFIGDIR"

echo "1. Generating shared accuracy datasets..."
if [[ ${#all_args[@]} -gt 0 ]]; then
  Rscript simulation/generate_accuracy_data.R "${all_args[@]}"
else
  Rscript simulation/generate_accuracy_data.R
fi

echo "2. GADGET accuracy..."
if [[ ${#all_args[@]} -gt 0 ]]; then
  Rscript simulation/benchmark_accuracy_gadget.R "${all_args[@]}"
else
  Rscript simulation/benchmark_accuracy_gadget.R
fi

echo "3. effector accuracy..."
if [[ ${#all_args[@]} -gt 0 ]]; then
  python3 -u simulation/benchmark_accuracy_effector.py "${all_args[@]}"
else
  python3 -u simulation/benchmark_accuracy_effector.py
fi

echo "4. Summarize and refresh paper figures..."
Rscript simulation/summarize_accuracy.R \
  --indir "$outdir" \
  --figdir "$figdir" \
  --paper-figdir "$paper_figdir"

echo "Done. Results in $outdir and figures in $figdir."
