#!/usr/bin/env Rscript
# Summarize ranger layout-sensitivity diagnostic timings.

args = commandArgs(trailingOnly = TRUE)
indir = "simulation/results/ranger_layout_sensitivity"

i = 1L
while (i <= length(args)) {
  if (args[i] == "--indir" && i < length(args)) {
    indir = args[i + 1L]; i = i + 2L
  } else {
    i = i + 1L
  }
}

library(data.table)
raw_path = file.path(indir, "ranger_layout_sensitivity.csv")
if (!file.exists(raw_path)) {
  stop("No diagnostic raw CSV found at ", raw_path)
}

dt = fread(raw_path)
dt[, predict_time_sec := as.numeric(predict_time_sec)]
summary = dt[, .(
  predict_median = stats::median(predict_time_sec),
  predict_q25 = stats::quantile(predict_time_sec, 0.25, names = FALSE),
  predict_q75 = stats::quantile(predict_time_sec, 0.75, names = FALSE),
  predict_mean = mean(predict_time_sec),
  predict_sd = stats::sd(predict_time_sec),
  n_rep = .N
), by = .(N, D, n_grid, feature, layout)]
summary[, rel_to_grid_major := predict_median / predict_median[layout == "grid_major"], by = .(N, D, n_grid, feature)]
setorder(summary, N, D, n_grid, feature, layout)

summary_path = file.path(indir, "ranger_layout_sensitivity_summary.csv")
fwrite(summary, summary_path)
message("Written: ", summary_path)
print(summary)
