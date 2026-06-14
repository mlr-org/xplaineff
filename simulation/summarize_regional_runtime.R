#!/usr/bin/env Rscript
# Summarize regional runtime benchmark results.

Sys.setenv(
  OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS", "1"),
  OMP_THREAD_LIMIT = Sys.getenv("OMP_THREAD_LIMIT", "1"),
  OMP_PROC_BIND = Sys.getenv("OMP_PROC_BIND", "FALSE"),
  KMP_INIT_AT_FORK = Sys.getenv("KMP_INIT_AT_FORK", "FALSE"),
  KMP_AFFINITY = Sys.getenv("KMP_AFFINITY", "disabled"),
  OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS", "1"),
  MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS", "1"),
  VECLIB_MAXIMUM_THREADS = Sys.getenv("VECLIB_MAXIMUM_THREADS", "1"),
  DATATABLE_NUM_THREADS = Sys.getenv("DATATABLE_NUM_THREADS", "1"),
  RCPP_PARALLEL_NUM_THREADS = Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "1")
)

args = commandArgs(trailingOnly = TRUE)
indir = "simulation/results/regional_runtime"
figdir = "simulation/results/paper_figures"

i = 1L
while (i <= length(args)) {
  if (args[i] == "--indir" && i < length(args)) {
    indir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--figdir" && i < length(args)) {
    figdir = args[i + 1L]; i = i + 2L
  } else {
    i = i + 1L
  }
}

dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(ggplot2)
library(scales)
setDTthreads(1L)

files = list.files(indir, pattern = "^regional_runtime_.*\\.csv$", full.names = TRUE)
if (!length(files)) {
  stop("No regional runtime CSV files found in ", indir)
}
dt = rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
if (!nrow(dt)) {
  stop("Regional runtime CSV files are empty.")
}
if (!("status" %in% names(dt))) dt[, status := "ok"]
dt[is.na(status) | status == "", status := "ok"]

for (col in c("precompute_time_sec", "split_time_sec", "total_time_sec")) {
  dt[, (col) := as.numeric(get(col))]
}
dt[, resolution := as.integer(resolution)]
dt[, n_split := as.integer(n_split)]

n_err = dt[status == "error", .N]
if (n_err > 0L) {
  message("Warning: ", n_err, " regional benchmark row(s) have status=error; see error_message.")
}

dt_ok = dt[status == "ok" & is.finite(total_time_sec)]
summary_dt = dt_ok[, .(
  precompute_median = stats::median(precompute_time_sec),
  precompute_q25 = stats::quantile(precompute_time_sec, 0.25, names = FALSE),
  precompute_q75 = stats::quantile(precompute_time_sec, 0.75, names = FALSE),
  split_median = stats::median(split_time_sec),
  split_q25 = stats::quantile(split_time_sec, 0.25, names = FALSE),
  split_q75 = stats::quantile(split_time_sec, 0.75, names = FALSE),
  total_median = stats::median(total_time_sec),
  total_q25 = stats::quantile(total_time_sec, 0.25, names = FALSE),
  total_q75 = stats::quantile(total_time_sec, 0.75, names = FALSE),
  n_rep = .N
), by = .(module, package, impl, effect, method, model_type, sub_experiment, N, D, resolution, n_split)]

summary_dt[, label := fcase(
  package == "gadget", "gadget",
  package == "effector", "effector",
  default = package
)]

write.csv(summary_dt, file.path(indir, "summary.csv"), row.names = FALSE)
message("Written: ", file.path(indir, "summary.csv"))

palette_values = c(
  "gadget" = "#1f77b4",
  "effector" = "#2ca02c"
)
shape_values = c(
  "gadget" = 16,
  "effector" = 17
)
x_offset_values = c(
  "gadget" = 0.985,
  "effector" = 1.015
)

plot_metric = function(data, metric, filename, title) {
  if (!nrow(data)) return(invisible(NULL))
  data = copy(data)
  median_col = sprintf("%s_median", metric)
  q25_col = sprintf("%s_q25", metric)
  q75_col = sprintf("%s_q75", metric)
  data[, time_median := get(median_col)]
  data[, time_q25 := get(q25_col)]
  data[, time_q75 := get(q75_col)]
  data[, effect_label := fifelse(effect == "pdp", "PDP", "ALE")]
  data[, model_label := fifelse(model_type == "rf", "RF model", "Toy model")]
  data[, panel := factor(paste(effect_label, model_label, sep = " / "), levels = c(
    "PDP / RF model",
    "PDP / Toy model",
    "ALE / RF model",
    "ALE / Toy model"
  ))]
  data[, sweep_label := fcase(
    sub_experiment == "vs_N", "Sample size N\nD = 20, resolution = 20, splits = 2",
    sub_experiment == "vs_D", "Feature dimension D\nN = 10,000, resolution = 20, splits = 2",
    sub_experiment == "vs_res", "Resolution\nN = 10,000, D = 20, splits = 2",
    sub_experiment == "vs_split", "Number of splits\nN = 10,000, D = 20, resolution = 20",
    default = sub_experiment
  )]
  data[, sweep_label := factor(sweep_label, levels = c(
    "Sample size N\nD = 20, resolution = 20, splits = 2",
    "Feature dimension D\nN = 10,000, resolution = 20, splits = 2",
    "Resolution\nN = 10,000, D = 20, splits = 2",
    "Number of splits\nN = 10,000, D = 20, resolution = 20"
  ))]
  data[, x_value := fcase(
    sub_experiment == "vs_N", as.numeric(N),
    sub_experiment == "vs_D", as.numeric(D),
    sub_experiment == "vs_res", as.numeric(resolution),
    sub_experiment == "vs_split", as.numeric(n_split),
    default = as.numeric(NA)
  )]
  data[, x_display := x_value * x_offset_values[label]]
  data[, time_median_plot := pmax(time_median, .Machine$double.eps)]
  data[, time_q25_plot := pmax(time_q25, .Machine$double.eps)]
  data[, time_q75_plot := pmax(time_q75, .Machine$double.eps)]

  p = ggplot(
    data,
    aes(x = x_display, y = time_median_plot, color = label, fill = label, group = label)
  ) +
    geom_ribbon(aes(ymin = time_q25_plot, ymax = time_q75_plot), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = label), size = 2.0, stroke = 0.2) +
    facet_grid(panel ~ sweep_label, scales = "free") +
    scale_x_log10(labels = comma, breaks = sort(unique(data$x_value))) +
    scale_y_log10(labels = label_number(accuracy = 0.1, trim = TRUE)) +
    scale_color_manual(values = palette_values, breaks = names(palette_values), drop = FALSE) +
    scale_fill_manual(values = alpha(palette_values, 0.25), breaks = names(palette_values), drop = FALSE, guide = "none") +
    scale_shape_manual(values = shape_values, breaks = names(shape_values), drop = FALSE, guide = "none") +
    labs(title = title, x = "Value of varied parameter", y = "Median runtime (seconds, log scale)", color = NULL) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.spacing.y = unit(1.4, "lines"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold", lineheight = 0.95)
    )

  ggsave(file.path(figdir, filename), p, width = 13, height = 12, dpi = 220)
  message("Written: ", file.path(figdir, filename))
}

plot_metric(summary_dt, "split", "regional_split_methods.png", "Regional split-search runtime")
plot_metric(summary_dt, "total", "regional_total_methods.png", "Regional total runtime")

message("Done.")
