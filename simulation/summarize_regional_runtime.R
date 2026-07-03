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
run_id = format(Sys.time(), "%Y%m%d_%H%M%S")
run_root = file.path("simulation/results/runtime_runs", run_id)
indir = file.path(run_root, "regional_runtime")
figdir = file.path(run_root, "paper_figures")
paper_figdir = ""
paper_filenames = "regional_split_methods_linear.png"

i = 1L
while (i <= length(args)) {
  if (args[i] == "--indir" && i < length(args)) {
    indir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--figdir" && i < length(args)) {
    figdir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--paper-figdir" && i < length(args)) {
    paper_figdir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--paper-filenames" && i < length(args)) {
    paper_filenames = trimws(strsplit(args[i + 1L], ",", fixed = TRUE)[[1L]])
    paper_filenames = paper_filenames[nzchar(paper_filenames)]
    i = i + 2L
  } else {
    i = i + 1L
  }
}

dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
if (nzchar(paper_figdir)) {
  dir.create(paper_figdir, showWarnings = FALSE, recursive = TRUE)
}

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
  package == "xplaineff", "xplaineff",
  package == "effector", "effector",
  default = package
)]

write.csv(summary_dt, file.path(indir, "summary.csv"), row.names = FALSE)
message("Written: ", file.path(indir, "summary.csv"))

palette_values = c(
  "xplaineff" = "#1f77b4",
  "effector" = "#2ca02c"
)
shape_values = c(
  "xplaineff" = 16,
  "effector" = 17
)
x_offset_values = c(
  "xplaineff" = 0.985,
  "effector" = 1.015
)

format_axis_number = function(x) {
  format(x, trim = TRUE, digits = 3L, scientific = FALSE, big.mark = ",")
}

save_plot = function(plot_obj, filename, width, height, dpi = 220) {
  out = file.path(figdir, filename)
  ggsave(out, plot_obj, width = width, height = height, dpi = dpi)
  message("Written: ", out)
  if (nzchar(paper_figdir) && filename %in% paper_filenames) {
    paper_out = file.path(paper_figdir, filename)
    ggsave(paper_out, plot_obj, width = width, height = height, dpi = dpi)
    message("Synced: ", paper_out)
  }
}

plot_metric = function(data, metric, filename, title, facet_layout = "grid", x_scale = "log10") {
  if (!nrow(data)) return(invisible(NULL))
  facet_layout = match.arg(facet_layout, c("grid", "wrap"))
  x_scale = match.arg(x_scale, c("log10", "continuous"))
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
  data[, sweep_short := fcase(
    sub_experiment == "vs_N", "vs n",
    sub_experiment == "vs_D", "vs p",
    sub_experiment == "vs_res", "PDP grid / ALE intervals",
    sub_experiment == "vs_split", "vs splits",
    default = sub_experiment
  )]
  data[, fixed_desc := fcase(
    sub_experiment == "vs_N", "p = 20, resolution = 20, splits = 2",
    sub_experiment == "vs_D", "n = 10,000, resolution = 20, splits = 2",
    sub_experiment == "vs_res", "n = 10,000, p = 20, splits = 2",
    sub_experiment == "vs_split", "n = 10,000, p = 20, resolution = 20",
    default = sub_experiment
  )]
  data[, sweep_label := fcase(
    sub_experiment == "vs_N", "Sample size n\np = 20, resolution = 20, splits = 2",
    sub_experiment == "vs_D", "Feature dimension p\nn = 10,000, resolution = 20, splits = 2",
    sub_experiment == "vs_res", "PDP grid / ALE intervals\nn = 10,000, p = 20, splits = 2",
    sub_experiment == "vs_split", "Number of splits\nn = 10,000, p = 20, resolution = 20",
    default = sub_experiment
  )]
  data[, sweep_label := factor(sweep_label, levels = c(
    "Sample size n\np = 20, resolution = 20, splits = 2",
    "Feature dimension p\nn = 10,000, resolution = 20, splits = 2",
    "PDP grid / ALE intervals\nn = 10,000, p = 20, splits = 2",
    "Number of splits\nn = 10,000, p = 20, resolution = 20"
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

  panel_levels = levels(data$panel)
  sweep_levels = c(
    "vs n\np = 20, resolution = 20, splits = 2",
    "vs p\nn = 10,000, resolution = 20, splits = 2",
    "PDP grid / ALE intervals\nn = 10,000, p = 20, splits = 2",
    "vs splits\nn = 10,000, p = 20, resolution = 20"
  )
  panel_title_levels = unlist(lapply(panel_levels, function(row) paste(row, sweep_levels, sep = " - ")))
  data[, panel_title := factor(paste(panel, paste(sweep_short, fixed_desc, sep = "\n"), sep = " - "),
    levels = panel_title_levels)]

  p = ggplot(
    data,
    aes(x = x_display, y = time_median_plot, color = label, fill = label, group = label)
  ) +
    geom_ribbon(aes(ymin = time_q25_plot, ymax = time_q75_plot), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = label), size = 2.0, stroke = 0.2)

  if (facet_layout == "wrap") {
    p = p +
      facet_wrap(~ panel_title, ncol = 4L, scales = "free")
  } else {
    p = p +
      facet_grid(panel ~ sweep_label, scales = "free")
  }

  if (x_scale == "continuous") {
    p = p +
      scale_x_continuous(labels = comma, breaks = sort(unique(data$x_value)))
  } else {
    p = p +
      scale_x_log10(labels = comma, breaks = sort(unique(data$x_value)))
  }

  p = p +
    scale_y_continuous(labels = format_axis_number) +
    scale_color_manual(values = palette_values, breaks = names(palette_values), drop = FALSE) +
    scale_fill_manual(
      values = alpha(palette_values, 0.25), breaks = names(palette_values), drop = FALSE, guide = "none"
    ) +
    scale_shape_manual(values = shape_values, breaks = names(shape_values), drop = FALSE, guide = "none") +
    labs(title = title, x = "Value of varied parameter", y = "Median runtime (seconds)", color = NULL) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.spacing.y = unit(1.4, "lines"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold", lineheight = 0.95)
    )

  save_plot(p, filename, width = 13, height = 12)
}

plot_metric(summary_dt, "split", "regional_split_methods_linear.png", "Regional split-search runtime",
  facet_layout = "wrap", x_scale = "log10")
plot_metric(summary_dt, "total", "regional_total_methods.png", "Regional total runtime")

message("Done.")
