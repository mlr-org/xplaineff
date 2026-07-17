#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
run_dir = ""
summary_file = ""
figdir = ""
paper_figdir = "paper/figures"
tag = "server_effector040_compare"
sync_paper = TRUE

usage = function() {
  cat(
    "Usage:\n",
    "  Rscript simulation/plot_regional_runtime_paperformat.R \\\n",
    "    --run-dir simulation/results/runtime_runs/<RUN_ID> \\\n",
    "    --tag server_effector040_compare\n\n",
    "Options:\n",
    "  --run-dir       Runtime run directory containing regional_runtime/.\n",
    "  --summary       Summary CSV to plot. Overrides --run-dir discovery.\n",
    "  --figdir        Output directory for PNGs. Defaults to <run-dir>/figures.\n",
    "  --paper-figdir  Paper figure sync directory. Defaults to paper/figures.\n",
    "  --tag           Output filename tag. Defaults to server_effector040_compare.\n",
    "  --no-paper-sync Do not sync PNGs to paper/figures.\n",
    "  --help          Show this message.\n",
    sep = ""
  )
}

i = 1L
while (i <= length(args)) {
  if (args[i] %in% c("--help", "-h")) {
    usage()
    quit(status = 0L)
  } else if (args[i] == "--run-dir" && i < length(args)) {
    run_dir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--summary" && i < length(args)) {
    summary_file = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--figdir" && i < length(args)) {
    figdir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--paper-figdir" && i < length(args)) {
    paper_figdir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--tag" && i < length(args)) {
    tag = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--no-paper-sync") {
    sync_paper = FALSE; i = i + 1L
  } else {
    stop("Unknown or incomplete argument: ", args[i], call. = FALSE)
  }
}

if (!nzchar(summary_file)) {
  if (!nzchar(run_dir)) {
    usage()
    stop("Please provide --run-dir or --summary.", call. = FALSE)
  }
  regional_dir = file.path(run_dir, "regional_runtime")
  summary_candidates = list.files(
    regional_dir,
    pattern = "^summary_.*\\.csv$|^summary\\.csv$",
    full.names = TRUE
  )
  if (!length(summary_candidates)) {
    stop("No summary CSV found in ", regional_dir, call. = FALSE)
  }
  summary_file = summary_candidates[grepl("^summary_", basename(summary_candidates))][1L]
  if (is.na(summary_file)) {
    summary_file = summary_candidates[1L]
  }
}
if (!nzchar(figdir)) {
  inferred_run_dir = if (nzchar(run_dir)) {
    run_dir
  } else if (basename(dirname(summary_file)) == "regional_runtime") {
    dirname(dirname(summary_file))
  } else {
    dirname(summary_file)
  }
  figdir = file.path(inferred_run_dir, "figures")
}

dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
if (sync_paper && nzchar(paper_figdir)) {
  dir.create(paper_figdir, showWarnings = FALSE, recursive = TRUE)
}

library(data.table)
library(ggplot2)
library(scales)
setDTthreads(1L)

summary_dt = fread(summary_file)
required_cols = c(
  "effect", "model_type", "sub_experiment", "N", "D", "resolution", "n_split", "package",
  "precompute_median", "precompute_q25", "precompute_q75",
  "split_median", "split_q25", "split_q75",
  "total_median", "total_q25", "total_q75"
)
missing_cols = setdiff(required_cols, names(summary_dt))
if (length(missing_cols)) {
  stop(
    "Summary CSV is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

palette_values = c(
  "effector" = "#2ca02c",
  "xplaineff" = "#1f77b4"
)
shape_values = c(
  "effector" = 17,
  "xplaineff" = 16
)
format_axis_number = function(x) {
  format(x, trim = TRUE, digits = 3L, scientific = FALSE, big.mark = ",")
}

format_x_breaks = function(limits) {
  max_limit = max(limits, na.rm = TRUE)
  if (max_limit <= 15) {
    c(2, 5, 8, 10)
  } else if (max_limit <= 60) {
    c(10, 20, 50)
  } else if (max_limit <= 150) {
    c(10, 20, 50, 100)
  } else {
    c(1000, 5000, 10000, 20000)
  }
}

save_plot = function(plot_obj, filename, width, height, dpi = 220L) {
  out = file.path(figdir, filename)
  ggsave(out, plot_obj, width = width, height = height, dpi = dpi)
  message("Written: ", out)
  if (sync_paper && nzchar(paper_figdir) && normalizePath(figdir) != normalizePath(paper_figdir)) {
    paper_out = file.path(paper_figdir, filename)
    ggsave(paper_out, plot_obj, width = width, height = height, dpi = dpi)
    message("Synced: ", paper_out)
  }
}

plot_metric = function(data, metric, filename, title, y_label, include_split_sweep, ncol, width, height) {
  data = copy(data)
  if (!include_split_sweep) {
    data = data[sub_experiment != "vs_split"]
  }
  if (!nrow(data)) {
    return(invisible(NULL))
  }

  median_col = sprintf("%s_median", metric)
  q25_col = sprintf("%s_q25", metric)
  q75_col = sprintf("%s_q75", metric)
  data[, time_median := get(median_col)]
  data[, time_q25 := get(q25_col)]
  data[, time_q75 := get(q75_col)]
  data[, display_label := factor(fifelse(package == "effector", "effector", "xplaineff"),
    levels = names(palette_values))]
  data[, effect_label := fifelse(effect == "pdp", "PDP", "ALE")]
  data[, model_label := fifelse(model_type == "rf", "bagged trees", "toy model")]
  data[, panel := factor(paste(effect_label, model_label, sep = " / "), levels = c(
    "PDP / bagged trees",
    "PDP / toy model",
    "ALE / bagged trees",
    "ALE / toy model"
  ))]
  data[, sweep_short := fcase(
    sub_experiment == "vs_N", "Sample size n",
    sub_experiment == "vs_D", "Feature dimension p",
    sub_experiment == "vs_res", "PDP grid / ALE intervals",
    sub_experiment == "vs_split", "Number of splits",
    default = sub_experiment
  )]
  data[, fixed_desc := fcase(
    sub_experiment == "vs_N" & include_split_sweep, "p = 20, K = 20, n_split = 2",
    sub_experiment == "vs_D" & include_split_sweep, "n = 10,000, K = 20, n_split = 2",
    sub_experiment == "vs_res" & include_split_sweep, "n = 10,000, p = 20, n_split = 2",
    sub_experiment == "vs_split" & include_split_sweep, "n = 10,000, p = 20, K = 20",
    sub_experiment == "vs_N", "p = 20, K = 20",
    sub_experiment == "vs_D", "n = 10,000, K = 20",
    sub_experiment == "vs_res", "n = 10,000, p = 20",
    default = sub_experiment
  )]
  data[, x_value := fcase(
    sub_experiment == "vs_N", as.numeric(N),
    sub_experiment == "vs_D", as.numeric(D),
    sub_experiment == "vs_res", as.numeric(resolution),
    sub_experiment == "vs_split", as.numeric(n_split),
    default = as.numeric(NA)
  )]

  sweep_levels = if (include_split_sweep) {
    c(
      "Sample size n\np = 20, K = 20, n_split = 2",
      "Feature dimension p\nn = 10,000, K = 20, n_split = 2",
      "PDP grid / ALE intervals\nn = 10,000, p = 20, n_split = 2",
      "Number of splits\nn = 10,000, p = 20, K = 20"
    )
  } else {
    c(
      "Sample size n\np = 20, K = 20",
      "Feature dimension p\nn = 10,000, K = 20",
      "PDP grid / ALE intervals\nn = 10,000, p = 20"
    )
  }
  panel_levels = levels(data$panel)
  panel_title_levels = unlist(lapply(panel_levels, function(row) paste(row, sweep_levels, sep = " - ")))
  data[, panel_title := factor(paste(panel, paste(sweep_short, fixed_desc, sep = "\n"), sep = " - "),
    levels = panel_title_levels)]

  p = ggplot(
    data,
    aes(x = x_value, y = time_median, color = display_label, fill = display_label, group = display_label)
  ) +
    geom_ribbon(aes(ymin = time_q25, ymax = time_q75), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = display_label), size = 2.0, stroke = 0.2) +
    facet_wrap(~ panel_title, ncol = ncol, scales = "free") +
    scale_x_continuous(labels = comma, breaks = format_x_breaks) +
    scale_y_continuous(labels = format_axis_number) +
    scale_color_manual(values = palette_values, breaks = names(palette_values), drop = FALSE) +
    scale_fill_manual(
      values = alpha(palette_values, 0.25), breaks = names(palette_values), drop = FALSE, guide = "none"
    ) +
    scale_shape_manual(values = shape_values, breaks = names(shape_values), drop = FALSE, guide = "none") +
    labs(title = title, x = "Value of varied parameter", y = y_label, color = NULL) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.spacing.y = unit(1.4, "lines"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = "bold", lineheight = 0.95)
    )

  save_plot(p, filename, width = width, height = height)
}

plot_metric(
  summary_dt,
  metric = "precompute",
  filename = sprintf("regional_precompute_runtime_%s.png", tag),
  title = "Regional global-effect precompute runtime",
  y_label = "Median global-effect precompute runtime (seconds)",
  include_split_sweep = FALSE,
  ncol = 3L,
  width = 11,
  height = 10.4
)
plot_metric(
  summary_dt,
  metric = "split",
  filename = sprintf("regional_split_runtime_%s.png", tag),
  title = "Regional split-search runtime",
  y_label = "Median split-search runtime (seconds)",
  include_split_sweep = TRUE,
  ncol = 4L,
  width = 13.6,
  height = 10.6
)
plot_metric(
  summary_dt,
  metric = "total",
  filename = sprintf("regional_total_runtime_%s.png", tag),
  title = "Regional total runtime",
  y_label = "Median total runtime (seconds)",
  include_split_sweep = TRUE,
  ncol = 4L,
  width = 13.6,
  height = 10.6
)
