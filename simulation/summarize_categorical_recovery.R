#!/usr/bin/env Rscript
# Aggregate the categorical recovery benchmark and write summary figures.
#
# Sweeps are summarized into one figure each:
#   - leakage:   the leakage-strength curve (the core result)
#   - K:         cardinality scaling
#   - N:         finite-sample convergence
#   - D:         noise-feature sensitivity
#   - slope_mag: response-side SNR sensitivity
#   - group_frac: signal-group-size sensitivity (optional)
# xplaineff ALE order_method variants (raw/random/mds/pca) are the focus, plus
# the oracle_partition method (exhaustive enumeration of binary partitions
# when K <= max_partition_K) which serves as an upper bound.

args = commandArgs(trailingOnly = TRUE)
indir = "simulation/results/categorical_recovery"
figdir = "simulation/results/paper_figures"
paper_figdir = ""  # not synced to the paper by default; set explicitly if wanted
sweeps = c("leakage", "K", "N", "D", "slope_mag", "group_frac")
ttest_baseline = "random"  # paired t-test baseline for mds/pca/raw

i = 1L
while (i <= length(args)) {
  if (args[i] == "--indir" && i < length(args)) { indir = args[i + 1L]; i = i + 2L }
  else if (args[i] == "--figdir" && i < length(args)) { figdir = args[i + 1L]; i = i + 2L }
  else if (args[i] == "--paper-figdir" && i < length(args)) { paper_figdir = args[i + 1L]; i = i + 2L }
  else if (args[i] == "--ttest-baseline" && i < length(args)) {
    ttest_baseline = args[i + 1L]; i = i + 2L
  }
  else if (args[i] == "--sweeps" && i < length(args)) {
    sweeps = strsplit(args[i + 1L], ",", fixed = TRUE)[[1L]]
    i = i + 2L
  }
  else { i = i + 1L }
}

dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
if (nzchar(paper_figdir)) dir.create(paper_figdir, showWarnings = FALSE, recursive = TRUE)

library(data.table)
library(ggplot2)

method_order_display_levels = c("raw", "random", "mds", "pca", "exhaustive search")
method_order_labels = c(
  "gadget_ale|raw" = "raw",
  "gadget_ale|random" = "random",
  "gadget_ale|mds" = "mds",
  "gadget_ale|pca" = "pca",
  "gadget_ale|oracle_partition" = "exhaustive search",
  "xplaineff_ale|raw" = "raw",
  "xplaineff_ale|random" = "random",
  "xplaineff_ale|mds" = "mds",
  "xplaineff_ale|pca" = "pca",
  "xplaineff_ale|oracle_partition" = "exhaustive search"
)
method_order_levels = names(method_order_labels)

# dgp_prefix is set lazily after the first sweep CSV is loaded; it disambiguates
# output filenames so that running the summarizer for binary-slope and group mechanisms into
# the same figure directory does not overwrite each other.
dgp_prefix = ""
plot_prefix = function(dgp_type) {
  if (identical(dgp_type, "digit")) "binary_slope"
  else dgp_type
}
prefixed = function(stem) {
  if (nzchar(dgp_prefix)) sprintf("categorical_recovery_%s_%s", dgp_prefix, stem)
  else sprintf("categorical_recovery_%s", stem)
}

save_plot = function(stem, plot_obj, width, height) {
  fname = paste0(prefixed(stem), ".png")
  ggsave(file.path(figdir, fname), plot_obj, width = width, height = height, dpi = 150)
  if (nzchar(paper_figdir)) {
    ggsave(file.path(paper_figdir, fname), plot_obj, width = width, height = height, dpi = 150)
  }
}

load_sweep = function(name) {
  fn = file.path(indir, sprintf("categorical_recovery_%s.csv", name))
  if (!file.exists(fn)) return(NULL)
  d = fread(fn)
  if (!"status" %in% names(d)) d[, status := "ok"]
  if (!"error_message" %in% names(d)) d[, error_message := NA_character_]
  if (!"dgp_type" %in% names(d)) d[, dgp_type := "digit"]
  if (!nzchar(dgp_prefix)) {
    found = unique(as.character(d$dgp_type))
    if (length(found) == 1L) dgp_prefix <<- plot_prefix(found)
  }
  d[, mo := factor(paste(method, order_method, sep = "|"), levels = method_order_levels)]
  d[, mo_label := factor(method_order_labels[as.character(mo)], levels = method_order_display_levels)]
  d
}

mean_or_na = function(x) {
  x = x[!is.na(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

aggregate_sweep = function(d, xvar) {
  agg_cols = list(
    exact = quote(mean_or_na(exact_recovery)),
    ari = quote(mean_or_na(ari)),
    node_acc = quote(mean_or_na(node_acc)),
    int_imp = quote(mean_or_na(int_imp)),
    order_correct = quote(mean_or_na(order_correct)),
    n_ok = quote(sum(status == "ok", na.rm = TRUE)),
    n_no_split = quote(sum(status == "no_split", na.rm = TRUE)),
    n_fit_fail = quote(sum(status == "fit_fail", na.rm = TRUE)),
    n_fail = quote(sum(!is.na(status) & status != "ok")),
    n_rep = quote(.N)
  )
  if ("elapsed_sec" %in% names(d)) {
    agg_cols$elapsed_sec = quote(mean_or_na(elapsed_sec))
  }
  d[, eval(as.call(c(quote(list), agg_cols))), by = c("mo", "mo_label", xvar)]
}

# Per-cell paired t-test of {ari, node_acc, exact_recovery, int_imp} for each
# non-baseline method against ttest_baseline (default "random"), pairing on the
# common seed within the cell. Returns a long data.table; cells where the paired
# differences have zero variance (e.g. degenerate leakage = 1) yield NA p-values.
paired_t_table = function(d, xvar, baseline = ttest_baseline) {
  metrics = c("exact_recovery", "ari", "node_acc", "int_imp")
  baseline_label = method_order_labels[paste0("xplaineff_ale|", baseline)]
  if (is.na(baseline_label)) {
    warning("paired_t_table: baseline '", baseline, "' not found in method labels; skipping")
    return(data.table())
  }
  d_base = d[order_method == baseline, .SD, .SDcols = c(xvar, "seed", metrics)]
  setnames(d_base, metrics, paste0(metrics, "_base"))
  d_other = d[order_method != baseline]
  merged = merge(d_other, d_base, by = c(xvar, "seed"), all = FALSE)
  out = list()
  cells = unique(merged[, c(xvar, "mo", "mo_label", "order_method"), with = FALSE])
  for (k in seq_len(nrow(cells))) {
    cell_filter = merged[
      get(xvar) == cells[[xvar]][k] & order_method == cells$order_method[k]
    ]
    for (m in metrics) {
      x = as.numeric(cell_filter[[m]])
      y = as.numeric(cell_filter[[paste0(m, "_base")]])
      pair_ok = is.finite(x) & is.finite(y)
      x = x[pair_ok]; y = y[pair_ok]
      tt = if (length(x) >= 2L && stats::sd(x - y) > 0) {
        stats::t.test(x, y, paired = TRUE, alternative = "two.sided")
      } else {
        NULL
      }
      out[[length(out) + 1L]] = data.table(
        sweep_axis = xvar,
        sweep_value = cells[[xvar]][k],
        order_method = cells$order_method[k],
        mo_label = cells$mo_label[k],
        baseline = baseline,
        metric = m,
        n_pairs = length(x),
        mean_diff = if (length(x)) mean(x - y) else NA_real_,
        t_stat = if (!is.null(tt)) as.numeric(tt$statistic) else NA_real_,
        p_value = if (!is.null(tt)) as.numeric(tt$p.value) else NA_real_,
        ci_low = if (!is.null(tt)) as.numeric(tt$conf.int[1L]) else NA_real_,
        ci_high = if (!is.null(tt)) as.numeric(tt$conf.int[2L]) else NA_real_
      )
    }
  }
  rbindlist(out, fill = TRUE)
}

# A sweep figure: exact recovery (solid) + ARI (dashed) vs the swept axis.
# subtitle_fixed: describes the non-swept axes held constant at base values.
sweep_plot = function(d, xvar, xlab, title, subtitle_fixed = "", logx = FALSE,
    x_breaks = NULL, x_labels = waiver(), x_text_angle = 0) {
  agg = aggregate_sweep(d, xvar)
  long = melt(agg, id.vars = c("mo", "mo_label", xvar),
    measure.vars = c("exact", "ari"),
    variable.name = "metric", value.name = "value")
  long[, metric := factor(metric, levels = c("exact", "ari"),
    labels = c("Exact recovery", "Adjusted Rand index"))]
  p = ggplot(long, aes(x = .data[[xvar]], y = value, color = mo_label, linetype = metric)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.0) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = title, subtitle = subtitle_fixed, x = xlab, y = "Recovery",
         color = NULL, linetype = NULL) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  if (logx) {
    p = p + scale_x_log10()
  } else if (!is.null(x_breaks)) {
    p = p + scale_x_continuous(
      breaks = x_breaks, labels = x_labels
    )
  }
  if (x_text_angle != 0) {
    p = p + theme(axis.text.x = element_text(angle = x_text_angle, hjust = 1, vjust = 1))
  }
  p
}

d_leakage = load_sweep("leakage")
d_K = load_sweep("K")
d_N = load_sweep("N")
d_D = load_sweep("D")
d_slope_mag = load_sweep("slope_mag")
d_group_frac = load_sweep("group_frac")

# write a combined long summary table
summ = rbindlist(list(
  if ("leakage" %in% sweeps && !is.null(d_leakage)) aggregate_sweep(d_leakage, "leakage")[, sweep := "leakage"],
  if ("K" %in% sweeps && !is.null(d_K)) aggregate_sweep(d_K, "K")[, sweep := "K"],
  if ("N" %in% sweeps && !is.null(d_N)) aggregate_sweep(d_N, "N")[, sweep := "N"],
  if ("D" %in% sweeps && !is.null(d_D)) aggregate_sweep(d_D, "D")[, sweep := "D"],
  if ("slope_mag" %in% sweeps && !is.null(d_slope_mag)) {
    aggregate_sweep(d_slope_mag, "slope_mag")[, sweep := "slope_mag"]
  },
  if ("group_frac" %in% sweeps && !is.null(d_group_frac)) {
    aggregate_sweep(d_group_frac, "group_frac")[, sweep := "group_frac"]
  }
), fill = TRUE)
fwrite(summ, file.path(indir, paste0(prefixed("summary"), ".csv")))
message("Written: ", file.path(indir, paste0(prefixed("summary"), ".csv")))

# write paired t-test table (mds/pca/raw vs random by default)
ttest_rows = list()
if ("leakage" %in% sweeps && !is.null(d_leakage)) {
  ttest_rows[[length(ttest_rows) + 1L]] = paired_t_table(d_leakage, "leakage")
}
if ("K" %in% sweeps && !is.null(d_K)) ttest_rows[[length(ttest_rows) + 1L]] = paired_t_table(d_K, "K")
if ("N" %in% sweeps && !is.null(d_N)) ttest_rows[[length(ttest_rows) + 1L]] = paired_t_table(d_N, "N")
if ("D" %in% sweeps && !is.null(d_D)) ttest_rows[[length(ttest_rows) + 1L]] = paired_t_table(d_D, "D")
if ("slope_mag" %in% sweeps && !is.null(d_slope_mag)) {
  ttest_rows[[length(ttest_rows) + 1L]] = paired_t_table(d_slope_mag, "slope_mag")
}
if ("group_frac" %in% sweeps && !is.null(d_group_frac)) {
  ttest_rows[[length(ttest_rows) + 1L]] = paired_t_table(d_group_frac, "group_frac")
}
ttest_tbl = rbindlist(ttest_rows, fill = TRUE)
if (nrow(ttest_tbl)) {
  fwrite(ttest_tbl, file.path(indir, paste0(prefixed("paired_ttest"), ".csv")))
  message("Written: ", file.path(indir, paste0(prefixed("paired_ttest"), ".csv")),
    "  (baseline = ", ttest_baseline, ")")
}

if ("leakage" %in% sweeps && !is.null(d_leakage)) {
  p = sweep_plot(d_leakage, "leakage",
    expression(paste("Leakage strength ", lambda["leak"], " (0 = none, 1 = deterministic)")),
    "Level-ordering diagnostic vs leakage strength",
    subtitle_fixed = expression(paste("Sweeping ", lambda["leak"], "; fixed M = 6, n = 2000, p = 10, ",
                                      beta[max], " = 4")),
    x_breaks = c(0, 0.05, 0.075, 0.10, 0.25, 0.50, 1.00),
    x_labels = c("0", "0.05", "0.075", "0.10", "0.25", "0.50", "1.00"),
    x_text_angle = 30)
  save_plot("leakage", p, width = 10, height = 5.5)
}
if ("K" %in% sweeps && !is.null(d_K)) {
  p = sweep_plot(d_K, "K", "Number of levels M",
    "Level-ordering diagnostic vs cardinality",
    subtitle_fixed = expression(paste("Sweeping M; fixed n = 2000, p = 10, ", lambda["leak"], " = 0.1, ",
                                      beta[max], " = 4")),
    x_breaks = c(6, 8, 12, 20))
  save_plot("K", p, width = 10, height = 5.5)
}
if ("N" %in% sweeps && !is.null(d_N)) {
  p = sweep_plot(d_N, "N", "Sample size n",
    "Level-ordering diagnostic vs sample size",
    subtitle_fixed = expression(paste("Sweeping n; fixed M = 6, p = 10, ", lambda["leak"], " = 0.1, ",
                                      beta[max], " = 4")),
    logx = TRUE)
  save_plot("N", p, width = 10, height = 5.5)
}
if ("D" %in% sweeps && !is.null(d_D)) {
  p = sweep_plot(d_D, "D", "Number of features p",
    "Level-ordering diagnostic vs noise features",
    subtitle_fixed = expression(paste("Sweeping p; fixed M = 6, n = 2000, ", lambda["leak"], " = 0.1, ",
                                      beta[max], " = 4")),
    logx = TRUE)
  save_plot("D", p, width = 10, height = 5.5)
}
if ("slope_mag" %in% sweeps && !is.null(d_slope_mag)) {
  p = sweep_plot(d_slope_mag, "slope_mag",
    expression(paste("Response-side slope magnitude ", beta[max])),
    "Level-ordering diagnostic vs slope size",
    subtitle_fixed = expression(paste("Sweeping ", beta[max], "; fixed M = 6, n = 2000, p = 10, ",
                                      lambda["leak"], " = 0.1")))
  save_plot("slope_mag", p, width = 10, height = 5.5)
}
if ("group_frac" %in% sweeps && !is.null(d_group_frac)) {
  p = sweep_plot(d_group_frac, "group_frac", "Signal-level fraction",
    "Level-ordering diagnostic vs signal-group size",
    subtitle_fixed = expression(paste("Sweeping signal fraction; fixed M = 6, n = 2000, p = 10, ",
                                      lambda["leak"], " = 0.1, ", beta[max], " = 4")))
  save_plot("group_frac", p, width = 10, height = 5.5)
}

message("Written figures to ", figdir)
