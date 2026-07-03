#!/usr/bin/env Rscript
# Diagnostic benchmark for ranger prediction sensitivity to PDP stacked-data layout.

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

if (!file.exists("DESCRIPTION") || readLines("DESCRIPTION", 1L) != "Package: xplaineff") {
  if (file.exists("../DESCRIPTION") && readLines("../DESCRIPTION", 1L) == "Package: xplaineff") {
    setwd("..")
  } else {
    stop("Run from xplaineff package root")
  }
}

if (!requireNamespace("ranger", quietly = TRUE)) {
  stop("Install ranger for this benchmark")
}

library(data.table)
setDTthreads(1L)

args = commandArgs(trailingOnly = TRUE)
run_id = format(Sys.time(), "%Y%m%d_%H%M%S")
run_root = file.path("simulation/results/runtime_runs", run_id)
datadir = "simulation/data/global_r_runtime"
outdir = file.path(run_root, "ranger_layout_sensitivity")
cells_arg = "10000:20:20,10000:20:50"
features_arg = "x1,x5"
reps = 30L
warmup = 2L
seed = 21L

parse_cells = function(x) {
  parts = trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  rows = lapply(parts[nzchar(parts)], function(part) {
    vals = as.integer(strsplit(part, ":", fixed = TRUE)[[1L]])
    if (length(vals) != 3L || anyNA(vals)) {
      stop("Cells must use N:D:G, e.g. 10000:20:20,10000:20:50")
    }
    data.frame(N = vals[1L], D = vals[2L], n_grid = vals[3L])
  })
  do.call(rbind, rows)
}

parse_chr_vec = function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])

i = 1L
while (i <= length(args)) {
  if (args[i] == "--datadir" && i < length(args)) {
    datadir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--outdir" && i < length(args)) {
    outdir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--cells" && i < length(args)) {
    cells_arg = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--features" && i < length(args)) {
    features_arg = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--reps" && i < length(args)) {
    reps = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--warmup" && i < length(args)) {
    warmup = as.integer(args[i + 1L]); i = i + 2L
  } else {
    i = i + 1L
  }
}

cells = parse_cells(cells_arg)
features = parse_chr_vec(features_arg)
features = features[nzchar(features)]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

load_data = function(N, D) {
  path = file.path(datadir, sprintf("benchmark_N%d_D%d_seed%d.csv", N, D, seed))
  if (!file.exists(path)) stop("Missing data file: ", path)
  as.data.frame(fread(path))
}

rf_config = list(
  num_trees = 100L,
  min_node_size = 1L,
  replace = TRUE,
  sample_fraction = 1.0,
  splitrule = "variance",
  respect_unordered_factors = "ignore",
  num_threads = 1L,
  seed = 21L
)

fit_rf = function(dat) {
  p = ncol(dat) - 1L
  ranger::ranger(
    y ~ .,
    data = dat,
    num.trees = rf_config$num_trees,
    mtry = p,
    min.node.size = rf_config$min_node_size,
    replace = rf_config$replace,
    sample.fraction = rf_config$sample_fraction,
    splitrule = rf_config$splitrule,
    respect.unordered.factors = rf_config$respect_unordered_factors,
    num.threads = rf_config$num_threads,
    seed = rf_config$seed
  )
}

grid_for_feature = function(x, n) {
  sort(unique(as.numeric(stats::quantile(x, probs = seq(0, 1, length.out = n), type = 7, na.rm = TRUE))))
}

build_grid_major = function(X, feature, grid) {
  n = nrow(X)
  out = as.data.frame(lapply(X, rep, times = length(grid)), optional = TRUE)
  out[[feature]] = rep(grid, each = n)
  row.names(out) = NULL
  out
}

build_observation_major = function(X, feature, grid) {
  n = nrow(X)
  out = X[rep(seq_len(n), each = length(grid)), , drop = FALSE]
  out[[feature]] = rep(grid, times = n)
  row.names(out) = NULL
  out
}

predict_once = function(model, newdata) {
  invisible(predict(model, data = newdata, num.threads = 1L)$predictions)
}

run_cell = function(N, D, n_grid) {
  message(sprintf("=== N=%d D=%d G=%d ===", N, D, n_grid))
  dat = load_data(N, D)
  X = dat[, setdiff(names(dat), "y"), drop = FALSE]
  model = fit_rf(dat)
  rows = list()
  for (feature in features) {
    if (!feature %in% names(X)) {
      stop("Feature ", feature, " is missing from the N=", N, ", D=", D, " data.")
    }
    grid = grid_for_feature(X[[feature]], n_grid)
    newdata = list(
      grid_major = build_grid_major(X, feature, grid),
      observation_major = build_observation_major(X, feature, grid)
    )
    for (layout in names(newdata)) {
      for (j in seq_len(warmup)) predict_once(model, newdata[[layout]])
    }
    set.seed(seed + N + D + n_grid)
    order_index = 0L
    for (r in seq_len(reps)) {
      for (layout in sample(names(newdata))) {
        order_index = order_index + 1L
        elapsed = as.numeric(system.time(predict_once(model, newdata[[layout]]))[["elapsed"]])
        rows[[length(rows) + 1L]] = data.frame(
          N = N,
          D = D,
          n_grid = length(grid),
          feature = feature,
          layout = layout,
          repetition = r,
          order_index = order_index,
          predict_time_sec = elapsed,
          nrow_newdata = nrow(newdata[[layout]]),
          ncol_newdata = ncol(newdata[[layout]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  rbindlist(rows)
}

raw = rbindlist(lapply(seq_len(nrow(cells)), function(i) run_cell(cells$N[i], cells$D[i], cells$n_grid[i])))
summary = raw[, .(
  predict_median = stats::median(predict_time_sec),
  predict_q25 = stats::quantile(predict_time_sec, 0.25, names = FALSE),
  predict_q75 = stats::quantile(predict_time_sec, 0.75, names = FALSE),
  predict_mean = mean(predict_time_sec),
  predict_sd = stats::sd(predict_time_sec),
  n_rep = .N
), by = .(N, D, n_grid, feature, layout)]
summary[, rel_to_grid_major := predict_median / predict_median[layout == "grid_major"], by = .(N, D, n_grid, feature)]
setorder(summary, N, D, n_grid, feature, layout)

raw_path = file.path(outdir, "ranger_layout_sensitivity.csv")
summary_path = file.path(outdir, "ranger_layout_sensitivity_summary.csv")
fwrite(raw, raw_path)
fwrite(summary, summary_path)
message("Written: ", raw_path)
message("Written: ", summary_path)
print(summary[, .(
  N, D, n_grid, feature, layout,
  predict_median = round(predict_median, 4L),
  rel_to_grid_major = round(rel_to_grid_major, 3L)
)])
