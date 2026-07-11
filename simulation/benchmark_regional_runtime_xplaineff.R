#!/usr/bin/env Rscript
# xplaineff regional PDP/ALE runtime benchmark.
# Measures the regional precompute, split-search, and total timings separately.

Sys.setenv(
  OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS", "1"),
  OMP_THREAD_LIMIT = Sys.getenv("OMP_THREAD_LIMIT", "1"),
  OMP_PROC_BIND = Sys.getenv("OMP_PROC_BIND", "FALSE"),
  KMP_INIT_AT_FORK = Sys.getenv("KMP_INIT_AT_FORK", "FALSE"),
  KMP_AFFINITY = Sys.getenv("KMP_AFFINITY", "disabled"),
  OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS", "1"),
  MKL_NUM_THREADS = Sys.getenv("MKL_NUM_THREADS", "1"),
  VECLIB_MAXIMUM_THREADS = Sys.getenv("VECLIB_MAXIMUM_THREADS", "1"),
  NUMEXPR_NUM_THREADS = Sys.getenv("NUMEXPR_NUM_THREADS", "1"),
  DATATABLE_NUM_THREADS = Sys.getenv("DATATABLE_NUM_THREADS", "1"),
  RCPP_PARALLEL_NUM_THREADS = Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "1")
)

parse_flag = function(x) {
  s = tolower(trimws(as.character(x)))
  if (s %in% c("true", "1", "yes", "y", "on")) return(TRUE)
  if (s %in% c("false", "0", "no", "n", "off")) return(FALSE)
  FALSE
}

load_xplaineff_for_benchmark = function() {
  use_load_all = parse_flag(Sys.getenv("XPLAINEFF_BENCH_LOAD_ALL", "true"))
  if (use_load_all) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Install pkgload or set XPLAINEFF_BENCH_LOAD_ALL=false to use the installed xplaineff package")
    }
    pkgload::load_all(".", quiet = TRUE)
    message("Loaded xplaineff from local source with pkgload::load_all().")
  } else if (!requireNamespace("xplaineff", quietly = TRUE)) {
    stop("Install xplaineff or set XPLAINEFF_BENCH_LOAD_ALL=true")
  } else {
    library(xplaineff)
    message("Loaded installed xplaineff package.")
  }
}

if (!file.exists("DESCRIPTION") || readLines("DESCRIPTION", 1L) != "Package: xplaineff") {
  if (file.exists("../DESCRIPTION") && readLines("../DESCRIPTION", 1L) == "Package: xplaineff") {
    setwd("..")
  } else {
    stop("Run from xplaineff package root")
  }
}

load_xplaineff_for_benchmark()

library(data.table)
setDTthreads(1L)

args = commandArgs(trailingOnly = TRUE)
run_id = format(Sys.time(), "%Y%m%d_%H%M%S")
run_root = file.path("simulation/results/runtime_runs", run_id)
datadir = "simulation/data/global_r_runtime"
outdir = file.path(run_root, "regional_runtime")
reps = 30L
N_vec = c(5000L, 10000L, 25000L, 50000L)
D_vec = c(10L, 20L, 50L, 100L)
fixed_N = 10000L
fixed_D = 20L
default_resolution = 20L
resolution_vec = c(10L, 20L, 50L)
default_n_split = 2L
n_split_vec = c(2L, 5L, 8L, 10L)
min_node_size = 50L
n_quantiles = 19L
fail_fast = FALSE
model_types = c("rf", "toy")
sub_experiments = c("vs_N", "vs_D", "vs_res", "vs_split")
output_suffix = ""

parse_int_vec = function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
parse_chr_vec = function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
parse_nullable_int = function(x) {
  x = tolower(trimws(as.character(x)))
  if (x %in% c("", "na", "null", "none", "all")) return(NULL)
  as.integer(x)
}

i = 1L
while (i <= length(args)) {
  if (args[i] == "--datadir" && i < length(args)) {
    datadir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--outdir" && i < length(args)) {
    outdir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--reps" && i < length(args)) {
    reps = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--N-vec" && i < length(args)) {
    N_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--D-vec" && i < length(args)) {
    D_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--fixed-N" && i < length(args)) {
    fixed_N = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--fixed-D" && i < length(args)) {
    fixed_D = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--resolution" && i < length(args)) {
    default_resolution = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--resolution-vec" && i < length(args)) {
    resolution_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-split" && i < length(args)) {
    default_n_split = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-split-vec" && i < length(args)) {
    n_split_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--min-node-size" && i < length(args)) {
    min_node_size = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-quantiles" && i < length(args)) {
    n_quantiles = parse_nullable_int(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--fail-fast" && i < length(args)) {
    fail_fast = parse_flag(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--models" && i < length(args)) {
    model_types = parse_chr_vec(args[i + 1L]); model_types = model_types[nzchar(model_types)]; i = i + 2L
  } else if (args[i] == "--sub-experiments" && i < length(args)) {
    sub_experiments = parse_chr_vec(args[i + 1L]); sub_experiments = sub_experiments[nzchar(sub_experiments)]
    i = i + 2L
  } else if (args[i] == "--output-suffix" && i < length(args)) {
    output_suffix = trimws(as.character(args[i + 1L])); i = i + 2L
  } else {
    i = i + 1L
  }
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

invalid_model_types = setdiff(model_types, c("rf", "toy"))
if (length(invalid_model_types)) {
  stop("Unsupported regional model type(s): ", paste(invalid_model_types, collapse = ", "), call. = FALSE)
}

load_data = function(N, D) {
  path = file.path(datadir, sprintf("benchmark_N%d_D%d_seed21.csv", N, D))
  if (!file.exists(path)) return(NULL)
  as.data.frame(fread(path))
}

toy_pred_fun = function(model, newdata) {
  x1 = newdata[["x1"]]
  x2 = newdata[["x2"]]
  x3 = newdata[["x3"]]
  5 * x1 + 5 * x2 + ifelse(x3 > 0, 10 * x1, 0) - ifelse(x3 > 0, 10 * x2, 0)
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
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Install ranger for the RF benchmark")
  }
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

rf_pred_fun = function(model, newdata) {
  as.numeric(predict(model, data = newdata, num.threads = 1L)$predictions)
}

fit_model = function(dat, model_type) {
  if (identical(model_type, "rf")) fit_rf(dat) else "toy"
}

predict_for_model = function(model_type) {
  if (identical(model_type, "rf")) rf_pred_fun else toy_pred_fun
}

xplaineff_predict_fun = function(model, pred_fun) {
  if (inherits(model, "ranger")) {
    return(NULL)
  }
  pred_fun
}

make_cells = function() {
  rows = list()
  make_cell = function(sub_experiment, N, D, resolution, n_split) {
    data.frame(
      sub_experiment = sub_experiment,
      N = as.integer(N),
      D = as.integer(D),
      resolution = as.integer(resolution),
      n_split = as.integer(n_split)
    )
  }
  if ("vs_N" %in% sub_experiments) {
    rows[[length(rows) + 1L]] = do.call(rbind, lapply(
      N_vec,
      function(N) make_cell("vs_N", N, fixed_D, default_resolution, default_n_split)
    ))
  }
  if ("vs_D" %in% sub_experiments) {
    rows[[length(rows) + 1L]] = do.call(rbind, lapply(
      D_vec,
      function(D) make_cell("vs_D", fixed_N, D, default_resolution, default_n_split)
    ))
  }
  if ("vs_res" %in% sub_experiments) {
    rows[[length(rows) + 1L]] = do.call(rbind, lapply(
      resolution_vec,
      function(resolution) make_cell("vs_res", fixed_N, fixed_D, resolution, default_n_split)
    ))
  }
  if ("vs_split" %in% sub_experiments) {
    rows[[length(rows) + 1L]] = do.call(rbind, lapply(
      n_split_vec,
      function(n_split) make_cell("vs_split", fixed_N, fixed_D, default_resolution, n_split)
    ))
  }
  unique(do.call(rbind, rows))
}

run_regional = function(effect, dat, model, pred_fun, resolution, n_split) {
  pred_fun = xplaineff_predict_fun(model, pred_fun)
  if (identical(effect, "pdp")) {
    strat = PdStrategy$new()
    tree = GadgetTree$new(
      strategy = strat,
      n_split = n_split,
      min_node_size = min_node_size,
      n_quantiles = n_quantiles
    )
    tree$fit(
      data = dat,
      target_feature_name = "y",
      model = model,
      predict_fun = pred_fun,
      n_grid = resolution,
      pd_engine = "cpp"
    )
  } else {
    strat = AleStrategy$new()
    tree = GadgetTree$new(
      strategy = strat,
      n_split = n_split,
      min_node_size = min_node_size,
      n_quantiles = n_quantiles
    )
    tree$fit(
      data = dat,
      target_feature_name = "y",
      model = model,
      n_intervals = resolution,
      predict_fun = pred_fun,
      order_method = "raw",
      ale_engine = "cpp"
    )
  }
  c(
    precompute = strat$fit_timing$global,
    split = strat$fit_timing$regional,
    total = strat$fit_timing$global + strat$fit_timing$regional
  )
}

record_row = function(rows, model_type, effect, cell, repetition, timing = NULL,
  status = "ok", error_message = NA_character_) {
  rows[[length(rows) + 1L]] = data.frame(
    module = "regional_runtime",
    package = "xplaineff",
    impl = "cpp",
    effect = effect,
    method = sprintf("regional_%s", effect),
    model_type = model_type,
    sub_experiment = cell$sub_experiment,
    N = cell$N,
    D = cell$D,
    resolution = cell$resolution,
    n_grid = if (identical(effect, "pdp")) cell$resolution else NA_integer_,
    n_intervals = if (identical(effect, "ale")) cell$resolution else NA_integer_,
    n_split = cell$n_split,
    n_quantiles = if (is.null(n_quantiles)) NA_integer_ else n_quantiles,
    n_candidates = if (is.null(n_quantiles)) NA_integer_ else n_quantiles,
    split_candidate_rule = if (is.null(n_quantiles)) "all_unique" else "quantile",
    repetition = repetition,
    precompute_time_sec = if (is.null(timing)) NA_real_ else timing[["precompute"]],
    split_time_sec = if (is.null(timing)) NA_real_ else timing[["split"]],
    total_time_sec = if (is.null(timing)) NA_real_ else timing[["total"]],
    status = status,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
  rows
}

run_model = function(model_type) {
  cells = make_cells()
  data_keys = unique(cells[, c("N", "D"), drop = FALSE])
  data_cache = list()
  model_cache = list()
  for (i in seq_len(nrow(data_keys))) {
    key = sprintf("N%d_D%d", data_keys$N[i], data_keys$D[i])
    dat = load_data(data_keys$N[i], data_keys$D[i])
    if (is.null(dat)) next
    data_cache[[key]] = dat
    model_cache[[key]] = fit_model(dat, model_type)
  }

  rows = list()
  pred_fun = predict_for_model(model_type)
  for (effect in c("pdp", "ale")) {
    for (i in seq_len(nrow(cells))) {
      cell = cells[i, , drop = FALSE]
      key = sprintf("N%d_D%d", cell$N, cell$D)
      if (!key %in% names(data_cache)) next
      dat = data_cache[[key]]
      model = model_cache[[key]]
      log_msg = sprintf(
        "[%s] xplaineff regional %s %s N=%d D=%d res=%d n_split=%d n_quantiles=%s",
        model_type, effect, cell$sub_experiment, cell$N, cell$D, cell$resolution, cell$n_split,
        if (is.null(n_quantiles)) "NULL" else as.character(n_quantiles)
      )
      message(log_msg, " | start")
      tryCatch(
        run_regional(effect, dat, model, pred_fun, cell$resolution, cell$n_split),
        error = function(e) {
          if (isTRUE(fail_fast)) stop(e)
          message(log_msg, " | warmup skipped/failed: ", conditionMessage(e))
          invisible(NULL)
        }
      )
      for (r in seq_len(reps)) {
        out = tryCatch(
          list(ok = TRUE, timing = run_regional(effect, dat, model, pred_fun, cell$resolution, cell$n_split),
            error = NA_character_),
          error = function(e) {
            if (isTRUE(fail_fast)) stop(e)
            list(ok = FALSE, timing = NULL, error = conditionMessage(e))
          }
        )
        rows = record_row(
          rows = rows,
          model_type = model_type,
          effect = effect,
          cell = cell,
          repetition = r,
          timing = out$timing,
          status = if (out$ok) "ok" else "error",
          error_message = out$error
        )
      }
      message(log_msg, " | done")
    }
  }
  rows
}

out_filename = function(stem) {
  if (nzchar(output_suffix)) sprintf("%s_%s.csv", stem, output_suffix) else sprintf("%s.csv", stem)
}

all_rows = list()
for (model_type in model_types) {
  all_rows = c(all_rows, run_model(model_type))
}

out = if (length(all_rows)) rbindlist(all_rows, fill = TRUE) else data.table()
fn = out_filename("regional_runtime_xplaineff")
fwrite(out, file.path(outdir, fn))
message("Written: ", fn)
