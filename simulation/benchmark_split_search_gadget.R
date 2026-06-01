#!/usr/bin/env Rscript
# GADGET split-search benchmark. Global precomputation is excluded from time_sec.

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

load_gadget_for_benchmark = function() {
  use_load_all = parse_flag(Sys.getenv("GADGET_BENCH_LOAD_ALL", "true"))
  if (use_load_all) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Install pkgload or set GADGET_BENCH_LOAD_ALL=false to use the installed gadget package")
    }
    pkgload::load_all(".", quiet = TRUE)
    message("Loaded gadget from local source with pkgload::load_all().")
  } else if (!requireNamespace("gadget", quietly = TRUE)) {
    stop("Install gadget or set GADGET_BENCH_LOAD_ALL=true")
  } else {
    library(gadget)
    message("Loaded installed gadget package.")
  }
}

if (!file.exists("DESCRIPTION") || readLines("DESCRIPTION", 1L) != "Package: gadget") {
  if (file.exists("../DESCRIPTION") && readLines("../DESCRIPTION", 1L) == "Package: gadget") {
    setwd("..")
  } else {
    stop("Run from GADGET package root")
  }
}

load_gadget_for_benchmark()

library(data.table)
setDTthreads(1L)

args = commandArgs(trailingOnly = TRUE)
datadir = "simulation/data/global_r_runtime"
outdir = "simulation/results/split_search_runtime"
reps = 20L
N_vec = c(500L, 1000L, 2500L, 5000L)
D = 10L
n_grid = 20L
n_intervals = 20L
fail_fast = FALSE

parse_int_vec = function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])

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
  } else if (args[i] == "--D" && i < length(args)) {
    D = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-grid" && i < length(args)) {
    n_grid = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-intervals" && i < length(args)) {
    n_intervals = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--fail-fast" && i < length(args)) {
    fail_fast = parse_flag(args[i + 1L]); i = i + 2L
  } else {
    i = i + 1L
  }
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

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

run_split = function(dat, method) {
  if (identical(method, "regional_pdp")) {
    strat = PdStrategy$new()
    tree = GadgetTree$new(strategy = strat, n_split = 2L, min_node_size = 50L)
    tree$fit(
      data = dat,
      target_feature_name = "y",
      model = "toy",
      predict_fun = toy_pred_fun,
      n_grid = n_grid,
      pd_engine = "cpp"
    )
  } else {
    strat = AleStrategy$new()
    tree = GadgetTree$new(strategy = strat, n_split = 2L, min_node_size = 50L)
    tree$fit(
      data = dat,
      target_feature_name = "y",
      model = "toy",
      predict_fun = toy_pred_fun,
      n_intervals = n_intervals,
      order_method = "raw",
      ale_engine = "cpp"
    )
  }
  strat$fit_timing$regional
}

record_row = function(rows, method, N, time_sec, repetition, status = "ok", error_message = NA_character_) {
  rows[[length(rows) + 1L]] = data.frame(
    module = "split_search",
    package = "gadget",
    impl = "split",
    method = method,
    model_type = "toy",
    N = N,
    D = D,
    n_grid = if (grepl("pdp", method)) n_grid else NA_integer_,
    n_intervals = if (grepl("ale", method)) n_intervals else NA_integer_,
    repetition = repetition,
    time_sec = time_sec,
    status = status,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
  rows
}

rows = list()
for (N in N_vec) {
  dat = load_data(N, D)
  if (is.null(dat)) next
  message(sprintf("=== gadget split: N=%d, D=%d ===", N, D))
  for (method in c("regional_pdp", "regional_ale")) {
    message("  ", method)
    tryCatch(
      run_split(dat, method),
      error = function(e) {
        if (isTRUE(fail_fast)) stop(e)
        message("    warmup skipped/failed: ", conditionMessage(e))
        invisible(NA_real_)
      }
    )
    for (r in seq_len(reps)) {
      out = tryCatch(
        list(ok = TRUE, time = run_split(dat, method), error = NA_character_),
        error = function(e) {
          if (isTRUE(fail_fast)) stop(e)
          list(ok = FALSE, time = NA_real_, error = conditionMessage(e))
        }
      )
      rows = record_row(
        rows = rows,
        method = method,
        N = N,
        time_sec = out$time,
        repetition = r,
        status = if (out$ok) "ok" else "error",
        error_message = out$error
      )
    }
  }
}

out = if (length(rows)) do.call(rbind, rows) else data.frame()
write.csv(out, file.path(outdir, "split_search_gadget_toy.csv"), row.names = FALSE)
message("Written: split_search_gadget_toy.csv")
