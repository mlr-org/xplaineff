#!/usr/bin/env Rscript
# Diagnose completeness of the global PDP/ALE runtime benchmark.

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

args = commandArgs(trailingOnly = TRUE)
run_id = format(Sys.time(), "%Y%m%d_%H%M%S")
run_root = file.path("simulation/results/runtime_runs", run_id)
indir = file.path(run_root, "global_r_runtime")
outdir = indir
reps = 20L
N_vec = c(1000L, 5000L, 10000L, 20000L)
D_vec = c(10L, 20L, 50L, 100L)
fixed_N = 10000L
fixed_D = 20L
default_n_grid = 20L
default_n_intervals = 20L
n_grid_vec = c(10L, 20L, 50L)
n_intervals_vec = c(10L, 20L, 50L)
model_types = c("rf", "toy")
sub_experiments = c("vs_N", "vs_D", "vs_res")
include_mlr3 = FALSE

parse_flag = function(x) {
  s = tolower(trimws(as.character(x)))
  if (s %in% c("true", "1", "yes", "y", "on")) return(TRUE)
  if (s %in% c("false", "0", "no", "n", "off")) return(FALSE)
  FALSE
}

parse_int_vec = function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
parse_chr_vec = function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1L]])

i = 1L
while (i <= length(args)) {
  if (args[i] == "--indir" && i < length(args)) {
    indir = args[i + 1L]; i = i + 2L
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
  } else if (args[i] == "--n-grid" && i < length(args)) {
    default_n_grid = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-intervals" && i < length(args)) {
    default_n_intervals = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-grid-vec" && i < length(args)) {
    n_grid_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--n-intervals-vec" && i < length(args)) {
    n_intervals_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--models" && i < length(args)) {
    model_types = parse_chr_vec(args[i + 1L])
    model_types = model_types[nzchar(model_types)]
    i = i + 2L
  } else if (args[i] == "--sub-experiments" && i < length(args)) {
    sub_experiments = parse_chr_vec(args[i + 1L])
    sub_experiments = sub_experiments[nzchar(sub_experiments)]
    i = i + 2L
  } else if (args[i] == "--include-mlr3" && i < length(args)) {
    include_mlr3 = parse_flag(args[i + 1L]); i = i + 2L
  } else {
    i = i + 1L
  }
}

if (isTRUE(include_mlr3) && !"mlr3_rf" %in% model_types) {
  model_types = c(model_types, "mlr3_rf")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

library(data.table)
source("simulation/global_runtime_io.R")

setDTthreads(1L)

method_specs = list(
  list(package = "xplaineff", impl = "r", method = "global_pdp", model_types = c("rf", "toy", "mlr3_rf")),
  list(package = "xplaineff", impl = "cpp", method = "global_pdp", model_types = c("rf", "toy", "mlr3_rf")),
  list(package = "xplaineff", impl = "r", method = "global_ale", model_types = c("rf", "toy", "mlr3_rf")),
  list(package = "xplaineff", impl = "cpp", method = "global_ale", model_types = c("rf", "toy", "mlr3_rf")),
  list(package = "pdp", impl = "default", method = "global_pdp", model_types = c("rf", "toy")),
  list(package = "iml", impl = "ice", method = "global_pdp", model_types = c("rf", "toy")),
  list(package = "iml", impl = "ale", method = "global_ale", model_types = c("rf", "toy")),
  list(package = "DALEX/ingredients", impl = "partial_dependence", method = "global_pdp", model_types = c("rf", "toy")),
  list(package = "DALEX/ingredients", impl = "accumulated_dependence", method = "global_ale", model_types = c("rf", "toy")),
  list(package = "ale", impl = "point_estimate", method = "global_ale", model_types = c("rf", "toy")),
  list(package = "effectplots", impl = "default", method = "global_pdp", model_types = c("rf", "toy")),
  list(package = "effectplots", impl = "default", method = "global_ale", model_types = c("rf", "toy"))
)

make_cells = function(method) {
  is_pdp = grepl("pdp", method)
  res_vec = if (is_pdp) n_grid_vec else n_intervals_vec
  make_cell = function(sub_experiment, N, D, res) {
    data.table(
      sub_experiment = sub_experiment,
      N = as.integer(N),
      D = as.integer(D),
      n_grid = if (is_pdp) as.integer(res) else NA_integer_,
      n_intervals = if (!is_pdp) as.integer(res) else NA_integer_
    )
  }

  cells = list()
  if ("vs_N" %in% sub_experiments) {
    cells[[length(cells) + 1L]] = rbindlist(lapply(N_vec, function(N) {
      make_cell("vs_N", N, fixed_D, if (is_pdp) default_n_grid else default_n_intervals)
    }))
  }
  if ("vs_D" %in% sub_experiments) {
    cells[[length(cells) + 1L]] = rbindlist(lapply(D_vec, function(D) {
      make_cell("vs_D", fixed_N, D, if (is_pdp) default_n_grid else default_n_intervals)
    }))
  }
  if ("vs_res" %in% sub_experiments) {
    cells[[length(cells) + 1L]] = rbindlist(lapply(res_vec, function(res) {
      make_cell("vs_res", fixed_N, fixed_D, res)
    }))
  }

  if (!length(cells)) return(data.table())
  rbindlist(cells)
}

expected_parts = list()
for (spec in method_specs) {
  applicable_models = intersect(model_types, spec$model_types)
  for (model_type in applicable_models) {
    cells = make_cells(spec$method)
    if (!nrow(cells)) next
    cells[, `:=`(
      module = "global_r",
      package = spec$package,
      impl = spec$impl,
      method = spec$method,
      model_type = model_type,
      expected_reps = reps
    )]
    expected_parts[[length(expected_parts) + 1L]] = cells
  }
}

expected = rbindlist(expected_parts, use.names = TRUE, fill = TRUE)
cell_key_cols = c(
  "module", "package", "impl", "method", "model_type", "sub_experiment", "N", "D", "n_grid", "n_intervals"
)
setcolorder(expected, c(cell_key_cols, "expected_reps"))
setorderv(expected, cell_key_cols)

dt = load_global_runtime_data(indir = indir, model_types = model_types, include_mlr3 = include_mlr3)
if (nrow(dt)) {
  if (!("error_message" %in% names(dt))) {
    dt[, error_message := NA_character_]
  }
  dt[, time_sec := as.numeric(time_sec)]
  collapse_values = function(x, n = 3L) {
    x = sort(unique(trimws(as.character(x[!is.na(x) & x != ""]))))
    if (!length(x)) return("")
    suffix = if (length(x) > n) " | ..." else ""
    paste0(paste(head(x, n), collapse = " | "), suffix)
  }
  observed = dt[, .(
    row_count = .N,
    observed_reps = uniqueN(repetition[!is.na(repetition)]),
    ok_reps = uniqueN(repetition[status == "ok" & is.finite(time_sec)]),
    error_reps = uniqueN(repetition[status == "error"]),
    skipped_reps = uniqueN(repetition[status == "skipped"]),
    other_statuses = collapse_values(setdiff(status, c("ok", "error", "skipped"))),
    source_files = collapse_values(source_file, n = 12L),
    error_messages = collapse_values(error_message[status == "error"])
  ), by = cell_key_cols]
} else {
  observed = data.table()
}

diagnostics = merge(expected, observed, by = cell_key_cols, all.x = TRUE)
count_cols = c("row_count", "observed_reps", "ok_reps", "error_reps", "skipped_reps")
for (col in count_cols) {
  diagnostics[is.na(get(col)), (col) := 0L]
  diagnostics[, (col) := as.integer(get(col))]
}
text_cols = c("other_statuses", "source_files", "error_messages")
for (col in text_cols) {
  diagnostics[is.na(get(col)), (col) := ""]
}

diagnostics[, missing_reps := pmax(expected_reps - observed_reps, 0L)]
diagnostics[, missing_ok_reps := pmax(expected_reps - ok_reps, 0L)]
diagnostics[, summary_point_status := fcase(
  ok_reps == expected_reps, "complete",
  ok_reps > 0L, "incomplete",
  default = "absent"
)]
diagnostics[, cell_status := fcase(
  observed_reps == 0L, "missing",
  ok_reps == expected_reps, "complete",
  error_reps > 0L & missing_reps > 0L, "partial_error",
  error_reps > 0L, "error",
  skipped_reps > 0L & missing_reps > 0L, "partial_skipped",
  skipped_reps > 0L, "skipped",
  missing_reps > 0L, "partial",
  default = "incomplete"
)]
setorderv(diagnostics, cell_key_cols)

summary_dt = diagnostics[, .(
  expected_cells = .N,
  complete_cells = sum(cell_status == "complete"),
  missing_cells = sum(cell_status == "missing"),
  partial_cells = sum(grepl("^partial", cell_status)),
  error_cells = sum(cell_status %in% c("error", "partial_error")),
  skipped_cells = sum(cell_status %in% c("skipped", "partial_skipped")),
  incomplete_summary_points = sum(summary_point_status == "incomplete"),
  absent_summary_points = sum(summary_point_status == "absent"),
  expected_reps = sum(expected_reps),
  ok_reps = sum(ok_reps),
  missing_reps = sum(missing_reps),
  error_reps = sum(error_reps),
  skipped_reps = sum(skipped_reps)
), by = .(model_type, sub_experiment, package, impl, method)]
setorderv(summary_dt, c("model_type", "sub_experiment", "package", "impl", "method"))

problem_dt = diagnostics[cell_status != "complete"]
missing_dt = diagnostics[observed_reps == 0L | missing_reps > 0L]

cells_file = file.path(outdir, "global_r_diagnostics_cells.csv")
summary_file = file.path(outdir, "global_r_diagnostics_summary.csv")
problem_file = file.path(outdir, "global_r_diagnostics_problems.csv")
missing_file = file.path(outdir, "global_r_diagnostics_missing_cells.csv")

fwrite(diagnostics, cells_file)
fwrite(summary_dt, summary_file)
fwrite(problem_dt, problem_file)
fwrite(missing_dt, missing_file)

message("Written: ", cells_file)
message("Written: ", summary_file)
message("Written: ", problem_file)
message("Written: ", missing_file)

overall = diagnostics[, .(
  expected_cells = .N,
  complete_cells = sum(cell_status == "complete"),
  problem_cells = sum(cell_status != "complete"),
  expected_reps = sum(expected_reps),
  ok_reps = sum(ok_reps),
  error_reps = sum(error_reps),
  skipped_reps = sum(skipped_reps),
  missing_reps = sum(missing_reps)
)]
print(overall)
