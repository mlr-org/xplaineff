#!/usr/bin/env Rscript
# Generate shared benchmark data for the runtime sweeps.
# Defaults match the current global and regional runtime benchmarks.

args = commandArgs(trailingOnly = TRUE)
outdir = "simulation/data/global_r_runtime"
seed = 21L
N_vec = c(5000L, 10000L, 25000L, 50000L)
D_vec = c(10L, 20L, 50L, 100L)
fixed_N = 10000L
fixed_D = 20L
full_grid = FALSE

parse_int_vec = function(x) {
  as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
}

parse_flag = function(x) {
  s = tolower(trimws(as.character(x)))
  if (s %in% c("true", "1", "yes", "y", "on")) return(TRUE)
  if (s %in% c("false", "0", "no", "n", "off")) return(FALSE)
  FALSE
}

i = 1L
while (i <= length(args)) {
  if (args[i] == "--outdir" && i < length(args)) {
    outdir = args[i + 1L]; i = i + 2L
  } else if (args[i] == "--seed" && i < length(args)) {
    seed = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--N-vec" && i < length(args)) {
    N_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--D-vec" && i < length(args)) {
    D_vec = parse_int_vec(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--fixed-N" && i < length(args)) {
    fixed_N = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--fixed-D" && i < length(args)) {
    fixed_D = as.integer(args[i + 1L]); i = i + 2L
  } else if (args[i] == "--full-grid" && i < length(args)) {
    full_grid = parse_flag(args[i + 1L]); i = i + 2L
  } else {
    i = i + 1L
  }
}

make_cells = function() {
  if (isTRUE(full_grid)) {
    return(unique(expand.grid(N = N_vec, D = D_vec)))
  }
  unique(rbind(
    data.frame(N = N_vec, D = fixed_D),
    data.frame(N = fixed_N, D = D_vec)
  ))
}

generate_one = function(N, D, seed, outdir) {
  set.seed(seed)
  X = matrix(runif(N * D, -1, 1), nrow = N, ncol = D)
  colnames(X) = paste0("x", seq_len(D))
  x1 = X[, 1L]
  x2 = X[, 2L]
  x3 = X[, 3L]
  y_det = 5 * x1 + 5 * x2 + ifelse(x3 > 0, 10 * x1, 0) - ifelse(x3 > 0, 10 * x2, 0)
  eps = rnorm(N, 0, 0.1 * stats::sd(y_det))
  dat = data.frame(X, y = y_det + eps)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  fname = file.path(outdir, sprintf("benchmark_N%d_D%d_seed%d.csv", N, D, seed))
  write.csv(dat, fname, row.names = FALSE)
  message("Written: ", fname)
}

cells = make_cells()
for (i in seq_len(nrow(cells))) {
  generate_one(cells$N[i], cells$D[i], seed, outdir)
}
