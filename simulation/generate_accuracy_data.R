#!/usr/bin/env Rscript
# Generate shared datasets for the structural recovery benchmark.
# Run from the package root.

args = commandArgs(trailingOnly = TRUE)
outdir = "simulation/data/accuracy"
n_seeds = 30L
N_vec = c(200L, 500L, 1000L, 5000L)
D_vec = c(5L, 10L, 20L)
variants = c("num_0", "num_04", "cat")
datadir_seen = FALSE
noise_scale = 0.3

parse_int_vec = function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
parse_chr_vec = function(x) strsplit(x, ",", fixed = TRUE)[[1L]]

i = 1L
while (i <= length(args)) {
  if (args[i] == "--datadir" && i < length(args)) {
    outdir = args[i + 1L]
    datadir_seen = TRUE
    i = i + 2L
  } else if (args[i] == "--outdir" && i < length(args) && !datadir_seen) {
    outdir = args[i + 1L]
    i = i + 2L
  } else if (args[i] == "--n-seeds" && i < length(args)) {
    n_seeds = as.integer(args[i + 1L])
    i = i + 2L
  } else if (args[i] == "--N-vec" && i < length(args)) {
    N_vec = parse_int_vec(args[i + 1L])
    i = i + 2L
  } else if (args[i] == "--D-vec" && i < length(args)) {
    D_vec = parse_int_vec(args[i + 1L])
    i = i + 2L
  } else if (args[i] == "--variants" && i < length(args)) {
    variants = parse_chr_vec(args[i + 1L])
    i = i + 2L
  } else {
    i = i + 1L
  }
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

variant_threshold = function(variant) {
  if (variant == "num_0") {
    return(0)
  }
  if (variant == "num_04") {
    return(0.4)
  }
  NA_real_
}

make_noise_feature = function(N, j) {
  if (j %% 3L == 1L) {
    return(round(stats::rnorm(N, mean = j / 2, sd = 1 + j / 10), 3))
  }
  if (j %% 3L == 2L) {
    return(round(stats::runif(N, min = -1, max = 1), 3))
  }
  stats::rbinom(N, size = 1L, prob = 0.35)
}

make_data = function(N, D, seed, variant) {
  set.seed(seed)

  dat = data.frame(
    x1 = round(stats::runif(N, -1, 1), 1),
    x2 = round(stats::runif(N, -1, 1), 3)
  )

  if (variant == "cat") {
    dat$x3 = factor(
      ifelse(stats::rbinom(N, size = 1L, prob = 0.5) == 1L, "1", "0"),
      levels = c("0", "1")
    )
  } else {
    dat$x3 = stats::runif(N, -1, 1)
  }

  if (D >= 4L) {
    for (j in 4:D) {
      dat[[paste0("x", j)]] = make_noise_feature(N, j)
    }
  }

  moderator_x1 = as.numeric(dat$x1 > 0)
  moderator_x3 = if (variant == "cat") {
    as.numeric(as.character(dat$x3) == "0")
  } else {
    as.numeric(dat$x3 <= variant_threshold(variant))
  }

  coef_x2 = -8 + 9 * moderator_x3 + 8 * moderator_x1
  y_det = 0.2 * dat$x1 + coef_x2 * dat$x2
  sd_eps = noise_scale * stats::sd(y_det)
  if (!is.finite(sd_eps) || sd_eps <= 0) {
    sd_eps = 0.01
  }

  dat$y = y_det + stats::rnorm(N, mean = 0, sd = sd_eps)
  dat
}

n_written = 0L
for (variant in variants) {
  for (N in N_vec) {
    for (D in D_vec) {
      for (s in seq_len(n_seeds)) {
        seed = 1000L + s
        dat = make_data(N, D, seed, variant)
        fn = file.path(outdir, sprintf("acc_N%d_D%d_%s_seed%d.csv", N, D, variant, seed))
        write.csv(dat, fn, row.names = FALSE)
        n_written = n_written + 1L
      }
    }
  }
}

message("Written ", n_written, " files to ", outdir)
