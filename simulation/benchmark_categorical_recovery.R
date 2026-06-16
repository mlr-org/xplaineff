#!/usr/bin/env Rscript
# Categorical split-feature recovery benchmark (GADGET-only).
#
# Question: can a GADGET-ALE tree recover a KNOWN level grouping of a categorical
# split feature, and how does order_method (raw/random/mds/pca) affect recovery?
#
# DGP: factor x3 with K levels.
#   - dgp_type = "digit": true grouping is a signal level set vs rest;
#     x2 slope flips sign across the two groups (heterogeneity lives in x2).
#   - dgp_type = "group": each level k has its own slope beta_k ~ U[-slope_mag, slope_mag];
#     no clean 2-partition. Signal levels are defined as the ceil(K/3) levels with the
#     largest |beta_k| so that exact_recovery / ARI metrics remain defined.
#   - In both DGPs, a covariate z leaks the binary "signal-vs-rest" indicator with
#     strength leakage in [0, 1]. z does NOT enter the model; it is a pure covariate.
#     leakage = 0 is the null case (mds/pca should collapse to random-order performance).
#   - factor level encoding order is shuffled per seed, so "raw" cannot exploit
#     a lexical coincidence.
#   - split_feature is forced to x3, isolating grouping recovery from the
#     feature-selection question (covered by the structural-recovery benchmark).
#
# Methods: raw / random / mds / pca + oracle_partition (exhaustive K-1 partition
#   enumeration, only enabled when K <= max_partition_K, default 12).
#
# Output: one row per (sweep cell, seed, method) written to a raw CSV; includes
#   elapsed_sec column for cost-vs-accuracy reporting.

if (!file.exists("DESCRIPTION") || readLines("DESCRIPTION", 1L) != "Package: gadget") {
  if (file.exists("../DESCRIPTION") && readLines("../DESCRIPTION", 1L) == "Package: gadget") {
    setwd("..")
  } else {
    stop("Run from the GADGET package root")
  }
}

suppressPackageStartupMessages({
  library(gadget)
  library(data.table)
})

args = commandArgs(trailingOnly = TRUE)

# ---- defaults ----------------------------------------------------------------
outdir = "simulation/results/categorical_recovery"
sweep = "leakage"        # one of: leakage, K, N, D, group_frac, slope_mag
n_seeds = 30L
cores = 1L
dgp_type = "digit"       # one of: digit, group
max_partition_K = 12L    # oracle_partition is only enumerated for K <= this

# base configuration (the non-swept axes are held at these values)
base_K = 6L
base_N = 2000L
base_D = 10L
base_leakage = 0.1
base_slope_mag = 4
group_frac = 1 / 3

# sweep grids
leakage_vec = c(0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.5, 1)
K_vec = c(6L, 8L, 12L, 20L)
N_vec = c(500L, 1000L, 2000L, 5000L)
D_vec = c(3L, 10L, 25L, 50L)
group_frac_vec = c(1 / 6, 1 / 3, 1 / 2)
slope_mag_vec = c(1, 2, 4, 8)

# fixed tree / DGP controls
foi_feature = "x2"
oracle_split_feature = "x3"
noise_scale = 0.3
n_intervals = 20L
n_split = 1L
min_node_size = 40L
impr_par = 0
order_methods = c("raw", "random", "mds", "pca")

parse_num_vec = function(x) as.numeric(strsplit(x, ",", fixed = TRUE)[[1L]])
parse_int_vec = function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])

i = 1L
while (i <= length(args)) {
  a = args[i]
  if (a == "--outdir" && i < length(args)) { outdir = args[i + 1L]; i = i + 2L }
  else if (a == "--sweep" && i < length(args)) { sweep = args[i + 1L]; i = i + 2L }
  else if (a == "--n-seeds" && i < length(args)) { n_seeds = as.integer(args[i + 1L]); i = i + 2L }
  else if (a == "--cores" && i < length(args)) { cores = as.integer(args[i + 1L]); i = i + 2L }
  else if (a == "--leakage-vec" && i < length(args)) { leakage_vec = parse_num_vec(args[i + 1L]); i = i + 2L }
  else if (a == "--K-vec" && i < length(args)) { K_vec = parse_int_vec(args[i + 1L]); i = i + 2L }
  else if (a == "--N-vec" && i < length(args)) { N_vec = parse_int_vec(args[i + 1L]); i = i + 2L }
  else if (a == "--D-vec" && i < length(args)) { D_vec = parse_int_vec(args[i + 1L]); i = i + 2L }
  else if (a == "--group-frac-vec" && i < length(args)) { group_frac_vec = parse_num_vec(args[i + 1L]); i = i + 2L }
  else if (a == "--slope-mag-vec" && i < length(args)) { slope_mag_vec = parse_num_vec(args[i + 1L]); i = i + 2L }
  else if (a == "--base-N" && i < length(args)) { base_N = as.integer(args[i + 1L]); i = i + 2L }
  else if (a == "--base-K" && i < length(args)) { base_K = as.integer(args[i + 1L]); i = i + 2L }
  else if (a == "--base-D" && i < length(args)) { base_D = as.integer(args[i + 1L]); i = i + 2L }
  else if (a == "--base-leakage" && i < length(args)) { base_leakage = as.numeric(args[i + 1L]); i = i + 2L }
  else if (a == "--base-slope-mag" && i < length(args)) { base_slope_mag = as.numeric(args[i + 1L]); i = i + 2L }
  else if (a == "--group-frac" && i < length(args)) { group_frac = as.numeric(args[i + 1L]); i = i + 2L }
  else if (a == "--dgp-type" && i < length(args)) { dgp_type = args[i + 1L]; i = i + 2L }
  else if (a == "--max-partition-K" && i < length(args)) { max_partition_K = as.integer(args[i + 1L]); i = i + 2L }
  else { i = i + 1L }
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

leakage_to_cohen_d = function(leakage) {
  ifelse(leakage >= 1, Inf, 2 * leakage / sqrt(1 - leakage^2))
}

# ---- DGP ---------------------------------------------------------------------
LEVEL_LABELS = LETTERS

make_noise_feature = function(N, j) {
  if (j %% 3L == 1L) return(round(stats::rnorm(N, 0, 1), 3))
  if (j %% 3L == 2L) return(round(stats::runif(N, -1, 1), 3))
  stats::rbinom(N, size = 1L, prob = 0.4)
}

# Returns: data.frame dat (with y), the signal level set, and the true
# binary group indicator g per row used by the leakage covariate z.
#
# In dgp_type = "digit": y = +/- slope_mag * x2, sign flipped between signal vs rest.
# In dgp_type = "group": y = beta_{x3} * x2 with per-level beta_k ~ U[-slope_mag, slope_mag];
#   "signal levels" = ceil(K/3) levels with the largest |beta_k|, used by metrics
#   and by z to define the binary leakage indicator g.
make_categorical_data = function(N, D, K, leakage, seed, signal_frac, slope_mag_value,
    dgp_type_value) {
  if (leakage < 0 || leakage > 1) {
    stop("leakage must be between 0 and 1")
  }
  if (!dgp_type_value %in% c("digit", "group")) {
    stop("dgp_type must be 'digit' or 'group'")
  }
  set.seed(seed)
  labels = LEVEL_LABELS[seq_len(K)]
  signal_level_size = max(1L, min(K - 1L, round(K * signal_frac)))
  level_order = sample(labels)  # shuffled encoding order -> "raw" is arbitrary

  x3_chr = sample(labels, N, replace = TRUE)
  x3 = factor(x3_chr, levels = level_order)

  if (dgp_type_value == "digit") {
    signal_levels = sample(labels, size = signal_level_size)
    g = as.integer(x3_chr %in% signal_levels)
    slope = ifelse(g == 1L, slope_mag_value, -slope_mag_value)
    betas = NULL
  } else {
    # Group DGP: each level gets an independent slope; signal levels are the ones
    # with the largest |beta_k| (so exact_recovery / ARI remain defined).
    betas = setNames(stats::runif(K, -slope_mag_value, slope_mag_value), labels)
    signal_levels = names(sort(abs(betas), decreasing = TRUE))[seq_len(signal_level_size)]
    g = as.integer(x3_chr %in% signal_levels)
    slope = betas[x3_chr]
  }

  z_noise = stats::rnorm(N)
  group_score = ifelse(g == 1L, 1, -1)
  z = leakage * group_score + sqrt(1 - leakage^2) * z_noise
  x2 = stats::runif(N, -1, 1)

  y_det = slope * x2
  sd_eps = noise_scale * stats::sd(y_det)
  if (!is.finite(sd_eps) || sd_eps <= 0) sd_eps = 0.01
  y = y_det + stats::rnorm(N, 0, sd_eps)

  dat = data.frame(x2 = x2, x3 = x3, z = z)
  # fill remaining columns with noise features to reach D total features
  # current feature count (excluding y): x2, x3, z = 3
  n_noise = max(0L, D - 3L)
  for (j in seq_len(n_noise)) {
    dat[[paste0("n", j)]] = make_noise_feature(N, j)
  }
  dat$y = y
  list(dat = dat, signal_levels = signal_levels, g = g, betas = betas)
}

# Oracle predictor.
#   digit DGP: basis (1, x2, x2 * g), g = 1[x3 in signal_levels].
#   group DGP: basis (1, x2, x2 * I[x3 = level]) for each level (level dummies on the
#     interaction). Both fit by OLS on the sampled data; isolates ordering recovery
#     from any model-fitting noise.
fit_oracle_predictor = function(dat, signal_levels, dgp_type_value) {
  if (dgp_type_value == "digit") {
    g = as.integer(as.character(dat$x3) %in% signal_levels)
    X = cbind(`(Intercept)` = 1, x2 = dat$x2, x2_g = dat$x2 * g)
    fit = stats::lm.fit(x = X, y = dat$y)
    coef = as.numeric(fit$coefficients)
    coef[is.na(coef)] = 0
    predict_fun = function(model, newdata) {
      gg = as.integer(as.character(newdata[["x3"]]) %in% model$signal_levels)
      model$coef[1L] + model$coef[2L] * newdata[["x2"]] + model$coef[3L] * (newdata[["x2"]] * gg)
    }
    return(list(model = list(coef = coef, signal_levels = signal_levels), predict_fun = predict_fun))
  }
  # group DGP: per-level interaction with x2
  lev = sort(unique(as.character(dat$x3)))
  Xint = sapply(lev, function(l) dat$x2 * (as.character(dat$x3) == l))
  X = cbind(`(Intercept)` = 1, x2 = dat$x2, Xint)
  colnames(X) = c("(Intercept)", "x2", paste0("x2_", lev))
  fit = stats::lm.fit(x = X, y = dat$y)
  coef = as.numeric(fit$coefficients)
  coef[is.na(coef)] = 0
  names(coef) = colnames(X)
  predict_fun = function(model, newdata) {
    pred = model$coef["(Intercept)"] + model$coef["x2"] * newdata[["x2"]]
    chr = as.character(newdata[["x3"]])
    for (l in model$lev) {
      cname = paste0("x2_", l)
      if (cname %in% names(model$coef)) {
        pred = pred + model$coef[cname] * newdata[["x2"]] * (chr == l)
      }
    }
    as.numeric(pred)
  }
  list(model = list(coef = coef, lev = lev, signal_levels = signal_levels), predict_fun = predict_fun)
}

# ---- metrics -----------------------------------------------------------------
# Adjusted Rand Index between two integer labelings of the same items.
adjusted_rand_index = function(a, b) {
  tab = table(a, b)
  sum_comb = function(x) sum(choose(x, 2))
  n = length(a)
  idx = sum_comb(as.vector(tab))
  a_i = sum_comb(rowSums(tab))
  b_j = sum_comb(colSums(tab))
  expected = a_i * b_j / choose(n, 2)
  maxv = (a_i + b_j) / 2
  if (maxv - expected == 0) return(1)
  (idx - expected) / (maxv - expected)
}

# Recover left/right level sets from the two root children via subset membership.
root_children_levels = function(tree, dat) {
  tl = tree$get_tree_list()
  if (length(tl) < 2L || length(tl[[2L]]) < 2L) {
    return(NULL)
  }
  lc = tl[[2L]][[1L]]
  rc = tl[[2L]][[2L]]
  if (is.null(lc$subset_idx) || is.null(rc$subset_idx)) return(NULL)
  list(
    left = sort(unique(as.character(dat$x3[lc$subset_idx]))),
    right = sort(unique(as.character(dat$x3[rc$subset_idx])))
  )
}

# Is the true signal-level grouping representable as a contiguous cut in the
# level order induced by the selected order_method?
order_correct = function(dat, signal_levels, order_method) {
  ord = tryCatch(
    gadget:::order_categorical_levels(
      x_cat = droplevels(dat$x3), data = dat, feature = "x3",
      target_feature_name = "y", order_method = order_method),
    error = function(e) NULL)
  if (is.null(ord)) return(NA)
  lev = levels(ord)
  m = length(signal_levels)
  pos = which(lev %in% signal_levels)
  setequal(pos, seq_len(m)) || setequal(pos, (length(lev) - m + 1L):length(lev))
}

record_row = function(method, order_method, K, N, D, leakage, slope_mag_value, dgp_type_value,
    seed, method_seed, signal_frac, dat, truth_g, signal_levels, tree, elapsed_sec = NA_real_) {
  si = tree$extract_split_info()
  root = si[si$node_type == "root" & !is.na(si$split_feature), ]
  int_imp = if (nrow(root)) as.numeric(root$int_imp[1L]) else NA_real_

  sides = root_children_levels(tree, dat)
  exact = NA; ari = NA_real_; node_acc = NA_real_
  if (!is.null(sides)) {
    all_levels = sort(unique(as.character(dat$x3)))
    rest = setdiff(all_levels, signal_levels)
    exact = (setequal(sides$left, signal_levels) && setequal(sides$right, rest)) ||
      (setequal(sides$left, rest) && setequal(sides$right, signal_levels))
    # level-level ARI: side label per level vs true group per level
    side_of_level = ifelse(all_levels %in% sides$left, 0L, 1L)
    grp_of_level = ifelse(all_levels %in% signal_levels, 1L, 0L)
    ari = adjusted_rand_index(side_of_level, grp_of_level)
    # observation node accuracy
    pred_left = as.character(dat$x3) %in% sides$left
    truth_left = truth_g == 1L
    node_acc = max(mean(pred_left == truth_left), mean(pred_left != truth_left))
  }

  data.table(
    method = method,
    order_method = order_method,
    sweep = sweep,
    dgp_type = dgp_type_value,
    K = K, N = N, D = D, leakage = leakage, cohen_d = leakage_to_cohen_d(leakage),
    slope_mag = slope_mag_value, seed = seed,
    method_seed = method_seed,
    group_frac = signal_frac,
    signal_level_size = length(signal_levels),
    status = if (!is.null(sides)) "ok" else "no_split",
    error_message = NA_character_,
    exact_recovery = exact,
    ari = ari,
    node_acc = node_acc,
    int_imp = int_imp,
    order_correct = if (order_method %in% c("raw", "random", "mds", "pca")) {
      order_correct(dat, signal_levels, order_method)
    } else {
      NA
    },
    elapsed_sec = elapsed_sec
  )
}

record_error_row = function(method, order_method, K, N, D, leakage, slope_mag_value, dgp_type_value,
    seed, method_seed, signal_frac, dat, signal_levels, error_message, elapsed_sec = NA_real_) {
  data.table(
    method = method,
    order_method = order_method,
    sweep = sweep,
    dgp_type = dgp_type_value,
    K = K, N = N, D = D, leakage = leakage, cohen_d = leakage_to_cohen_d(leakage),
    slope_mag = slope_mag_value, seed = seed,
    method_seed = method_seed,
    group_frac = signal_frac,
    signal_level_size = length(signal_levels),
    status = "fit_fail",
    error_message = error_message,
    exact_recovery = NA,
    ari = NA_real_,
    node_acc = NA_real_,
    int_imp = NA_real_,
    order_correct = if (order_method %in% c("raw", "random", "mds", "pca")) {
      order_correct(dat, signal_levels, order_method)
    } else {
      NA
    },
    elapsed_sec = elapsed_sec
  )
}

run_gadget_ale = function(dat, bundle, order_method) {
  tree = GadgetTree$new(strategy = AleStrategy$new(), n_split = n_split,
    impr_par = impr_par, min_node_size = min_node_size)
  tree$fit(data = dat, target_feature_name = "y", model = bundle$model,
    n_intervals = n_intervals, feature_set = foi_feature,
    split_feature = oracle_split_feature, predict_fun = bundle$predict_fun,
    order_method = order_method, ale_engine = "cpp")
  tree
}

# ---- oracle exhaustive partition ---------------------------------------------
# Enumerates every binary partition of the K factor levels (2^(K-1) - 1 candidates),
# fits a fresh ALE tree forced to that ordering by inserting the candidate level set
# as a contiguous prefix, and picks the partition with the largest int_imp on the
# root split. This gives raw/random/mds/pca an absolute upper bound. Only invoked
# when K <= max_partition_K because the candidate count grows as 2^(K-1).
enumerate_partitions = function(labels) {
  K = length(labels)
  # Each unordered binary partition {S, labels \ S} is represented uniquely by the
  # subset S that contains labels[1] and is a non-empty proper subset of labels.
  # That subset is determined by an arbitrary subset of labels[-1] (size 0..K-2),
  # which gives exactly 2^(K-1) - 1 partitions (we exclude the full set, which has
  # an empty complement).
  rest = labels[-1L]
  out = vector("list", 2L^(K - 1L) - 1L)
  idx = 1L
  for (m in seq_len(K - 1L)) {  # m = size of S, ranges 1..K-1
    cmb = utils::combn(length(rest), m - 1L)
    for (k in seq_len(ncol(cmb))) {
      out[[idx]] = c(labels[1L], rest[cmb[, k]])
      idx = idx + 1L
    }
  }
  out
}

run_oracle_partition = function(dat, bundle, signal_levels) {
  labels = sort(unique(as.character(dat$x3)))
  partitions = enumerate_partitions(labels)
  partitions = partitions[!sapply(partitions, is.null)]
  best_tree = NULL; best_int_imp = -Inf
  for (sub in partitions) {
    rest = setdiff(labels, sub)
    if (length(rest) == 0L) next
    new_order = c(sub, rest)
    dat_local = dat
    dat_local$x3 = factor(as.character(dat_local$x3), levels = new_order, ordered = TRUE)
    tree_try = tryCatch(run_gadget_ale(dat_local, bundle, order_method = "raw"),
      error = function(e) NULL)
    if (is.null(tree_try)) next
    si = tree_try$extract_split_info()
    root = si[si$node_type == "root" & !is.na(si$split_feature), ]
    if (nrow(root) == 0L) next
    int_imp_val = as.numeric(root$int_imp[1L])
    if (is.finite(int_imp_val) && int_imp_val > best_int_imp) {
      best_int_imp = int_imp_val
      best_tree = tree_try
    }
  }
  best_tree
}

# ---- build the list of cells for the requested sweep -------------------------
cells = list()
if (sweep == "leakage") {
  for (l in leakage_vec) {
    cells[[length(cells) + 1L]] = list(K = base_K, N = base_N, D = base_D, leakage = l,
      group_frac = group_frac, slope_mag = base_slope_mag)
  }
} else if (sweep == "K") {
  for (k in K_vec) {
    cells[[length(cells) + 1L]] = list(K = k, N = base_N, D = base_D, leakage = base_leakage,
      group_frac = group_frac, slope_mag = base_slope_mag)
  }
} else if (sweep == "N") {
  for (nn in N_vec) {
    cells[[length(cells) + 1L]] = list(K = base_K, N = nn, D = base_D, leakage = base_leakage,
      group_frac = group_frac, slope_mag = base_slope_mag)
  }
} else if (sweep == "D") {
  for (dd in D_vec) {
    cells[[length(cells) + 1L]] = list(K = base_K, N = base_N, D = dd, leakage = base_leakage,
      group_frac = group_frac, slope_mag = base_slope_mag)
  }
} else if (sweep == "group_frac") {
  for (gf in group_frac_vec) {
    cells[[length(cells) + 1L]] = list(K = base_K, N = base_N, D = base_D, leakage = base_leakage,
      group_frac = gf, slope_mag = base_slope_mag)
  }
} else if (sweep == "slope_mag") {
  for (sm in slope_mag_vec) {
    cells[[length(cells) + 1L]] = list(K = base_K, N = base_N, D = base_D, leakage = base_leakage,
      group_frac = group_frac, slope_mag = sm)
  }
} else {
  stop("Unknown --sweep: ", sweep)
}

run_seed = function(cell, seed) {
  gen = make_categorical_data(cell$N, cell$D, cell$K, cell$leakage, seed,
    cell$group_frac, cell$slope_mag, dgp_type)
  bundle = fit_oracle_predictor(gen$dat, gen$signal_levels, dgp_type)
  rows = list()
  # Heuristic ordering methods
  for (om in order_methods) {
    method_seed = seed * 100L + match(om, order_methods)
    set.seed(method_seed)
    timing = system.time({
      result = tryCatch(
        list(tree = run_gadget_ale(gen$dat, bundle, om), error_message = NA_character_),
        error = function(e) list(tree = NULL, error_message = conditionMessage(e))
      )
    })
    elapsed = as.numeric(timing["elapsed"])
    rows[[length(rows) + 1L]] = if (!is.null(result$tree)) {
      record_row("gadget_ale", om, cell$K, cell$N, cell$D,
        cell$leakage, cell$slope_mag, dgp_type, seed, method_seed, cell$group_frac,
        gen$dat, gen$g, gen$signal_levels, result$tree, elapsed_sec = elapsed)
    } else {
      record_error_row("gadget_ale", om, cell$K, cell$N, cell$D,
        cell$leakage, cell$slope_mag, dgp_type, seed, method_seed, cell$group_frac,
        gen$dat, gen$signal_levels, result$error_message, elapsed_sec = elapsed)
    }
  }
  # Oracle exhaustive partition (only when K small enough)
  if (cell$K <= max_partition_K) {
    method_seed = seed * 100L + length(order_methods) + 1L
    set.seed(method_seed)
    timing = system.time({
      tree_oracle = tryCatch(run_oracle_partition(gen$dat, bundle, gen$signal_levels),
        error = function(e) NULL)
    })
    elapsed = as.numeric(timing["elapsed"])
    rows[[length(rows) + 1L]] = if (!is.null(tree_oracle)) {
      record_row("gadget_ale", "oracle_partition", cell$K, cell$N, cell$D,
        cell$leakage, cell$slope_mag, dgp_type, seed, method_seed, cell$group_frac,
        gen$dat, gen$g, gen$signal_levels, tree_oracle, elapsed_sec = elapsed)
    } else {
      record_error_row("gadget_ale", "oracle_partition", cell$K, cell$N, cell$D,
        cell$leakage, cell$slope_mag, dgp_type, seed, method_seed, cell$group_frac,
        gen$dat, gen$signal_levels,
        error_message = "oracle_partition failed", elapsed_sec = elapsed)
    }
  }
  rbindlist(rows, fill = TRUE)
}

all_rows = list()
for (cell in cells) {
  message(sprintf("[sweep=%s dgp=%s] K=%d N=%d D=%d leakage=%g group_frac=%g slope_mag=%g  (%d seeds, %d cores)",
    sweep, dgp_type, cell$K, cell$N, cell$D, cell$leakage, cell$group_frac, cell$slope_mag,
    n_seeds, cores))
  seeds = 1000L + seq_len(n_seeds)
  per_seed = if (cores > 1L) {
    parallel::mclapply(seeds, function(s) run_seed(cell, s), mc.cores = cores)
  } else {
    lapply(seeds, function(s) run_seed(cell, s))
  }
  all_rows[[length(all_rows) + 1L]] = rbindlist(per_seed, fill = TRUE)
}

out = rbindlist(all_rows, fill = TRUE)
fout = file.path(outdir, sprintf("categorical_recovery_%s.csv", sweep))
fwrite(out, fout)
message("Written: ", fout, "  (", nrow(out), " rows)")
