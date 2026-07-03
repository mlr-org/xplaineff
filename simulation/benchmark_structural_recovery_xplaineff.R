#!/usr/bin/env Rscript
# Accuracy benchmark (xplaineff): recover the oracle first split for the x2 regional effect.

if (!file.exists("DESCRIPTION") || readLines("DESCRIPTION", 1L) != "Package: xplaineff") {
  if (file.exists("../DESCRIPTION") && readLines("../DESCRIPTION", 1L) == "Package: xplaineff") {
    setwd("..")
  } else {
    stop("Run from the xplaineff package root")
  }
}

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(xplaineff)
  }
  library(data.table)
})

args = commandArgs(trailingOnly = TRUE)
outdir = "simulation/results/structural_recovery"
datadir = "simulation/data/structural_recovery"
n_seeds = 30L
N_vec = c(200L, 500L, 1000L, 5000L)
D_vec = c(5L, 10L, 20L)
variants = c("num_0", "num_04", "cat")

foi_feature = "x2"
oracle_split_feature = "x3"
n_grid = 20L
n_intervals = 20L
n_split = 1L
min_node_size = 40L
impr_par = 0
n_quantiles = NULL

parse_int_vec = function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
parse_chr_vec = function(x) strsplit(x, ",", fixed = TRUE)[[1L]]

i = 1L
while (i <= length(args)) {
  if (args[i] == "--outdir" && i < length(args)) {
    outdir = args[i + 1L]
    i = i + 2L
  } else if (args[i] == "--datadir" && i < length(args)) {
    datadir = args[i + 1L]
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

load_accuracy_csv = function(datadir, N, D, variant, seed) {
  fn = file.path(datadir, sprintf("acc_N%d_D%d_%s_seed%d.csv", N, D, variant, seed))
  if (!file.exists(fn)) {
    stop("Missing dataset ", fn, ". Run: Rscript simulation/generate_structural_recovery_data.R", call. = FALSE)
  }
  dat = as.data.frame(data.table::fread(fn))
  if (variant == "cat" && !is.factor(dat$x3)) {
    dat$x3 = factor(as.character(dat$x3), levels = c("0", "1"))
  }
  dat
}

toy_pred_fun_factory = function(variant) {
  function(model, newdata) {
    moderator_x1 = as.numeric(newdata[["x1"]] > 0)
    moderator_x3 = if (variant == "cat") {
      as.numeric(as.character(newdata[["x3"]]) == "0")
    } else {
      as.numeric(newdata[["x3"]] <= variant_threshold(variant))
    }
    coef = model$coef
    coef[[1L]] +
      coef[[2L]] * newdata[["x1"]] +
      coef[[3L]] * newdata[["x2"]] +
      coef[[4L]] * (newdata[["x2"]] * moderator_x3) +
      coef[[5L]] * (newdata[["x2"]] * moderator_x1)
  }
}

fit_oracle_predictor = function(dat, variant) {
  moderator_x1 = as.numeric(dat$x1 > 0)
  moderator_x3 = if (variant == "cat") {
    as.numeric(as.character(dat$x3) == "0")
  } else {
    as.numeric(dat$x3 <= variant_threshold(variant))
  }

  X = cbind(
    `(Intercept)` = 1,
    x1 = dat$x1,
    x2 = dat$x2,
    x2_x3 = dat$x2 * moderator_x3,
    x2_x1 = dat$x2 * moderator_x1
  )
  fit = lm.fit(x = X, y = dat$y)
  coef = as.numeric(fit$coefficients)
  coef[is.na(coef)] = 0
  list(
    model = list(coef = coef, variant = variant),
    predict_fun = toy_pred_fun_factory(variant)
  )
}

true_left_mask = function(dat, variant) {
  if (variant == "cat") {
    return(as.character(dat$x3) == "0")
  }
  dat$x3 <= variant_threshold(variant)
}

root_split_info = function(split_info) {
  root_row = split_info[split_info$node_type == "root" & !is.na(split_info$split_feature), ]
  if (nrow(root_row) == 0L) {
    return(list(feature = NA_character_, value = NA))
  }
  list(
    feature = as.character(root_row$split_feature[1L]),
    value = root_row$split_value[1L]
  )
}

predicted_left_mask = function(dat, split_feature, split_value, method_label) {
  if (is.na(split_feature) || is.na(split_value)) {
    return(rep(NA, nrow(dat)))
  }

  z = dat[[split_feature]]
  if (is.factor(z)) {
    if (identical(method_label, "xplaineff_ale")) {
      return(xplaineff:::ordered_categorical_left_mask(z, split_value))
    }
    return(z == split_value)
  }

  z <= as.numeric(split_value)
}

record_xplaineff = function(method_label, variant, N, D, seed, dat, tree) {
  split_info = tree$extract_split_info()
  root_split = root_split_info(split_info)
  split_feature = root_split$feature
  split_value = root_split$value
  oracle_hit = !is.na(split_feature) && identical(split_feature, oracle_split_feature)
  node_acc = if (!is.na(split_feature) && !is.na(split_value)) {
    pred_left = predicted_left_mask(dat, split_feature, split_value, method_label)
    truth_left = true_left_mask(dat, variant)
    max(mean(pred_left == truth_left), mean((!pred_left) == truth_left))
  } else {
    NA_real_
  }

  data.table(
    package = "xplaineff",
    method = method_label,
    variant = variant,
    N = N,
    D = D,
    seed = seed,
    foi_feature = foi_feature,
    oracle_split_feature = oracle_split_feature,
    selected_split_feature = if (is.na(split_feature)) NA_character_ else split_feature,
    selected_split_value = if (is.na(split_value)) NA_character_ else as.character(split_value),
    split_feat_correct = oracle_hit,
    split_pt_error = if (variant == "cat" || !oracle_hit || is.na(split_value)) {
      NA_real_
    } else {
      abs(as.numeric(split_value) - variant_threshold(variant))
    },
    node_acc = node_acc
  )
}

run_xplaineff_pdp = function(dat, model, pred) {
  split_features = setdiff(colnames(dat), c("y", foi_feature))
  effect = xplaineff:::calculate_pd(
    model = model,
    data = dat,
    target_feature_name = "y",
    feature_set = foi_feature,
    predict_fun = pred,
    n_grid = n_grid,
    pd_engine = "cpp"
  )
  tree = GadgetTree$new(
    strategy = PdStrategy$new(),
    n_split = n_split,
    impr_par = impr_par,
    min_node_size = min_node_size,
    n_quantiles = n_quantiles
  )
  tree$fit(
    data = dat,
    target_feature_name = "y",
    effect = effect,
    feature_set = foi_feature,
    split_feature = split_features
  )
  tree
}

run_xplaineff_ale = function(dat, model, pred) {
  split_features = setdiff(colnames(dat), c("y", foi_feature))
  tree = GadgetTree$new(
    strategy = AleStrategy$new(),
    n_split = n_split,
    impr_par = impr_par,
    min_node_size = min_node_size,
    n_quantiles = n_quantiles
  )
  tree$fit(
    data = dat,
    target_feature_name = "y",
    model = model,
    n_intervals = n_intervals,
    feature_set = foi_feature,
    split_feature = split_features,
    predict_fun = pred,
    order_method = "raw",
    ale_engine = "cpp"
  )
  tree
}

rows = list()
for (variant in variants) {
  for (N in N_vec) {
    for (D in D_vec) {
      for (s in seq_len(n_seeds)) {
        seed = 1000L + s
        dat = load_accuracy_csv(datadir, N, D, variant, seed)
        pred_bundle = fit_oracle_predictor(dat, variant)

        tree_pdp = run_xplaineff_pdp(dat, pred_bundle$model, pred_bundle$predict_fun)
        rows[[length(rows) + 1L]] = record_xplaineff("xplaineff_pdp", variant, N, D, seed, dat, tree_pdp)

        tree_ale = run_xplaineff_ale(dat, pred_bundle$model, pred_bundle$predict_fun)
        rows[[length(rows) + 1L]] = record_xplaineff("xplaineff_ale", variant, N, D, seed, dat, tree_ale)
      }
    }
  }
}

out = rbindlist(rows, fill = TRUE)
fout = file.path(outdir, "structural_recovery_xplaineff.csv")
fwrite(out, fout)
message("Written: ", fout)
