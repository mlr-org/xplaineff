### Shared imports and helpers for the 26-07-13 selective-early-stopping examples.

library(xplaineff)
library(iml)
library(mlr3)
library(mlr3learners)
library(ranger)

# TODO: Support the usage of a simple, closed-formula function instead of a model

# Train a ranger model on some given data (target "y") and report its hold-out fit.
# `test_data_size` rows are held out (as an absolute count) purely to report an out-of-sample
# fit; the model is trained on the remaining rows.
# Split out from make_effect() so that variation studies which change only the ICE settings
# (grid size, which rows the curves cover) or only the tree settings (tau) can reuse the model
# instead of refitting an identical forest.
train_model = function(dat, test_data_size = 200) {
  test_idx = sample(nrow(dat), size = min(test_data_size, nrow(dat) / 2))
  train_data = dat[-test_idx, ]
  test_data = dat[test_idx, ]
  task = TaskRegr$new("syn", backend = train_data, target = "y")
  learner = lrn("regr.ranger")
  learner$train(task)
  prediction_test = learner$predict_newdata(test_data)$response
  rmse = sqrt(mean((test_data$y - prediction_test)^2))
  r_squared = 1 - sum((test_data$y - prediction_test)^2) /
    sum((test_data$y - mean(test_data$y))^2)
  cat(sprintf("\n\nTrained ranger model with n_test = %d and n_train = %d; performance on hold-out test data: RMSE = %.3f, R^2 = %.3f\n\n",
    nrow(test_data), nrow(train_data), rmse, r_squared))
  list(learner = learner, train_data = train_data, test_data = test_data,
    train_idx = setdiff(seq_len(nrow(dat)), test_idx), test_idx = sort(test_idx),
    rmse = rmse, r_squared = r_squared)
}

# Restrict an ICE effect to a subset of observations, keeping the grid fixed.
# This is exact: the raw ICE curve of observation i depends only on that observation and the
# model, and GADGET mean-centers per node anyway, so no full-data centering is baked in.
# Preferred over recomputing the curves on the subset, which would additionally re-derive the
# grid from the subset's feature range and thus confound "which rows" with "which grid".
subset_effect_rows = function(effect, rows) {
  parts = xplaineff:::mean_center_ice(effect, mean_center = FALSE)
  structure(
    list(Y = lapply(parts$Y, function(mat) mat[rows, , drop = FALSE]), grid = parts$grid),
    class = "xplaineff_pd_matrix"
  )
}

# ICE feature effects for an already trained model.
# `effect_data` selects the rows the curves are computed on: "all" keeps effects and GADGET tree
# on the same data (the default everywhere else), "train"/"test" separate them, which is how one
# sees whether the regional effects rest on structure the model overfitted on.
# The rows actually used come back as an attribute, because the tree must be fitted on exactly
# those rows -- the ICE matrices carry one row per observation.
ice_effects = function(model_fit, dat, grid.size = 50, effect_data = "all") {
  features = setdiff(colnames(dat), "y")
  effect_rows = switch(effect_data,
    all = dat, train = model_fit$train_data, test = model_fit$test_data)
  predictor = Predictor$new(model_fit$learner,
    data = effect_rows[, features], y = effect_rows$y)
  effects = FeatureEffects$new(predictor, grid.size = grid.size, method = "ice")
  attr(effects, "effect_rows") = effect_rows
  effects
}

# Train a model and return its ICE feature effects in one step (the common case).
make_effect = function(dat, grid.size = 50, test_data_size = 200, effect_data = "all") {
  ice_effects(train_model(dat, test_data_size), dat, grid.size = grid.size,
    effect_data = effect_data)
}

# Fit a PD GADGET tree, optionally with selective early stopping.
fit_tree = function(dat, effect, gadget_improvements = NULL,
  tau = NULL, n_split = 4, impr_par = 0.05, min_node_size = 30, verbose = 0) {
  tree = GadgetTree$new(strategy = PdStrategy$new(), n_split = n_split,
    impr_par = impr_par, min_node_size = min_node_size)
  tree$fit(data = dat, target_feature_name = "y", effect = effect,
    gadget_improvements = gadget_improvements,
    gadget_impr_args = if (is.null(tau)) NULL else list(tau = tau), verbose = verbose)
  tree
}

# Walk the tree, printing each node's still-interacting features (vecb_remaining_features).
walk_remaining = function(node, prefix = "root") {
  if (is.null(node)) return(invisible())
  remaining_features = node$vecb_remaining_features
  # This is null in case selective early stopping is disabled.
  tag = if (is.null(remaining_features)) "NULL (disabled)" else
    paste(names(remaining_features)[remaining_features], collapse = ",")
  cat(sprintf("  %-9s id=%-2d depth=%d  remaining: %s\n", prefix, node$id, node$depth, tag))
  if (!is.null(node$children)) {
    walk_remaining(node$children$left_child, paste0(prefix, ".L"))
    walk_remaining(node$children$right_child, paste0(prefix, ".R"))
  }
}

split_cols = c("depth", "id", "node_type", "split_feature", "split_value", "int_imp")
show_tree = function(tree) print(tree$extract_split_info()[, split_cols])

# Selective early stopping methods (Section 5.1), Methods 1-4.
early_stopping_methods = c("plain_risk", "risk_reduction", "interaction_fraction", "interaction_fraction_total")
