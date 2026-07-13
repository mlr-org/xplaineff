### Shared imports and helpers for the 26-07-13 selective-early-stopping examples.
### Source this from example_CB_CC_26-07-13_first_test.R and
### example_CB_CC_26-07-13_first_benchmark_noise.R.

library(xplaineff)
library(iml)
library(mlr3)
library(mlr3learners)
library(ranger)

# Train a ranger model on `dat` (target "y") and return its ICE feature effects.
make_effect = function(dat) {
  feats = setdiff(colnames(dat), "y")
  task = TaskRegr$new("syn", backend = dat, target = "y")
  learner = lrn("regr.ranger")
  learner$train(task)
  predictor = Predictor$new(learner, data = dat[, feats], y = dat$y)
  FeatureEffects$new(predictor, grid.size = 50, method = "ice")
}

# Fit a PD GADGET tree, optionally with selective early stopping.
fit_tree = function(dat, effect, gadget_improvements = NULL, tau = NULL,
  n_split = 4, impr_par = 0.05, min_node_size = 30) {
  tr = GadgetTree$new(strategy = PdStrategy$new(), n_split = n_split,
    impr_par = impr_par, min_node_size = min_node_size)
  tr$fit(data = dat, target_feature_name = "y", effect = effect,
    gadget_improvements = gadget_improvements,
    gadget_impr_args = if (is.null(tau)) NULL else list(tau = tau))
  tr
}

# Walk the tree, printing each node's still-interacting features (vecb_remaining_features).
walk_remaining = function(node, prefix = "root") {
  if (is.null(node)) return(invisible())
  rf = node$vecb_remaining_features
  tag = if (is.null(rf)) "NULL (disabled)" else paste(names(rf)[rf], collapse = ",")
  cat(sprintf("  %-9s id=%-2d depth=%d  remaining: %s\n", prefix, node$id, node$depth, tag))
  if (!is.null(node$children)) {
    walk_remaining(node$children$left_child, paste0(prefix, ".L"))
    walk_remaining(node$children$right_child, paste0(prefix, ".R"))
  }
}

split_cols = c("depth", "id", "node_type", "split_feature", "split_value", "int_imp")
show_tree = function(tr) print(tr$extract_split_info()[, split_cols])
