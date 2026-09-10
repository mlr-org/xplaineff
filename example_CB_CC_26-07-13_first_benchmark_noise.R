### Benchmark: how does selective early stopping shrink the GADGET tree as we vary
### the observation noise and the early-stopping threshold tau?
### 2D grid (noise level x tau) -> number of nodes.

source("example_CB_CC_26-07-13_input.R")

cat("\nTook ~ 8 min on 23rd Aug 26\n")

noise_levels = c(0.001, 0.1, 0.2, 0.5, 2.0)
tau_values = c(0.005, 0.05, 0.5, 0.95)

# Column 1: disabled baseline; then every (method, tau) combination.
configs = c(
  list(list(improvement_method = NULL, tau = NULL, label = "disabled")),
  unlist(lapply(early_stopping_methods, function(method) {
    lapply(tau_values, function(tau) list(improvement_method = method, tau = tau,
      label = paste0(method, "|tau=", tau)))
  }), recursive = FALSE)
)

n_nodes = matrix(NA_integer_, nrow = length(noise_levels), ncol = length(configs),
  dimnames = list(noise = paste0("sd=", noise_levels),
    config = vapply(configs, function(cfg) cfg$label, character(1))))
root_remaining = n_nodes  # same shape, holds the root's kept-feature summary
storage.mode(root_remaining) = "character"

# Set to TRUE to additionally keep the full extract_split_info() of every single run
# (tree structure, remaining features per node, the per-feature early stopping statistic
# early_stopping_stat_*, int_imp / int_imp_remaining and int_imp_* per feature). Stored in
# `detailed_results[[noise]][[config]]` and printed at the end.
detailed = FALSE
detailed_results = list()

set.seed(42)
n = 2000
dat3 = data.frame(
  x1 = runif(n, -1, 1), x2 = runif(n, -1, 1), x3 = runif(n, -1, 1),
  x4 = runif(n, -1, 1), x5 = runif(n, -1, 1)
)
# Two independent interactions (x1:x3, x2:x4) plus main effects; x5 is pure noise.
signal = ifelse(dat3$x3 > 0, 3 * dat3$x1, -3 * dat3$x1) +
  ifelse(dat3$x4 > 0, 2 * dat3$x2, -2 * dat3$x2) + dat3$x3 -
  dat3$x2 + 0.5 * dat3$x1 + 2 * dat3$x4
for (i in seq_along(noise_levels)) {
  dat3$y = signal + rnorm(n, sd = noise_levels[i])
  effect = make_effect(dat3)
  if (detailed) detailed_results[[rownames(n_nodes)[i]]] = list()
  for (j in seq_along(configs)) {
    tree = fit_tree(dat3, effect, configs[[j]]$improvement_method, tau = configs[[j]]$tau)
    split_info = tree$extract_split_info()
    n_nodes[i, j] = nrow(split_info)
    remaining_features = tree$root$vecb_remaining_features
    root_remaining[i, j] = if (is.null(remaining_features)) "all" else
      paste(names(remaining_features)[remaining_features], collapse = ",")
    if (detailed) detailed_results[[rownames(n_nodes)[i]]][[configs[[j]]$label]] = split_info
  }
}

cat("\n===== number of nodes (rows = noise sd, cols = tau) =====\n")
print(n_nodes)
cat("\n===== features kept at the root =====\n")
print(root_remaining)

# Full per-run trees, only when requested.
if (detailed) {
  detail_cols = function(split_info) {
    keep = c("depth", "id", "n_obs", "split_feature", "split_value", "remaining_features",
      "node_objective", "node_objective_remaining", "int_imp", "int_imp_remaining",
      grep("^(int_imp_|early_stopping_stat_)", colnames(split_info), value = TRUE))
    split_info[, intersect(keep, colnames(split_info)), drop = FALSE]
  }
  for (noise_label in names(detailed_results)) {
    for (config_label in names(detailed_results[[noise_label]])) {
      cat(sprintf("\n----- %s | %s -----\n", noise_label, config_label))
      print(detail_cols(detailed_results[[noise_label]][[config_label]]), row.names = FALSE)
    }
  }
}

# One heatmap per method over the (noise x tau) grid. The disabled baseline is repeated as the
# leftmost column of every plot, so each method can be read against it without flipping back.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  config_method = vapply(configs, function(cfg)
    if (is.null(cfg$improvement_method)) "disabled" else cfg$improvement_method, character(1))
  config_tau = vapply(configs, function(cfg)
    if (is.null(cfg$tau)) NA_real_ else cfg$tau, numeric(1))

  results = data.frame(
    noise = rep(rownames(n_nodes), times = length(configs)),
    method = rep(config_method, each = nrow(n_nodes)),
    tau = rep(config_tau, each = nrow(n_nodes)),
    n_nodes = as.vector(n_nodes),
    stringsAsFactors = FALSE
  )
  baseline = results[results$method == "disabled", ]
  tau_levels = c("off (baseline)", format(tau_values))

  for (method in early_stopping_methods) {
    df = rbind(baseline, results[results$method == method, ])
    df$tau_label = factor(ifelse(is.na(df$tau), "off (baseline)", format(df$tau)), levels = tau_levels)
    df$noise = factor(df$noise, levels = rownames(n_nodes))
    p = ggplot(df, aes(x = tau_label, y = noise, fill = n_nodes)) +
      geom_tile(color = "white") +
      geom_text(aes(label = n_nodes)) +
      scale_fill_gradient(low = "#f7fbff", high = "#08519c") +
      labs(title = paste("GADGET tree size (# nodes) under selective early stopping:", method),
        x = "early-stopping threshold tau", y = "observation noise (sd)") +
      theme_minimal()
    print(p)
  }
}
