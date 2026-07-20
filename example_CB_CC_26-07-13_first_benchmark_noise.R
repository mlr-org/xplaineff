### Benchmark: how does selective early stopping shrink the GADGET tree as we vary
### the observation noise and the early-stopping threshold tau?
### 2D grid (noise level x tau) -> number of nodes.

source("example_CB_CC_26-07-13_input.R")

noise_levels = c(0.001, 0.1, 0.2, 0.5, 2.0)
tau_values = c(0.005, 0.05, 0.5, 0.95)
# Selective early stopping methods (Section 5.1), Methods 1-4.
methods = c("plain_risk", "risk_reduction", "interaction_fraction", "interaction_fraction_total")

# Column 1: disabled baseline; then every (method, tau) combination.
configs = c(
  list(list(gi = NULL, tau = NULL, label = "disabled")),
  unlist(lapply(methods, function(m) {
    lapply(tau_values, function(tt) list(gi = m, tau = tt, label = paste0(m, "|tau=", tt)))
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
  eff = make_effect(dat3)
  if (detailed) detailed_results[[rownames(n_nodes)[i]]] = list()
  for (j in seq_along(configs)) {
    tr = fit_tree(dat3, eff, configs[[j]]$gi, tau = configs[[j]]$tau)
    si = tr$extract_split_info()
    n_nodes[i, j] = nrow(si)
    rf = tr$root$vecb_remaining_features
    root_remaining[i, j] = if (is.null(rf)) "all" else paste(names(rf)[rf], collapse = ",")
    if (detailed) detailed_results[[rownames(n_nodes)[i]]][[configs[[j]]$label]] = si
  }
}

cat("\n===== number of nodes (rows = noise sd, cols = tau) =====\n")
print(n_nodes)
cat("\n===== features kept at the root =====\n")
print(root_remaining)

# Full per-run trees, only when requested.
if (detailed) {
  detail_cols = function(si) {
    keep = c("depth", "id", "n_obs", "split_feature", "split_value", "remaining_features",
      "node_objective", "node_objective_remaining", "int_imp", "int_imp_remaining",
      grep("^(int_imp_|early_stopping_stat_)", colnames(si), value = TRUE))
    si[, intersect(keep, colnames(si)), drop = FALSE]
  }
  for (noise_label in names(detailed_results)) {
    for (config_label in names(detailed_results[[noise_label]])) {
      cat(sprintf("\n----- %s | %s -----\n", noise_label, config_label))
      print(detail_cols(detailed_results[[noise_label]][[config_label]]), row.names = FALSE)
    }
  }
}

# Optional heatmap of tree size across the (noise, method x tau) grid.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  df = expand.grid(noise = rownames(n_nodes), config = colnames(n_nodes),
    stringsAsFactors = FALSE)
  df$n_nodes = as.vector(n_nodes)
  df$config = factor(df$config, levels = colnames(n_nodes))
  p = ggplot(df, aes(x = config, y = noise, fill = n_nodes)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_nodes)) +
    scale_fill_gradient(low = "#f7fbff", high = "#08519c") +
    labs(title = "GADGET tree size under selective early stopping",
      x = "early-stopping method and threshold", y = "observation noise (sd)", fill = "# nodes") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  print(p)
}
