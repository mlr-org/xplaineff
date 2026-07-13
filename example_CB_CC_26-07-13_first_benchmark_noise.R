### Benchmark: how does selective early stopping shrink the GADGET tree as we vary
### the observation noise and the early-stopping threshold tau?
### 2D grid (noise level x tau) -> number of nodes.

source("example_CB_CC_26-07-13_input.R")

noise_levels = c(0.001, 0.1, 0.2, 0.5, 2.0)
tau_values = c(0.005, 0.05, 0.5)

n_nodes = matrix(NA_integer_, nrow = length(noise_levels), ncol = length(tau_values) + 1,
  dimnames = list(noise = paste0("sd=", noise_levels),
    tau = c("disabled", paste0("tau=", tau_values))))
root_remaining = n_nodes  # same shape, holds the root's kept-feature summary
storage.mode(root_remaining) = "character"

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
# Column 1: disabled baseline; remaining columns: plain_risk for each tau.
configs = c(list(list(gi = NULL, tau = NULL)),
  lapply(tau_values, function(tt) list(gi = "plain_risk", tau = tt)))

for (i in seq_along(noise_levels)) {
  dat3$y = signal + rnorm(n, sd = noise_levels[i])
  eff = make_effect(dat3)
  for (j in seq_along(configs)) {
    tr = fit_tree(dat3, eff, configs[[j]]$gi, tau = configs[[j]]$tau)
    n_nodes[i, j] = nrow(tr$extract_split_info())
    rf = tr$root$vecb_remaining_features
    root_remaining[i, j] = if (is.null(rf)) "all" else paste(names(rf)[rf], collapse = ",")
  }
}

cat("\n===== number of nodes (rows = noise sd, cols = tau) =====\n")
print(n_nodes)
cat("\n===== features kept at the root =====\n")
print(root_remaining)

# Optional heatmap of tree size across the (noise, tau) grid.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  df = expand.grid(noise = rownames(n_nodes), tau = colnames(n_nodes),
    stringsAsFactors = FALSE)
  df$n_nodes = as.vector(n_nodes)
  df$tau = factor(df$tau, levels = colnames(n_nodes))
  p = ggplot(df, aes(x = tau, y = noise, fill = n_nodes)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_nodes)) +
    scale_fill_gradient(low = "#f7fbff", high = "#08519c") +
    labs(title = "GADGET tree size under selective early stopping",
      x = "early-stopping threshold", y = "observation noise (sd)", fill = "# nodes") +
    theme_minimal()
  print(p)
}
